function mdOut = ladd_enrich_metadata(metadataCsv, varargin)
% LADD_ENRICH_METADATA  Fill known_u_ppm / known_th_ppm / known_sm_ppm from
% a standard lookup table, then write an enriched metadata CSV ready for
% the public parent-isotope reducers.
%
% The public parent reducers call this helper automatically when
% ReferenceLookupFile is supplied. It can also be run separately to inspect
% or save the enriched metadata. Existing non-blank values are never
% overwritten.
%
% USAGE
%   % Enrich in place (returns table, no file written unless saveAs given)
%   md = ladd_enrich_metadata('my_metadata.csv', ...
%       'lookupTable','reference_material_lookup.csv');
%
%   % Write an enriched CSV for review or later reuse
%   md = ladd_enrich_metadata('my_metadata.csv', ...
%       'lookupTable','reference_material_lookup.csv', ...
%       'saveAs','my_metadata_enriched.csv');
%   out = ladd_reduce_zircon('/run/folder/', ...
%       'my_metadata_enriched.csv','he_pit_volumes.csv', ...
%       'uth_pit_volumes.csv','BridgeStandardName','ReferenceMaterial', ...
%       'OutputFile','parents_reduced.csv');
%
%   % The lookup table is always user supplied. The public package does not
%   % embed laboratory reference-material concentrations.
%
% INPUTS
%   metadataCsv  : path to your metadata CSV, OR an already-loaded table
%
% OPTIONAL NAME-VALUE PAIRS
%   'lookupTable' : REQUIRED path to the user's standard lookup CSV
%   'saveAs'      : write enriched table to this path (default '' = no write)
%   'verbose'     : true/false — print matching summary (default true)
%   'runOrder'    : optional cell array of filename prefixes in chronological run order,
%                   e.g. {'RunA', 'RunB'}. Use this when two or more
%                   independent runs share a folder and their files have
%                   overlapping index numbers (_1…N each). Prefixes not listed
%                   are appended alphabetically after those specified.
%                   If omitted, the metadata row order is preserved.
%                   Example:
%                     md = ladd_enrich_metadata('meta.csv', ...
%                            'runOrder', {'RunA','RunB'}, ...
%                            'saveAs',   'meta_enriched.csv');
%
% LOOKUP TABLE FORMAT
%   stdname        — standard name (case-insensitive normalized match)
%   known_u_ppm       — U concentration in ppm
%   known_th_ppm      — Th concentration in ppm
%   known_sm_ppm      — Sm concentration on the basis documented by the user;
%                       pass the matching SmReferenceBasis to apatite reduction
%   known_u_1sd_ppm   — optional 1SD uncertainty on known_u_ppm
%   known_th_1sd_ppm  — optional 1SD uncertainty on known_th_ppm
%   known_sm_1sd_ppm  — optional 1SD uncertainty on known_sm_ppm
%   sm_reference_basis — recommended provenance column; not interpreted here
%   reference      — citation (informational only, ignored by code)
%
% MATCHING LOGIC
%   For each metadata row with a non-empty stdname, the code normalises both
%   the metadata stdname and each lookup table entry (lowercase, strip spaces,
%   underscores, and hyphens). Exact normalized matches are preferred. A
%   unique substring match is accepted for decorated run labels; ambiguous
%   matches stop with an error instead of silently choosing one standard.
%   Existing non-blank known_* values in the metadata are NEVER overwritten.
%
% EXAMPLES OF MATCHING
%   metadata stdname     lookup stdname    matches?
%   'NIST612'            'NIST612'         yes
%   'G_NIST612'          'NIST612'         yes  (contains 'nist612')
%   'nist612_g01'        'NIST612'         yes
%   'ReferenceMaterial'  'ReferenceMaterial' yes
%   'REF_A1'             'ReferenceMaterial' NO (use exact stdname in metadata)
%
% TIP: Use consistent stdname values in your metadata that contain the
% lookup table key as a substring. Prefer the exact lookup-table name.

% ── Parse arguments ───────────────────────────────────────────────────────
p = inputParser;
p.addParameter('lookupTable', '', @(x)ischar(x)||isstring(x));
p.addParameter('saveAs',      '', @(x)ischar(x)||isstring(x));
p.addParameter('verbose',     true, @islogical);
p.addParameter('runOrder',    {},   @(x) iscell(x) || isstring(x));
p.parse(varargin{:});
lookupPath = p.Results.lookupTable;
saveAs     = char(string(p.Results.saveAs));
verbose    = p.Results.verbose;
runOrder   = cellstr(p.Results.runOrder);   % normalise to cell array of char

% ── Require a user-supplied lookup table ─────────────────────────────────
if isempty(lookupPath)
    error('ladd_enrich_metadata:LookupRequired', ...
        ['The public package does not embed laboratory reference-material ' ...
         'concentrations. Supply an explicit lookupTable CSV.']);
end

% ── Read metadata ─────────────────────────────────────────────────────────
if ischar(metadataCsv) || isstring(metadataCsv)
    md = readtable(char(metadataCsv), 'VariableNamingRule','preserve', ...
                   'TextType','string', 'Delimiter',',');
elseif istable(metadataCsv)
    md = metadataCsv;
else
    error('metadataCsv must be a file path string or a table.');
end

% Lowercase all column names (reduce_core expects lowercase)
md.Properties.VariableNames = lower(md.Properties.VariableNames);

% ── Natural sort by run order then trailing index ─────────────────────────
% Problem: when two independent runs (e.g. RunA_1…175 and RunB_1…97) share a
% folder, sorting on the trailing index alone interleaves them as
% RunA_1, RunB_1, RunA_2, RunB_2 … which scrambles drift correction and
% bridge-standard bracketing across independent runs.
%
% Run order is determined by:
%   1. 'runOrder' parameter if supplied — explicit prefix order, e.g.
%      {'RunA', 'RunB'}. Any prefixes not listed are appended last
%      in alphabetical order. Use this when filesystem timestamps are
%      unreliable (e.g. Dropbox re-upload resets mtimes).
%   2. Existing metadata row order otherwise. LADDReduce never invents a
%      chronological order from alphabetical filenames.
%
% Within each prefix group, files are always sorted by trailing index (natural
% numeric: _1, _2 … _9, _10, _11, not lexicographic).
if ismember('file', md.Properties.VariableNames) && ~isempty(runOrder)
    fnames   = cellstr(md.file);
    prefixes = cell(height(md), 1);
    idx_nums = zeros(height(md), 1);

    for ii = 1:height(md)
        [~, bare, ext] = fileparts(fnames{ii});
        fname_bare = [bare ext];
        tok = regexp(fname_bare, '^(.+)_(\d+)\.csv$', 'tokens');
        if ~isempty(tok)
            prefixes{ii} = regexprep(tok{1}{1}, '_+$', '');  % strip any trailing underscores
            idx_nums(ii) = str2double(tok{1}{2});
        else
            prefixes{ii} = fname_bare;
            idx_nums(ii) = Inf;
        end
    end

    unique_pfx = unique(prefixes);

    % Validate supplied prefixes
    unknown_pfx = setdiff(runOrder, unique_pfx);
    if ~isempty(unknown_pfx)
        warning('ladd_enrich_metadata: runOrder contains prefix(es) not found in metadata: %s', ...
            strjoin(unknown_pfx, ', '));
    end
    % Build ordered list: supplied prefixes first, then any remainder alphabetically
    remainder = setdiff(unique_pfx, runOrder);
    ordered_pfx = [runOrder(:); remainder(:)];

    % Assign a group rank to each row based on ordered_pfx
    pfx_rank = containers.Map(ordered_pfx, num2cell(1:numel(ordered_pfx)));
    rank_key = zeros(height(md), 1);
    for ii = 1:height(md)
        rank_key(ii) = pfx_rank(prefixes{ii});
    end

    [~, sort_ord] = sortrows([rank_key, idx_nums]);
    md = md(sort_ord, :);

    fprintf('  Sorted to natural run order (user-specified runOrder, then file index)\n');
    fprintf('  Run order: %s\n', strjoin(ordered_pfx, ' → '));
end

% Ensure concentration and optional 1SD columns exist.
knownColumns = {'known_u_ppm','known_th_ppm','known_sm_ppm', ...
    'known_u_1sd_ppm','known_th_1sd_ppm','known_sm_1sd_ppm'};
for col = knownColumns
    if ~ismember(col{1}, md.Properties.VariableNames)
        md.(col{1}) = NaN(height(md), 1);
    end
end

% Ensure stdname column exists
if ~ismember('stdname', md.Properties.VariableNames)
    md.stdname = strings(height(md), 1);
end

% Convert known_* columns to double if they came in as string
for col = knownColumns
    v = md.(col{1});
    if isstring(v) || iscell(v)
        v = str2double(string(v));
        md.(col{1}) = v;
    end
end

% ── Read lookup table ─────────────────────────────────────────────────────
lk = readtable(lookupPath, 'VariableNamingRule','preserve', ...
               'TextType','string', 'Delimiter',',');
lk.Properties.VariableNames = lower(lk.Properties.VariableNames);
assert(ismember('stdname',lk.Properties.VariableNames), ...
    'ladd_enrich_metadata: lookup table must contain stdname.');

% Ensure lookup has concentration and optional 1SD columns.
for col = knownColumns
    if ~ismember(col{1}, lk.Properties.VariableNames)
        lk.(col{1}) = NaN(height(lk), 1);
    end
end
% Convert lookup known_* to double
for col = knownColumns
    v = lk.(col{1});
    if isstring(v) || iscell(v)
        v = str2double(string(v));
        lk.(col{1}) = v;
    end
end

% Normalised lookup keys for matching
lk_norm = arrayfun(@norm_id, lk.stdname, 'UniformOutput', false);
[uniqueLookupKeys, ~, lookupGroup] = unique(lk_norm);
duplicateKeys = uniqueLookupKeys(accumarray(lookupGroup(:),1) > 1);
assert(isempty(duplicateKeys), ...
    'ladd_enrich_metadata: lookup table contains duplicate normalized stdname values: %s', ...
    strjoin(duplicateKeys, ', '));

% ── Fill in missing known values ──────────────────────────────────────────
nFilled  = 0;
nSkipped = 0;   % already had values
nNoMatch = 0;

for i = 1:height(md)
    % Guard against <missing> string values (MATLAB's NaN equivalent for strings)
    raw_sn = md.stdname(i);
    if ismissing(raw_sn) || strlength(raw_sn) == 0
        continue;
    end
    sn = strtrim(char(raw_sn));
    if isempty(sn), continue; end

    sn_norm = norm_id(sn);

    % Prefer an exact normalized match. A unique substring match supports
    % decorated labels such as run-specific prefixes/suffixes without making
    % lookup-table row order scientifically consequential.
    matchIdx = find(strcmp(sn_norm, lk_norm));
    if isempty(matchIdx)
        candidates = find(cellfun(@(key) contains(sn_norm,key) || ...
            contains(key,sn_norm), lk_norm));
        if numel(candidates) > 1
            candidateNames = string(lk.stdname(candidates));
            error('ladd_enrich_metadata:AmbiguousStandardMatch', ...
                ['Metadata stdname "%s" matches multiple lookup entries: %s. ' ...
                 'Use an exact, unique standard name.'], ...
                sn, strjoin(candidateNames, ', '));
        end
        matchIdx = candidates;
    end

    if isempty(matchIdx)
        nNoMatch = nNoMatch + 1;
        continue;
    end

    % Fill only blank/NaN cells — never overwrite existing values
    filled_any = false;
    for col = knownColumns
        existing = md.(col{1})(i);
        lk_val   = lk.(col{1})(matchIdx);
        if (isnan(existing) || existing == 0) && isfinite(lk_val)
            md.(col{1})(i) = lk_val;
            filled_any = true;
        else
            nSkipped = nSkipped + 1;
        end
    end
    if filled_any, nFilled = nFilled + 1; end
end

if verbose
    fprintf('ladd_enrich_metadata: %d rows enriched from lookup table\n', nFilled);
    % List non-empty stdnames that had no lookup match
    has_stdname = ~ismissing(md.stdname) & strlength(md.stdname) > 0;
    no_u        = has_stdname & ~isfinite(md.known_u_ppm);
    unmatched   = unique(md.stdname(no_u));
    unmatched   = unmatched(~ismissing(unmatched));
    if ~isempty(unmatched)
        fprintf('  Stdnames with no lookup match (add them to the supplied lookup table if needed):\n');
        for k = 1:numel(unmatched)
            fprintf('    "%s"\n', unmatched(k));
        end
    end
    % Print summary of what was filled
    hasCal = isfinite(md.known_u_ppm);
    fprintf('  Rows with known_u_ppm: %d / %d\n', sum(hasCal), height(md));
end

mdOut = md;

% ── Write enriched CSV if requested ───────────────────────────────────────
if ~isempty(saveAs)
    writetable(mdOut, saveAs);
    fprintf('  Wrote enriched metadata: %s\n', saveAs);
end

end % ── END MAIN ──────────────────────────────────────────────────────────


function s = norm_id(id)
% Normalise a standard name for fuzzy matching:
% lowercase, strip spaces/underscores/hyphens
    s = lower(strtrim(char(id)));
    s = regexprep(s, '[\s_\-]+', '');
end
