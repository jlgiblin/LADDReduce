function T = ladd_add_grainid(uthTable, metadataCsv)
% LADD_ADD_GRAINID  Re-join GrainID from metadata CSV onto reduce_core output.
%
% reduce_core does not carry the GrainID column through to its output, but
% it does keep the 'file' column which matches the 'file' column in your
% metadata CSV. This function uses that to join GrainID back on.
%
% USAGE
%   % After running reduce_core / ladd_reduce_zircon:
%   out = ladd_reduce_zircon(folder, metaEnr, 'saveAs','UTh_out.csv');
%   out = ladd_add_grainid(out, metaEnr);
%
%   % Or if loading from saved CSV:
%   out = readtable('UTh_out.csv', 'VariableNamingRule','preserve');
%   out = ladd_add_grainid(out, 'metadata_enriched.csv');
%
% INPUTS
%   uthTable    : reduce_core output table (or path to saved CSV)
%   metadataCsv : path to metadata CSV (enriched or original — just needs
%                 'file' and 'grainid' columns)
%
% OUTPUT
%   T : same table as input with 'GrainID' column inserted after 'file'

% ── Read inputs ───────────────────────────────────────────────────────────
if ischar(uthTable) || isstring(uthTable)
    T = readtable(char(uthTable), 'VariableNamingRule','preserve', ...
        'TextType','string');
else
    T = uthTable;
end
T.Properties.VariableNames = lower(T.Properties.VariableNames);

% Accept either a file path or an already-loaded table
if istable(metadataCsv)
    md = metadataCsv;
    md.Properties.VariableNames = lower(md.Properties.VariableNames);
else
    md = readtable(char(metadataCsv), 'VariableNamingRule','preserve', ...
        'TextType','string', 'Delimiter',',', 'MissingRule','fill');
    md.Properties.VariableNames = lower(md.Properties.VariableNames);
end

% ── Validate ──────────────────────────────────────────────────────────────
assert(ismember('file', T.Properties.VariableNames), ...
    'ladd_add_grainid: UTh table must have a ''file'' column (reduce_core output).');
assert(ismember('file', md.Properties.VariableNames), ...
    'ladd_add_grainid: metadata CSV must have a ''file'' column.');
assert(ismember('grainid', md.Properties.VariableNames), ...
    'ladd_add_grainid: metadata CSV must have a ''grainid'' column.');

% ── Build file->grainid lookup from metadata ──────────────────────────────
% Strip any path from filenames so bare filename is the key
strip_path = @(f) char(regexp(string(f), '[^/\\]+$', 'match', 'once'));

md_files   = cellfun(strip_path, cellstr(md.file), 'UniformOutput', false);
md_grainid = cellstr(md.grainid);

% A duplicate filename makes every downstream metadata and GrainID join
% ambiguous. Stop instead of silently taking the first matching row.
md_keys = lower(string(md_files));
[uniqueKeys,~,keyGroup] = unique(md_keys);
duplicateKeys = uniqueKeys(accumarray(keyGroup,1) > 1);
assert(isempty(duplicateKeys), ...
    'ladd_add_grainid: duplicate metadata filename(s): %s', ...
    strjoin(duplicateKeys, ', '));

% ── Join GrainID onto output table ────────────────────────────────────────
T_files = cellfun(strip_path, cellstr(T.file), 'UniformOutput', false);

GrainID = repmat({''}, height(T), 1);
nMatched = 0;
for i = 1:height(T)
    idx = find(strcmp(T_files{i}, md_files), 1);
    if ~isempty(idx)
        GrainID{i} = char(md_grainid{idx});
        nMatched = nMatched + 1;
    end
end

nMissing = height(T) - nMatched;
fprintf('ladd_add_grainid: %d matched, %d unmatched\n', nMatched, nMissing);
if nMissing > 0
    miss = T.file(cellfun(@isempty, GrainID));
    fprintf('  Files with no GrainID match:\n');
    for k = 1:min(5, numel(miss))
        fprintf('    %s\n', miss{k});
    end
end

% Insert GrainID as second column (after 'file')
GrainID_col = string(GrainID);
T = addvars(T, GrainID_col, 'After','file', 'NewVariableNames','grainid');

fprintf('  GrainID column added to table (%d rows)\n', height(T));
end
