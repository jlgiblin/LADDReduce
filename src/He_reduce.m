function outtbl = He_reduce(inputFile, pitVolFile, airStdFile, typeMapFile, cfg, varargin)
% HE_REDUCE  Reduce raw 4He mass-spec data to atoms/g per grain.
%
% USAGE
%   cfg = he_default_cfg('zircon');
%   out = He_reduce('helium_input.csv', 'he_pit_volumes.csv', ...
%                   'air_calibration.csv', 'sample_types.csv', cfg, ...
%                   'saveAs','helium_reduced.csv');
%
% INPUTS
%   inputFile   : He raw file. Format auto-detected:
%                   'prep' — your clean prep CSV (has SampleName + GrainID cols)
%                   'old'  — ArArCalc/Helix-MC Excel export
%                   'new'  — Pychron CSV export
%   pitVolFile  : pit-volume CSV  — GrainID | PitVol_um3 | PV1SD_um3
%   airStdFile  : air standard CSV — AirID | FourHeAir | FourHeAir1SD
%                 FourHeAir = atoms/cps from tank depletion master sheet
%   typeMapFile : type map CSV  — SampleName | RunScript
%                 RunScript: 1=air  2=blank  3=unknown  4=He mineral std
%   cfg         : struct from he_default_cfg('zircon' or 'apatite')
%
% OPTIONAL NAME-VALUE PAIRS
%   'saveAs'        — output CSV filename (default '' = no save)
%   'minBlankFrac'  — blank/signal review threshold for HIGH_BLANK (default
%                     1.0). At the default, a row is flagged only when the
%                     session median blank exceeds that row's He signal.
%                     This flag never excludes or changes a value.
%   'skipWindow'    — unknowns on each side of an explicitly excluded std
%                     to flag POSSIBLE_SKIP (default 3). Statistical
%                     STD_REVIEW alone does not trigger this flag.
%   'stdOutlierMAD' — MAD multiplier used to flag
%                     STD_REVIEW (default 5). Review flags NEVER exclude a
%                     standard from calibration by themselves.
%   'stdExcludeIDs' — explicit GrainIDs to mark STD_EXCLUDE (default empty).
%                     Use only for independently documented problems such as
%                     a laser misfire, wrong grain, or confirmed sequence skip.
%
% OUTPUT COLUMNS
%   GrainID, SampleName, RunScript (3=unknown 4=mineralStd),
%   He4_cps, He4_1SD, He4_1SDpct,
%   PitVol_um3, PV1SD_um3, Wt_g, Wt1SD_g,
%   AirStdID, FourHeAir, FourHeAir1SD,
%   He4Unk_atoms, He4Unk1SD_atoms, He4Unk1SDpct,
%   He4Unk_atoms_g, He4Unk1SD_atoms_g, He4Unk1SDpct_g,
%   flags
%
% FLAGS
%   HIGH_BLANK    — median absolute blank/signal fraction > minBlankFrac;
%                   review only, never an automatic exclusion
%   NO_PITVOL     — no pit volume match found
%   NO_AIRSTD     — no FourHeAir available for nearest air shot
%   ZERO_HE       — He4_cps <= 0
%   NOISY_HE      — He4_1SDpct > 5%
%   STD_REVIEW    — mineral standard > stdOutlierMAD*MAD from run median;
%                   diagnostic only, never an automatic exclusion
%   STD_EXCLUDE   — explicitly excluded by stdExcludeIDs
%   POSSIBLE_SKIP — unknown within skipWindow rows of an explicitly
%                   excluded standard; investigate against run notes /
%                   laser log

% ── Parse optional arguments ──────────────────────────────────────────────
p = inputParser;
p.addParameter('saveAs',        '',   @ischar);
p.addParameter('minBlankFrac',  1.00, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('skipWindow',    3,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('stdOutlierMAD', 5,    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('stdExcludeIDs', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.parse(varargin{:});
saveAs        = p.Results.saveAs;
minBlankFrac  = p.Results.minBlankFrac;
skipWindow    = p.Results.skipWindow;
stdOutlierMAD = p.Results.stdOutlierMAD;
stdExcludeIDs = p.Results.stdExcludeIDs;

% ── Validate cfg ──────────────────────────────────────────────────────────
for fld = {'mineral','density'}
    assert(isfield(cfg, fld{1}), ...
        'cfg.%s is required. Call he_default_cfg(''zircon'') to build cfg.', fld{1});
end
rho          = cfg.density;
hasRenameMap = isfield(cfg,'renameMap') && ~isempty(cfg.renameMap);

% ── Detect format and read raw He file ───────────────────────────────────
[~,~,ext] = fileparts(inputFile);
if any(strcmpi(ext, {'.xlsx','.xls'}))
    fmt = 'old';
    raw = read_old_format(inputFile);
else
    % Peek at headers to distinguish prep CSV from raw Pychron CSV
    fid       = fopen(inputFile, 'r', 'n', 'UTF-8');
    firstline = fgetl(fid);
    fclose(fid);
    firstline    = regexprep(firstline, '^\xEF\xBB\xBF', '');  % strip BOM
    hdrs         = lower(strtrim(strsplit(firstline, ',')));
    if any(strcmp(hdrs,'samplename')) && any(strcmp(hdrs,'grainid'))
        fmt = 'prep';
        raw = read_prep_format(inputFile);
    else
        fmt = 'new';
        raw = read_new_format(inputFile);
    end
end
fprintf('\nHe_reduce: %d rows from %s (%s format)\n', height(raw), inputFile, fmt);

% ── Read type map and classify rows ──────────────────────────────────────
typeMap = read_typemap(typeMapFile);
raw     = classify_rows(raw, typeMap);

isAir  = raw.RunScript == 1;
isBlnk = raw.RunScript == 2;
isUnk  = raw.RunScript == 3;
isStd  = raw.RunScript == 4;
isUnkn = raw.RunScript == 0;
apply  = isUnk | isStd;

fprintf('  Air=%d  Blank=%d  Unknown=%d  MineralStd=%d  Unclassified=%d\n', ...
    sum(isAir), sum(isBlnk), sum(isUnk), sum(isStd), sum(isUnkn));

if any(isUnkn)
    unk_names = unique(raw.SampleName(isUnkn));
    fprintf('  WARNING: unclassified sample names (add to typeMapFile):\n');
    for k = 1:numel(unk_names)
        fprintf('    "%s"\n', unk_names{k});
    end
end

% ── Read air standard calibration table ──────────────────────────────────
airstd = read_airstd(airStdFile);
fprintf('  Air std table: %d entries from %s\n', height(airstd), airStdFile);

% ── Attach FourHeAir values to air rows in raw ───────────────────────────
raw_airnorm = cellfun(@he_normID, raw.AirID, 'UniformOutput', false);
as_norm     = cellfun(@he_normID, airstd.AirID, 'UniformOutput', false);

raw.FourHeAir    = NaN(height(raw), 1);
raw.FourHeAir1SD = NaN(height(raw), 1);
nAirMatched = 0;
for i = 1:height(raw)
    if ~isAir(i), continue; end
    idx = find(strcmp(raw_airnorm{i}, as_norm), 1);
    if ~isempty(idx)
        raw.FourHeAir(i)    = airstd.FourHeAir(idx);
        raw.FourHeAir1SD(i) = airstd.FourHeAir1SD(idx);
        nAirMatched         = nAirMatched + 1;
    end
end
nAirMissing = sum(isAir) - nAirMatched;
fprintf('  Air std matching: %d matched, %d missing FourHeAir entry\n', ...
    nAirMatched, nAirMissing);
if nAirMissing > 0
    missingAir = raw.AirID(isAir & isnan(raw.FourHeAir));
    for k = 1:numel(missingAir)
        fprintf('    "%s" not found in airStdFile\n', missingAir{k});
    end
end

% ── Nearest-air-shot bracketing ───────────────────────────────────────────
% Each unknown/standard is assigned FourHeAir from the closest valid air row.
airRows = find(isAir & ~isnan(raw.FourHeAir));

raw.AirStdID      = repmat({''}, height(raw), 1);
raw.AirStdRowUsed = NaN(height(raw), 1);

if isempty(airRows)
    warning('He_reduce: no air rows with valid FourHeAir found. All grains will get NO_AIRSTD flag.');
else
    for i = 1:height(raw)
        if ~apply(i), continue; end
        [~, nearest_k]       = min(abs(airRows - i));
        nearest_row          = airRows(nearest_k);
        raw.AirStdID{i}      = raw.AirID{nearest_row};
        raw.AirStdRowUsed(i) = nearest_row;
        raw.FourHeAir(i)     = raw.FourHeAir(nearest_row);
        raw.FourHeAir1SD(i)  = raw.FourHeAir1SD(nearest_row);
    end
end

% ── Apply rename map ─────────────────────────────────────────────────────
if hasRenameMap
    rm = cfg.renameMap;
    for k = 1:size(rm,1)
        hit = strcmp(raw.GrainID, rm{k,1});
        if any(hit)
            raw.GrainID(hit) = {rm{k,2}};
            fprintf('  Renamed "%s" -> "%s" (%d rows)\n', rm{k,1}, rm{k,2}, sum(hit));
        end
    end
end

% ── Match pit volumes ─────────────────────────────────────────────────────
pv       = read_pit_volumes(pitVolFile);
raw_norm = cellfun(@he_normID, raw.GrainID, 'UniformOutput', false);
pv_norm  = cellfun(@he_normID, pv.GrainID,  'UniformOutput', false);

raw.PitVol_um3 = NaN(height(raw), 1);
raw.PV1SD_um3  = NaN(height(raw), 1);
for i = 1:height(raw)
    if ~apply(i), continue; end
    idx = find(strcmp(raw_norm{i}, pv_norm), 1);
    if ~isempty(idx)
        raw.PitVol_um3(i) = pv.PitVol_um3(idx);
        raw.PV1SD_um3(i)  = pv.PV1SD_um3(idx);
    end
end

nPVmatched = sum(~isnan(raw.PitVol_um3) & apply);
nPVmissing = sum( isnan(raw.PitVol_um3) & apply);
fprintf('  Pit volumes: %d matched, %d missing\n', nPVmatched, nPVmissing);
if nPVmissing > 0
    miss = raw.GrainID(isnan(raw.PitVol_um3) & apply);
    for k = 1:min(numel(miss),15)
        fprintf('    "%s"  ->  norm="%s"\n', miss{k}, he_normID(miss{k}));
    end
    if numel(miss)>15
        fprintf('    ...and %d more. Use cfg.renameMap to fix.\n', numel(miss)-15);
    end
end

% ── Weight from pit volume ────────────────────────────────────────────────
% Wt (g) = PitVol (um3) * 1e-12 (cm3/um3) * density (g/cm3)
raw.Wt_g    = raw.PitVol_um3 .* 1e-12 .* rho;
raw.Wt1SD_g = (raw.PV1SD_um3 ./ raw.PitVol_um3) .* raw.Wt_g;

% ── He4 percent 1SD ──────────────────────────────────────────────────────
raw.He4_1SDpct = (raw.He4_1SD ./ abs(raw.He4_cps)) * 100;

% ── cps -> atoms via FourHeAir ────────────────────────────────────────────
% He4_atoms = FourHeAir * He4_cps
% Quadrature error: cps measurement uncertainty + FourHeAir calibration uncertainty
raw.He4Unk_atoms    = NaN(height(raw), 1);
raw.He4Unk1SD_atoms = NaN(height(raw), 1);

hasAir = apply & ~isnan(raw.FourHeAir);
raw.He4Unk_atoms(hasAir)    = raw.FourHeAir(hasAir) .* raw.He4_cps(hasAir);
raw.He4Unk1SD_atoms(hasAir) = sqrt( ...
    (raw.He4_1SD(hasAir)    .* raw.FourHeAir(hasAir)   ).^2 + ...
    (raw.He4_cps(hasAir)    .* raw.FourHeAir1SD(hasAir) ).^2 );

raw.He4Unk1SDpct = (raw.He4Unk1SD_atoms ./ raw.He4Unk_atoms) * 100;

% ── atoms -> atoms/g ─────────────────────────────────────────────────────
raw.He4Unk_atoms_g    = raw.He4Unk_atoms ./ raw.Wt_g;
raw.He4Unk1SD_atoms_g = sqrt( ...
    (raw.He4Unk1SD_atoms ./ raw.Wt_g                   ).^2 + ...
    (raw.He4Unk_atoms    .* raw.Wt1SD_g ./ raw.Wt_g.^2 ).^2 );
raw.He4Unk1SDpct_g = (raw.He4Unk1SD_atoms_g ./ raw.He4Unk_atoms_g) * 100;

% ── QC flags ─────────────────────────────────────────────────────────────
raw.flags = strings(height(raw), 1);

% BLANK/SIGNAL DIAGNOSTIC + HIGH_BLANK REVIEW FLAG
% Save the numeric ratio for every analyzed grain so the threshold is
% transparent and can be re-screened later without rerunning reduction.
raw.BlankMedian_cps = NaN(height(raw), 1);
raw.BlankSignalFraction = NaN(height(raw), 1);
raw.BlankReviewFractionThreshold = repmat(minBlankFrac, height(raw), 1);
if any(isBlnk)
    blankMed  = median(raw.He4_cps(isBlnk), 'omitnan');
    raw.BlankMedian_cps(:) = blankMed;
    raw.BlankSignalFraction(apply) = abs(blankMed) ./ abs(raw.He4_cps(apply));
    highBlank = apply & isfinite(raw.BlankSignalFraction) & ...
                raw.BlankSignalFraction > minBlankFrac;
    raw.flags(highBlank) = strtrim(raw.flags(highBlank) + " HIGH_BLANK");
end

% NO_PITVOL
noPV = apply & isnan(raw.PitVol_um3);
raw.flags(noPV) = strtrim(raw.flags(noPV) + " NO_PITVOL");

% NO_AIRSTD
noAir = apply & isnan(raw.FourHeAir);
raw.flags(noAir) = strtrim(raw.flags(noAir) + " NO_AIRSTD");

% ZERO_HE
zeroHe = apply & raw.He4_cps <= 0;
raw.flags(zeroHe) = strtrim(raw.flags(zeroHe) + " ZERO_HE");

% NOISY_HE: He4 1SD > 5%
noisyHe = apply & raw.He4_1SDpct > 5;
raw.flags(noisyHe) = strtrim(raw.flags(noisyHe) + " NOISY_HE");

% STANDARD REVIEW + EXPLICIT EXCLUSION
% The robust-deviation test is a SCREENING diagnostic only. A low/high
% standard can be scientifically real, especially while diagnosing a
% systematic He bias, so statistical unusualness must not silently become an
% exclusion. Only stdExcludeIDs creates STD_EXCLUDE.
raw.StdDeviationMAD       = NaN(height(raw), 1);
raw.StdReviewMADThreshold = repmat(stdOutlierMAD, height(raw), 1);
raw.StdExplicitExclude    = false(height(raw), 1);
reviewIdx                 = [];
if any(isStd) && sum(isfinite(raw.He4Unk_atoms_g(isStd))) >= 3
    std_vals   = raw.He4Unk_atoms_g(isStd);
    std_median = median(std_vals, 'omitnan');
    std_MAD    = median(abs(std_vals - std_median), 'omitnan');
    std_thresh = stdOutlierMAD * std_MAD;

    stdRows    = find(isStd);
    if std_MAD > 0
        raw.StdDeviationMAD(stdRows) = ...
            abs(raw.He4Unk_atoms_g(stdRows) - std_median) ./ std_MAD;
        reviewIdx = stdRows(abs(raw.He4Unk_atoms_g(stdRows) - std_median) > std_thresh);
    end

    if ~isempty(reviewIdx)
        for k = 1:numel(reviewIdx)
            oRow = reviewIdx(k);
            raw.flags(oRow) = strtrim(raw.flags(oRow) + " STD_REVIEW");
            dev = raw.StdDeviationMAD(oRow);
            fprintf('  STD_REVIEW (included unless explicitly excluded): %s  atoms_g=%.3e  median=%.3e  dev=%.1f x MAD\n', ...
                raw.GrainID{oRow}, raw.He4Unk_atoms_g(oRow), std_median, dev);
        end
    else
        fprintf('  Mineral std QC: all %d stds within %.0f x MAD of median (%.3e).\n', ...
            numel(stdRows), stdOutlierMAD, std_median);
    end
end

excludeList = to_cellstr(stdExcludeIDs);
excludeNorm = cellfun(@he_normID, excludeList, 'UniformOutput', false);
rawNormNow = cellfun(@he_normID, raw.GrainID, 'UniformOutput', false);
explicitExclude = isStd & ismember(rawNormNow, excludeNorm);
raw.StdExplicitExclude(explicitExclude) = true;
raw.flags(explicitExclude) = strtrim(raw.flags(explicitExclude) + " STD_EXCLUDE");
if ~isempty(excludeNorm)
    unmatched = excludeNorm(~ismember(excludeNorm, rawNormNow(isStd)));
    if ~isempty(unmatched)
        warning('He_reduce:stdExcludeIDsNotFound', ...
            'stdExcludeIDs not found among mineral standards: %s', strjoin(unmatched, ', '));
    end
end
if any(explicitExclude)
    fprintf('  Explicitly excluded mineral standards (%d): %s\n', ...
        sum(explicitExclude), strjoin(raw.GrainID(explicitExclude), ', '));
end

% Flag nearby unknowns only when a standard has been explicitly excluded.
% A statistical STD_REVIEW may simply reflect natural standard scatter and
% is not evidence that adjacent unknowns were skipped or misidentified.
% This remains review-only and never alters values.
checkRows = find(explicitExclude);
rowIdx = (1:height(raw)).';
for k = 1:numel(checkRows)
    nearRows = find(isUnk & abs(rowIdx - checkRows(k)) <= skipWindow);
    if ~isempty(nearRows)
        raw.flags(nearRows) = strtrim(raw.flags(nearRows) + " POSSIBLE_SKIP");
        fprintf('    POSSIBLE_SKIP near %s: %s\n', ...
            raw.GrainID{checkRows(k)}, strjoin(raw.GrainID(nearRows), ', '));
    end
end

% ── Build output table ────────────────────────────────────────────────────
wantCols = {'GrainID','SampleName','RunScript', ...
            'He4_cps','He4_1SD','He4_1SDpct', ...
            'PitVol_um3','PV1SD_um3','Wt_g','Wt1SD_g', ...
            'AirStdID','FourHeAir','FourHeAir1SD', ...
            'He4Unk_atoms','He4Unk1SD_atoms','He4Unk1SDpct', ...
            'He4Unk_atoms_g','He4Unk1SD_atoms_g','He4Unk1SDpct_g', ...
            'BlankMedian_cps','BlankSignalFraction','BlankReviewFractionThreshold', ...
            'StdDeviationMAD','StdReviewMADThreshold','StdExplicitExclude', ...
            'flags'};
haveCols = intersect(wantCols, raw.Properties.VariableNames, 'stable');
outtbl   = raw(apply, haveCols);

nFlagged = sum(strlength(outtbl.flags) > 0);
fprintf('He_reduce: output %d rows (%d unknowns, %d mineral stds, %d flagged)\n\n', ...
    height(outtbl), sum(isUnk), sum(isStd), nFlagged);

if ~isempty(saveAs)
    writetable(outtbl, saveAs);
    fprintf('Wrote: %s\n', saveAs);
end

end % ── END MAIN ──────────────────────────────────────────────────────────


%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function raw = read_prep_format(f)
% Standardised prep CSV: RunID | SampleName | GrainID | He4_cps | He4_1SD
% SFT format: air shot number is in GrainID column when SampleName='Air'
    T          = readtable(f, 'VariableNamingRule','preserve', 'Delimiter',',');
    grainID    = to_cellstr(col_str(T, {'GrainID','grainid','Grain_ID'}));
    sampleName = to_cellstr(col_str(T, {'SampleName','samplename','Sample','sample'}));
    runID      = to_cellstr(col_str(T, {'RunID','runid','Run_ID','run_id'}));
    % For SFT format the air-shot number lives in GrainID. Other prep files
    % Some exports leave GrainID blank and store the identifier in RunID.
    % Prefer a populated air-row GrainID, otherwise retain RunID.
    isAirRow = strcmpi(sampleName, 'air') | strcmpi(runID, 'a-01-s');
    airID    = runID;
    airGrainID = isAirRow & ~cellfun(@isempty, grainID);
    airID(airGrainID) = grainID(airGrainID);
    empty      = cellfun(@isempty, grainID);
    grainID(empty) = runID(empty);
    he4    = col_or_nan(T, {'He4_cps','He4cps','he4_cps','He4_','He4'});
    he4_er = col_or_nan(T, {'He4_1SD','He4_1sd','He4_Er','He4_er','He4_error'});
    rs     = zeros(height(T), 1);
    raw    = table(grainID, sampleName, airID, rs, he4, he4_er, ...
        'VariableNames', {'GrainID','SampleName','AirID','RunScript','He4_cps','He4_1SD'});
end


function raw = read_old_format(f)
% ArArCalc / Helix-MC Excel export.
% GrainID from Comment column; SampleName from Sample column; AirID from Run_ID.
    T       = readtable(f, 'VariableNamingRule','preserve', 'ReadVariableNames',true);
    runID   = to_cellstr(col_str(T, {'Run_ID','RunID'}));
    grainID = to_cellstr(col_str(T, {'Comment','comment','AnalysisName','analysisname'}));
    empty   = cellfun(@isempty, grainID);
    grainID(empty) = runID(empty);
    sampleName = to_cellstr(col_str(T, {'Sample','SampleName','sample'}));
    airID      = runID;
    he4        = col_or_nan(T, {'He4_','He4','He4_bl_corrected'});
    he4_er     = col_or_nan(T, {'He4_Er','He4_Er_','He4_1SD','He4_error'});
    rs         = zeros(height(T), 1);
    raw        = table(grainID, sampleName, airID, rs, he4, he4_er, ...
        'VariableNames', {'GrainID','SampleName','AirID','RunScript','He4_cps','He4_1SD'});
end


function raw = read_new_format(f)
% Pychron CSV export.
% GrainID from comment column; SampleName from sample column.
    T          = readtable(f, 'VariableNamingRule','preserve', 'Delimiter',',');
    cn         = T.Properties.VariableNames;
    comment    = to_cellstr(col_str(T, {'comment','Comment'}));
    identifier = to_cellstr(col_str(T, {'identifier','Identifier','aliquot'}));
    grainID    = comment;
    grainID(cellfun(@isempty, grainID)) = identifier(cellfun(@isempty, grainID));
    sampleName = to_cellstr(col_str(T, {'sample','Sample'}));
    airID      = comment;
    atype      = lower(to_cellstr(col_str(T, {'analysis_type','analysistype','type'})));
    rs         = zeros(height(T), 1);
    rs(contains(atype,'air') & ~contains(atype,'blank')) = 1;
    rs(contains(atype,'blank'))   = 2;
    rs(contains(atype,'unknown')) = 3;
    he4    = col_or_nan(T, {'He4_bl_corrected','He4_ic_decay_corrected','He4_'});
    erCols = find(strcmpi(cn,'error'));
    if ~isempty(erCols)
        he4_er = double(T{:, erCols(1)});
    else
        he4_er = col_or_nan(T, {'He4_error','He4_1SD','error'});
    end
    raw = table(grainID, sampleName, airID, rs, he4, he4_er, ...
        'VariableNames', {'GrainID','SampleName','AirID','RunScript','He4_cps','He4_1SD'});
end


function raw = classify_rows(raw, typeMap)
% Assign RunScript from typeMap by SampleName. Unmatched rows get 0.
    raw_norm = cellfun(@he_normID, raw.SampleName, 'UniformOutput', false);
    tm_norm  = cellfun(@he_normID, typeMap.SampleName, 'UniformOutput', false);
    raw.RunScript = zeros(height(raw), 1);
    for i = 1:height(raw)
        idx = find(strcmp(raw_norm{i}, tm_norm), 1);
        if ~isempty(idx)
            raw.RunScript(i) = typeMap.RunScript(idx);
        end
    end
end


function typeMap = read_typemap(f)
    T       = readtable(f, 'VariableNamingRule','preserve');
    cn_norm = lower(regexprep(T.Properties.VariableNames,'[^a-zA-Z0-9]+',''));
    sIdx    = find(strcmp(cn_norm,'samplename'), 1);
    if isempty(sIdx), sIdx = 1; end
    rIdx    = find(strcmp(cn_norm,'runscript'), 1);
    if isempty(rIdx), rIdx = 2; end
    sampleName = to_cellstr(T{:, sIdx});
    rs         = double(T{:, rIdx});
    typeMap    = table(sampleName, rs, 'VariableNames', {'SampleName','RunScript'});
    typeMap    = typeMap(~cellfun(@isempty, typeMap.SampleName), :);
    fprintf('  Type map: %d entries from %s\n', height(typeMap), f);
end


function airstd = read_airstd(f)
    T      = readtable(f, 'VariableNamingRule','preserve');
    cn     = T.Properties.VariableNames;
    airID  = to_cellstr(T{:,1});
    fha    = col_or_nan(T, {'FourHeAir','fourheair','FourHe_Air','He4Air','AirCal'});
    if all(isnan(fha)) && width(T) >= 2
        fha = double(T{:,2});
        fprintf('  AirStd: using column 2 (%s) as FourHeAir\n', cn{2});
    end
    fha1sd = col_or_nan(T, {'FourHeAir1SD','FourHe_Air1SD','He4Air1SD','AirCal1SD'});
    if all(isnan(fha1sd)) && width(T) >= 3
        fha1sd = double(T{:,3});
        fprintf('  AirStd: using column 3 (%s) as FourHeAir1SD\n', cn{3});
    end
    airstd = table(airID, fha, fha1sd, ...
        'VariableNames', {'AirID','FourHeAir','FourHeAir1SD'});
    airstd = airstd(~cellfun(@isempty, airstd.AirID), :);
end


function pv = read_pit_volumes(f)
    T       = readtable(f, 'VariableNamingRule','preserve');
    cn      = T.Properties.VariableNames;
    grainID = to_cellstr(col_str(T, {'GrainID','AnalysisName','Name','ID','Sample',cn{1}}));
    vol     = col_or_nan(T, {'PitVol_um3','PitVolume_um3','MIV_um3','MeanIntVol','Volume_um3','PitVol'});
    if all(isnan(vol)) && width(T) >= 2
        vol = double(T{:,2});
        fprintf('  PitVol: using column 2 (%s)\n', cn{2});
    end
    pv1sd   = col_or_nan(T, {'PV1SD_um3','PV1SD','PropSD_um3','MIVpropSD'});
    if all(isnan(pv1sd)) && width(T) >= 3
        pv1sd = double(T{:,end});
        fprintf('  PV1SD: using last column (%s)\n', cn{end});
    end
    pv = table(grainID, vol, pv1sd, ...
        'VariableNames', {'GrainID','PitVol_um3','PV1SD_um3'});
    pv = pv(~cellfun(@isempty, pv.GrainID), :);
    fprintf('  Pit volumes: %d grains from %s\n', height(pv), f);
end


function id_norm = he_normID(id)
    s = lower(strtrim(char(id)));
    s = regexprep(s, '\s+',  '-');
    s = regexprep(s, '[_]+', '-');
    s = regexprep(s, '-+',   '-');
    s = regexprep(s, '-ft$', '');
    id_norm = strtrim(s);
end


function v = col_or_nan(T, aliases)
    v       = NaN(height(T), 1);
    cn_norm = lower(regexprep(T.Properties.VariableNames,'[^a-zA-Z0-9]+',''));
    for a = 1:numel(aliases)
        key = lower(regexprep(char(aliases{a}),'[^a-zA-Z0-9]+',''));
        idx = find(strcmp(cn_norm, key), 1);
        if ~isempty(idx)
            c = T{:,idx};
            if ~isnumeric(c), c = str2double(c); end
            v = c; return;
        end
    end
end


function v = col_str(T, aliases)
    v       = repmat({''},height(T),1);
    cn_norm = lower(regexprep(T.Properties.VariableNames,'[^a-zA-Z0-9]+',''));
    for a = 1:numel(aliases)
        key = lower(regexprep(char(aliases{a}),'[^a-zA-Z0-9]+',''));
        idx = find(strcmp(cn_norm, key), 1);
        if ~isempty(idx)
            v = T{:,idx}; return;
        end
    end
end


function c = to_cellstr(v)
    if iscell(v)
        c = cellfun(@(x) strtrim(char(string(x))), v, 'UniformOutput', false);
    elseif isstring(v),   c = strtrim(cellstr(v));
    elseif isnumeric(v),  c = arrayfun(@num2str, v, 'UniformOutput', false);
    elseif isdatetime(v), c = cellstr(datestr(v));
    else,                 c = strtrim(cellstr(string(v)));
    end
    c(strcmp(c,'NaN'))       = {''};
    c(strcmp(c,'<missing>')) = {''};
end
