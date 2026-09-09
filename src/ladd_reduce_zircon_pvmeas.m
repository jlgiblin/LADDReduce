function outtbl = ladd_reduce_zircon_pvmeas(folder, metadataCsv, hePitVolCsv, ...
    uthPitVolSource, varargin)
% LADD_REDUCE_ZIRCON_PVMEAS  U-Th reduction for zircon with measured He
% and U-Th pit volumes. uthPitVolSource can be either a measured pit-volume
% CSV or one explicitly declared session-average U-Th pit volume (um3).
%
% This is the zircon counterpart of ladd_reduce_apatite_pvmeas_hybrid.
% The calibration and review policies are intentionally parallel:
%   1. NIST612 drift correction is run first by reduce_core_pv.
%   2. The NIST612-anchored U and Th values are retained as QC columns.
%   3. A caller-selected zircon bridge standard calibrates unknowns using
%      the same PV-normalized formula for U and Th:
%
%        X_ppm_unknown =
%          (X_cps_unknown / PV_unknown) /
%          (X_cps_bridge  / PV_bridge) * known_X_ppm_bridge
%
%   4. Bridge-standard Hampel results are review flags only. No calibration
%      row is removed unless it is named in excludeBridgeIDs.
%   5. The bridge result is checked against the independent NIST612 result.
%
% Zircon-specific differences from the apatite hybrid:
%   - density = 4.65 g cm^-3
%   - Si internal standard (29Si)
%   - no Sm measurement or Sm production term
%
% IMPORTANT STANDARD POLICY
%   bridgeStandardName is required. The function never silently assumes a
%   bridge identity, and it never substitutes a NIST-derived bridge concentration
%   when known_u_ppm or known_th_ppm is missing.
%
%   Metadata rows used by reduce_core_pv as NIST612 anchors must identify
%   actual NIST612 glass. Inconsistent type/stdname combinations stop the
%   run for correction rather than being silently repaired. At least one
%   genuine NIST612 row is required unless the caller explicitly selects
%   the documented bridge-only path.
%
% USAGE
%   out = ladd_reduce_zircon_pvmeas(folder, metadataCsv, ...
%             hePitVolCsv, uthPitVolCsv, ...
%             'bridgeStandardName','ReferenceMaterial', ...
%             'anchorMode','median', ...
%             'saveAs','UTh_out.csv');
%   out = ladd_reduce_zircon_pvmeas(folder, metadataCsv, ...
%             hePitVolCsv, 74200, 'uthAverage1sd',8885, ...
%             'bridgeStandardName','ReferenceMaterial', 'saveAs','UTh_out.csv');
%
% OPTIONAL NAME-VALUE PAIRS
%   'saveAs'             : output CSV filename (default '' = no save)
%   'bridgeStandardName' : REQUIRED exact stdname used as the zircon bridge
%   'anchorMode'         : 'following' | 'median' | 'single'
%                          (default 'following')
%   'singleBridgeID'     : exact bridge GrainID required for 'single' mode
%   'windowOverridesCsv' : optional file,t0,t1 table for controlled/manual windows
%   'parentScalarInterp' : 'linear' | 'pchip' (default 'linear')
%   'kMAD'               : NIST review threshold (default 6)
%   'nistHalfWin'        : NIST review half-window (default 12)
%   'autoExcludeNistReviews' : must remain false; review rows are retained
%   'excludeNistFiles'   : documented bad NIST filenames to exclude
%   'kMAD_std'           : bridge review threshold (default 5)
%   'halfWinStd'         : bridge review half-window (default 12)
%   'autoExcludeBridgeReviews' : must remain false; review rows are retained
%   'excludeBridgeIDs'   : documented bad bridge GrainIDs to exclude
%   'nistCheckTolPct'    : bridge-vs-NIST warning threshold (default 25)
%   'allowBridgeOnly'    : false by default. Set true only for a documented
%                          session with no genuine NIST612 glass. The named
%                          bridge is then used as the sole parent standard;
%                          no independent NIST comparison is reported.
%   'uthAverage1sd'      : required positive 1SD (um3) when uthPitVolSource
%                          is a numeric average; ignored for measured CSVs
%
% anchorMode='following' reproduces the historical next-standard bracket.
% anchorMode='median' uses one session-wide median
% known_ppm/(cps/PV) factor and is retained as a sensitivity option.

% ── Parse arguments ──────────────────────────────────────────────────────
p = inputParser;
p.addParameter('saveAs',             '',       @ischar);
p.addParameter('bridgeStandardName', '',       @(x)ischar(x)||isstring(x));
p.addParameter('anchorMode',         'following',@(x)ischar(x)||isstring(x));
p.addParameter('singleBridgeID',     '',       @(x)ischar(x)||isstring(x));
p.addParameter('windowOverridesCsv', '',       @(x)ischar(x)||isstring(x));
p.addParameter('parentScalarInterp', 'linear', @(x)ischar(x)||isstring(x));
p.addParameter('kMAD',               6,        @(x)isnumeric(x)&&isscalar(x));
p.addParameter('nistHalfWin',        12,       @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('autoExcludeNistReviews', false, ...
    @(x)islogical(x)||isnumeric(x));
p.addParameter('excludeNistFiles', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.addParameter('kMAD_std',           5,        @(x)isnumeric(x)&&isscalar(x));
p.addParameter('halfWinStd',         12,       @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('autoExcludeBridgeReviews', false, ...
    @(x)islogical(x)||isnumeric(x));
p.addParameter('excludeBridgeIDs', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.addParameter('nistCheckTolPct',    25,       @(x)isnumeric(x)&&isscalar(x));
p.addParameter('allowBridgeOnly',    false,    @(x)islogical(x)||isnumeric(x));
p.addParameter('uthAverage1sd',      NaN,      @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});

saveAs     = p.Results.saveAs;
bridgeName = char(string(p.Results.bridgeStandardName));
anchorMode = lower(char(string(p.Results.anchorMode)));
singleBridgeID = char(string(p.Results.singleBridgeID));
windowOverridesCsv = char(string(p.Results.windowOverridesCsv));
interpMth  = char(string(p.Results.parentScalarInterp));
kMAD       = p.Results.kMAD;
nistHalfWin= p.Results.nistHalfWin;
kMAD_std   = p.Results.kMAD_std;
halfWinStd = p.Results.halfWinStd;
nistTolPct = p.Results.nistCheckTolPct;
allowBridgeOnly = logical(p.Results.allowBridgeOnly);
uthAverage1sd = p.Results.uthAverage1sd;

useAverageUthPit = isnumeric(uthPitVolSource);
if useAverageUthPit
    assert(isscalar(uthPitVolSource) && isfinite(uthPitVolSource) && ...
        uthPitVolSource > 0, ...
        'Numeric uthPitVolSource must be one positive average volume in um3.');
    assert(isfinite(uthAverage1sd) && uthAverage1sd > 0, ...
        ['uthAverage1sd must be explicitly supplied and positive when using ' ...
         'an average U-Th pit volume.']);
end

autoExcludeNistReviews = logical(p.Results.autoExcludeNistReviews);
autoExcludeBridgeReviews = logical(p.Results.autoExcludeBridgeReviews);
excludeNistFiles = p.Results.excludeNistFiles;
excludeBridgeIDs = p.Results.excludeBridgeIDs;

assert(~isempty(strtrim(bridgeName)), ...
    ['ladd_reduce_zircon_pvmeas: bridgeStandardName is required. ' ...
     'Pass the exact stdname (for example, ''ReferenceMaterial'').']);
assert(any(strcmp(anchorMode, {'following','median','single'})), ...
    ['ladd_reduce_zircon_pvmeas: anchorMode must be ''following'', ', ...
     '''median'', or ''single''.']);
if strcmp(anchorMode,'single')
    assert(~isempty(strtrim(singleBridgeID)), ...
        'ladd_reduce_zircon_pvmeas: singleBridgeID is required for single mode.');
end
assert(~autoExcludeNistReviews, ...
    ['ladd_reduce_zircon_pvmeas: autoExcludeNistReviews=true is disabled. ' ...
     'Use excludeNistFiles for independently documented bad analyses.']);
assert(~autoExcludeBridgeReviews, ...
    ['ladd_reduce_zircon_pvmeas: autoExcludeBridgeReviews=true is disabled. ' ...
     'Use excludeBridgeIDs for independently documented bad analyses.']);

fprintf('  Zircon bridge standard: "%s" (exact stdname; blank-stdname GrainID fallback)\n', bridgeName);
fprintf('  Anchor mode: "%s"\n', anchorMode);
fprintf('  NIST and bridge Hampel policy: review only; no automatic exclusion\n');

NA          = 6.02214076e23;
MW_U        = 238.02891;
MW_Th       = 232.03806;
rho_zircon  = 4.65;
cm3_per_um3 = 1e-12;

% ── Step 0: Preflight metadata and protect the NIST-only pass ────────────
[coreMetadataCsv, tempMetadataCsv, nTrueNist, nMetadataReclassified, ...
    nBridgeCoreRows] = prepare_core_metadata( ...
        metadataCsv, bridgeName, allowBridgeOnly);
cleanupCoreMetadata = onCleanup(@() delete_if_exists(tempMetadataCsv));

hasIndependentNist = nTrueNist > 0;
if hasIndependentNist
    calibrationPath = "NIST612_THEN_BRIDGE";
else
    calibrationPath = "BRIDGE_ONLY_NO_NIST";
end

fprintf('  Genuine NIST612 rows available to reduce_core_pv: %d\n', nTrueNist);
if ~hasIndependentNist
    fprintf(['  Bridge-only session: %d "%s" rows are used as the sole ', ...
        'parent calibration. No independent NIST612 check is available.\n'], ...
        nBridgeCoreRows, bridgeName);
end

% ── Step 1: Extract and calibrate parent signals ─────────────────────────
if hasIndependentNist
    fprintf('ladd_reduce_zircon_pvmeas: running NIST612 calibration extraction...\n');
else
    fprintf('ladd_reduce_zircon_pvmeas: running bridge-only signal extraction...\n');
end

cfg = struct( ...
    'mineral',          'zircon',   ...
    'density',          rho_zircon, ...
    'mode',             'nois',     ...
    'IS_element',       'Si',       ...
    'IS_ppm_std',       3.10e5,     ...
    'IS_ppm_unk',       1.52e5,     ...
    'pick',             struct( ...
        'U',  "238U",  ...
        'Th', "232Th", ...
        'IS', "29Si"), ...
    'includeSm',        false,      ...
    'matrixScalarName', "NONE"      ...
);

outtbl = reduce_core_pv(folder, coreMetadataCsv, cfg, ...
    'parentScalarInterp', interpMth, ...
    'kMAD', kMAD, ...
    'nistHalfWin', nistHalfWin, ...
    'showNistSummary', hasIndependentNist, ...
    'autoExcludeNistReviews', autoExcludeNistReviews, ...
    'excludeNistFiles', excludeNistFiles, ...
    'windowOverridesCsv',windowOverridesCsv);
clear cleanupCoreMetadata
fprintf('  reduce_core_pv complete: %d rows\n', height(outtbl));

cn0 = lower(outtbl.Properties.VariableNames);
col_u0  = find(strcmp(cn0, 'u_ppm'),  1);
col_th0 = find(strcmp(cn0, 'th_ppm'), 1);
col_u_se0  = find(strcmp(cn0, 'u_ppm_se'),  1);
col_th_se0 = find(strcmp(cn0, 'th_ppm_se'), 1);
assert(~isempty(col_u0) && ~isempty(col_th0), ...
    'ladd_reduce_zircon_pvmeas: reduce_core_pv did not return U_ppm and Th_ppm.');
assert(~isempty(col_u_se0) && ~isempty(col_th_se0), ...
    'ladd_reduce_zircon_pvmeas: reduce_core_pv did not return parent concentration uncertainties.');
outtbl.u_ppm_core_reference  = double(outtbl{:, col_u0});
outtbl.th_ppm_core_reference = double(outtbl{:, col_th0});
outtbl.u_ppm_core_reference_se  = double(outtbl{:, col_u_se0});
outtbl.th_ppm_core_reference_se = double(outtbl{:, col_th_se0});
if hasIndependentNist
    outtbl.u_ppm_nist612  = outtbl.u_ppm_core_reference;
    outtbl.th_ppm_nist612 = outtbl.th_ppm_core_reference;
    outtbl.u_ppm_nist612_se  = outtbl.u_ppm_core_reference_se;
    outtbl.th_ppm_nist612_se = outtbl.th_ppm_core_reference_se;
else
    outtbl.u_ppm_nist612  = NaN(height(outtbl),1);
    outtbl.th_ppm_nist612 = NaN(height(outtbl),1);
    outtbl.u_ppm_nist612_se  = NaN(height(outtbl),1);
    outtbl.th_ppm_nist612_se = NaN(height(outtbl),1);
end

% ── Step 1b: Add GrainID and known concentrations ────────────────────────
outtbl = ladd_add_grainid(outtbl, metadataCsv);

mdJoin = readtable(metadataCsv, 'VariableNamingRule','preserve', ...
    'TextType','string', 'Delimiter',',');
mdJoin.Properties.VariableNames = lower(mdJoin.Properties.VariableNames);
assert(all(ismember({'file','known_u_ppm','known_th_ppm'}, ...
    mdJoin.Properties.VariableNames)), ...
    'ladd_reduce_zircon_pvmeas: enriched metadata lacks known U/Th columns.');

stripPath = @(f) char(regexp(string(f), '[^/\\]+$', 'match', 'once'));
mdFiles   = cellfun(stripPath, cellstr(mdJoin.file), 'UniformOutput', false);
outFiles  = cellfun(stripPath, cellstr(outtbl.file), 'UniformOutput', false);

knownU  = NaN(height(outtbl),1);
knownTh = NaN(height(outtbl),1);
sourceType = strings(height(outtbl),1);
sourceStdName = strings(height(outtbl),1);
for ii = 1:height(outtbl)
    jj = find(strcmp(outFiles{ii}, mdFiles), 1);
    if isempty(jj), continue; end
    knownU(ii)  = double(mdJoin.known_u_ppm(jj));
    knownTh(ii) = double(mdJoin.known_th_ppm(jj));
    sourceType(ii) = string(mdJoin.type(jj));
    sourceStdName(ii) = string(mdJoin.stdname(jj));
end
outtbl.known_u_ppm  = knownU;
outtbl.known_th_ppm = knownTh;
% A bridge-only core temporarily presents the bridge rows to reduce_core_pv
% as calibration anchors. Restore the exact source identities before any
% matching, reporting, or output is performed.
outtbl.type = sourceType;
outtbl.stdname = sourceStdName;
if ~hasIndependentNist
    % reduce_core_pv needs a temporary calibration identity to extract and
    % integrate the raw CPS consistently. Do not expose those internal rows
    % as real NIST memberships in a bridge-only output.
    outtbl.nist_u_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    outtbl.nist_th_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    outtbl.nist_sm_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    nistKCol = find(strcmpi(outtbl.Properties.VariableNames, ...
        'nist_review_kmad'),1);
    if ~isempty(nistKCol), outtbl{:,nistKCol} = NaN; end
    outtbl.nist_review_half_window(:) = NaN;
    outtbl.nist_auto_exclude_reviews(:) = false;
end
fprintf('  Joined known U/Th concentrations: %d/%d rows\n', ...
    sum(isfinite(knownU) & isfinite(knownTh)), height(outtbl));

% ── Step 2 (pre): Read and assign measured pit volumes ───────────────────
fprintf('ladd_reduce_zircon_pvmeas: reading pit volumes...\n');
hpv = read_pv_csv(hePitVolCsv,  'He');
if useAverageUthPit
    fprintf('  Using declared average UTh pit volume: %.6g +/- %.6g um3 (1SD)\n', ...
        uthPitVolSource, uthAverage1sd);
else
    upv = read_pv_csv(uthPitVolSource, 'UTh');
end

% GrainIDs are expected to be standardized before reduction. Matching is
% case-insensitive and ignores punctuation only; historical aliases are not
% inferred here.
gidsOut = cellfun(@cal_norm_id, cellstr(outtbl.grainid), 'UniformOutput', false);
gidsHpv = cellfun(@cal_norm_id, hpv.ids, 'UniformOutput', false);
if ~useAverageUthPit
    gidsUpv = cellfun(@cal_norm_id, upv.ids, 'UniformOutput', false);
end

nR = height(outtbl);
PV_He      = NaN(nR,1);
PV_He_1sd  = NaN(nR,1);
PV_UTh     = NaN(nR,1);
PV_UTh_1sd = NaN(nR,1);

isUnknown = strcmpi(string(outtbl.type), 'Unknown');
stdNameText = strtrim(string(outtbl.stdname));
grainText   = string(outtbl.grainid);
isExactBridge = strcmpi(stdNameText, strtrim(string(bridgeName)));
isBlankStd = ismissing(stdNameText) | strlength(stdNameText) == 0;
isGrainFallback = isBlankStd & contains(grainText, bridgeName, 'IgnoreCase', true);
isBridge = isExactBridge | isGrainFallback;

nBridgeMatched = sum(isBridge);
fprintf('  Rows matching bridge standard "%s": %d\n', bridgeName, nBridgeMatched);
if nBridgeMatched < 1
    error(['ladd_reduce_zircon_pvmeas: no rows matched bridgeStandardName="%s". ' ...
           'Check spelling against metadata stdname.'], bridgeName);
end

for ii = 1:nR
    if isUnknown(ii)
        ih = find(strcmp(gidsOut{ii}, gidsHpv), 1);
        if ~isempty(ih)
            PV_He(ii) = hpv.vols(ih);
            PV_He_1sd(ii) = hpv.sds(ih);
        end
        if useAverageUthPit
            PV_UTh(ii) = uthPitVolSource;
            PV_UTh_1sd(ii) = uthAverage1sd;
        else
            iu = find(strcmp(gidsOut{ii}, gidsUpv), 1);
            if ~isempty(iu)
                PV_UTh(ii) = upv.vols(iu);
                PV_UTh_1sd(ii) = upv.sds(iu);
            end
        end
    elseif isBridge(ii)
        if useAverageUthPit
            PV_UTh(ii) = uthPitVolSource;
            PV_UTh_1sd(ii) = uthAverage1sd;
        else
            iu = find(strcmp(gidsOut{ii}, gidsUpv), 1);
            if ~isempty(iu)
                PV_UTh(ii) = upv.vols(iu);
                PV_UTh_1sd(ii) = upv.sds(iu);
            end
        end
    end
end

PV_used = NaN(nR,1);
PV_used_1sd = NaN(nR,1);
hasBoth = isUnknown & isfinite(PV_UTh) & isfinite(PV_He);
hasUThOnly = isUnknown & isfinite(PV_UTh) & ~isfinite(PV_He);
hasNeither = isUnknown & ~isfinite(PV_UTh) & ...
    strlength(string(outtbl.grainid)) > 0;

PV_used(hasBoth) = PV_UTh(hasBoth) - PV_He(hasBoth);
PV_used_1sd(hasBoth) = sqrt(PV_UTh_1sd(hasBoth).^2 + PV_He_1sd(hasBoth).^2);
if any(hasUThOnly)
    fprintf('  WARNING: %d unknowns missing He PV; using full UTh PV (NO_HEPV)\n', ...
        sum(hasUThOnly));
    PV_used(hasUThOnly) = PV_UTh(hasUThOnly);
    PV_used_1sd(hasUThOnly) = PV_UTh_1sd(hasUThOnly);
end

badPV = hasBoth & PV_used <= 0;
if any(badPV)
    fprintf('  WARNING: %d unknowns have PV_UTh <= PV_He (BAD_PV)\n', sum(badPV));
    PV_used(badPV) = NaN;
end

PV_bridge = PV_UTh;
nUnknown = sum(isUnknown & strlength(string(outtbl.grainid)) > 0);
fprintf('  He PV matched:     %d/%d unknowns\n', sum(isfinite(PV_He) & isUnknown), nUnknown);
fprintf('  UTh PV matched:    %d/%d unknowns\n', sum(isfinite(PV_UTh) & isUnknown), nUnknown);
fprintf('  Bridge UTh matched:%d/%d rows\n', sum(isfinite(PV_bridge) & isBridge), nBridgeMatched);

% ── Step 2: Review and apply the zircon bridge calibration ───────────────
fprintf('ladd_reduce_zircon_pvmeas: applying reviewed bridge-standard calibration...\n');

cn = lower(outtbl.Properties.VariableNames);
colCpsU    = require_col(cn, 'cpsu');
colCpsTh   = require_col(cn, 'cpsth');
colCpsUSe  = require_col(cn, 'cpsu_se');
colCpsThSe = require_col(cn, 'cpsth_se');
colUPpm    = require_col(cn, 'u_ppm');
colThPpm   = require_col(cn, 'th_ppm');
colUPpmSe  = require_col(cn, 'u_ppm_se');
colThPpmSe = require_col(cn, 'th_ppm_se');

cpsU    = double(outtbl{:, colCpsU});
cpsTh   = double(outtbl{:, colCpsTh});
cpsUSe  = double(outtbl{:, colCpsUSe});
cpsThSe = double(outtbl{:, colCpsThSe});

isEligible = isBridge & ...
    isfinite(cpsU) & cpsU > 0 & ...
    isfinite(cpsTh) & cpsTh > 0 & ...
    isfinite(knownU) & knownU > 0 & ...
    isfinite(knownTh) & knownTh > 0 & ...
    isfinite(PV_bridge) & PV_bridge > 0;

nEligible = sum(isEligible);
fprintf('  Bridge rows eligible before review: %d\n', nEligible);
if nEligible < 1
    print_bridge_diagnostics(isBridge, cpsU, cpsTh, knownU, knownTh, ...
        PV_bridge, bridgeName);
    error(['ladd_reduce_zircon_pvmeas: no usable bridge rows. ' ...
           'See the diagnostic breakdown above.']);
end

seqAll = (1:nR)';
bridgeIdxAll = seqAll(isEligible);
impliedU  = knownU(bridgeIdxAll) ./ ...
    (cpsU(bridgeIdxAll) ./ PV_bridge(bridgeIdxAll));
impliedTh = knownTh(bridgeIdxAll) ./ ...
    (cpsTh(bridgeIdxAll) ./ PV_bridge(bridgeIdxAll));

[~, keepUIdx, nReviewU] = hampel_filter(impliedU, bridgeIdxAll, ...
    kMAD_std, halfWinStd);
[~, keepThIdx, nReviewTh] = hampel_filter(impliedTh, bridgeIdxAll, ...
    kMAD_std, halfWinStd);
keepHampel = intersect(keepUIdx, keepThIdx);
reviewedIdx = setdiff(bridgeIdxAll, keepHampel);

bridgeIDs = cellstr(string(outtbl.grainid));
excludeBridgeCell = cellstr(string(excludeBridgeIDs));
excludeBridgeCell = excludeBridgeCell(~cellfun(@isempty, excludeBridgeCell));
bridgeNorm = cellfun(@cal_norm_id, bridgeIDs, 'UniformOutput', false);
excludeBridgeNorm = cellfun(@cal_norm_id, excludeBridgeCell, 'UniformOutput', false);
manualExcludedIdx = bridgeIdxAll( ...
    ismember(bridgeNorm(bridgeIdxAll), excludeBridgeNorm));

if strcmp(anchorMode,'single')
    singleNorm = cal_norm_id(singleBridgeID);
    keepIdx = bridgeIdxAll(ismember(bridgeNorm(bridgeIdxAll),singleNorm));
    assert(numel(keepIdx) == 1, ...
        'ladd_reduce_zircon_pvmeas: singleBridgeID must match exactly one eligible bridge row.');
else
    keepIdx = setdiff(bridgeIdxAll, manualExcludedIdx);
end
if ~strcmp(anchorMode,'single') && numel(keepIdx) < 3
    error(['ladd_reduce_zircon_pvmeas: fewer than three bridge anchors remain ' ...
           'after explicit exclusions. Review excludeBridgeIDs.']);
end

fprintf('  Bridge Hampel review counts: U=%d, Th=%d, combined rows=%d\n', ...
    nReviewU, nReviewTh, numel(reviewedIdx));
fprintf('  Explicitly excluded bridge rows: %d\n', numel(manualExcludedIdx));
if ~isempty(reviewedIdx)
    fprintf('  Bridge rows flagged for review but retained: %s\n', ...
        strjoin(bridgeIDs(reviewedIdx), ', '));
end
if ~isempty(manualExcludedIdx)
    fprintf('  Bridge rows explicitly excluded: %s\n', ...
        strjoin(bridgeIDs(manualExcludedIdx), ', '));
end
if ~isempty(excludeBridgeNorm)
    unmatched = excludeBridgeNorm( ...
        ~ismember(excludeBridgeNorm, bridgeNorm(bridgeIdxAll)));
    if ~isempty(unmatched)
        warning('ladd_reduce_zircon_pvmeas:excludeBridgeIDsNotFound', ...
            'excludeBridgeIDs not found among eligible bridge rows: %s', ...
            strjoin(unmatched, ', '));
    end
end

bridgeSeqs = sort(keepIdx);
fprintf('  Bridge rows used as anchors: %d\n', numel(bridgeSeqs));

% Independent NIST612-vs-known bridge diagnostic, when the session actually
% contains NIST612 glass. A bridge-only session must not manufacture this QA
% comparison by comparing the bridge against itself.
if hasIndependentNist
    qaU  = outtbl.u_ppm_nist612(bridgeSeqs) ./ knownU(bridgeSeqs);
    qaTh = outtbl.th_ppm_nist612(bridgeSeqs) ./ knownTh(bridgeSeqs);
    medQaU  = median(qaU,  'omitnan');
    medQaTh = median(qaTh, 'omitnan');
    medNistBridgeU  = median(outtbl.u_ppm_nist612(bridgeSeqs), 'omitnan');
    medNistBridgeTh = median(outtbl.th_ppm_nist612(bridgeSeqs), 'omitnan');
    relNistBridgeU  = relative_robust_scatter(outtbl.u_ppm_nist612(bridgeSeqs));
    relNistBridgeTh = relative_robust_scatter(outtbl.th_ppm_nist612(bridgeSeqs));

    fprintf('  QA: bridge NIST612-anchored ppm / known ppm: U=%.2f, Th=%.2f\n', ...
        medQaU, medQaTh);
    nistMismatch = abs(medQaU - 1)*100 > nistTolPct || ...
                   abs(medQaTh - 1)*100 > nistTolPct;
    if nistMismatch
        fprintf(['  QA WARNING: bridge and NIST612 calibration differ by more than %.0f%%. ' ...
                 'Review matrix effects, NIST drift, bridge identity, and reference values.\n'], ...
                 nistTolPct);
    end
else
    medQaU = NaN;
    medQaTh = NaN;
    medNistBridgeU = NaN;
    medNistBridgeTh = NaN;
    relNistBridgeU = NaN;
    relNistBridgeTh = NaN;
    nistMismatch = false;
    fprintf('  QA: independent NIST612 comparison unavailable (bridge-only session)\n');
end

% Session-wide pooled factors (always computed for provenance).
KUi  = knownU(bridgeSeqs) ./ (cpsU(bridgeSeqs) ./ PV_bridge(bridgeSeqs));
KThi = knownTh(bridgeSeqs) ./ (cpsTh(bridgeSeqs) ./ PV_bridge(bridgeSeqs));
KUMed  = median(KUi,  'omitnan');
KThMed = median(KThi, 'omitnan');
relKU  = relative_robust_scatter(KUi);
relKTh = relative_robust_scatter(KThi);

if any(strcmp(anchorMode, {'median','single'}))
    fprintf('  anchorMode=%s pooled factors from %d bridge rows:\n', anchorMode, numel(bridgeSeqs));
    fprintf('    U:  %.5g (relative scatter %.1f%%)\n', KUMed, relKU*100);
    fprintf('    Th: %.5g (relative scatter %.1f%%)\n', KThMed, relKTh*100);
end

assignedBridge = NaN(nR,1);
usedFallback = false(nR,1);
if strcmp(anchorMode, 'following')
    for ii = 1:nR
        if ~isUnknown(ii), continue; end
        following = bridgeSeqs(bridgeSeqs > ii);
        if ~isempty(following)
            assignedBridge(ii) = following(1);
        else
            [~, jj] = min(abs(bridgeSeqs - ii));
            assignedBridge(ii) = bridgeSeqs(jj);
            usedFallback(ii) = true;
        end
    end
    if any(usedFallback & isUnknown)
        fprintf('  %d unknowns used last-bridge fallback (no following bridge)\n', ...
            sum(usedFallback & isUnknown));
    end
end
if strcmp(anchorMode,'single')
    assignedBridge(isUnknown) = bridgeSeqs(1);
end

uPpmNew    = double(outtbl{:, colUPpm});
thPpmNew   = double(outtbl{:, colThPpm});
uPpmSeNew  = double(outtbl{:, colUPpmSe});
thPpmSeNew = double(outtbl{:, colThPpmSe});

for ii = 1:nR
    if ~isUnknown(ii), continue; end
    if strcmp(anchorMode, 'following') && isnan(assignedBridge(ii)), continue; end
    pvUnknown = PV_used(ii);
    if ~(isfinite(pvUnknown) && pvUnknown > 0), continue; end

    pvUnknownRse = 0;
    if isfinite(PV_used_1sd(ii)) && PV_used_1sd(ii) > 0
        pvUnknownRse = PV_used_1sd(ii) / pvUnknown;
    end

    if any(strcmp(anchorMode, {'median','single'}))
        uPpmNew(ii) = (cpsU(ii) / pvUnknown) * KUMed;
        thPpmNew(ii)= (cpsTh(ii) / pvUnknown) * KThMed;
        uPpmSeNew(ii) = abs(uPpmNew(ii)) * sqrt( ...
            (cpsUSe(ii) / max(abs(cpsU(ii)), eps))^2 + ...
            relKU^2 + pvUnknownRse^2);
        thPpmSeNew(ii) = abs(thPpmNew(ii)) * sqrt( ...
            (cpsThSe(ii) / max(abs(cpsTh(ii)), eps))^2 + ...
            relKTh^2 + pvUnknownRse^2);
        continue;
    end

    bi = assignedBridge(ii);
    pvBridge = PV_bridge(bi);
    if ~(isfinite(pvBridge) && pvBridge > 0)
        error('Assigned bridge row %d (%s) has no matched UTh pit volume.', ...
            bi, string(outtbl.grainid(bi)));
    end
    pvBridgeRse = 0;
    if isfinite(PV_UTh_1sd(bi)) && PV_UTh_1sd(bi) > 0
        pvBridgeRse = PV_UTh_1sd(bi) / pvBridge;
    end

    uPpmNew(ii) = (cpsU(ii) / pvUnknown) / ...
        (cpsU(bi) / pvBridge) * knownU(bi);
    thPpmNew(ii) = (cpsTh(ii) / pvUnknown) / ...
        (cpsTh(bi) / pvBridge) * knownTh(bi);

    uPpmSeNew(ii) = abs(uPpmNew(ii)) * sqrt( ...
        (cpsUSe(ii) / max(abs(cpsU(ii)), eps))^2 + ...
        (cpsUSe(bi) / max(abs(cpsU(bi)), eps))^2 + ...
        pvUnknownRse^2 + pvBridgeRse^2);
    thPpmSeNew(ii) = abs(thPpmNew(ii)) * sqrt( ...
        (cpsThSe(ii) / max(abs(cpsTh(ii)), eps))^2 + ...
        (cpsThSe(bi) / max(abs(cpsTh(bi)), eps))^2 + ...
        pvUnknownRse^2 + pvBridgeRse^2);
end

outtbl{:, colUPpm}    = uPpmNew;
outtbl{:, colThPpm}   = thPpmNew;
outtbl{:, colUPpmSe}  = uPpmSeNew;
outtbl{:, colThPpmSe} = thPpmSeNew;
outtbl.bridge_seq = assignedBridge;
fprintf('  Bridge-standard calibration complete\n');

% ── Step 3: Store pit-volume columns ─────────────────────────────────────
outtbl.pv_he_um3      = PV_He;
outtbl.pv_he_1sd      = PV_He_1sd;
outtbl.pv_uth_um3     = PV_UTh;
outtbl.pv_uth_1sd     = PV_UTh_1sd;
outtbl.pv_used_um3    = PV_used;
outtbl.pv_used_1sd    = PV_used_1sd;

% ── Step 4: Convert final U/Th ppm to atoms ──────────────────────────────
fprintf('ladd_reduce_zircon_pvmeas: converting final ppm to atoms...\n');
grainMass = PV_used .* rho_zircon .* cm3_per_um3;

uAtoms  = (uPpmNew  / 1e6) / MW_U  .* NA .* grainMass;
thAtoms = (thPpmNew / 1e6) / MW_Th .* NA .* grainMass;

uPpmRelSe  = uPpmSeNew  ./ max(abs(uPpmNew), eps);
thPpmRelSe = thPpmSeNew ./ max(abs(thPpmNew), eps);
pvRelSe = PV_used_1sd ./ max(PV_used, eps);

uAtomsSe  = abs(uAtoms)  .* sqrt(uPpmRelSe.^2  + pvRelSe.^2);
thAtomsSe = abs(thAtoms) .* sqrt(thPpmRelSe.^2 + pvRelSe.^2);
uAtomsG    = uAtoms ./ grainMass;
thAtomsG   = thAtoms ./ grainMass;
uAtomsGSe  = uAtomsSe ./ grainMass;
thAtomsGSe = thAtomsSe ./ grainMass;

notUnknown = ~isUnknown;
uAtoms(notUnknown) = NaN;
thAtoms(notUnknown) = NaN;
uAtomsSe(notUnknown) = NaN;
thAtomsSe(notUnknown) = NaN;
uAtomsG(notUnknown) = NaN;
thAtomsG(notUnknown) = NaN;
uAtomsGSe(notUnknown) = NaN;
thAtomsGSe(notUnknown) = NaN;

outtbl.u_atoms = uAtoms;
outtbl.th_atoms = thAtoms;
outtbl.u_atoms_se = uAtomsSe;
outtbl.th_atoms_se = thAtomsSe;
outtbl.u_atoms_g = uAtomsG;
outtbl.th_atoms_g = thAtomsG;
outtbl.u_atoms_g_se = uAtomsGSe;
outtbl.th_atoms_g_se = thAtomsGSe;

% Same parent-production bookkeeping as the apatite hybrid, without Sm.
outtbl = add_parent_production(outtbl);

% ── Step 5: Flags and calibration provenance ─────────────────────────────
if ~ismember('flags', lower(outtbl.Properties.VariableNames))
    outtbl.flags = strings(nR,1);
end
flags = string(outtbl.flags);
noPit = ~isfinite(PV_used) & isUnknown & ...
    strlength(string(outtbl.grainid)) > 0;
flags(noPit) = strtrim(flags(noPit) + " NO_PITVOL");
flags(hasUThOnly) = strtrim(flags(hasUThOnly) + " NO_HEPV");
flags(badPV) = strtrim(flags(badPV) + " BAD_PV");
flags(hasNeither) = strtrim(flags(hasNeither) + " NO_UTHPV");
noParentSignal = isUnknown & ...
    (~isfinite(cpsU) | cpsU <= 0) & (~isfinite(cpsTh) | cpsTh <= 0);
flags(noParentSignal) = strtrim(flags(noParentSignal) + " NO_PARENT_SIGNAL");
flags(usedFallback & isUnknown) = strtrim( ...
    flags(usedFallback & isUnknown) + " BRIDGE_NEAREST_FALLBACK");
if ~isempty(reviewedIdx)
    flags(reviewedIdx) = strtrim(flags(reviewedIdx) + " STD_HAMPEL_REVIEW");
end
if ~isempty(manualExcludedIdx)
    flags(manualExcludedIdx) = strtrim( ...
        flags(manualExcludedIdx) + " STD_CAL_EXCLUDED_MANUAL");
end
if nistMismatch
    flags(bridgeIdxAll) = strtrim( ...
        flags(bridgeIdxAll) + " BRIDGE_NIST_MISMATCH_REVIEW");
end
outtbl.flags = flags;

bridgeRole = repmat("NOT_BRIDGE_STANDARD", nR,1);
bridgeRole(bridgeIdxAll) = "ANCHOR_USED";
bridgeRole(intersect(reviewedIdx, keepIdx)) = "ANCHOR_REVIEW_INCLUDED";
bridgeRole(manualExcludedIdx) = "EXCLUDED_MANUAL";
outtbl.bridge_anchor_role = bridgeRole;
outtbl.bridge_standard_used = repmat(string(bridgeName), nR,1);
outtbl.bridge_anchor_mode = repmat(string(anchorMode), nR,1);
outtbl.bridge_kMAD_std = repmat(kMAD_std, nR,1);
outtbl.bridge_halfWinStd = repmat(halfWinStd, nR,1);
outtbl.bridge_autoExcludeReviews = repmat(autoExcludeBridgeReviews, nR,1);
outtbl.bridge_nAnchorsUsed = repmat(numel(bridgeSeqs), nR,1);
outtbl.bridge_pooled_k_u = repmat(KUMed, nR,1);
outtbl.bridge_pooled_k_th = repmat(KThMed, nR,1);
outtbl.bridge_pooled_relscatter_u = repmat(relKU, nR,1);
outtbl.bridge_pooled_relscatter_th = repmat(relKTh, nR,1);
if hasIndependentNist
    nistReviewKUsed = kMAD;
    nistReviewHalfWinUsed = nistHalfWin;
else
    nistReviewKUsed = NaN;
    nistReviewHalfWinUsed = NaN;
end
outtbl.nist_review_kMAD_used = repmat(nistReviewKUsed, nR,1);
outtbl.nist_review_halfWin_used = repmat(nistReviewHalfWinUsed, nR,1);
outtbl.nist_autoExcludeReviews_used = repmat(autoExcludeNistReviews, nR,1);
outtbl.bridge_nist_median_u_ppm = repmat(medNistBridgeU, nR,1);
outtbl.bridge_nist_median_th_ppm = repmat(medNistBridgeTh, nR,1);
outtbl.bridge_nist_ratio_u = repmat(medQaU, nR,1);
outtbl.bridge_nist_ratio_th = repmat(medQaTh, nR,1);
outtbl.bridge_nist_relscatter_u = repmat(relNistBridgeU, nR,1);
outtbl.bridge_nist_relscatter_th = repmat(relNistBridgeTh, nR,1);
outtbl.bridge_nist_mismatch_review = repmat(nistMismatch, nR,1);
outtbl.parent_calibration_path = repmat(calibrationPath, nR,1);
outtbl.independent_nist_check_available = repmat(hasIndependentNist, nR,1);
outtbl.bridge_rows_used_for_core_extraction = repmat(nBridgeCoreRows, nR,1);
if useAverageUthPit
    uthPitMode = "SESSION_AVERAGE";
    uthPitSourceText = sprintf('%.15g +/- %.15g um3 (1SD)', ...
        uthPitVolSource, uthAverage1sd);
else
    uthPitMode = "MEASURED_PER_ANALYSIS";
    uthPitSourceText = char(string(uthPitVolSource));
end
outtbl.uth_pit_volume_mode = repmat(uthPitMode, nR,1);
outtbl.uth_pit_volume_source = repmat(string(uthPitSourceText), nR,1);
outtbl.metadata_rows_reclassified = repmat(nMetadataReclassified, nR,1);

nCalibrated = sum(isfinite(uAtoms) & isUnknown);
fprintf('\n-- Zircon reduction summary -------------------------------\n');
fprintf('  Total analyses:          %d\n', nR);
fprintf('  Unknown grains:          %d\n', nUnknown);
fprintf('  Successfully calibrated: %d (%.0f%%)\n', ...
    nCalibrated, 100*nCalibrated/max(nUnknown,1));
fprintf('-----------------------------------------------------------\n\n');

if ~isempty(saveAs)
    writetable(outtbl, saveAs);
    fprintf('Wrote: %s\n', saveAs);
end

end % main function


function [coreMetadataCsv, tempMetadataCsv, nTrueNist, nRelabeled, ...
    nBridgeCoreRows] = prepare_core_metadata( ...
        metadataCsv, bridgeName, allowBridgeOnly)
md = readtable(metadataCsv, 'VariableNamingRule','preserve', ...
    'TextType','string', 'Delimiter',',');
md.Properties.VariableNames = lower(md.Properties.VariableNames);
assert(all(ismember({'file','type'}, md.Properties.VariableNames)), ...
    'ladd_reduce_zircon_pvmeas: metadata must include file and type.');
if ~ismember('stdname', md.Properties.VariableNames)
    md.stdname = strings(height(md),1);
end

typeText = strtrim(string(md.type));
stdText = strtrim(string(md.stdname));
isTypeNist612 = strcmpi(typeText, 'NIST612');
isBlankStd = ismissing(stdText) | strlength(stdText) == 0;
isTrueNist = isTypeNist612 & (strcmpi(stdText, 'NIST612') | isBlankStd);
isInconsistentNist = isTypeNist612 & ~isTrueNist;
nTrueNist = sum(isTrueNist);
nRelabeled = 0;
nBridgeCoreRows = 0;

if any(isInconsistentNist)
    badNames = unique(stdText(isInconsistentNist));
    error(['ladd_reduce_zircon_pvmeas: inconsistent NIST metadata. Rows with ' ...
           'type=NIST612 must identify NIST612 glass, but found stdname(s): %s. ' ...
           'Correct the metadata before reduction.'], strjoin(badNames, ', '));
end

if nTrueNist < 1
    if ~allowBridgeOnly
        error(['ladd_reduce_zircon_pvmeas: no genuine NIST612 glass rows were found. ' ...
               'The function will not treat a zircon bridge as NIST glass unless ' ...
               'allowBridgeOnly=true is explicitly documented for the session.']);
    end
    isBridgeCore = strcmpi(stdText, strtrim(string(bridgeName)));
    nBridgeCoreRows = sum(isBridgeCore);
    assert(nBridgeCoreRows >= 3, ...
        ['ladd_reduce_zircon_pvmeas: bridge-only mode requires at least three ' ...
         'rows whose stdname exactly matches "%s"; found %d.'], ...
        bridgeName, nBridgeCoreRows);
    md.type(isBridgeCore) = "NIST612";
    md.stdname(isBridgeCore) = "NIST612";
end

coreMetadataCsv = metadataCsv;
tempMetadataCsv = '';
if nBridgeCoreRows > 0
    tempMetadataCsv = [tempname, '.csv'];
    writetable(md, tempMetadataCsv);
    coreMetadataCsv = tempMetadataCsv;
end
end


function delete_if_exists(pathText)
if ~isempty(pathText) && isfile(pathText)
    delete(pathText);
end
end


function pv = read_pv_csv(csvPath, label)
t = readtable(csvPath, 'VariableNamingRule','preserve', 'TextType','string');
t.Properties.VariableNames = lower(t.Properties.VariableNames);
cn = lower(regexprep(t.Properties.VariableNames, '[^a-zA-Z0-9]', ''));
gidCol = find_col(cn, {'grainid','id','name'});
volCol = find_col(cn, {'pitvol_um3','pitvolumeum3','pitvol','volume'});
sdCol  = find_col(cn, {'pv1sd_um3','pv1sd','propsd','sdum3','sd_um3'});
pv.ids = cellstr(string(t{:, gidCol}));
pv.vols = double(t{:, volCol});
pv.sds = double(t{:, sdCol});
fprintf('  Read %s pit-volume CSV: %d rows from %s\n', label, height(t), csvPath);
end


function idx = find_col(cn, aliases)
for kk = 1:numel(aliases)
    key = regexprep(char(aliases{kk}), '[^a-zA-Z0-9]', '');
    hit = find(strcmpi(cn, key), 1);
    if ~isempty(hit)
        idx = hit;
        return;
    end
end
error('ladd_reduce_zircon_pvmeas: cannot find column matching [%s]', ...
    strjoin(aliases, '/'));
end


function idx = require_col(cn, name)
idx = find(strcmpi(cn, name), 1);
if isempty(idx)
    error('ladd_reduce_zircon_pvmeas: required output column "%s" is missing.', name);
end
end


function idNorm = cal_norm_id(id)
idNorm = lower(strtrim(char(string(id))));
idNorm = regexprep(idNorm, '[^a-z0-9]+', '');
% Treat zero-padded and unpadded trailing grain numbers as the same ID
% (for example DS_25_z04 and DS_25_z4).
idNorm = regexprep(idNorm, '([a-z])0+([0-9]+)$', '$1$2');
end


function [KClean, seqClean, nReview] = hampel_filter(K, seq, kMAD, halfWin)
if nargin < 3, kMAD = 5; end
if nargin < 4, halfWin = 12; end
n = numel(K);
keep = true(n,1);
for ii = 1:n
    lo = max(1, ii-halfWin);
    hi = min(n, ii+halfWin);
    neighbours = K([lo:ii-1, ii+1:hi]);
    neighbours = neighbours(isfinite(neighbours));
    if numel(neighbours) < 2, continue; end
    med = median(neighbours, 'omitnan');
    sigma = 1.4826 * mad(neighbours, 1);
    if sigma > 0 && abs(K(ii)-med) > kMAD*sigma
        keep(ii) = false;
    end
end
nReview = sum(~keep);
if sum(keep) < 3
    warning('hampel_filter: review would leave fewer than three anchors; retaining all.');
    keep(:) = true;
    nReview = 0;
end
KClean = K(keep);
seqClean = seq(keep);
end


function rel = relative_robust_scatter(values)
values = values(isfinite(values));
if isempty(values)
    rel = NaN;
    return;
end
med = median(values, 'omitnan');
rel = 1.4826 * mad(values, 1) / max(eps, abs(med));
end


function print_bridge_diagnostics(isBridge, cpsU, cpsTh, knownU, knownTh, ...
    pvBridge, bridgeName)
fprintf('  Diagnostic breakdown for %d rows matching "%s":\n', ...
    sum(isBridge), bridgeName);
fprintf('    cpsu > 0:          %d\n', sum(isBridge & isfinite(cpsU) & cpsU > 0));
fprintf('    cpsth > 0:         %d\n', sum(isBridge & isfinite(cpsTh) & cpsTh > 0));
fprintf('    known_u_ppm > 0:   %d\n', sum(isBridge & isfinite(knownU) & knownU > 0));
fprintf('    known_th_ppm > 0:  %d\n', sum(isBridge & isfinite(knownTh) & knownTh > 0));
fprintf('    bridge PV > 0:     %d\n', ...
    sum(isBridge & isfinite(pvBridge) & pvBridge > 0));
end


function outtbl = add_parent_production(outtbl)
f238 = 0.992742;
f235 = 0.007204;
lam238 = 1.55125e-10;
lam235 = 9.8485e-10;
lam232 = 4.9475e-11;

U = grab(outtbl, 'u_atoms_g');
USe = grab(outtbl, 'u_atoms_g_se');
Th = grab(outtbl, 'th_atoms_g');
ThSe = grab(outtbl, 'th_atoms_g_se');

N238 = f238 .* U;
N238Se = f238 .* USe;
N235 = f235 .* U;
N235Se = f235 .* USe;
N232 = Th;
N232Se = ThSe;

production = 8*lam238.*N238 + 7*lam235.*N235 + 6*lam232.*N232;
productionSe = sqrt((8*lam238.*N238Se).^2 + ...
    (7*lam235.*N235Se).^2 + (6*lam232.*N232Se).^2);

outtbl.n238_atoms_g = N238;
outtbl.n238_atoms_g_se = N238Se;
outtbl.n235_atoms_g = N235;
outtbl.n235_atoms_g_se = N235Se;
outtbl.n232_atoms_g = N232;
outtbl.n232_atoms_g_se = N232Se;
outtbl.parentprod_atoms_g_yr = production;
outtbl.parentprod_atoms_g_yr_se = productionSe;
end


function values = grab(tbl, name)
cn = tbl.Properties.VariableNames;
idx = find(strcmpi(cn, name), 1);
if isempty(idx)
    values = zeros(height(tbl),1);
else
    values = double(tbl{:, idx});
    values(~isfinite(values)) = 0;
end
end
