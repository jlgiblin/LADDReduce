function outtbl = ladd_reduce_apatite_pvmeas_hybrid(folder, metadataCsv, hePitVolCsv, ...
    uthPitVolSource, varargin)
% LADD_REDUCE_APATITE_PVMEAS_HYBRID  UThSm reduction for apatite, combining
% the excel-match bracket arithmetic with the Hampel-filtered/NIST-tracked
% robustness of ladd_reduce_apatite_pvmeas.
%
% WHAT THIS IS
%   Core bridge-standard calibration formula:
%       X_ppm_unknown = (X_cps_unknown/PV_unknown) / ...
%                       (X_cps_bridge/PV_bridge) * known_X_ppm_bridge
%   applied identically to U, Th, and Sm. NIST-anchored results are retained
%   as an independent QA comparison rather than substituted into the
%   mineral-bridge calculation.
%
%   Robustness layered in from ladd_reduce_apatite_pvmeas:
%     1. When NIST612 was analyzed, its Hampel-reviewed drift correction runs
%        first (via reduce_core_pv, mode='nois') and its ppm values are retained as
%        u_ppm_nist612 / th_ppm_nist612 / sm_ppm_nist612 — a visible,
%        trackable QC cross-check — but are NOT used as the bracket
%        multiplier, preserving the excel-match arithmetic. NIST review
%        rows remain included by default; only explicit filename exclusions
%        remove them.
%     2. Bridge-standard anchors receive a Hampel REVIEW before assignment.
%        A replicate whose implied ppm (known_ppm / (cps/PV)) deviates from
%        local neighbours is flagged STD_HAMPEL_REVIEW but remains included
%        by default. Only excludeBridgeIDs removes a documented bad analysis.
%     3. A QA print (and optional flag) compares each bridge-standard
%        replicate's NIST-anchored ppm to its declared known ppm — an
%        early-warning signal for a standard-identity or drift problem.
%
%   Parent-anchor selection is ALSO a call-time switch via 'anchorMode':
%     'nearest' (default) — original excel-match behavior: each unknown is
%        calibrated against the single nearest-following included
%        standard replicate. This is most appropriate when the bridge
%        standard is demonstrably homogeneous at the analyzed scale.
%     'median'  — each unknown is calibrated against a SINGLE session-wide
%        pooled factor: the median of known_ppm/(cps/PV) across every
%        included anchor, applied identically to U, Th, and Sm.
%        This pooled calibration can be more stable than a single-analysis
%        bracket when the bridge standard has real replicate scatter. It is
%        applied symmetrically to all three
%        elements instead of singling out Th, and built from the SAME
%        anchor pool/formula as 'nearest' rather than a separate mechanism.
%        Appropriate when the selected bridge standard has meaningful
%        replicate-to-replicate spread.
%     'nist_following' — direct NIST612 bracketing:
%        each unknown is calibrated against the next following NIST612
%        analysis using the same cps/PV ratio equation. The terminal block
%        falls back to the nearest NIST612. This bypasses the mineral bridge
%        standard while retaining the same raw-signal and pit-volume path.
%     'nist_interpolated' — retain reduce_core_pv's time-interpolated
%        NIST612 U, Th, and Sm concentrations as the final parent
%        calibration. This is the continuous-drift NIST control.
%
%   Bridge standard is a CALL-TIME switch, not a metadata edit.
%   'bridgeStandardName' is required and picks the exact stdname that acts
%   as the calibration anchor for a given run.
%   To select a bridge, pass its exact metadata name, for example
%   'bridgeStandardName','ReferenceMaterial'. The selected standard must
%   already have its own
%   known_u_ppm/known_th_ppm/known_sm_ppm populated under its real
%   stdname in the metadata.
%
% USAGE
%   out = ladd_reduce_apatite_pvmeas_hybrid(folder, metadataCsv, ...
%             hePitVolCsv, uthPitVolCsv, ...
%             'bridgeStandardName','ReferenceMaterial', ...
%             'smReferenceBasis','total', 'saveAs', 'UThSm_out.csv');
%   out = ladd_reduce_apatite_pvmeas_hybrid(folder, metadataCsv, ...
%             hePitVolCsv, uthPitVolCsv, ...
%             'bridgeStandardName', 'ReferenceMaterial', ...
%             'smReferenceBasis','total', 'anchorMode', 'median');
%
% OPTIONAL NAME-VALUE PAIRS
%   'saveAs'             : output CSV filename (default '' = no save)
%   'bridgeStandardName' : REQUIRED exact stdname used as the calibration
%                          bridge standard. GrainID substring matching is
%                          retained only as a fallback when stdname is blank.
%   'anchorMode'         : 'nearest' | 'median' | 'nist_following' |
%                          'nist_interpolated'
%                          (default 'nearest') — see above
%   'smReferenceBasis'   : REQUIRED 'total' or '147isotope'. This declares
%                          whether known_sm_ppm for the selected bridge is
%                          elemental Sm ppm or isotope-specific 147Sm ppm.
%   'parentScalarInterp' : 'linear'|'pchip' (default 'linear') — passed to reduce_core_pv
%   'kMAD'               : MAD multiplier for ablation-window plateau detection (default 6)
%                          and loose NIST-anchor review
%   'nistHalfWin'        : half-window for NIST-anchor review (default 12)
%   'autoExcludeNistReviews' : must remain false under the preserve-and-review
%                          policy; NIST review rows remain in calibration
%   'excludeNistFiles'   : explicit NIST measurement filenames to exclude
%                          (default empty), for documented misfires only
%   'kMAD_std'            : MAD multiplier for bridge-standard review (default 5)
%   'halfWinStd'          : bridge review half-window (default 12)
%   'autoExcludeBridgeReviews' : retained only for call compatibility and
%                          must remain false. Review flags cannot auto-exclude.
%   'excludeBridgeIDs'    : explicit bridge-standard GrainIDs to exclude
%                          (default empty). Use only for independently
%                          documented analytical problems.
%   Review thresholds do not automatically exclude bridge analyses. Validate
%   these settings for the reference material and acquisition design in use.
%   'SmMaxPlaus'         : max plausible Sm ppm, for NIST-QC column only (default 7500)
%   'SmMinCps'            : min Sm cps, for NIST-QC column only (default 150)
%   'SmRelSEmax'          : max relative SE on Sm, for NIST-QC column only (default 0.40)
%   'nistCheckTolPct'     : QA warn threshold for NIST-vs-known bridge mismatch, percent (default 25)
%   'allowBridgeOnly'     : false by default. Set true for a documented
%                          bridge-calibrated session with no NIST612 glass.
%                          Direct NIST anchor modes still require NIST612.
%   'uthAverage1sd'       : required positive 1SD (um3) when uthPitVolSource
%                          is a numeric session average

% ── Parse arguments ────────────────────────────────────────────────────────
p = inputParser;
p.addParameter('saveAs',             '',       @ischar);
p.addParameter('bridgeStandardName', '',       @(x)ischar(x)||isstring(x));
p.addParameter('anchorMode',         'nearest',@(x)ischar(x)||isstring(x));
p.addParameter('smReferenceBasis',   '',       @(x)ischar(x)||isstring(x));
p.addParameter('parentScalarInterp', 'linear', @(x)ischar(x)||isstring(x));
p.addParameter('kMAD',               6,        @(x)isnumeric(x)&&isscalar(x));
p.addParameter('nistHalfWin',        12,       @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('autoExcludeNistReviews', false, ...
    @(x)islogical(x)||isnumeric(x));
p.addParameter('excludeNistFiles', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.addParameter('kMAD_std',           5,        @(x)isnumeric(x)&&isscalar(x));
p.addParameter('halfWinStd',         12,       @(x)isnumeric(x)&&isscalar(x));
p.addParameter('autoExcludeBridgeReviews', false, ...
    @(x)islogical(x)||isnumeric(x));
p.addParameter('excludeBridgeIDs', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.addParameter('SmMaxPlaus',         7500,     @(x)isnumeric(x)&&isscalar(x));
p.addParameter('SmMinCps',           150,      @(x)isnumeric(x)&&isscalar(x));
p.addParameter('SmRelSEmax',         0.40,     @(x)isnumeric(x)&&isscalar(x));
p.addParameter('nistCheckTolPct',    25,       @(x)isnumeric(x)&&isscalar(x));
p.addParameter('allowBridgeOnly',    false,    @(x)islogical(x)||isnumeric(x));
p.addParameter('uthAverage1sd',      NaN,      @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});

saveAs      = p.Results.saveAs;
bridgeName  = char(string(p.Results.bridgeStandardName));
anchorMode  = lower(char(string(p.Results.anchorMode)));
smReferenceBasis = lower(char(string(p.Results.smReferenceBasis)));
allowBridgeOnly = logical(p.Results.allowBridgeOnly);
uthAverage1sd = p.Results.uthAverage1sd;
useAverageUthPit = isnumeric(uthPitVolSource);
if useAverageUthPit
    assert(isscalar(uthPitVolSource) && isfinite(uthPitVolSource) && ...
        uthPitVolSource > 0, ...
        'Numeric uthPitVolSource must be one positive average volume in um3.');
    assert(isfinite(uthAverage1sd) && uthAverage1sd > 0, ...
        ['uthAverage1sd must be explicitly supplied and positive when using ', ...
         'an average U-Th pit volume.']);
end
assert(any(strcmp(anchorMode,{'nearest','median','nist_following','nist_interpolated'})), ...
    ['ladd_reduce_apatite_pvmeas_hybrid: anchorMode must be ''nearest'', ', ...
     '''median'', ''nist_following'', or ''nist_interpolated''.']);
usesDirectNist = any(strcmp(anchorMode,{'nist_following','nist_interpolated'}));
if usesDirectNist
    % Direct NIST control: the calibration anchor is NIST612 itself, not
    % the mineral bridge named by bridgeStandardName.
    bridgeName = 'NIST612';
else
    assert(~isempty(strtrim(bridgeName)), ...
        ['ladd_reduce_apatite_pvmeas_hybrid: bridgeStandardName is required. ' ...
         'Pass the full stdname explicitly (for example, ''ReferenceMaterial'').']);
end
assert(any(strcmp(smReferenceBasis,{'total','147isotope'})), ...
    ['ladd_reduce_apatite_pvmeas_hybrid: smReferenceBasis is required and must be ' ...
     '''total'' or ''147isotope''. Check the source of known_sm_ppm for the selected bridge.']);
interpMth   = p.Results.parentScalarInterp;
kMAD        = p.Results.kMAD;
nistHalfWin = p.Results.nistHalfWin;
autoExcludeNistReviews = logical(p.Results.autoExcludeNistReviews);
excludeNistFiles = p.Results.excludeNistFiles;
assert(~autoExcludeNistReviews, ...
    ['ladd_reduce_apatite_pvmeas_hybrid: autoExcludeNistReviews=true is disabled under ' ...
     'the preserve-and-review policy. Name documented bad NIST files with excludeNistFiles instead.']);
kMAD_std    = p.Results.kMAD_std;
halfWinStd  = p.Results.halfWinStd;
autoExcludeBridgeReviews = logical(p.Results.autoExcludeBridgeReviews);
excludeBridgeIDs = p.Results.excludeBridgeIDs;
assert(~autoExcludeBridgeReviews, ...
    ['ladd_reduce_apatite_pvmeas_hybrid: autoExcludeBridgeReviews=true is disabled under ' ...
     'the preserve-and-review policy. Name documented bad anchors with excludeBridgeIDs instead.']);
nistTolPct  = p.Results.nistCheckTolPct;
fprintf('  Parent anchor standard for this run: "%s" (exact stdname; GrainID fallback)\n', bridgeName);
fprintf('  Parent anchor mode: "%s"\n', anchorMode);
fprintf('  Bridge known_sm_ppm basis: "%s"\n', smReferenceBasis);
fprintf('  Bridge review policy: Hampel review only; automatic exclusion=%d\n', ...
    autoExcludeBridgeReviews);
fprintf('  NIST review policy: Hampel review only; automatic exclusion=%d\n', ...
    autoExcludeNistReviews);

NA           = 6.02214076e23;
MW_U         = 238.02891;
MW_Th        = 232.03806;
if strcmp(smReferenceBasis, 'total')
    MW_Sm = 150.36;
else
    MW_Sm = 147;
end
rho_apatite  = 3.19;
cm3_per_um3  = 1e-12;

% ── Step 1: Run reduce_core_pv (NIST612 Hampel-filtered drift correction) ──
[coreMetadataCsv, tempMetadataCsv, nTrueNist, nBridgeCoreRows] = ...
    prepare_core_metadata(metadataCsv, bridgeName, allowBridgeOnly, usesDirectNist);
cleanupCoreMetadata = onCleanup(@() delete_if_exists(tempMetadataCsv));
hasIndependentNist = nTrueNist > 0;
if usesDirectNist
    calibrationPath = "NIST612_DIRECT";
elseif hasIndependentNist
    calibrationPath = "NIST612_THEN_BRIDGE";
else
    calibrationPath = "BRIDGE_ONLY_NO_NIST";
end
fprintf('  Genuine NIST612 rows available: %d\n', nTrueNist);
if ~hasIndependentNist
    fprintf(['  Bridge-only session: %d "%s" rows supply the parent ', ...
        'calibration; no independent NIST612 comparison is available.\n'], ...
        nBridgeCoreRows, bridgeName);
end
fprintf('ladd_reduce_apatite_pvmeas_hybrid: extracting parent signals...\n');

cfg = struct( ...
    'mineral',             'apatite',  ...
    'density',             3.19,       ...
    'mode',                'nois',     ...
    'IS_element',          'Ca',       ...
    'IS_ppm_std',          0,          ...
    'IS_ppm_unk',          3.99e5,     ...
    'pick',                struct( ...
        'U',  "238U",  ...
        'Th', "232Th", ...
        'IS', "44Ca",  ...
        'Sm', "147Sm"), ...
    'includeSm',           true,       ...
    'Sm_isotope',          '147Sm',    ...
    'Sm_known_is_isotope', false,      ... % Excel-match: keep 147Sm ppm as isotope ppm
    'matrixScalarName',    "NONE",     ...  % U/Th/Sm remain symmetric
    'SmMaxPlaus',          p.Results.SmMaxPlaus, ...
    'SmMinCps',            p.Results.SmMinCps,   ...
    'SmRelSEmax',          p.Results.SmRelSEmax  ...
);

outtbl = reduce_core_pv(folder, coreMetadataCsv, cfg, ...
    'parentScalarInterp', interpMth, ...
    'kMAD', kMAD, ...
    'nistHalfWin', nistHalfWin, ...
    'showNistSummary', hasIndependentNist, ...
    'autoExcludeNistReviews', autoExcludeNistReviews, ...
    'excludeNistFiles', excludeNistFiles);
clear cleanupCoreMetadata
fprintf('  reduce_core_pv complete: %d rows\n', height(outtbl));

% Stash the NIST612-anchored ppm as QC columns BEFORE anything overwrites them.
% These are never used as the bracket multiplier — they are a visible,
% trackable cross-check against the bridge-standard-based calibration below.
% NOTE: reduce_core_pv returns capitalized names (U_ppm/Th_ppm/Sm_ppm), not
% lowercase — look them up case-insensitively rather than using dot access,
% same pattern used everywhere else in this file.
cn0     = lower(outtbl.Properties.VariableNames);
col_u0  = find(strcmp(cn0,'u_ppm'),  1);
col_th0 = find(strcmp(cn0,'th_ppm'), 1);
col_sm0 = find(strcmp(cn0,'sm_ppm'), 1);
if hasIndependentNist
    outtbl.u_ppm_nist612  = double(outtbl{:,col_u0});
    outtbl.th_ppm_nist612 = double(outtbl{:,col_th0});
    if ~isempty(col_sm0)
        outtbl.sm_ppm_nist612 = double(outtbl{:,col_sm0});
    else
        outtbl.sm_ppm_nist612 = NaN(height(outtbl),1);
    end
else
    outtbl.u_ppm_nist612 = NaN(height(outtbl),1);
    outtbl.th_ppm_nist612 = NaN(height(outtbl),1);
    outtbl.sm_ppm_nist612 = NaN(height(outtbl),1);
end

% ── Step 1b: Add GrainID ──────────────────────────────────────────────────
outtbl = ladd_add_grainid(outtbl, metadataCsv);

% ── Step 1c: Join known concentrations from metadata ──────────────────────
md_join = readtable(metadataCsv, 'VariableNamingRule','preserve', ...
    'TextType','string', 'Delimiter',',');
md_join.Properties.VariableNames = lower(md_join.Properties.VariableNames);
assert(all(ismember({'file','type','stdname','known_u_ppm','known_th_ppm'}, ...
    md_join.Properties.VariableNames)), ...
    'ladd_reduce_apatite_pvmeas_hybrid: metadata lacks required columns.');
for optionalColumn = {'known_u_1sd_ppm','known_th_1sd_ppm','known_sm_1sd_ppm'}
    if ~ismember(optionalColumn{1},md_join.Properties.VariableNames)
        md_join.(optionalColumn{1}) = NaN(height(md_join),1);
    end
end
strip_path = @(f) char(regexp(string(f), '[^/\\]+$', 'match', 'once'));
md_files_j  = cellfun(strip_path, cellstr(md_join.file), 'UniformOutput', false);
out_files_j = cellfun(strip_path, cellstr(outtbl.file),  'UniformOutput', false);

ku_join=NaN(height(outtbl),1); kth_join=NaN(height(outtbl),1); ksm_join=NaN(height(outtbl),1);
ku1sd_join=NaN(height(outtbl),1); kth1sd_join=NaN(height(outtbl),1); ksm1sd_join=NaN(height(outtbl),1);
sourceType = strings(height(outtbl),1); sourceStdName = strings(height(outtbl),1);
for ii = 1:height(outtbl)
    idx_j = find(strcmp(out_files_j{ii}, md_files_j), 1);
    if ~isempty(idx_j)
        ku_join(ii)  = double(md_join.known_u_ppm(idx_j));
        kth_join(ii) = double(md_join.known_th_ppm(idx_j));
        ku1sd_join(ii)  = double(md_join.known_u_1sd_ppm(idx_j));
        kth1sd_join(ii) = double(md_join.known_th_1sd_ppm(idx_j));
        if ismember('known_sm_ppm', md_join.Properties.VariableNames)
            ksm_join(ii) = double(md_join.known_sm_ppm(idx_j));
        end
        ksm1sd_join(ii) = double(md_join.known_sm_1sd_ppm(idx_j));
        sourceType(ii) = string(md_join.type(idx_j));
        sourceStdName(ii) = string(md_join.stdname(idx_j));
    end
end
outtbl.known_u_ppm  = ku_join;
outtbl.known_th_ppm = kth_join;
outtbl.known_sm_ppm = ksm_join;
outtbl.known_u_1sd_ppm  = ku1sd_join;
outtbl.known_th_1sd_ppm = kth1sd_join;
outtbl.known_sm_1sd_ppm = ksm1sd_join;
outtbl.type = sourceType;
outtbl.stdname = sourceStdName;
if ~hasIndependentNist
    outtbl.nist_u_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    outtbl.nist_th_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    outtbl.nist_sm_anchor_role(:) = "NOT_AVAILABLE_BRIDGE_ONLY";
    outtbl.nist_review_kmad(:) = NaN;
    outtbl.nist_review_half_window(:) = NaN;
    outtbl.nist_auto_exclude_reviews(:) = false;
end
fprintf('  Joined known ppm: %d/%d rows matched\n', sum(isfinite(ku_join)), height(outtbl));

% ── Step 2 (pre): Read and assign pit volumes ─────────────────────────────
fprintf('ladd_reduce_apatite_pvmeas_hybrid: reading pit volumes...\n');
hpv = read_pv_csv(hePitVolCsv,  'He');
if useAverageUthPit
    fprintf('  Using declared average UTh pit volume: %.6g +/- %.6g um3 (1SD)\n', ...
        uthPitVolSource, uthAverage1sd);
else
    upv = read_pv_csv(uthPitVolSource, 'UTh');
end

norm_id  = @(s) lower(regexprep(strtrim(char(string(s))), '[\s_\-]+', '-'));
gids_out = cellfun(@(x) norm_id(x), cellstr(outtbl.grainid), 'UniformOutput', false);
gids_hpv = cellfun(@(x) norm_id(x), hpv.ids,                'UniformOutput', false);
if ~useAverageUthPit
    gids_upv = cellfun(@(x) norm_id(x), upv.ids, 'UniformOutput', false);
end

nR     = height(outtbl);
PV_He  = NaN(nR,1); PV_He_1sd  = NaN(nR,1);
PV_UTh = NaN(nR,1); PV_UTh_1sd = NaN(nR,1);

is_unknown_pre = strcmpi(string(outtbl.type), 'Unknown');

for i = 1:nR
    if ~is_unknown_pre(i), continue; end
    ih = find(strcmp(gids_out{i}, gids_hpv), 1);
    if ~isempty(ih), PV_He(i)  = hpv.vols(ih); PV_He_1sd(i)  = hpv.sds(ih);  end
    if useAverageUthPit
        PV_UTh(i) = uthPitVolSource;
        PV_UTh_1sd(i) = uthAverage1sd;
    else
        iu = find(strcmp(gids_out{i}, gids_upv), 1);
        if ~isempty(iu), PV_UTh(i) = upv.vols(iu); PV_UTh_1sd(i) = upv.sds(iu); end
    end
end

% UTh pit volumes for bridge-standard rows (no He pit — use full UTh PV)
% Match the declared stdname exactly. A GrainID substring is accepted only
% where stdname is blank, for older metadata files that did not populate it.
stdNameText = strtrim(string(outtbl.stdname));
grainText   = string(outtbl.grainid);
isExactStd  = strcmpi(stdNameText, strtrim(string(bridgeName)));
isBlankStd  = ismissing(stdNameText) | strlength(stdNameText) == 0;
isGrainFallback = isBlankStd & contains(grainText, bridgeName, 'IgnoreCase', true);
is_bridge_pre = isExactStd | isGrainFallback;
n_bridge_rows = sum(is_bridge_pre);
fprintf('  Rows matching bridge standard "%s": %d\n', bridgeName, n_bridge_rows);
if n_bridge_rows < 1
    error(['ladd_reduce_apatite_pvmeas_hybrid: no rows matched bridgeStandardName="%s". ', ...
           'Check spelling against the stdname column in your metadata CSV.'], bridgeName);
end
for i = 1:nR
    if ~is_bridge_pre(i), continue; end
    if useAverageUthPit
        PV_UTh(i) = uthPitVolSource;
        PV_UTh_1sd(i) = uthAverage1sd;
    else
        iu = find(strcmp(gids_out{i}, gids_upv), 1);
        if ~isempty(iu), PV_UTh(i) = upv.vols(iu); PV_UTh_1sd(i) = upv.sds(iu); end
    end
end

PV_used     = NaN(nR,1); PV_used_1sd = NaN(nR,1);
has_both     = isfinite(PV_UTh) & isfinite(PV_He)  & is_unknown_pre;
has_uth_only = isfinite(PV_UTh) & ~isfinite(PV_He) & is_unknown_pre;
has_neither  = ~isfinite(PV_UTh) & is_unknown_pre & strlength(string(outtbl.grainid)) > 0;

PV_used(has_both)     = PV_UTh(has_both) - PV_He(has_both);
PV_used_1sd(has_both) = sqrt(PV_UTh_1sd(has_both).^2 + PV_He_1sd(has_both).^2);
if any(has_uth_only)
    fprintf('  WARNING: %d unknowns missing He PV — using full UTh PV (NO_HEPV)\n', sum(has_uth_only));
    PV_used(has_uth_only)     = PV_UTh(has_uth_only);
    PV_used_1sd(has_uth_only) = PV_UTh_1sd(has_uth_only);
end
bad_res = has_both & PV_used <= 0;
if any(bad_res)
    fprintf('  WARNING: %d grains PV_UTh <= PV_He — flagged BAD_PV\n', sum(bad_res));
    PV_used(bad_res) = NaN;
end

PV_used_bridge = PV_UTh;   % NaN for non-bridge rows; filled for bridge-standard rows above

n_unk = sum(is_unknown_pre & strlength(string(outtbl.grainid)) > 0);
fprintf('  He PV matched:  %d/%d unknowns\n', sum(isfinite(PV_He)  & is_unknown_pre), n_unk);
fprintf('  UTh PV matched: %d/%d unknowns\n', sum(isfinite(PV_UTh) & is_unknown_pre), n_unk);

% ── Step 2: Hampel-reviewed, excel-style bridge-standard bracketing ───────
fprintf('ladd_reduce_apatite_pvmeas_hybrid: applying Hampel-reviewed bridge-standard bracketing...\n');

is_bridge = is_bridge_pre;
is_unknown  = is_unknown_pre;
seq_all     = (1:height(outtbl))';
cn_tbl      = lower(outtbl.Properties.VariableNames);

col_cpsu     = find(strcmp(cn_tbl,'cpsu'),      1);
col_cpsth    = find(strcmp(cn_tbl,'cpsth'),     1);
col_cpssm    = find(strcmp(cn_tbl,'cpssm'),     1);
col_cpsu_se  = find(strcmp(cn_tbl,'cpsu_se'),   1);
col_cpsth_se = find(strcmp(cn_tbl,'cpsth_se'),  1);
col_cpssm_se = find(strcmp(cn_tbl,'cpssm_se'),  1);

col_uppm     = find(strcmp(cn_tbl,'u_ppm'),     1);
col_thppm    = find(strcmp(cn_tbl,'th_ppm'),    1);
col_smppm    = find(strcmp(cn_tbl,'sm_ppm'),    1);
col_uppmse   = find(strcmp(cn_tbl,'u_ppm_se'),  1);
col_thppmse  = find(strcmp(cn_tbl,'th_ppm_se'), 1);
col_smppmse  = find(strcmp(cn_tbl,'sm_ppm_se'), 1);

col_ku       = find(strcmp(cn_tbl,'known_u_ppm'),  1);
col_kth      = find(strcmp(cn_tbl,'known_th_ppm'), 1);
col_ksm      = find(strcmp(cn_tbl,'known_sm_ppm'), 1);

hasSm = ~isempty(col_cpssm) && ~isempty(col_smppm) && ~isempty(col_ksm);

cpsu_all     = double(outtbl{:,col_cpsu});
cpsth_all    = double(outtbl{:,col_cpsth});
cpsu_se_all  = double(outtbl{:,col_cpsu_se});
cpsth_se_all = double(outtbl{:,col_cpsth_se});

if hasSm
    cpssm_all = double(outtbl{:,col_cpssm});
    if ~isempty(col_cpssm_se)
        cpssm_se_all = double(outtbl{:,col_cpssm_se});
    else
        cpssm_se_all = zeros(height(outtbl),1);
    end
end

known_u_all  = double(outtbl{:,col_ku});
known_th_all = double(outtbl{:,col_kth});

% ---- Initial anchor pool: U/Th eligibility only ---------------------------
% Sm eligibility is intentionally NOT required here. A standard missing only
% its known_sm_ppm reference (e.g. a lookup-table gap) should not disqualify
% it from anchoring U/Th — Sm is handled as its own, independent check below,
% and any unknown bracketed to a no-Sm-reference anchor simply keeps its
% NIST612-anchored Sm value (sm_ppm_nist612) instead of being bracket-corrected.
is_valid_bridge = is_bridge & ...
               isfinite(cpsu_all)  & cpsu_all  > 0 & ...
               isfinite(cpsth_all) & cpsth_all > 0 & ...
               isfinite(known_u_all)  & known_u_all  > 0 & ...
               isfinite(known_th_all) & known_th_all > 0 & ...
               isfinite(PV_used_bridge)  & PV_used_bridge  > 0;

if hasSm
    known_sm_all = double(outtbl{:,col_ksm});
    n_sm_missing = sum(is_valid_bridge & ...
        ~(isfinite(cpssm_all) & cpssm_all > 0 & isfinite(known_sm_all) & known_sm_all > 0));
    if n_sm_missing > 0
        fprintf(['  NOTE: %d bridge-standard analyses have valid U/Th references but no usable ', ...
                 'known_sm_ppm/cpssm — they will still anchor U/Th calibration; Sm for any unknown ', ...
                 'assigned to them falls back to the NIST612-anchored value instead.\n'], n_sm_missing);
    end
else
    known_sm_all = NaN(height(outtbl),1);
end

n_pre_filter = sum(is_valid_bridge);
fprintf('  Bridge-standard analyses eligible before Hampel filter: %d\n', n_pre_filter);
if n_pre_filter < 1
    % Break down exactly which criterion is failing, among the rows that
    % matched the bridge-standard name, so this doesn't require a second
    % round of guessing.
    n_matched   = sum(is_bridge);
    n_cpsu_ok   = sum(is_bridge & isfinite(cpsu_all)  & cpsu_all  > 0);
    n_cpsth_ok  = sum(is_bridge & isfinite(cpsth_all) & cpsth_all > 0);
    n_ku_ok     = sum(is_bridge & isfinite(known_u_all)  & known_u_all  > 0);
    n_kth_ok    = sum(is_bridge & isfinite(known_th_all) & known_th_all > 0);
    n_pv_ok     = sum(is_bridge & isfinite(PV_used_bridge)  & PV_used_bridge  > 0);
    fprintf(['  Diagnostic breakdown for %d rows matching bridge standard "%s":\n', ...
             '    cpsu > 0:            %d\n', ...
             '    cpsth > 0:           %d\n', ...
             '    known_u_ppm > 0:     %d\n', ...
             '    known_th_ppm > 0:    %d\n', ...
             '    bridge PV > 0:       %d\n'], ...
        n_matched, bridgeName, n_cpsu_ok, n_cpsth_ok, n_ku_ok, n_kth_ok, n_pv_ok);
    if n_pv_ok == 0
        fprintf(['  -> Likely cause: none of these rows'' grainids matched an entry in your ', ...
                 'UTh pit-volume CSV. Check grainid spelling/formatting there against the ', ...
                 'bridge-standard rows'' grainid in your metadata.\n']);
    elseif n_ku_ok == 0 || n_kth_ok == 0
        fprintf(['  -> Likely cause: known_u_ppm/known_th_ppm never got joined for these rows — ', ...
                 'check that ladd_enrich_metadata actually matched this stdname against your ', ...
                 'lookup table (see its "Rows with known_u_ppm: X/Y" console output).\n']);
    elseif n_cpsu_ok == 0 || n_cpsth_ok == 0
        fprintf('  -> Likely cause: no usable cps signal on these rows (check for NOWIN/NO_SEGMENT flags).\n');
    end
    error(['ladd_reduce_apatite_pvmeas_hybrid: no usable bridge-standard rows. ', ...
           'See diagnostic breakdown above.']);
end

% ---- Hampel-filter the anchor pool on PV-normalized implied ppm ----------
% K*_X = known_X_ppm / (cps_X / PV) for each candidate anchor. This is
% element-specific — deviations here reflect a mismeasured pit volume, a
% laser hiccup, or genuine standard heterogeneity, not real drift (that's
% NIST612's job). U and Th are required for every anchor (they were already
% required to reach this point). Sm is filtered only among anchors that
% actually HAVE a usable known_sm_ppm/cpssm — an anchor lacking an Sm
% reference is neither rewarded nor penalized on Sm, it simply won't be
% used for Sm bracketing downstream (see NOTE above).
bridge_idx_all = seq_all(is_valid_bridge);

implied_U  = known_u_all(bridge_idx_all)  ./ (cpsu_all(bridge_idx_all)  ./ PV_used_bridge(bridge_idx_all));
implied_Th = known_th_all(bridge_idx_all) ./ (cpsth_all(bridge_idx_all) ./ PV_used_bridge(bridge_idx_all));

[~, keepU_idx,  nRejU]  = hampel_filter(implied_U,  bridge_idx_all, kMAD_std, halfWinStd);
[~, keepTh_idx, nRejTh] = hampel_filter(implied_Th, bridge_idx_all, kMAD_std, halfWinStd);
keep_idx_hampel = intersect(keepU_idx, keepTh_idx);

if hasSm
    has_sm_ref = isfinite(cpssm_all(bridge_idx_all)) & cpssm_all(bridge_idx_all) > 0 & ...
                 isfinite(known_sm_all(bridge_idx_all)) & known_sm_all(bridge_idx_all) > 0;
    sm_idx_subset = bridge_idx_all(has_sm_ref);
    if numel(sm_idx_subset) >= 1
        implied_Sm = known_sm_all(sm_idx_subset) ./ (cpssm_all(sm_idx_subset) ./ PV_used_bridge(sm_idx_subset));
        [~, keepSm_idx, nRejSm] = hampel_filter(implied_Sm, sm_idx_subset, kMAD_std, halfWinStd);
        failedSm_idx = setdiff(sm_idx_subset, keepSm_idx);
        % Only drop anchors that HAD an Sm reference and failed the filter —
        % anchors with no Sm reference at all stay eligible for U/Th.
        keep_idx_hampel = setdiff(keep_idx_hampel, failedSm_idx);
    else
        nRejSm = 0;
    end
else
    nRejSm = 0;
end

fprintf('  Bridge-standard Hampel rejects — U: %d, Th: %d, Sm: %d\n', nRejU, nRejTh, nRejSm);

reviewed_idx = setdiff(bridge_idx_all, keep_idx_hampel);

bridgeIDsCell = cellstr(string(outtbl.grainid));
excludeBridgeCell = cellstr(string(excludeBridgeIDs));
excludeBridgeCell = excludeBridgeCell(~cellfun(@isempty, excludeBridgeCell));
bridgeNorm = cellfun(@cal_norm_id, bridgeIDsCell, 'UniformOutput', false);
excludeBridgeNorm = cellfun(@cal_norm_id, excludeBridgeCell, 'UniformOutput', false);
manual_excluded_idx = bridge_idx_all(ismember(bridgeNorm(bridge_idx_all), excludeBridgeNorm));

if autoExcludeBridgeReviews
    auto_excluded_idx = reviewed_idx;
    keep_idx = setdiff(keep_idx_hampel, manual_excluded_idx);
else
    auto_excluded_idx = [];
    keep_idx = setdiff(bridge_idx_all, manual_excluded_idx);
end

if numel(keep_idx) < 3
    keep_without_auto = setdiff(bridge_idx_all, manual_excluded_idx);
    if autoExcludeBridgeReviews && numel(keep_without_auto) >= 3
        warning(['ladd_reduce_apatite_pvmeas_hybrid: <3 bridge-standard anchors remain after ', ...
                 'automatic Hampel exclusion; restoring reviewed anchors while retaining ', ...
                 'explicit exclusions.']);
        keep_idx = keep_without_auto;
        auto_excluded_idx = [];
    else
        error(['ladd_reduce_apatite_pvmeas_hybrid: fewer than 3 bridge anchors remain after ', ...
               'explicit exclusions. Review excludeBridgeIDs; explicit exclusions are never ', ...
               'silently restored.']);
    end
end

rejected_idx = union(manual_excluded_idx, auto_excluded_idx);
fprintf('  Hampel review rows: %d; explicitly excluded: %d; automatically excluded: %d\n', ...
    numel(reviewed_idx), numel(manual_excluded_idx), numel(auto_excluded_idx));
if ~isempty(reviewed_idx)
    fprintf('  Bridge standards flagged for review: %s\n', ...
        strjoin(bridgeIDsCell(reviewed_idx), ', '));
end
if ~isempty(manual_excluded_idx)
    fprintf('  Bridge standards explicitly excluded: %s\n', ...
        strjoin(bridgeIDsCell(manual_excluded_idx), ', '));
end
if ~isempty(excludeBridgeNorm)
    unmatchedBridge = excludeBridgeNorm(~ismember(excludeBridgeNorm, bridgeNorm(bridge_idx_all)));
    if ~isempty(unmatchedBridge)
        warning('ladd_reduce_apatite_pvmeas_hybrid:excludeBridgeIDsNotFound', ...
            'excludeBridgeIDs not found in eligible bridge standards: %s', ...
            strjoin(unmatchedBridge, ', '));
    end
end

is_valid_bridge(:) = false;
is_valid_bridge(keep_idx) = true;

bridge_seqs = seq_all(is_valid_bridge);
fprintf('  Bridge-standard analyses used as anchors: %d\n', numel(bridge_seqs));

% ---- QA: NIST612-anchored ppm vs declared known ppm for the bridge std ---
% Purely diagnostic — does not affect bridge calibration. It is available
% only when actual NIST612 glass was included in the session.
if hasIndependentNist
    qa_u  = outtbl.u_ppm_nist612(bridge_seqs)  ./ known_u_all(bridge_seqs);
    qa_th = outtbl.th_ppm_nist612(bridge_seqs) ./ known_th_all(bridge_seqs);
    med_qa_u  = median(qa_u,  'omitnan');
    med_qa_th = median(qa_th, 'omitnan');
    med_nist_bridge_u  = median(outtbl.u_ppm_nist612(bridge_seqs), 'omitnan');
    med_nist_bridge_th = median(outtbl.th_ppm_nist612(bridge_seqs), 'omitnan');
    med_nist_bridge_sm = median(outtbl.sm_ppm_nist612(bridge_seqs), 'omitnan');
    qaRobStd = @(x) 1.4826*mad(x(isfinite(x)),1);
    rel_nist_bridge_u  = qaRobStd(outtbl.u_ppm_nist612(bridge_seqs))  / max(eps,abs(med_nist_bridge_u));
    rel_nist_bridge_th = qaRobStd(outtbl.th_ppm_nist612(bridge_seqs)) / max(eps,abs(med_nist_bridge_th));
    rel_nist_bridge_sm = qaRobStd(outtbl.sm_ppm_nist612(bridge_seqs)) / max(eps,abs(med_nist_bridge_sm));
    med_qa_sm = NaN;
    if strcmp(smReferenceBasis,'total') && any(isfinite(known_sm_all(bridge_seqs)) & known_sm_all(bridge_seqs)>0)
        qa_sm = outtbl.sm_ppm_nist612(bridge_seqs) ./ known_sm_all(bridge_seqs);
        med_qa_sm = median(qa_sm, 'omitnan');
    end
    fprintf('  QA — bridge standard NIST612-anchored ppm / known ppm: U=%.2f  Th=%.2f (1.00 = perfect agreement)\n', ...
        med_qa_u, med_qa_th);
    nistMismatch = abs(med_qa_u-1)*100 > nistTolPct || abs(med_qa_th-1)*100 > nistTolPct;
    if nistMismatch
        fprintf(['  QA WARNING: bridge standard disagrees with NIST612-anchored calibration by >%.0f%%. ', ...
                 'This does not affect the bridge result, but is worth investigating ', ...
                 '(matrix mismatch, NIST612 drift, or bridge-standard identity).\n'], nistTolPct);
    end
else
    med_qa_u = NaN; med_qa_th = NaN; med_qa_sm = NaN;
    med_nist_bridge_u = NaN; med_nist_bridge_th = NaN; med_nist_bridge_sm = NaN;
    rel_nist_bridge_u = NaN; rel_nist_bridge_th = NaN; rel_nist_bridge_sm = NaN;
    nistMismatch = false;
    fprintf('  QA — independent NIST612 comparison unavailable (bridge-only session)\n');
end

% ---- Pooled (session-wide) calibration factors, for anchorMode='median' --
% K_X = known_X_ppm / (cps_X/PV), pooled by median across every Hampel-
% included anchor rather than picking the single nearest one. This is the
% A robust multi-analysis factor can be more stable than a single-analysis
% bracket when the bridge standard has real replicate-to-replicate scatter.
% It is applied symmetrically to U, Th, and Sm.
robstd = @(x) 1.4826*mad(x,1);
K_U_i  = known_u_all(bridge_seqs)  ./ (cpsu_all(bridge_seqs)  ./ PV_used_bridge(bridge_seqs));
K_Th_i = known_th_all(bridge_seqs) ./ (cpsth_all(bridge_seqs) ./ PV_used_bridge(bridge_seqs));
K_U_med  = median(K_U_i,  'omitnan'); relK_U  = robstd(K_U_i)  / max(eps, abs(K_U_med));
K_Th_med = median(K_Th_i, 'omitnan'); relK_Th = robstd(K_Th_i) / max(eps, abs(K_Th_med));
relRef_U  = median_reference_rse(known_u_all(bridge_seqs), ku1sd_join(bridge_seqs));
relRef_Th = median_reference_rse(known_th_all(bridge_seqs), kth1sd_join(bridge_seqs));
K_Sm_med = NaN; relK_Sm = NaN;
relRef_Sm = 0;
if hasSm
    has_sm_ref_final = isfinite(cpssm_all(bridge_seqs)) & cpssm_all(bridge_seqs) > 0 & ...
                        isfinite(known_sm_all(bridge_seqs)) & known_sm_all(bridge_seqs) > 0;
    sm_anchor_idx = bridge_seqs(has_sm_ref_final);
    if ~isempty(sm_anchor_idx)
        K_Sm_i   = known_sm_all(sm_anchor_idx) ./ (cpssm_all(sm_anchor_idx) ./ PV_used_bridge(sm_anchor_idx));
        K_Sm_med = median(K_Sm_i, 'omitnan'); relK_Sm = robstd(K_Sm_i) / max(eps, abs(K_Sm_med));
        relRef_Sm = median_reference_rse(known_sm_all(sm_anchor_idx), ksm1sd_join(sm_anchor_idx));
    end
end
if strcmpi(anchorMode, 'median')
    fprintf('  anchorMode=median — pooled K (known_ppm / [cps/PV]) from %d anchors:\n', numel(bridge_seqs));
    fprintf('    U:  %.4g (relative scatter %.1f%%)\n', K_U_med,  relK_U*100);
    fprintf('    Th: %.4g (relative scatter %.1f%%)\n', K_Th_med, relK_Th*100);
    if hasSm && isfinite(K_Sm_med)
        fprintf('    Sm: %.4g (relative scatter %.1f%%, from %d anchors with a usable Sm reference)\n', ...
            K_Sm_med, relK_Sm*100, numel(sm_anchor_idx));
    else
        fprintf('    Sm: no anchors with a usable known_sm_ppm — Sm stays on the NIST612-anchored value.\n');
    end
end

% ---- Per-unknown nearest-following bracket assignment --------------------
% Computed regardless of anchorMode, so bridge_seq is always populated for
% record-keeping. It drives the actual ppm calculation for both the mineral-
% standard 'nearest' mode and the direct-NIST 'nist_following' mode.
assigned_bridge  = NaN(height(outtbl),1);
used_fallback = false(height(outtbl),1);
for i = 1:height(outtbl)
    if ~is_unknown(i), continue; end
    fwd = bridge_seqs(bridge_seqs > i);
    if ~isempty(fwd)
        assigned_bridge(i) = fwd(1);
    else
        [~,ni] = min(abs(bridge_seqs-i));
        assigned_bridge(i) = bridge_seqs(ni);
        used_fallback(i) = true;
    end
end
if ~strcmp(anchorMode, 'nist_interpolated') && sum(used_fallback & is_unknown) > 0
    fprintf('  %d unknowns used nearest-bridge-standard fallback\n', sum(used_fallback & is_unknown));
end
if strcmp(anchorMode, 'nist_interpolated')
    % No single discrete bracket is used in this mode. Keep these audit
    % fields empty rather than implying that the assigned row controlled the
    % interpolated calibration.
    assigned_bridge(:) = NaN;
    used_fallback(:) = false;
end

u_ppm_new    = double(outtbl{:,col_uppm});
th_ppm_new   = double(outtbl{:,col_thppm});
u_ppm_se_new = double(outtbl{:,col_uppmse});
th_ppm_se_new= double(outtbl{:,col_thppmse});

if hasSm
    sm_ppm_new = double(outtbl{:,col_smppm});
    if ~isempty(col_smppmse)
        sm_ppm_se_new = double(outtbl{:,col_smppmse});
    else
        sm_ppm_se_new = NaN(height(outtbl),1);
    end
end
smBridgeApplied = false(height(outtbl),1);
uReferenceRel1sdUsed = NaN(height(outtbl),1);
thReferenceRel1sdUsed = NaN(height(outtbl),1);
smReferenceRel1sdUsed = NaN(height(outtbl),1);

for i = 1:height(outtbl)
    if ~is_unknown(i) || isnan(assigned_bridge(i)), continue; end
    mi = assigned_bridge(i);

    pv_grain = PV_used(i);
    if ~(isfinite(pv_grain) && pv_grain > 0)
        continue;  % later flagged as NO_PITVOL/BAD_PV
    end

    if strcmp(anchorMode, 'nist_interpolated')
        % u_ppm_new/th_ppm_new/sm_ppm_new were initialized from
        % reduce_core_pv's NIST612 time-interpolated calibration. Retain
        % those values and their original analytical uncertainties.
        if hasSm && isfinite(sm_ppm_new(i))
            smBridgeApplied(i) = true;
        end
        isNist = strcmpi(string(outtbl.type),'NIST612');
        uReferenceRel1sdUsed(i) = median_reference_rse( ...
            known_u_all(isNist), ku1sd_join(isNist));
        thReferenceRel1sdUsed(i) = median_reference_rse( ...
            known_th_all(isNist), kth1sd_join(isNist));
        smReferenceRel1sdUsed(i) = median_reference_rse( ...
            known_sm_all(isNist), ksm1sd_join(isNist));
        continue;
    end

    pv_grain_rse = 0;
    if isfinite(PV_used_1sd(i)) && PV_used_1sd(i) > 0
        pv_grain_rse = PV_used_1sd(i) / pv_grain;
    end

    if strcmpi(anchorMode, 'median')
        % ---- Pooled-median calibration: session-wide K, not one anchor ---
        u_ppm_new(i) = (cpsu_all(i) / pv_grain) * K_U_med;
        u_rse = sqrt( (cpsu_se_all(i) / max(abs(cpsu_all(i)), eps))^2 + ...
            relK_U^2 + relRef_U^2 + pv_grain_rse^2 );
        u_ppm_se_new(i) = abs(u_ppm_new(i)) * u_rse;

        th_ppm_new(i) = (cpsth_all(i) / pv_grain) * K_Th_med;
        th_rse = sqrt( (cpsth_se_all(i) / max(abs(cpsth_all(i)), eps))^2 + ...
            relK_Th^2 + relRef_Th^2 + pv_grain_rse^2 );
        th_ppm_se_new(i) = abs(th_ppm_new(i)) * th_rse;
        uReferenceRel1sdUsed(i) = relRef_U;
        thReferenceRel1sdUsed(i) = relRef_Th;

        if hasSm && isfinite(K_Sm_med) && isfinite(cpssm_all(i)) && cpssm_all(i) > 0
            sm_ppm_new(i) = (cpssm_all(i) / pv_grain) * K_Sm_med;
            sm_rse = sqrt( (cpssm_se_all(i) / max(abs(cpssm_all(i)), eps))^2 + ...
                relK_Sm^2 + relRef_Sm^2 + pv_grain_rse^2 );
            sm_ppm_se_new(i) = abs(sm_ppm_new(i)) * sm_rse;
            smBridgeApplied(i) = true;
            smReferenceRel1sdUsed(i) = relRef_Sm;
        end
        continue;
    end

    % ---- anchorMode = 'nearest' (default): original per-analysis bracket -
    pv_mi = PV_used_bridge(mi);
    if ~(isfinite(pv_mi) && pv_mi > 0)
        error('Assigned bridge-standard row %d (%s) has no matched UTh pit volume.', ...
              mi, string(outtbl.grainid(mi)));
    end
    pv_mer_rse = 0;
    if isfinite(PV_UTh_1sd(mi)) && PV_UTh_1sd(mi) > 0
        pv_mer_rse = PV_UTh_1sd(mi) / pv_mi;
    end

    % U — identical structure for U/Th/Sm, per excel-match convention
    u_ppm_new(i) = (cpsu_all(i) / pv_grain) / (cpsu_all(mi) / pv_mi) * known_u_all(mi);
    uRefRse = scalar_reference_rse(known_u_all(mi), ku1sd_join(mi));
    u_rse = sqrt( ...
        (cpsu_se_all(i)  / max(abs(cpsu_all(i)),  eps))^2 + ...
        (cpsu_se_all(mi) / max(abs(cpsu_all(mi)), eps))^2 + ...
        pv_grain_rse^2 + pv_mer_rse^2 + uRefRse^2 );
    u_ppm_se_new(i) = abs(u_ppm_new(i)) * u_rse;
    uReferenceRel1sdUsed(i) = uRefRse;

    % Th
    th_ppm_new(i) = (cpsth_all(i) / pv_grain) / (cpsth_all(mi) / pv_mi) * known_th_all(mi);
    thRefRse = scalar_reference_rse(known_th_all(mi), kth1sd_join(mi));
    th_rse = sqrt( ...
        (cpsth_se_all(i)  / max(abs(cpsth_all(i)),  eps))^2 + ...
        (cpsth_se_all(mi) / max(abs(cpsth_all(mi)), eps))^2 + ...
        pv_grain_rse^2 + pv_mer_rse^2 + thRefRse^2 );
    th_ppm_se_new(i) = abs(th_ppm_new(i)) * th_rse;
    thReferenceRel1sdUsed(i) = thRefRse;

    % Sm — same formula as U and Th (no separate convention)
    if hasSm && isfinite(cpssm_all(i)) && cpssm_all(i) > 0 && ...
               isfinite(cpssm_all(mi)) && cpssm_all(mi) > 0 && ...
               isfinite(known_sm_all(mi)) && known_sm_all(mi) > 0
        sm_ppm_new(i) = (cpssm_all(i) / pv_grain) / (cpssm_all(mi) / pv_mi) * known_sm_all(mi);
        smRefRse = scalar_reference_rse(known_sm_all(mi), ksm1sd_join(mi));
        sm_rse = sqrt( ...
            (cpssm_se_all(i)  / max(abs(cpssm_all(i)),  eps))^2 + ...
            (cpssm_se_all(mi) / max(abs(cpssm_all(mi)), eps))^2 + ...
            pv_grain_rse^2 + pv_mer_rse^2 + smRefRse^2 );
        sm_ppm_se_new(i) = abs(sm_ppm_new(i)) * sm_rse;
        smBridgeApplied(i) = true;
        smReferenceRel1sdUsed(i) = smRefRse;
    end
end

outtbl{:,col_uppm}    = u_ppm_new;
outtbl{:,col_thppm}   = th_ppm_new;
outtbl{:,col_uppmse}  = u_ppm_se_new;
outtbl{:,col_thppmse} = th_ppm_se_new;
if hasSm
    outtbl{:,col_smppm} = sm_ppm_new;
    if ~isempty(col_smppmse)
        outtbl{:,col_smppmse} = sm_ppm_se_new;
    end
end
smCalSource = repmat("NOT_APPLICABLE", height(outtbl), 1);
if strcmp(anchorMode, 'nist_interpolated')
    hasNistSm = is_unknown & isfinite(outtbl.sm_ppm_nist612);
    smCalSource(hasNistSm) = "NIST612_INTERPOLATED";
    smCalSource(is_unknown & ~hasNistSm) = "NOT_CALCULATED";
else
    needsSmFallback = is_unknown & ~smBridgeApplied;
    hasNistSm = needsSmFallback & hasIndependentNist & ...
        isfinite(outtbl.sm_ppm_nist612);
    smCalSource(hasNistSm) = "NIST612_FALLBACK";
    smCalSource(needsSmFallback & ~hasNistSm) = "NOT_CALCULATED";
    if strcmp(anchorMode, 'median')
        smCalSource(smBridgeApplied) = "BRIDGE_MEDIAN";
    elseif strcmp(anchorMode, 'nist_following')
        smCalSource(smBridgeApplied) = "NIST612_FOLLOWING";
    else
        smCalSource(smBridgeApplied) = "BRIDGE_NEAREST";
    end
end
outtbl.sm_calibration_source = smCalSource;
outtbl.bridge_seq = assigned_bridge;
outtbl.bridge_anchor_fallback = used_fallback;
outtbl.u_reference_rel_1sd_used = uReferenceRel1sdUsed;
outtbl.th_reference_rel_1sd_used = thReferenceRel1sdUsed;
outtbl.sm_reference_rel_1sd_used = smReferenceRel1sdUsed;
referenceUncertaintyStatus = repmat("NOT_APPLICABLE",height(outtbl),1);
hasReferenceUncertainty = is_unknown & ...
    ((isfinite(uReferenceRel1sdUsed) & uReferenceRel1sdUsed > 0) | ...
     (isfinite(thReferenceRel1sdUsed) & thReferenceRel1sdUsed > 0) | ...
     (isfinite(smReferenceRel1sdUsed) & smReferenceRel1sdUsed > 0));
referenceUncertaintyStatus(is_unknown) = "NOT_SUPPLIED_ASSUMED_ZERO";
referenceUncertaintyStatus(hasReferenceUncertainty) = "PROPAGATED";
outtbl.reference_uncertainty_status = referenceUncertaintyStatus;

if strcmp(anchorMode, 'nist_interpolated')
    fprintf('  NIST612 time-interpolated parent calibration retained\n');
else
    fprintf('  Bridge-standard bracketing complete\n');
end

% ── Step 3: Store pit volume columns ─────────────────────────────────────
outtbl.pv_he_um3    = PV_He;    outtbl.pv_he_1sd    = PV_He_1sd;
outtbl.pv_uth_um3   = PV_UTh;   outtbl.pv_uth_1sd   = PV_UTh_1sd;
outtbl.pv_used_um3  = PV_used;  outtbl.pv_used_1sd  = PV_used_1sd;
if useAverageUthPit
    uthPitMode = "SESSION_AVERAGE";
    uthPitSourceText = sprintf('%.15g +/- %.15g um3 (1SD)', ...
        uthPitVolSource, uthAverage1sd);
else
    uthPitMode = "MEASURED_PER_ANALYSIS";
    uthPitSourceText = char(string(uthPitVolSource));
end
outtbl.uth_pit_volume_mode = repmat(uthPitMode, nR, 1);
outtbl.uth_pit_volume_source = repmat(string(uthPitSourceText), nR, 1);
% ── Step 4: Convert ppm -> atoms ─────────────────────────────────────────
fprintf('ladd_reduce_apatite_pvmeas_hybrid: converting ppm to atoms...\n');
cn_out = lower(outtbl.Properties.VariableNames);
assert(ismember('u_ppm', cn_out) && ismember('th_ppm', cn_out), ...
    'ladd_reduce_apatite_pvmeas_hybrid: u_ppm or th_ppm missing');

grain_mass = PV_used .* rho_apatite .* cm3_per_um3;

u_atoms  = (outtbl.u_ppm  / 1e6) / MW_U  .* NA .* grain_mass;
th_atoms = (outtbl.th_ppm / 1e6) / MW_Th .* NA .* grain_mass;
sm_atoms = zeros(nR,1); sm_atoms_se = zeros(nR,1);
if ismember('sm_ppm', cn_out)
    sm_atoms = (outtbl.sm_ppm / 1e6) / MW_Sm .* NA .* grain_mass;
    if ismember('sm_ppm_se', cn_out)
        sm_ppm_se_rel = outtbl.sm_ppm_se ./ max(abs(outtbl.sm_ppm), eps);
        sm_atoms_se   = abs(sm_atoms) .* sm_ppm_se_rel;
    end
end

ppm_rel_se_u  = outtbl.u_ppm_se  ./ max(abs(outtbl.u_ppm),  eps);
ppm_rel_se_th = outtbl.th_ppm_se ./ max(abs(outtbl.th_ppm), eps);

% The concentration uncertainties already include the applicable pit-volume
% terms. Do not add them again while constructing the atoms/g values used
% by the age calculation.
u_atoms_se  = abs(u_atoms)  .* ppm_rel_se_u;
th_atoms_se = abs(th_atoms) .* ppm_rel_se_th;

u_atoms_g    = u_atoms    ./ grain_mass;  th_atoms_g    = th_atoms    ./ grain_mass;
u_atoms_g_se = (outtbl.u_ppm_se  / 1e6) / MW_U  .* NA;
th_atoms_g_se= (outtbl.th_ppm_se / 1e6) / MW_Th .* NA;
sm_atoms_g   = sm_atoms ./ grain_mass;
sm_atoms_g_se= (outtbl.sm_ppm_se / 1e6) / MW_Sm .* NA;

mask_nan = ~is_unknown;
u_atoms(mask_nan)=NaN;      th_atoms(mask_nan)=NaN;     sm_atoms(mask_nan)=NaN;
u_atoms_se(mask_nan)=NaN;   th_atoms_se(mask_nan)=NaN;  sm_atoms_se(mask_nan)=NaN;
u_atoms_g(mask_nan)=NaN;    th_atoms_g(mask_nan)=NaN;   sm_atoms_g(mask_nan)=NaN;
u_atoms_g_se(mask_nan)=NaN; th_atoms_g_se(mask_nan)=NaN;sm_atoms_g_se(mask_nan)=NaN;

outtbl.u_atoms=u_atoms;       outtbl.th_atoms=th_atoms;      outtbl.sm_atoms=sm_atoms;
outtbl.u_atoms_se=u_atoms_se; outtbl.th_atoms_se=th_atoms_se;outtbl.sm_atoms_se=sm_atoms_se;
outtbl.u_atoms_g=u_atoms_g;   outtbl.th_atoms_g=th_atoms_g;  outtbl.sm_atoms_g=sm_atoms_g;
outtbl.u_atoms_g_se=u_atoms_g_se; outtbl.th_atoms_g_se=th_atoms_g_se;
outtbl.sm_atoms_g_se=sm_atoms_g_se;


% ── Step 4b: Re-evaluate Sm QC flags AFTER bridge correction ─────────────
% reduce_core_pv's SM_UNSTABLE/SM_WEAK/SM_OUTLIER gate was evaluated against
% the NIST612-anchored sm_ppm/sm_ppm_se, BEFORE Step 2 above overwrote those
% columns with the bridge-standard-corrected values. That NIST612-only
% estimate can be noisier than the mineral-bridge result because glass and
% apatite have different Sm response. Re-evaluate the review flags after the
% final bridge calibration so they describe the reported values.
fprintf('ladd_reduce_apatite_pvmeas_hybrid: re-evaluating Sm QC flags post-bridge...\n');

SmMaxPlaus_post = p.Results.SmMaxPlaus;
SmMinCps_post   = p.Results.SmMinCps;
SmRelSEmax_post = p.Results.SmRelSEmax;

if ~ismember('flags', outtbl.Properties.VariableNames)
    outtbl.flags = strings(height(outtbl),1);
end
nStaleUnstable = sum(contains(outtbl.flags, "SM_UNSTABLE"));

% strip the stale pre-bridge Sm flags before re-adding fresh ones
outtbl.flags = strtrim(erase(outtbl.flags, ["SM_UNSTABLE","SM_WEAK","SM_OUTLIER"]));

relSE_post    = outtbl.sm_ppm_se ./ max(abs(outtbl.sm_ppm), eps);
weakCps_post  = ~isfinite(outtbl.cpssm) | outtbl.cpssm < SmMinCps_post;
outlier_post  = isfinite(outtbl.sm_ppm) & outtbl.sm_ppm > SmMaxPlaus_post;
unstable_post = isfinite(relSE_post)    & relSE_post    > SmRelSEmax_post;
reviewMask_post = outlier_post | unstable_post | weakCps_post;

addFlag_post = strings(height(outtbl),1);
addFlag_post(outlier_post)  = strtrim(addFlag_post(outlier_post)  + " " + "SM_OUTLIER");
addFlag_post(unstable_post) = strtrim(addFlag_post(unstable_post) + " " + "SM_UNSTABLE");
addFlag_post(weakCps_post)  = strtrim(addFlag_post(weakCps_post)  + " " + "SM_WEAK");
outtbl.flags = strtrim(outtbl.flags + " " + addFlag_post);

% Review only: preserve every finite bridge-calibrated Sm value. These
% thresholds identify measurements worth inspecting; they are not evidence
% of a laser misfire and therefore do not silently remove Sm from an age.
outtbl.sm_review_only = reviewMask_post;
outtbl.sm_review_max_ppm = repmat(p.Results.SmMaxPlaus, height(outtbl), 1);
outtbl.sm_review_min_cps = repmat(p.Results.SmMinCps, height(outtbl), 1);
outtbl.sm_review_max_relse = repmat(p.Results.SmRelSEmax, height(outtbl), 1);

fprintf('  Sm QC post-bridge (review only; values preserved): %d/%d unstable, %d weak, %d outlier  (was %d/%d flagged SM_UNSTABLE pre-bridge)\n', ...
    sum(unstable_post), height(outtbl), sum(weakCps_post), sum(outlier_post), nStaleUnstable, height(outtbl));


% ── Step 5: Parent production columns ────────────────────────────────────
outtbl = add_parent_production(outtbl, smReferenceBasis);

% ── Step 6: Flags ────────────────────────────────────────────────────────
if ~ismember('flags', lower(outtbl.Properties.VariableNames))
    outtbl.flags = strings(nR,1);
end
fl = string(outtbl.flags);
fl(~isfinite(PV_used) & is_unknown & strlength(string(outtbl.grainid))>0) = ...
    strtrim(fl(~isfinite(PV_used) & is_unknown & strlength(string(outtbl.grainid))>0) + " NO_PITVOL");
fl(has_uth_only) = strtrim(fl(has_uth_only) + " NO_HEPV");
fl(bad_res)      = strtrim(fl(bad_res)      + " BAD_PV");
fl(has_neither)  = strtrim(fl(has_neither)  + " NO_UTHPV");
smFallback = is_unknown & outtbl.sm_calibration_source == "NIST612_FALLBACK";
fl(smFallback) = strtrim(fl(smFallback) + " SM_NIST612_FALLBACK");
smNotCalculated = is_unknown & outtbl.sm_calibration_source == "NOT_CALCULATED";
fl(smNotCalculated) = strtrim(fl(smNotCalculated) + " SM_NOT_CALCULATED");
if ~isempty(reviewed_idx)
    fl(reviewed_idx) = strtrim(fl(reviewed_idx) + " STD_HAMPEL_REVIEW");
end
if ~isempty(manual_excluded_idx)
    fl(manual_excluded_idx) = strtrim(fl(manual_excluded_idx) + " STD_CAL_EXCLUDED_MANUAL");
end
if ~isempty(auto_excluded_idx)
    fl(auto_excluded_idx) = strtrim(fl(auto_excluded_idx) + " STD_CAL_EXCLUDED_AUTO");
end
if nistMismatch
    fl(bridge_idx_all) = strtrim(fl(bridge_idx_all) + " BRIDGE_NIST_MISMATCH_REVIEW");
end
outtbl.flags = fl;

bridgeRole = repmat("NOT_BRIDGE_STANDARD", nR, 1);
bridgeRole(bridge_idx_all) = "ANCHOR_USED";
bridgeRole(intersect(reviewed_idx, keep_idx)) = "ANCHOR_REVIEW_INCLUDED";
bridgeRole(manual_excluded_idx) = "EXCLUDED_MANUAL";
bridgeRole(auto_excluded_idx) = "EXCLUDED_AUTO";
outtbl.bridge_anchor_role = bridgeRole;
outtbl.bridge_standard_used = repmat(string(bridgeName), nR, 1);
outtbl.bridge_anchor_mode = repmat(string(anchorMode), nR, 1);
if strcmp(anchorMode,'nist_interpolated')
    parentCalibrationMode = "NIST612_INTERPOLATED";
elseif strcmp(anchorMode,'nist_following')
    parentCalibrationMode = "NIST612_FOLLOWING";
elseif strcmp(anchorMode,'median')
    parentCalibrationMode = "BRIDGE_MEDIAN";
else
    parentCalibrationMode = "BRIDGE_NEAREST";
end
outtbl.parent_calibration_mode = repmat(parentCalibrationMode, nR, 1);
outtbl.primary_calibration_path = repmat(calibrationPath, nR, 1);
outtbl.independent_nist_check_available = repmat(hasIndependentNist, nR, 1);
outtbl.sm_reference_basis = repmat(string(smReferenceBasis), nR, 1);
outtbl.bridge_kMAD_std = repmat(kMAD_std, nR, 1);
outtbl.bridge_halfWinStd = repmat(halfWinStd, nR, 1);
outtbl.bridge_autoExcludeReviews = repmat(autoExcludeBridgeReviews, nR, 1);
outtbl.bridge_nAnchorsUsed = repmat(numel(bridge_seqs), nR, 1);
if hasIndependentNist
    nistReviewKUsed = kMAD;
    nistReviewHalfWinUsed = nistHalfWin;
else
    nistReviewKUsed = NaN;
    nistReviewHalfWinUsed = NaN;
end
outtbl.nist_review_kMAD_used = repmat(nistReviewKUsed, nR, 1);
outtbl.nist_review_halfWin_used = repmat(nistReviewHalfWinUsed, nR, 1);
outtbl.nist_autoExcludeReviews_used = repmat(autoExcludeNistReviews, nR, 1);
outtbl.bridge_nist_median_u_ppm = repmat(med_nist_bridge_u, nR, 1);
outtbl.bridge_nist_median_th_ppm = repmat(med_nist_bridge_th, nR, 1);
outtbl.bridge_nist_median_sm_ppm = repmat(med_nist_bridge_sm, nR, 1);
outtbl.bridge_nist_ratio_u = repmat(med_qa_u, nR, 1);
outtbl.bridge_nist_ratio_th = repmat(med_qa_th, nR, 1);
outtbl.bridge_nist_ratio_sm = repmat(med_qa_sm, nR, 1);
outtbl.bridge_nist_relscatter_u = repmat(rel_nist_bridge_u, nR, 1);
outtbl.bridge_nist_relscatter_th = repmat(rel_nist_bridge_th, nR, 1);
outtbl.bridge_nist_relscatter_sm = repmat(rel_nist_bridge_sm, nR, 1);
outtbl.bridge_nist_mismatch_review = repmat(nistMismatch, nR, 1);

n_cal = sum(isfinite(u_atoms) & is_unknown);
fprintf('\n── Reduction summary ────────────────────────────────────\n');
fprintf('  Total analyses:          %d\n', nR);
fprintf('  Unknown grains:          %d\n', n_unk);
fprintf('  Successfully calibrated: %d  (%.0f%%)\n', n_cal, 100*n_cal/max(n_unk,1));
fprintf('\n');

if ~isempty(saveAs)
    writetable(outtbl, saveAs);
    fprintf('Wrote: %s\n', saveAs);
end

end % ── END MAIN ──────────────────────────────────────────────────────────


function [coreMetadataCsv, tempMetadataCsv, nTrueNist, nBridgeCoreRows] = ...
    prepare_core_metadata(metadataCsv, bridgeName, allowBridgeOnly, usesDirectNist)
md = readtable(metadataCsv, 'VariableNamingRule','preserve', ...
    'TextType','string', 'Delimiter',',');
md.Properties.VariableNames = lower(md.Properties.VariableNames);
assert(all(ismember({'file','type'}, md.Properties.VariableNames)), ...
    'ladd_reduce_apatite_pvmeas_hybrid: metadata must include file and type.');
if ~ismember('stdname', md.Properties.VariableNames)
    md.stdname = strings(height(md),1);
end

typeText = strtrim(string(md.type));
stdText = strtrim(string(md.stdname));
isTypeNist612 = strcmpi(typeText, 'NIST612');
isBlankStd = ismissing(stdText) | strlength(stdText) == 0;
isTrueNist = isTypeNist612 & (strcmpi(stdText, 'NIST612') | isBlankStd);
isInconsistentNist = isTypeNist612 & ~isTrueNist;
if any(isInconsistentNist)
    badNames = unique(stdText(isInconsistentNist));
    error(['ladd_reduce_apatite_pvmeas_hybrid: rows with type=NIST612 must ', ...
        'identify actual NIST612 glass, but found stdname(s): %s.'], ...
        strjoin(badNames, ', '));
end
nTrueNist = sum(isTrueNist);
nBridgeCoreRows = 0;

if nTrueNist < 1
    if usesDirectNist
        error(['ladd_reduce_apatite_pvmeas_hybrid: direct NIST AnchorMode ', ...
            'requires genuine NIST612 rows. Use a bridge mode when NIST612 ', ...
            'was not run.']);
    end
    if ~allowBridgeOnly
        error(['ladd_reduce_apatite_pvmeas_hybrid: no genuine NIST612 rows ', ...
            'were found. Set AllowBridgeOnly=true to use the explicitly ', ...
            'named mineral reference material as the sole calibration.']);
    end
    isBridgeCore = strcmpi(stdText, strtrim(string(bridgeName)));
    nBridgeCoreRows = sum(isBridgeCore);
    assert(nBridgeCoreRows >= 3, ...
        ['ladd_reduce_apatite_pvmeas_hybrid: bridge-only mode requires at ', ...
         'least three rows whose stdname exactly matches "%s"; found %d.'], ...
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
    cn = lower(regexprep(t.Properties.VariableNames, '[^a-zA-Z0-9]',''));
    gid_col = find_col(cn, {'grainid','id','name'});
    vol_col = find_col(cn, {'pitvol_um3','pitvolumeum3','pitvol','volume'});
    sd_col  = find_col(cn, {'pv1sd_um3','pv1sd','propsd','sdum3','sd_um3'});
    pv.ids  = lower(strtrim(cellstr(t{:,gid_col})));
    pv.vols = double(t{:,vol_col});
    pv.sds  = double(t{:,sd_col});
    fprintf('  Read %s pit vol CSV: %d rows from %s\n', label, height(t), csvPath);
end


function rse = scalar_reference_rse(knownValue, known1sd)
if isfinite(knownValue) && knownValue > 0 && ...
        isfinite(known1sd) && known1sd >= 0
    rse = known1sd / knownValue;
else
    rse = 0;
end
end


function rse = median_reference_rse(knownValues, known1sdValues)
valid = isfinite(knownValues) & knownValues > 0 & ...
    isfinite(known1sdValues) & known1sdValues >= 0;
if any(valid)
    rse = median(known1sdValues(valid) ./ knownValues(valid), 'omitnan');
else
    rse = 0;
end
end


function outtbl = add_parent_production(outtbl, smReferenceBasis)
f238=0.992742; f235=0.007204; f147Sm=0.1499;
lam238=1.55125e-10; lam235=9.8485e-10; lam232=4.9475e-11; lam147=6.54e-12;
U_ag=grab(outtbl,'u_atoms_g');   U_ag_se=grab(outtbl,'u_atoms_g_se');
Th_ag=grab(outtbl,'th_atoms_g'); Th_ag_se=grab(outtbl,'th_atoms_g_se');
Sm_ag=grab(outtbl,'sm_atoms_g'); Sm_ag_se=grab(outtbl,'sm_atoms_g_se');
N238=f238.*U_ag; N238_se=f238.*U_ag_se;
N235=f235.*U_ag; N235_se=f235.*U_ag_se;
N232=Th_ag;       N232_se=Th_ag_se;

if strcmp(smReferenceBasis, 'total')
    N147=f147Sm.*Sm_ag; N147_se=f147Sm.*Sm_ag_se;
else
    N147=Sm_ag;          N147_se=Sm_ag_se;
end

R    = 8*lam238.*N238 + 7*lam235.*N235 + 6*lam232.*N232 + lam147.*N147;
% 238U and 235U share one total-U measurement and one uncertainty term.
uProductionCoefficient = 8*lam238*f238 + 7*lam235*f235;
R_se = sqrt((uProductionCoefficient.*U_ag_se).^2 + ...
            (6*lam232.*N232_se).^2+(lam147.*N147_se).^2);
outtbl.n238_atoms_g=N238; outtbl.n238_atoms_g_se=N238_se;
outtbl.n235_atoms_g=N235; outtbl.n235_atoms_g_se=N235_se;
outtbl.n232_atoms_g=N232; outtbl.n232_atoms_g_se=N232_se;
outtbl.n147_atoms_g=N147; outtbl.n147_atoms_g_se=N147_se;
outtbl.parentprod_atoms_g_yr=R; outtbl.parentprod_atoms_g_yr_se=R_se;
end


function v = grab(T, name)
    cn = lower(T.Properties.VariableNames);
    idx = find(strcmp(cn, lower(name)), 1);
    if ~isempty(idx), v=double(T{:,idx}); v(~isfinite(v))=0;
    else, v=zeros(height(T),1); end
end


function idx = find_col(cn, aliases)
    idx = [];
    for k = 1:numel(aliases)
        key = lower(regexprep(char(aliases{k}), '[^a-zA-Z0-9]',''));
        hit = find(strcmp(cn, key), 1);
        if ~isempty(hit), idx=hit; return; end
    end
    error('ladd_reduce_apatite_pvmeas_hybrid: cannot find column matching [%s]', strjoin(aliases,'/'));
end


function id_norm = cal_norm_id(id)
    id_norm = lower(strtrim(char(string(id))));
    id_norm = regexprep(id_norm, '[^a-z0-9]+', '');
end


function [K_clean, seq_clean, n_rej] = hampel_filter(K, seq, kMAD, half_win)
if nargin<3, kMAD=3; end; if nargin<4, half_win=5; end
n=numel(K); keep=true(n,1); robstd=@(x)1.4826*mad(x,1);
for ii=1:n
    lo=max(1,ii-half_win); hi=min(n,ii+half_win);
    nb=K([lo:ii-1,ii+1:hi]); if numel(nb)<2, continue; end
    m=median(nb,'omitnan'); s=robstd(nb);
    if s>0 && abs(K(ii)-m)>kMAD*s, keep(ii)=false; end
end
n_rej=sum(~keep);
if sum(keep)<3, warning('hampel_filter: <3 anchors; keeping all.'); keep(:)=true; n_rej=0; end
K_clean=K(keep); seq_clean=seq(keep);
end
