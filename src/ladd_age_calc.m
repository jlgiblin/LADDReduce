function ages = ladd_age_calc(heFile, uthFile, varargin)
% LADD_AGE_CALC  Compute (U-Th)/He ages from He_reduce and UTh_reduce outputs.
%
% Joins He and U-Th-Sm data by GrainID, then solves the (U-Th)/He age
% equation iteratively for each grain. No alpha-ejection correction (FT)
% is applied — this is correct for LADD sub-grain pit analyses.
%
% USAGE
%   ages = ladd_age_calc('helium_reduced.csv', 'parents_reduced.csv');
%   ages = ladd_age_calc(heFile, uthFile, 'saveAs','ladd_ages.csv');
%   ages = ladd_age_calc(heFile, uthFile, 'mineral','zircon');
%   ages = ladd_age_calc(heFile, uthFile, ...
%       'excludeGrainIDs',"SampleA-grain-012",'saveAs','ladd_ages.csv');
%
% INPUTS
%   heFile  : output CSV from He_reduce  (must have GrainID, He4Unk_atoms_g,
%             He4Unk1SD_atoms_g, RunScript columns)
%   uthFile : output CSV from reduce_core/ladd_reduce_* (must have GrainID
%             or file column, U_atoms_g, U_atoms_g_se, Th_atoms_g,
%             Th_atoms_g_se, and Sm_atoms_g / Sm_atoms_g_se for apatite)
%
% OPTIONAL NAME-VALUE PAIRS
%   'saveAs'    — output CSV filename (default '' = no save)
%   'mineral'   — 'zircon' | 'apatite' (default: auto-detect from Sm column)
%   'maxIter'   — max iterations for solver (default 10)
%   'tolPct'    — convergence tolerance in % He difference (default 0.001)
%   'verbose'   — print matching and age summary (default true)
%
% AGE EQUATION (Reiners et al.)
%   4He = 8*N238*(exp(lam238*t)-1) + 7*N235*(exp(lam235*t)-1)
%        + 6*N232*(exp(lam232*t)-1) [+ N147*(exp(lam147*t)-1) for apatite]
%
%   Where N238 = f238 * U_atoms_g  (f238 = 0.992742)
%         N235 = f235 * U_atoms_g  (f235 = 0.007204)
%         N232 = Th_atoms_g
%         N147 = explicit n147_atoms_g from the reducer (apatite only)
%
% SOLVER: iterative Newton-like approach (matches spreadsheet method)
%   1. Compute linear age as first guess: t0 = He / P  where P = parent prod rate
%   2. Compute expected He at t0, compare to measured
%   3. Adjust t by the fractional He difference
%   4. Repeat until |delta_He%| < tolPct or maxIter reached
%
% ERROR PROPAGATION
%   Age uncertainty from quadrature propagation of He, U, Th (and Sm)
%   uncertainties through the age equation, evaluated at the converged age.
%
% OUTPUT COLUMNS
%   GrainID         — grain identifier
%   SampleName      — from He file
%   RunScript        — from He file (3=unknown, 4=He mineral std)
%   He4_atoms_g     — measured He atoms/g
%   He4_1SD_atoms_g — 1SD on He
%   U_atoms_g       — U atoms/g
%   U_atoms_g_se    — 1SD on U
%   Th_atoms_g      — Th atoms/g
%   Th_atoms_g_se   — 1SD on Th
%   Sm_atoms_g      — Sm atoms/g (apatite only, NaN for zircon)
%   Sm_atoms_g_se   — 1SD on Sm (apatite only)
%   ThU             — Th/U molar ratio
%   Age_Ma          — iterative (U-Th)/He age (Ma)
%   Age_1SD_Ma      — 1SD on age (Ma)
%   Age_2SD_Ma      — 2SD on age (Ma)
%   Age_1SDpct      — %1SD on age
%   LinearAge_Ma    — linear approximation age (Ma) — fast check
%   nIter           — iterations to convergence
%   converged       — true if solver converged within maxIter
%   ManualExclude   — true only for an explicitly named exclusion; the row
%                     is retained and its reported age fields are NaN
%   UThMatchCount   — number of U-Th rows sharing the normalized GrainID
%   UThSelectedRow  — selected row number in the filtered U-Th table
%   UThSelectedFile — source filename for the selected U-Th row, when present
%   UThMatchDecision— UNIQUE_MATCH, DUPLICATE_AUTO_RESOLVED,
%                     DUPLICATE_REVIEW_REQUIRED, DUPLICATE_NO_VALID_MATCH,
%                     or NO_MATCH
%   flags           — inherited flags from He and U-Th files
%
% IMPLEMENTATION NOTE: he_flags / uth_flags are read as MATLAB `string` arrays.
% A truly-blank cell read via readtable(...,'TextType','string') becomes
% the special value <missing>, NOT "". String concatenation propagates
% <missing> ("abc" + <missing> == <missing>), so any row whose UTh (or He)
% flags cell was blank would silently poison combined_flags for that row —
% and every flag appended after that point (NO_CONVERGE, NEG_AGE,
% NO_UTH_MATCH, LOW_PARENT) was lost too, since they're built by appending
% onto the already-<missing> string. Fixed by coercing <missing> to "" at
% the point each flags string is read/combined.

% ── Constants ─────────────────────────────────────────────────────────────
lam238 = 1.55125e-10;   % yr^-1
lam235 = 9.8485e-10;    % yr^-1
lam232 = 4.9475e-11;    % yr^-1
lam147 = 6.54e-12;      % yr^-1
f238   = 0.992742;      % 238U natural abundance
f235   = 0.007204;      % 235U natural abundance
f147   = 0.1499;        % 147Sm natural abundance
yr2Ma  = 1e-6;          % years to Ma

% ── Parse arguments ────────────────────────────────────────────────────────
p = inputParser;
p.addParameter('saveAs',  '',        @ischar);
p.addParameter('mineral', 'auto',    @(s)ischar(s)||isstring(s));
p.addParameter('maxIter', 10,        @(x)isnumeric(x)&&isscalar(x));
p.addParameter('tolPct',  0.001,     @(x)isnumeric(x)&&isscalar(x));
p.addParameter('verbose', true,      @islogical);
p.addParameter('excludeGrainIDs', strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.parse(varargin{:});
saveAs  = p.Results.saveAs;
mineral = lower(char(p.Results.mineral));
maxIter = p.Results.maxIter;
tolPct  = p.Results.tolPct;
verbose = p.Results.verbose;
excludeGrainIDs = p.Results.excludeGrainIDs;

% ── Read He file ──────────────────────────────────────────────────────────
He = readtable(heFile, 'VariableNamingRule','preserve', 'TextType','string');
He.Properties.VariableNames = lower_unique_names(He.Properties.VariableNames);
fprintf('ladd_age_calc: He file: %d rows from %s\n', height(He), heFile);

% Filter to unknowns only (RunScript=3); mineral stds (4) kept for QC
if ismember('runscript', He.Properties.VariableNames)
    rs = double(He.runscript);
    He = He(rs==3 | rs==4, :);
    fprintf('  After filtering unknowns+mineralStds: %d rows\n', height(He));
end

% ── Read U-Th file ────────────────────────────────────────────────────────
UTh = readtable(uthFile, 'VariableNamingRule','preserve', 'TextType','string');
UTh.Properties.VariableNames = lower_unique_names(UTh.Properties.VariableNames);
fprintf('  UTh file: %d rows from %s\n', height(UTh), uthFile);

% Filter UTh to unknowns only
if ismember('type', UTh.Properties.VariableNames)
    UTh = UTh(strcmpi(UTh.type,'Unknown') | strcmpi(UTh.type,'Sample'), :);
    fprintf('  After filtering UTh unknowns: %d rows\n', height(UTh));
end

% ── Detect mineral type ───────────────────────────────────────────────────
if strcmp(mineral,'auto')
    % First try: presence of a populated sm_atoms / sm_atoms_g column.
    % Column names are already lowercased above.
    hasSm = (ismember('sm_atoms', UTh.Properties.VariableNames) && ...
             any(isfinite(str2double(string(UTh.sm_atoms))))) || ...
            (ismember('sm_atoms_g', UTh.Properties.VariableNames) && ...
             any(isfinite(str2double(string(UTh.sm_atoms_g)))));

    if hasSm
        mineral = 'apatite';
    else
        % Sm column absent or all-NaN (e.g. all grains below SmMinCps).
        % Fall back to GrainID suffix: _ap## = apatite, _z## = zircon.
        all_ids = [cellstr(UTh.grainid); cellstr(He.grainid)];
        n_ap = sum(~cellfun(@isempty, regexp(all_ids, '_ap\d', 'match')));
        n_z  = sum(~cellfun(@isempty, regexp(all_ids, '_z\d',  'match')));
        if n_ap > n_z
            mineral = 'apatite';
        elseif n_z > n_ap
            mineral = 'zircon';
        else
            % Can't determine from GrainIDs either — check for sm_atoms
            % column existence regardless of values (column present = apatite)
            if ismember('sm_atoms', UTh.Properties.VariableNames) || ...
               ismember('sm_atoms_g', UTh.Properties.VariableNames)
                mineral = 'apatite';
            else
                mineral = 'zircon';
            end
        end
    end
    fprintf('  Mineral auto-detected: %s\n', mineral);
end
includeSm = strcmp(mineral,'apatite');

% ── Resolve GrainID in UTh file ──────────────────────────────────────────
% reduce_core output uses 'file' column not 'grainid'
% Try grainid first, then fall back to file-based matching
if ~ismember('grainid', UTh.Properties.VariableNames)
    if ismember('file', UTh.Properties.VariableNames)
        % Extract GrainID from filename: strip path, extension, and run prefix
        UTh.grainid = UTh.file;
        fprintf('  Note: UTh file has no grainid column — using file column for matching\n');
        fprintf('  Add a GrainID column to your metadata CSV for reliable matching\n');
    else
        error('ladd_age_calc: UTh file has no grainid or file column for matching.');
    end
end

% ── Match grains ──────────────────────────────────────────────────────────
he_ids  = cellfun(@norm_id, cellstr(He.grainid),  'UniformOutput', false);
uth_ids = cellfun(@norm_id, cellstr(UTh.grainid), 'UniformOutput', false);

He.match_idx       = zeros(height(He), 1);
uthMatchCount      = zeros(height(He), 1);
uthSelectedRow     = NaN(height(He), 1);
uthSelectedFile    = strings(height(He), 1);
uthMatchDecision   = repmat("NO_MATCH", height(He), 1);
hasParentProd      = any(strcmpi(UTh.Properties.VariableNames, 'parentprod_atoms_g_yr'));

for i = 1:height(He)
    idxAll = find(strcmp(he_ids{i}, uth_ids));
    uthMatchCount(i) = numel(idxAll);

    if isscalar(idxAll)
        selected = idxAll(1);
        uthMatchDecision(i) = "UNIQUE_MATCH";
    elseif numel(idxAll) > 1
        % Never silently take the first duplicate. Auto-resolve only when
        % exactly one candidate has physically usable parent measurements.
        physicallyUsable = false(numel(idxAll), 1);
        for k = 1:numel(idxAll)
            row = idxAll(k);
            uNow  = safe_num(UTh, row, {'u_atoms_g'});
            thNow = safe_num(UTh, row, {'th_atoms_g'});
            physicallyUsable(k) = isfinite(uNow) && uNow > 0 && ...
                                  isfinite(thNow) && thNow >= 0;
            if hasParentProd
                pNow = safe_num(UTh, row, {'parentprod_atoms_g_yr'});
                physicallyUsable(k) = physicallyUsable(k) && ...
                                      isfinite(pNow) && pNow > 0;
            end
        end

        usableIdx = idxAll(physicallyUsable);
        if isscalar(usableIdx)
            selected = usableIdx(1);
            uthMatchDecision(i) = "DUPLICATE_AUTO_RESOLVED";
        elseif isempty(usableIdx)
            selected = [];
            uthMatchDecision(i) = "DUPLICATE_NO_VALID_MATCH";
        else
            selected = [];
            uthMatchDecision(i) = "DUPLICATE_REVIEW_REQUIRED";
        end
    else
        selected = [];
    end

    if ~isempty(selected)
        He.match_idx(i) = selected;
        uthSelectedRow(i) = selected;
        if ismember('file', UTh.Properties.VariableNames)
            selectedFile = string(UTh.file(selected));
            if ~ismissing(selectedFile), uthSelectedFile(i) = selectedFile; end
        end
    end
end

nMatched = sum(He.match_idx > 0);
nMissing = sum(uthMatchCount == 0);
nDuplicateResolved = sum(uthMatchDecision == "DUPLICATE_AUTO_RESOLVED");
nDuplicateReview = sum(startsWith(uthMatchDecision, "DUPLICATE_") & ...
                       uthMatchDecision ~= "DUPLICATE_AUTO_RESOLVED");
fprintf('  GrainID matching: %d matched, %d unmatched, %d duplicate auto-resolved, %d duplicate unresolved\n', ...
    nMatched, nMissing, nDuplicateResolved, nDuplicateReview);
if nDuplicateResolved > 0
    resolvedRows = find(uthMatchDecision == "DUPLICATE_AUTO_RESOLVED");
    for k = 1:numel(resolvedRows)
        i = resolvedRows(k);
        fprintf('  DUPLICATE_AUTO_RESOLVED: %s -> UTh row %d (%s)\n', ...
            string(He.grainid(i)), uthSelectedRow(i), uthSelectedFile(i));
    end
end
if nDuplicateReview > 0
    reviewRows = find(startsWith(uthMatchDecision, "DUPLICATE_") & ...
                      uthMatchDecision ~= "DUPLICATE_AUTO_RESOLVED");
    for k = 1:numel(reviewRows)
        i = reviewRows(k);
        fprintf('  %s: %s (%d candidate rows; age withheld)\n', ...
            uthMatchDecision(i), string(He.grainid(i)), uthMatchCount(i));
    end
end
if nMissing > 0 && verbose
    miss_he = He.grainid(uthMatchCount == 0);
    fprintf('  He grains with no UTh match:\n');
    for k = 1:min(10, numel(miss_he))
        fprintf('    "%s"\n', miss_he{k});
    end
    if numel(miss_he) > 10
        fprintf('    ...and %d more\n', numel(miss_he)-10);
    end
    % Also show a sample of UTh GrainIDs so mismatches are easy to diagnose
    fprintf('  UTh file GrainID sample (first 5):\n');
    for k = 1:min(5, numel(uth_ids))
        fprintf('    "%s"\n', UTh.grainid{k});
    end
    fprintf('  He  file GrainID sample (first 5):\n');
    for k = 1:min(5, numel(he_ids))
        fprintf('    "%s"\n', He.grainid{k});
    end
end

% ── Extract numeric values ────────────────────────────────────────────────
get_col = @(T, names) get_numeric(T, names);

He4     = get_col(He,  {'he4unk_atoms_g','he4_atoms_g'});
He4_se  = get_col(He,  {'he4unk1sd_atoms_g','he4_1sd_atoms_g','he4unk1sd_atoms_g'});

% Build joined arrays
n = height(He);
U_ag    = NaN(n,1);  U_se   = NaN(n,1);
Th_ag   = NaN(n,1);  Th_se  = NaN(n,1);
Sm_ag   = NaN(n,1);  Sm_se  = NaN(n,1);
N147_ag = NaN(n,1);  N147_se_ag = NaN(n,1);
uth_flags = strings(n,1);

for i = 1:n
    j = He.match_idx(i);
    if j == 0, continue; end
    U_ag(i)  = safe_num(UTh, j, {'u_atoms_g'});
    U_se(i)  = safe_num(UTh, j, {'u_atoms_g_se'});
    Th_ag(i) = safe_num(UTh, j, {'th_atoms_g'});
    Th_se(i) = safe_num(UTh, j, {'th_atoms_g_se'});
    if includeSm
        % Ages require atoms/g, never absolute atoms. Prefer the reducer's
        % explicit 147Sm columns so the isotope fraction is not applied
        % twice. The fallback supports older total-Sm output files.
        Sm_ag(i) = safe_num(UTh, j, {'sm_atoms_g'});
        Sm_se(i) = safe_num(UTh, j, {'sm_atoms_g_se'});
        N147_ag(i) = safe_num(UTh, j, {'n147_atoms_g'});
        N147_se_ag(i) = safe_num(UTh, j, {'n147_atoms_g_se'});
        if ~isfinite(N147_ag(i)) && isfinite(Sm_ag(i))
            N147_ag(i) = f147 .* Sm_ag(i);
        end
        if ~isfinite(N147_se_ag(i)) && isfinite(Sm_se(i))
            N147_se_ag(i) = f147 .* Sm_se(i);
        end
    end
    if ismember('flags', UTh.Properties.VariableNames)
        % Guard against <missing> string (blank cell) poisoning downstream
        % concatenation -- see the implementation note above.
        v = string(UTh.flags(j));
        if ismissing(v), v = ""; end
        uth_flags(i) = v;
    end
    if uthMatchDecision(i) == "DUPLICATE_AUTO_RESOLVED"
        uth_flags(i) = strtrim(uth_flags(i) + " DUPLICATE_AUTO_RESOLVED");
    end
end

unresolvedDuplicate = startsWith(uthMatchDecision, "DUPLICATE_") & ...
                      uthMatchDecision ~= "DUPLICATE_AUTO_RESOLVED";
uth_flags(unresolvedDuplicate) = strtrim(uth_flags(unresolvedDuplicate) + " " + ...
                                        uthMatchDecision(unresolvedDuplicate));

% ── Compute parent atom counts ────────────────────────────────────────────
N238 = f238 .* U_ag;
N235 = f235 .* U_ag;
N232 = Th_ag;
N147 = N147_ag;           N147_se = N147_se_ag;
N147(~includeSm | ~isfinite(N147)) = 0;
N147_se(~includeSm | ~isfinite(N147_se)) = 0;

% Parent production rate (atoms/g/yr) at t=0 (linear approximation denom)
% P = 8*lam238*N238 + 7*lam235*N235 + 6*lam232*N232 + lam147*N147
P = 8*lam238.*N238 + 7*lam235.*N235 + 6*lam232.*N232 + lam147.*N147;

% ── Linear age (first approximation) ─────────────────────────────────────
% t_linear = He4 / P  (valid for small ages; good starting guess for all)
t_linear = He4 ./ max(P, eps);   % in years
LinearAge_Ma = t_linear .* yr2Ma;

% ── Iterative solver ──────────────────────────────────────────────────────
Age_yr    = t_linear;             % start from linear age
nIter_out = zeros(n,1);
converged = false(n,1);

for i = 1:n
    % Skip grains with no He, no parent atoms (unmatched UTh), or zero production
    if ~isfinite(He4(i)) || ~isfinite(U_ag(i)) || ~isfinite(Th_ag(i))
        Age_yr(i) = NaN;
        continue;
    end
    if He4(i) <= 0 || U_ag(i) <= 0 || P(i) <= 0
        Age_yr(i) = NaN;
        continue;
    end
    % Skip grains where parent production is implausibly low
    % (linear age > 5000 Ma — flags near-zero U+Th measurements)
    if He4(i) / max(P(i), eps) > 5000 * 1e6   % 5000 Ma in years
        Age_yr(i) = NaN;
        continue;   % LOW_PARENT flag applied after combined_flags is built
    end

    t = t_linear(i);
    for iter = 1:maxIter
        % Expected He at current age estimate
        He_exp = 8*N238(i)*(exp(lam238*t)-1) + ...
                 7*N235(i)*(exp(lam235*t)-1) + ...
                 6*N232(i)*(exp(lam232*t)-1) + ...
                 N147(i) *(exp(lam147*t)-1);

        delta_He     = He_exp - He4(i);
        delta_He_pct = delta_He / He4(i) * 100;

        % Age adjustment: dHe/dt at current t
        dHe_dt = 8*lam238*N238(i)*exp(lam238*t) + ...
                 7*lam235*N235(i)*exp(lam235*t) + ...
                 6*lam232*N232(i)*exp(lam232*t) + ...
                 lam147 *N147(i) *exp(lam147*t);

        t = t - delta_He / max(dHe_dt, eps);
        t = max(t, 0);   % age can't be negative

        nIter_out(i) = iter;
        if abs(delta_He_pct) < tolPct
            converged(i) = true;
            break;
        end
    end
    Age_yr(i) = t;
end

Age_Ma = Age_yr .* yr2Ma;

% ── Uncertainty propagation ───────────────────────────────────────────────
% dAge/dHe  = 1 / (dHe_prod/dt) evaluated at converged age
% dAge/dN238 = -8*(exp(lam238*t)-1) / (dHe_prod/dt)
% etc.

Age_1SD_Ma = NaN(n,1);

for i = 1:n
    if ~converged(i) && ~isfinite(Age_yr(i)), continue; end
    t = Age_yr(i);

    dHe_dt = 8*lam238*N238(i)*exp(lam238*t) + ...
             7*lam235*N235(i)*exp(lam235*t) + ...
             6*lam232*N232(i)*exp(lam232*t) + ...
             lam147 *N147(i) *exp(lam147*t);

    if dHe_dt <= 0, continue; end

    % Partial derivatives of age w.r.t. each measured quantity
    dT_dHe   =  1 / dHe_dt;
    dT_dN238 = -(8*(exp(lam238*t)-1) + 7*(exp(lam235*t)-1)) / dHe_dt;  % U feeds both
    dT_dN232 = -(6*(exp(lam232*t)-1)) / dHe_dt;
    dT_dN147 = -(exp(lam147*t)-1) / dHe_dt;

    % Combined U partial derivative (238U and 235U both scale with U_atoms_g)
    dT_dU = -(8*f238*(exp(lam238*t)-1) + 7*f235*(exp(lam235*t)-1)) / dHe_dt;
    dT_dTh= dT_dN232;
    var_age = (dT_dHe  * He4_se(i)).^2 + ...
              (dT_dU   * U_se(i)  ).^2 + ...
              (dT_dTh  * Th_se(i) ).^2;

    if includeSm && isfinite(N147_se(i)) && N147_se(i) > 0
        var_age = var_age + (dT_dN147 * N147_se(i)).^2;
    end

    Age_1SD_Ma(i) = sqrt(var_age) * yr2Ma;
end

Age_2SD_Ma  = 2 .* Age_1SD_Ma;
Age_1SDpct  = Age_1SD_Ma ./ max(Age_Ma, eps) * 100;
ThU         = Th_ag ./ max(U_ag, eps) .* (232/238);  % molar Th/U

% ── Combine flags ─────────────────────────────────────────────────────────
he_flags = strings(n,1);
if ismember('flags', He.Properties.VariableNames)
    % Guard against <missing> string (blank cell) poisoning downstream
    % concatenation -- see the implementation note above.
    he_flags = string(He.flags);
    he_flags(ismissing(he_flags)) = "";
end
combined_flags = strtrim(he_flags + " " + uth_flags);
% Belt-and-suspenders: if anything still slipped through as <missing>
% (e.g. a future column added without the same guard), don't let it
% silently swallow the NO_CONVERGE/NEG_AGE/NO_UTH_MATCH/LOW_PARENT flags
% appended below.
combined_flags(ismissing(combined_flags)) = "";
combined_flags(~converged & isfinite(Age_Ma)) = ...
    strtrim(combined_flags(~converged & isfinite(Age_Ma)) + " NO_CONVERGE");
combined_flags(Age_Ma < 0) = strtrim(combined_flags(Age_Ma < 0) + " NEG_AGE");
% Flag grains with no UTh match — age is meaningless for these
noUTh = ~isfinite(U_ag) | ~isfinite(Th_ag);
combined_flags(noUTh) = strtrim(combined_flags(noUTh) + " NO_UTH_MATCH");
Age_Ma(noUTh)       = NaN;
Age_1SD_Ma(noUTh)   = NaN;
Age_2SD_Ma(noUTh)   = NaN;
Age_1SDpct(noUTh)   = NaN;
LinearAge_Ma(noUTh) = NaN;

% Flag and NaN out LOW_PARENT grains (near-zero U+Th, linear age > 5000 Ma)
lowP = isfinite(He4) & isfinite(P) & P > 0 & (He4 ./ P) > 5000e6;
combined_flags(lowP) = strtrim(combined_flags(lowP) + " LOW_PARENT");
Age_Ma(lowP)       = NaN;
Age_1SD_Ma(lowP)   = NaN;
Age_2SD_Ma(lowP)   = NaN;
Age_1SDpct(lowP)   = NaN;
LinearAge_Ma(lowP) = NaN;

% Flag and NaN out grains with invalid production (non-positive He4, U, or P) —
% these fall through both noUTh (finite but non-positive, not non-finite) and
% lowP (which requires P>0). Without this, Age_yr keeps its pre-loop t_linear
% value from He4./max(P,eps), and a negative/zero P divides through the eps
% floor into an astronomically large, meaningless finite age.
invalidParent = ~noUTh & (He4 <= 0 | U_ag <= 0 | P <= 0);
combined_flags(invalidParent) = strtrim(combined_flags(invalidParent) + " INVALID_PARENT");
Age_Ma(invalidParent)       = NaN;
Age_1SD_Ma(invalidParent)   = NaN;
Age_2SD_Ma(invalidParent)   = NaN;
Age_1SDpct(invalidParent)   = NaN;
LinearAge_Ma(invalidParent) = NaN;

% Explicit manual exclusions are preserved as rows with a reason-bearing
% flag. Nothing is deleted, and statistical review flags never reach here.
excludeCell = cellstr(string(excludeGrainIDs));
excludeCell = excludeCell(~cellfun(@isempty, excludeCell));
ageNorm = cellfun(@norm_id, cellstr(string(He.grainid)), 'UniformOutput', false);
excludeNorm = cellfun(@norm_id, excludeCell, 'UniformOutput', false);
manualExclude = ismember(ageNorm, excludeNorm);
if ~isempty(excludeNorm)
    unmatched = excludeNorm(~ismember(excludeNorm, ageNorm));
    if ~isempty(unmatched)
        warning('ladd_age_calc:excludeGrainIDsNotFound', ...
            'excludeGrainIDs not found in joined He rows: %s', strjoin(unmatched, ', '));
    end
end
combined_flags(manualExclude) = strtrim(combined_flags(manualExclude) + " MANUAL_EXCLUDE");
Age_Ma(manualExclude)       = NaN;
Age_1SD_Ma(manualExclude)   = NaN;
Age_2SD_Ma(manualExclude)   = NaN;
Age_1SDpct(manualExclude)   = NaN;
LinearAge_Ma(manualExclude) = NaN;

% ── Build output table ────────────────────────────────────────────────────
GrainID      = He.grainid;
SampleName   = strings(n,1);
if ismember('samplename', He.Properties.VariableNames)
    SampleName = string(He.samplename);
end
RunScript_out = zeros(n,1);
if ismember('runscript', He.Properties.VariableNames)
    RunScript_out = double(He.runscript);
end

ages = table(GrainID, SampleName, RunScript_out, ...
    He4, He4_se, U_ag, U_se, Th_ag, Th_se, Sm_ag, Sm_se, ...
    ThU, Age_Ma, Age_1SD_Ma, Age_2SD_Ma, Age_1SDpct, LinearAge_Ma, ...
    nIter_out, converged, manualExclude, uthMatchCount, uthSelectedRow, ...
    uthSelectedFile, uthMatchDecision, combined_flags, ...
    'VariableNames', { ...
    'GrainID','SampleName','RunScript', ...
    'He4_atoms_g','He4_1SD_atoms_g','U_atoms_g','U_1SD_atoms_g', ...
    'Th_atoms_g','Th_1SD_atoms_g','Sm_atoms_g','Sm_1SD_atoms_g', ...
    'ThU','Age_Ma','Age_1SD_Ma','Age_2SD_Ma','Age_1SDpct', ...
    'LinearAge_Ma','nIter','converged','ManualExclude', ...
    'UThMatchCount','UThSelectedRow','UThSelectedFile','UThMatchDecision', ...
    'flags'});

% ── Print summary ─────────────────────────────────────────────────────────
if verbose
    unk  = ages(ages.RunScript == 3, :);
    good = unk(isfinite(unk.Age_Ma) & unk.Age_Ma > 0, :);
    fprintf('\n── Age Summary (%s) ──────────────────────────────\n', mineral);
    fprintf('  Total grains:    %d\n', height(unk));
    fprintf('  With valid age:  %d\n', height(good));
    fprintf('  Converged:       %d\n', sum(unk.converged));
    if ~isempty(good)
        fprintf('  Age range:       %.1f – %.1f Ma\n', ...
            min(good.Age_Ma), max(good.Age_Ma));
        fprintf('  Median age:      %.1f Ma\n', median(good.Age_Ma,'omitnan'));
        fprintf('  Median 1SD:      %.1f Ma (%.1f%%)\n', ...
            median(good.Age_1SD_Ma,'omitnan'), ...
            median(good.Age_1SDpct,'omitnan'));
    end
    flagged = sum(strlength(strtrim(unk.flags)) > 0);
    if flagged > 0
        fprintf('  Flagged grains:  %d\n', flagged);
    end
    fprintf('\n');
end

% ── Save ──────────────────────────────────────────────────────────────────
if ~isempty(saveAs)
    writetable(ages, saveAs);
    fprintf('Wrote: %s\n', saveAs);
end

end % ── END MAIN ──────────────────────────────────────────────────────────


%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function v = get_numeric(T, names)
    v = NaN(height(T), 1);
    cn = T.Properties.VariableNames;
    for k = 1:numel(names)
        idx = find(strcmpi(cn, names{k}), 1);
        if ~isempty(idx)
            c = T{:,idx};
            if ~isnumeric(c), c = str2double(string(c)); end
            v = c; return;
        end
    end
end

function v = safe_num(T, row, names)
    v = NaN;
    cn = T.Properties.VariableNames;
    for k = 1:numel(names)
        idx = find(strcmpi(cn, names{k}), 1);
        if ~isempty(idx)
            c = T{row, idx};
            if iscell(c), c = c{1}; end
            if isstring(c) || ischar(c), c = str2double(c); end
            v = double(c); return;
        end
    end
end

function s = norm_id(id)
    s = lower(strtrim(char(string(id))));
    s = regexprep(s, '[\s_\-]+', '-');
    s = regexprep(s, '-ft$', '');
    s = regexprep(s, '-+', '-');
end

function names = lower_unique_names(namesIn)
    % Some older audit columns differ only by capitalization (for example,
    % kMAD versus kmad). Lowercase them without creating an invalid table.
    names = matlab.lang.makeUniqueStrings(lower(string(namesIn)),{}, ...
        namelengthmax);
    names = cellstr(names);
end
