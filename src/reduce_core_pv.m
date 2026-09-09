function outtbl = reduce_core_pv(folder, metadataCsv, cfg, varargin)
% REDUCE_CORE_PV — generic LA-ICP-MS reducer for U-Th(-Sm)
% Integration is median signal-plateau CPS minus median pre-ablation
% background CPS. Calibrate via NIST612 drift
% interpolation. Intended for use with the
% pit-volume atom count path in ladd_reduce_zircon_pv.
% Works for zircon (no-IS) and apatite (Ca-IS or no-IS, +Sm).
% Configure with CFG fields (see "CFG fields" section below).
%
% Hampel review applied to NIST K* anchor time series before interpolation.
% By default, reviewed NIST analyses remain in the calibration. Automatic
% statistical exclusion is disabled unless explicitly requested; documented
% bad NIST files can be named with excludeNistFiles.

% -------- CFG fields (required/typical) ----------------------------------
% cfg.mineral            : 'zircon' | 'apatite' (used for defaults/messages)
% cfg.density            : g/cm^3 (zircon ~4.65; apatite ~3.19) [optional in this core]
% cfg.mode               : 'si'|'nois' for zircon; 'ca'|'nois' for apatite
% cfg.IS_element         : 'Si'|'Ca' (ignored for 'nois')
% cfg.IS_ppm_std         : IS ppm in calibration standard (e.g., NIST612 Si~3.10e5). 0 if unknown.
% cfg.IS_ppm_unk         : IS ppm in unknown (zircon Si~1.52e5; apatite Ca~3.99e5)
% cfg.pick.U             : channel alias for U (e.g., "238U")
% cfg.pick.Th            : channel alias for Th (e.g., "232Th")
% cfg.pick.IS            : channel alias for internal standard (e.g., "29Si" or "43Ca")
% cfg.pick.Sm (optional) : channel alias for Sm isotope (e.g., "147Sm") for apatite
% cfg.includeSm          : true/false (apatite true; zircon false)
% cfg.Sm_isotope         : '147Sm' if includeSm=true (affects conversion to total Sm)
% cfg.matrixScalarName   : exact reference-material name, or 'NONE'
% ------------------------------------------------------------------------

% ------------------------- Parameters ------------------------------------
p = inputParser;
p.addParameter('kMAD',6);
p.addParameter('nistHalfWin',12,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('autoExcludeNistReviews',false,@(x)islogical(x)||isnumeric(x));
p.addParameter('showNistSummary',true,@(x)islogical(x)||isnumeric(x));
p.addParameter('excludeNistFiles',strings(0,1), ...
    @(x)ischar(x)||isstring(x)||iscell(x));
p.addParameter('minIS',1e5);
p.addParameter('parentScalarInterp','linear');  % linear closer to discrete bracketing
p.addParameter('saveAs','',@ischar);
p.addParameter('SmMaxPlaus', 7500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SmMinCps',   150,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SmRelSEmax', 0.40, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('windowOverridesCsv','',@(x)ischar(x)||isstring(x));


p.parse(varargin{:});
kMAD        = p.Results.kMAD;
nistHalfWin = p.Results.nistHalfWin;
autoExcludeNistReviews = logical(p.Results.autoExcludeNistReviews);
showNistSummary = logical(p.Results.showNistSummary);
excludeNistFiles = string(p.Results.excludeNistFiles);
minIS       = p.Results.minIS;
gInterp     = lower(p.Results.parentScalarInterp);
saveAs      = p.Results.saveAs;
SmMaxPlaus  = p.Results.SmMaxPlaus;
SmMinCps    = p.Results.SmMinCps;
SmRelSEmax  = p.Results.SmRelSEmax;
windowOverridesCsv = string(p.Results.windowOverridesCsv);


% ------------------------- Constants -------------------------------------
NA      = 6.02214076e23;
MW_U    = 238.02891;  MW_Th = 232.03806;  MW_Sm = 150.36;
lam238  = 1.55125e-10; lam235 = 9.8485e-10; lam232 = 4.9475e-11;
lam147  = 6.54e-12;   % 147Sm (yr^-1)
f238    = 0.992742;   f235 = 0.007204;
f147Sm  = 0.1499;     % natural abundance fraction of 147Sm

% ------------------------- Read metadata ---------------------------------
md = readtable(metadataCsv,'TextType','string','ReadVariableNames',true, ...
               'VariableNamingRule','preserve','Delimiter',',');
md.Properties.VariableNames = lower(md.Properties.VariableNames);

needCols = {'file','type'};
assert(all(ismember(needCols, md.Properties.VariableNames)), ...
       'metadata must have columns: file, type');

% stdname should be STRING; known_* should be numeric
if ~ismember('stdname', md.Properties.VariableNames)
    md.stdname = strings(height(md),1);
elseif ~isstring(md.stdname)
    md.stdname = string(md.stdname);
end
if ~ismember('known_u_ppm', md.Properties.VariableNames),  md.known_u_ppm  = nan(height(md),1); end
if ~ismember('known_th_ppm', md.Properties.VariableNames), md.known_th_ppm = nan(height(md),1); end
if ~ismember('known_sm_ppm', md.Properties.VariableNames), md.known_sm_ppm = nan(height(md),1); end

% Optional audit-only/manual window table: file,t0,t1. The default remains
% fully automatic. This is used for controlled comparisons with historical
% manually clipped integrations without altering raw acquisition files.
windowOverrides = table();
if strlength(windowOverridesCsv) > 0
    assert(isfile(windowOverridesCsv), ...
        'reduce_core_pv: window override file not found: %s', windowOverridesCsv);
    windowOverrides = readtable(windowOverridesCsv,'TextType','string', ...
        'VariableNamingRule','preserve','Delimiter',',');
    windowOverrides.Properties.VariableNames = lower(windowOverrides.Properties.VariableNames);
    assert(all(ismember({'file','t0','t1'},windowOverrides.Properties.VariableNames)), ...
        'reduce_core_pv: window override CSV must contain file,t0,t1.');
    windowOverrides.file = lower(regexprep(strtrim(string(windowOverrides.file)), ...
        '^.*[\\/]', ''));
    windowOverrides.t0 = double(windowOverrides.t0);
    windowOverrides.t1 = double(windowOverrides.t1);
end

% ---------------------- Per-analysis extraction --------------------------
rows = []; seqCounter = 0;
for i = 1:height(md)
  f = fullfile(folder, md.file{i});
  if ~isfile(f), warning('Missing file: %s', f); continue; end

  try
    [t, M, colmap] = read_icpms_csv(f);
  catch ME
    warning('Read failed %s: %s', f, ME.message); continue
  end

  u  = pickCol(M,colmap, aliasSet(cfg.pick.U));
  th = pickCol(M,colmap, aliasSet(cfg.pick.Th));
  IS = [];
  if ~strcmpi(cfg.mode,'nois'), IS = pickCol(M,colmap, aliasSet(cfg.pick.IS)); end
  sm = []; if isfield(cfg,'includeSm') && cfg.includeSm && isfield(cfg.pick,'Sm')
    sm = pickCol(M,colmap, aliasSet(cfg.pick.Sm));
  end

  need = ~isempty(u) && ~isempty(th);
  if strcmpi(cfg.mode,'si') || strcmpi(cfg.mode,'ca'), need = need && ~isempty(IS); end
  if ~need, warning('Missing required channels in %s', f); continue; end

  u=u(:); th=th(:); if ~isempty(IS), IS=IS(:); else, IS=zeros(size(u)); end
  if ~isempty(sm), sm=sm(:); else, sm=zeros(size(u)); end

  % Composite for window
  comp = smoothdata( (IS>0).*IS + u + th + sm, 'movmedian', max(5, round(numel(u)/200)) );
  [i0,i1,notes] = robust_window(comp, kMAD);
  if ~isempty(windowOverrides)
      fileNorm = lower(regexprep(strtrim(string(md.file(i))), '^.*[\\/]', ''));
      hitOverride = find(windowOverrides.file == fileNorm);
      assert(numel(hitOverride) <= 1, ...
          'reduce_core_pv: duplicate window overrides for %s.', fileNorm);
      if numel(hitOverride) == 1
          requestedT0 = windowOverrides.t0(hitOverride);
          requestedT1 = windowOverrides.t1(hitOverride);
          idxOverride = find(t >= requestedT0 & t <= requestedT1);
          assert(numel(idxOverride) >= 2, ...
              'reduce_core_pv: override window for %s contains fewer than two rows.', fileNorm);
          i0 = idxOverride(1);
          i1 = idxOverride(end);
          notes = "WINDOW_OVERRIDE";
      end
  end
    if isnan(i0) || isnan(i1) || i1 <= i0
        warning('No ablation window found in %s', f);
        seqCounter = seqCounter + 1;  % still advance sequence
        rows = [rows; makeRow(md(i,:), f, t, i0, i1, ...
                      NaN, NaN, NaN, NaN, NaN, ...      % u, th, IS, RU, RTh
                      NaN, NaN, ...                     % Upp, Thpp
                      "NOWIN", notes, seqCounter, ...   % flag, notes, seq
                      NaN, NaN, NaN, NaN, ...           % cpsU_se, cpsTh_se, cpsIS_se, cpsSm_se
                      NaN)];                             % med_sm
        continue
    end

  pre = 1:max(1,i0-1); win = i0:i1;

  % Qtegra exports each time slice in counts per second (CPS). Preserve that
  % unit by taking the robust signal plateau minus the robust pre-ablation
  % background. Summing the rows would multiply CPS by the number of time
  % slices and make otherwise identical integrations depend on window length.
  med = @(x) median(x(win)) - median(x(pre));
  med_u  = med(u);
  med_th = med(th);
  med_is = med(IS);
  med_sm = med(sm);

  % Robust standard-error estimate used by the validated historical
  % reductions.
  robstd = @(x) 1.4826 * mad(x, 1);
  sig_w  = @(x) robstd(x(win));
  sig_b  = @(x) robstd(x(pre));
  Neff   = max(1, numel(win));
  Neffb  = max(1, numel(pre));
  se_u   = sqrt((sig_w(u).^2)  / Neff + (sig_b(u).^2)  / Neffb);
  se_th  = sqrt((sig_w(th).^2) / Neff + (sig_b(th).^2) / Neffb);
  se_is  = sqrt((sig_w(IS).^2) / Neff + (sig_b(IS).^2) / Neffb);
  se_sm  = sqrt((sig_w(sm).^2) / Neff + (sig_b(sm).^2) / Neffb);

  % Ratios if IS mode
  if ~strcmpi(cfg.mode,'nois')
    RU   = med_u  / max(med_is, eps);
    RTh  = med_th / max(med_is, eps);
    RSm  = NaN; if ~isempty(sm), RSm = med_sm / max(med_is,eps); end
  else
    RU=NaN; RTh=NaN; RSm=NaN;
  end

    seqCounter = seqCounter + 1;
    rows = [rows; makeRow(md(i,:), f, t, i0, i1, ...
                med_u, med_th, med_is, RU, RTh, ...
                NaN, NaN, ...                  % Upp, Thpp (filled later)
                "", notes, seqCounter, ...
                se_u, se_th, se_is, se_sm, ... % all four SEs
                med_sm)];                       % pass the Sm median cps


end

if isempty(rows), outtbl = table(); return; end
raw = struct2table(rows);
raw.integration_mode = repmat("median_cps",height(raw),1);
if strlength(windowOverridesCsv) > 0
    raw.window_override_source = repmat(windowOverridesCsv,height(raw),1);
else
    raw.window_override_source = repmat("",height(raw),1);
end

% NIST calibration membership is saved row-by-row for audit. A review flag
% never changes membership under the default preserve-and-review policy.
raw.nist_u_anchor_role  = repmat("NOT_NIST", height(raw), 1);
raw.nist_th_anchor_role = repmat("NOT_NIST", height(raw), 1);
raw.nist_sm_anchor_role = repmat("NOT_NIST", height(raw), 1);
raw.nist_review_kMAD = repmat(kMAD, height(raw), 1);
raw.nist_review_half_window = repmat(nistHalfWin, height(raw), 1);
raw.nist_auto_exclude_reviews = repmat(autoExcludeNistReviews, height(raw), 1);

% Explicit NIST exclusions use measurement filenames because NIST rows do
% not always carry a stable GrainID. Full paths and basenames are accepted.
rawFileNorm = lower(regexprep(string(raw.file), '^.*[\\/]', ''));
excludeNistNorm = lower(regexprep(strtrim(excludeNistFiles(:)), '^.*[\\/]', ''));
excludeNistNorm = excludeNistNorm(strlength(excludeNistNorm) > 0);



% ---------------------- Calibration selection flags ----------------------
hasStdCol = ismember('stdname', raw.Properties.VariableNames);
is612_type= strcmpi(raw.type,'NIST612');
is612_std = false(height(raw),1); if hasStdCol, is612_std = strcmpi(raw.stdname,'NIST612'); end
isCAL     = is612_type | is612_std;

% ---------------------- Mode: IS vs no-IS -------------------------------
switch lower(cfg.mode)
  case {'si','ca'}    % internal-standard path (zircon: Si, apatite: Ca)
    okU  = raw.cpssi>minIS & isCAL & isfinite(raw.R_U)  & raw.R_U>0;
    okTh = raw.cpssi>minIS & isCAL & isfinite(raw.R_Th) & raw.R_Th>0;

    KU_i  = raw.known_u_ppm(okU)  ./ raw.R_U(okU);
    KTh_i = raw.known_th_ppm(okTh)./ raw.R_Th(okTh);
    KU    = median(KU_i,'omitnan');     KTh = median(KTh_i,'omitnan');

    % calibration scatter
    robstd = @(x) 1.4826*mad(x,1);
    sigKU = robstd(KU_i);  sigKTh = robstd(KTh_i);

    % ratio SE from cps SEs
    Upos = max(raw.cpsu,eps); Thpos = max(raw.cpsth,eps); ISpos = max(raw.cpssi,eps);
    RU_se  = abs(raw.R_U ) .* sqrt( (raw.cpsu_se ./Upos).^2 + (raw.cpssi_se./ISpos).^2 );
    RTh_se = abs(raw.R_Th) .* sqrt( (raw.cpsth_se./Thpos).^2 + (raw.cpssi_se./ISpos).^2 );

    conc_corr = cfg.IS_ppm_unk / max(cfg.IS_ppm_std, eps);  % if std unknown (0), this will zero out → prefer nois
    if cfg.IS_ppm_std==0
      warning('%s-IS requested but IS_ppm_std==0; consider mode="nois".', upper(cfg.IS_element));
      conc_corr = 1; % avoid zeroing; you can switch to 'nois' instead
    end

    raw.U_ppm  = KU  .* raw.R_U  * conc_corr;
    raw.Th_ppm = KTh .* raw.R_Th * conc_corr;
    raw.U_ppm_se  = sqrt( (raw.R_U *conc_corr*sigKU ).^2 + (KU *conc_corr*RU_se ).^2 );
    raw.Th_ppm_se = sqrt( (raw.R_Th*conc_corr*sigKTh).^2 + (KTh*conc_corr*RTh_se).^2 );

% ----- Sm (optional; IS mode; generic anchors) -----
if isfield(cfg,'includeSm') && cfg.includeSm && ismember('R_Sm', raw.Properties.VariableNames)
    hasKnownSm = ismember('known_sm_ppm', raw.Properties.VariableNames);

    okSmAnch = (raw.cpssi > minIS) & isfinite(raw.R_Sm) & (raw.R_Sm > 0);
    if hasKnownSm
        okSmAnch = okSmAnch & isfinite(raw.known_sm_ppm) & (raw.known_sm_ppm > 0);

    end

    if any(okSmAnch)
        % Per-anchor K for Sm (ppm per ratio unit)
        KSm_i = raw.known_sm_ppm(okSmAnch) ./ raw.R_Sm(okSmAnch);

        % ---- Robust trim of anchor outliers (±4 MAD) ----
        robstd = @(x) 1.4826 * mad(x,1);
        medKSm = median(KSm_i, 'omitnan');
        sigKSm = robstd(KSm_i);
        keepSm = abs(KSm_i - medKSm) <= 4 * sigKSm;
        if any(~keepSm)
            fprintf('Removed %d Sm outlier anchors (IS mode)\n', nnz(~keepSm));
            KSm_i = KSm_i(keepSm);
            medKSm = median(KSm_i, 'omitnan');
            sigKSm = robstd(KSm_i);
        end
        % -----------------------------------------------

        KSm    = medKSm;          % central calibration factor
        sigKSm = sigKSm;          % absolute scatter of KSm across anchors

        % Ratio SE from CPS SEs (if cpssm present; otherwise leave NaN)
        ISpos = max(raw.cpssi, eps);
        if ismember('cpssm', raw.Properties.VariableNames) && ismember('cpssm_se', raw.Properties.VariableNames)
            Smpos  = max(raw.cpssm, eps);
            RSm_se = abs(raw.R_Sm) .* sqrt( (raw.cpssm_se./Smpos).^2 + (raw.cpssi_se./ISpos).^2 );
        else
            RSm_se = NaN(height(raw),1);
        end

% Concentration correction (IS ppm unknown in std → conc_corr=1)
conc_corr = cfg.IS_ppm_unk / max(cfg.IS_ppm_std, eps);
if cfg.IS_ppm_std == 0, conc_corr = 1; end

% Apply to all rows
raw.Sm_ppm = KSm .* raw.R_Sm * conc_corr;

% NOTE: raw.R_Sm is a vector; keep element-wise multiply with conc_corr
%       and multiply by scalar sigKSm. RSm_se is a vector.
raw.Sm_ppm_se = sqrt( ((raw.R_Sm .* conc_corr) * sigKSm).^2 + ((KSm * conc_corr) .* RSm_se).^2 );

    end

end
% ---------------- Sm QC: keep outliers, only NaN weak/noisy ---------------
if isfield(cfg,'includeSm') && cfg.includeSm && ismember('Sm_ppm', raw.Properties.VariableNames)
    % thresholds (or override via Name-Value at call)
    if ~exist('SmMaxPlaus','var'),  SmMaxPlaus  = 7500; end  % ppm
    if ~exist('SmMinCps','var'),    SmMinCps    = 150;  end  % cps
    if ~exist('SmRelSEmax','var'),  SmRelSEmax  = 0.40; end  % unitless

    SmWeak  = raw.cpssm < SmMinCps;

    SmRelSE = zeros(height(raw),1);
    if ismember('cpssm_se', raw.Properties.VariableNames)
        SmRelSE = raw.cpssm_se ./ max(raw.cpssm, eps);
    end
    SmNoisy = SmRelSE > SmRelSEmax;

    SmOut   = raw.Sm_ppm > SmMaxPlaus;   % outliers get flagged, NOT NaN’d

    % Only NaN truly unusable values:
    kill = SmWeak | SmNoisy;
    raw.Sm_ppm(kill)    = NaN;
    raw.Sm_ppm_se(kill) = NaN;

    % Append flags
    addflag = strings(height(raw),1);
    addflag(SmWeak) = strtrim(addflag(SmWeak) + " " + "SM_WEAK");
    addflag(SmNoisy)= strtrim(addflag(SmNoisy)+ " " + "SM_UNSTABLE");
    addflag(SmOut  )= strtrim(addflag(SmOut  )+ " " + "SM_OUTLIER");

    if ~ismember('flags', raw.Properties.VariableNames), raw.flags = strings(height(raw),1); end
    needs = strlength(addflag)>0;
    raw.flags(needs) = strtrim(raw.flags(needs) + " " + addflag(needs));
end
% --------------------------------------------------------------------------



% ---- convert 147Sm -> total Sm if metadata/anchors were isotope-specific
if isfield(cfg,'Sm_known_is_isotope') && cfg.Sm_known_is_isotope
    f147Sm = 0.1499; % natural abundance
    raw.Sm_ppm = raw.Sm_ppm ./ f147Sm;
    if ismember('Sm_ppm_se', raw.Properties.VariableNames) && ~all(isnan(raw.Sm_ppm_se))
        raw.Sm_ppm_se = raw.Sm_ppm_se ./ f147Sm;
    end
end

% gentle hint if IS is weak
badIS = (raw.cpssi < minIS) | ~isfinite(raw.R_U) | ~isfinite(raw.R_Th);
raw.flags(badIS) = strtrim(raw.flags(badIS) + " " + "NOIS_FALLBACK_RECOMMENDED");

  case 'nois'        % direct cps->ppm via NIST612 anchors (with seq drift interp)
    seq = raw.seq; is612 = strcmpi(raw.type,'NIST612');
    raw.nist_u_anchor_role(is612) = "UNUSABLE_NIST";
    raw.nist_th_anchor_role(is612) = "UNUSABLE_NIST";
    raw.nist_sm_anchor_role(is612) = "UNUSABLE_NIST";
    explicitNistExclude = is612 & ismember(rawFileNorm, excludeNistNorm);
    okU   = is612 & ~explicitNistExclude & isfinite(raw.cpsu)  & raw.cpsu>0 & isfinite(raw.known_u_ppm);
    okTh  = is612 & ~explicitNistExclude & isfinite(raw.cpsth) & raw.cpsth>0 & isfinite(raw.known_th_ppm);

    KUstar_i  = raw.known_u_ppm(okU)   ./ raw.cpsu(okU);
    KThstar_i = raw.known_th_ppm(okTh) ./ raw.cpsth(okTh);
    seqU = seq(okU); seqTh = seq(okTh);
    idxU = find(okU); idxTh = find(okTh);
    robstd = @(x) 1.4826*mad(x,1);

    % --- Loose Hampel REVIEW on NIST K* time series ---------------------
    % Reviews are retained by default. A NIST analysis is excluded only
    % when named in excludeNistFiles or when the caller deliberately opts
    % into autoExcludeNistReviews=true.
    [KUstar_i, seqU, nU_review, reviewU] = hampel_filter( ...
        KUstar_i, seqU, kMAD, nistHalfWin, autoExcludeNistReviews);
    [KThstar_i, seqTh, nTh_review, reviewTh] = hampel_filter( ...
        KThstar_i, seqTh, kMAD, nistHalfWin, autoExcludeNistReviews);

    raw.nist_u_anchor_role(idxU) = "USED";
    raw.nist_th_anchor_role(idxTh) = "USED";
    if any(reviewU)
        raw.flags(idxU(reviewU)) = strtrim(raw.flags(idxU(reviewU)) + " NIST_U_REVIEW");
        if autoExcludeNistReviews
            raw.nist_u_anchor_role(idxU(reviewU)) = "EXCLUDED_REVIEW";
        else
            raw.nist_u_anchor_role(idxU(reviewU)) = "USED_REVIEW_INCLUDED";
        end
    end
    if any(reviewTh)
        raw.flags(idxTh(reviewTh)) = strtrim(raw.flags(idxTh(reviewTh)) + " NIST_TH_REVIEW");
        if autoExcludeNistReviews
            raw.nist_th_anchor_role(idxTh(reviewTh)) = "EXCLUDED_REVIEW";
        else
            raw.nist_th_anchor_role(idxTh(reviewTh)) = "USED_REVIEW_INCLUDED";
        end
    end
    raw.nist_u_anchor_role(explicitNistExclude) = "EXCLUDED_EXPLICIT";
    raw.nist_th_anchor_role(explicitNistExclude) = "EXCLUDED_EXPLICIT";
    raw.nist_sm_anchor_role(explicitNistExclude) = "EXCLUDED_EXPLICIT";
    raw.flags(explicitNistExclude) = strtrim(raw.flags(explicitNistExclude) + " NIST_EXCLUDE");
    if nU_review > 0
        fprintf('NIST Hampel review: %d U anchors (%s)\n', nU_review, ...
            ternary_text(autoExcludeNistReviews, 'excluded', 'retained'));
    end
    if nTh_review > 0
        fprintf('NIST Hampel review: %d Th anchors (%s)\n', nTh_review, ...
            ternary_text(autoExcludeNistReviews, 'excluded', 'retained'));
    end
    if any(explicitNistExclude)
        fprintf('NIST explicit exclusions: %d file(s)\n', sum(explicitNistExclude));
    end
    % -----------------------------------------------------

    medKU = median(KUstar_i,'omitnan'); medKTh = median(KThstar_i,'omitnan');
    relKU = robstd(KUstar_i)/max(eps,medKU);   relKTh = robstd(KThstar_i)/max(eps,medKTh);

    if numel(seqU)>=2,  KUrow  = interp1(seqU, KUstar_i,  seq, 'pchip','extrap'); else, KUrow  = repmat(medKU, size(seq)); end
    if numel(seqTh)>=2, KThrow = interp1(seqTh,KThstar_i, seq, 'pchip','extrap'); else, KThrow = repmat(medKTh,size(seq)); end
    sigKUrow  = relKU  * KUrow;  sigKThrow = relKTh * KThrow;

    raw.U_ppm  = KUrow  .* raw.cpsu;
    raw.Th_ppm = KThrow .* raw.cpsth;
    raw.U_ppm_se  = sqrt( (KUrow  .* raw.cpsu_se ).^2 + (raw.cpsu  .* sigKUrow ).^2 );
    raw.Th_ppm_se = sqrt( (KThrow .* raw.cpsth_se).^2 + (raw.cpsth .* sigKThrow).^2 );

% ----- Sm (optional; NOIS path; NIST612 anchors on cpssm) -----
if isfield(cfg,'includeSm') && cfg.includeSm && ismember('cpssm', raw.Properties.VariableNames)
    hasKnownSm = ismember('known_sm_ppm', raw.Properties.VariableNames);
    if hasKnownSm
        % Keep the saved "NIST612" comparison genuinely NIST-only. The
        % hybrid reducer performs the separate mineral-bridge step later.
        okSmAnch = is612 & ~explicitNistExclude & ...
                   isfinite(raw.cpssm) & raw.cpssm>0 & ...
                   isfinite(raw.known_sm_ppm) & raw.known_sm_ppm>0;

        if any(okSmAnch)
            KSmstar_i = raw.known_sm_ppm(okSmAnch) ./ raw.cpssm(okSmAnch);
            seqSm     = raw.seq(okSmAnch);
            idxSm     = find(okSmAnch);

            [KSmstar_i, seqSm, nSm_review, reviewSm] = hampel_filter( ...
                KSmstar_i, seqSm, kMAD, nistHalfWin, autoExcludeNistReviews);
            raw.nist_sm_anchor_role(idxSm) = "USED";
            if any(reviewSm)
                raw.flags(idxSm(reviewSm)) = strtrim(raw.flags(idxSm(reviewSm)) + " NIST_SM_REVIEW");
                if autoExcludeNistReviews
                    raw.nist_sm_anchor_role(idxSm(reviewSm)) = "EXCLUDED_REVIEW";
                else
                    raw.nist_sm_anchor_role(idxSm(reviewSm)) = "USED_REVIEW_INCLUDED";
                end
            end
            if nSm_review > 0
                fprintf('NIST Hampel review: %d Sm anchors (%s)\n', nSm_review, ...
                    ternary_text(autoExcludeNistReviews, 'excluded', 'retained'));
            end
            medKSm  = median(KSmstar_i,'omitnan');
            sigKSm  = robstd(KSmstar_i);
            relKSm = sigKSm / max(eps, medKSm);

            % Drift interpolation across full sequence
            if numel(seqSm) >= 2
                KSmrow = interp1(seqSm, KSmstar_i, raw.seq, 'pchip', 'extrap');
            else
                KSmrow = repmat(medKSm, size(raw.seq));
            end
            sigKSmrow = relKSm * KSmrow;

            % Apply to all rows
            raw.Sm_ppm = KSmrow .* raw.cpssm;
            if ismember('cpssm_se', raw.Properties.VariableNames)
                raw.Sm_ppm_se = sqrt( (KSmrow .* raw.cpssm_se).^2 + (raw.cpssm .* sigKSmrow).^2 );
            else
                raw.Sm_ppm_se = NaN(height(raw),1);
            end
        else
            raw.Sm_ppm    = NaN(height(raw),1);
            raw.Sm_ppm_se = NaN(height(raw),1);
        end
    else
        raw.Sm_ppm    = NaN(height(raw),1);
        raw.Sm_ppm_se = NaN(height(raw),1);
    end
else
    raw.Sm_ppm    = NaN(height(raw),1);
    raw.Sm_ppm_se = NaN(height(raw),1);
end

% ---------------- Sm QC: keep outliers, only NaN weak/noisy ---------------
if isfield(cfg,'includeSm') && cfg.includeSm && ismember('Sm_ppm', raw.Properties.VariableNames)
    % thresholds (or override via Name-Value at call)
    if ~exist('SmMaxPlaus','var'),  SmMaxPlaus  = 7500; end  % ppm
    if ~exist('SmMinCps','var'),    SmMinCps    = 150;  end  % cps
    if ~exist('SmRelSEmax','var'),  SmRelSEmax  = 0.40; end  % unitless

    SmWeak  = raw.cpssm < SmMinCps;

    SmRelSE = zeros(height(raw),1);
    if ismember('cpssm_se', raw.Properties.VariableNames)
        SmRelSE = raw.cpssm_se ./ max(raw.cpssm, eps);
    end
    SmNoisy = SmRelSE > SmRelSEmax;

    SmOut   = raw.Sm_ppm > SmMaxPlaus;   % outliers get flagged, NOT NaN’d

    % Only NaN truly unusable values:
    kill = SmWeak | SmNoisy;
    raw.Sm_ppm(kill)    = NaN;
    raw.Sm_ppm_se(kill) = NaN;

    % Append flags
    addflag = strings(height(raw),1);
    addflag(SmWeak) = strtrim(addflag(SmWeak) + " " + "SM_WEAK");
    addflag(SmNoisy)= strtrim(addflag(SmNoisy)+ " " + "SM_UNSTABLE");
    addflag(SmOut  )= strtrim(addflag(SmOut  )+ " " + "SM_OUTLIER");

    if ~ismember('flags', raw.Properties.VariableNames), raw.flags = strings(height(raw),1); end
    needs = strlength(addflag)>0;
    raw.flags(needs) = strtrim(raw.flags(needs) + " " + addflag(needs));
end
% --------------------------------------------------------------------------


% ---- convert 147Sm -> total Sm if metadata/anchors were isotope-specific
if isfield(cfg,'Sm_known_is_isotope') && cfg.Sm_known_is_isotope
    f147Sm = 0.1499; % natural abundance
    raw.Sm_ppm    = raw.Sm_ppm    ./ f147Sm;
    if ismember('Sm_ppm_se', raw.Properties.VariableNames) && ~all(isnan(raw.Sm_ppm_se))
        raw.Sm_ppm_se = raw.Sm_ppm_se ./ f147Sm;
    end
end
end

% ---------------------- Matrix parent scalar (generic) --------------------
matrixScalarName = strtrim(string(cfg.matrixScalarName));
matrixScalarEnabled = isscalar(matrixScalarName) && ~ismissing(matrixScalarName) && ...
    strlength(matrixScalarName) > 0 && ~strcmpi(matrixScalarName, "NONE");
if matrixScalarEnabled
    isMat = hasStdCol & contains(lower(string(raw.stdname)), lower(matrixScalarName));
    hasKnownU  = isfinite(raw.known_u_ppm);
    hasKnownTh = isfinite(raw.known_th_ppm);
    if ismember('known_sm_ppm', raw.Properties.VariableNames)
        hasKnownSm = isfinite(raw.known_sm_ppm);        % logical vector, per-row
    else
        hasKnownSm = false(height(raw),1);              % no column → all false
    end


    % anchors must have measured ppm and known ppm
    ok = isMat & isfinite(raw.U_ppm) & isfinite(raw.Th_ppm) & ...
         isfinite(raw.known_u_ppm) & isfinite(raw.known_th_ppm) & ...
         raw.known_u_ppm > 0 & raw.known_th_ppm >= 0;

    hasSmTerm = isfield(cfg,'includeSm') && cfg.includeSm && ...
                ismember('Sm_ppm', raw.Properties.VariableNames) && any(isfinite(raw.Sm_ppm));

    if hasSmTerm
        ok = ok & isfinite(raw.Sm_ppm) & hasKnownSm & raw.known_sm_ppm >= 0;
    end


    if nnz(ok) >= 1
        % measured & known parents (atoms/g)
        NA = 6.02214076e23; MW_U = 238.02891; MW_Th = 232.03806; MW_Sm = 150.36;
        f238 = 0.992742; f235 = 0.007204; f147Sm = 0.1499;
        lam238 = 1.55125e-10; lam235 = 9.8485e-10; lam232 = 4.9475e-11; lam147 = 6.54e-12;

        Um  = raw.U_ppm(ok)      * 1e-6 * NA / MW_U;   Uk  = raw.known_u_ppm(ok)  * 1e-6 * NA / MW_U;
        Thm = raw.Th_ppm(ok)     * 1e-6 * NA / MW_Th;  Thk = raw.known_th_ppm(ok) * 1e-6 * NA / MW_Th;
        N238m = f238*Um; N235m = f235*Um; NThm = Thm;
        N238k = f238*Uk; N235k = f235*Uk; NThk = Thk;

        if hasSmTerm
            Smm = raw.Sm_ppm(ok)      * 1e-6 * NA / MW_Sm;
            Smk = raw.known_sm_ppm(ok)* 1e-6 * NA / MW_Sm;
            N147m = f147Sm*Smm; N147k = f147Sm*Smk;
        else
            N147m = 0; N147k = 0;
        end

        Rm = 8*lam238.*N238m + 7*lam235.*N235m + 6*lam232.*NThm + 1*lam147.*N147m;
        Rk = 8*lam238.*N238k + 7*lam235.*N235k + 6*lam232.*NThk + 1*lam147.*N147k;

        g_i = Rk ./ Rm;
        seqAll = raw.seq; seqOk = raw.seq(ok);
        if numel(g_i) >= 2
            g_seq = interp1(seqOk, g_i, seqAll, lower(p.Results.parentScalarInterp), 'extrap');
        else
            g_seq = repmat(median(g_i, 'omitnan'), size(seqAll));
        end

        % >>> Apply ONLY to non-NIST rows <<<
        is_nist   = strcmpi(raw.type, 'NIST612');
        applyMask = ~is_nist;

        % ppm
        raw.U_ppm(applyMask)  = raw.U_ppm(applyMask)  .* g_seq(applyMask);
        raw.Th_ppm(applyMask) = raw.Th_ppm(applyMask) .* g_seq(applyMask);

        % ppm SE
        if ismember('U_ppm_se',  raw.Properties.VariableNames) && ~all(isnan(raw.U_ppm_se))
            raw.U_ppm_se(applyMask)  = raw.U_ppm_se(applyMask)  .* g_seq(applyMask);
        end
        if ismember('Th_ppm_se', raw.Properties.VariableNames) && ~all(isnan(raw.Th_ppm_se))
            raw.Th_ppm_se(applyMask) = raw.Th_ppm_se(applyMask) .* g_seq(applyMask);
        end

        % Sm (optional)
        if hasSmTerm
            if ismember('Sm_ppm', raw.Properties.VariableNames)
                raw.Sm_ppm(applyMask) = raw.Sm_ppm(applyMask) .* g_seq(applyMask);
            end
            if ismember('Sm_ppm_se', raw.Properties.VariableNames) && ~all(isnan(raw.Sm_ppm_se))
                raw.Sm_ppm_se(applyMask) = raw.Sm_ppm_se(applyMask) .* g_seq(applyMask);
            end
        end

        fprintf('%s matrix scalar: %.0f anchors; median g = %.3f\n', ...
            char(matrixScalarName), sum(ok), double(median(g_i,'omitnan')));
    else
        fprintf('%s matrix scalar: no valid anchors; skipped.\n', matrixScalarName);
    end
end  % closes: if matrixScalarEnabled

% ---------------------- Global ppm -> atoms/g ----------------------------
toAtoms = @(ppm,MW) ppm * 1e-6 * NA / MW;

if ismember('U_ppm', raw.Properties.VariableNames)
    raw.U_atoms_g = toAtoms(raw.U_ppm,  MW_U);
    if ismember('U_ppm_se', raw.Properties.VariableNames)
        raw.U_atoms_g_se = toAtoms(raw.U_ppm_se, MW_U);
    end
end
if ismember('Th_ppm', raw.Properties.VariableNames)
    raw.Th_atoms_g = toAtoms(raw.Th_ppm, MW_Th);
    if ismember('Th_ppm_se', raw.Properties.VariableNames)
        raw.Th_atoms_g_se = toAtoms(raw.Th_ppm_se, MW_Th);
    end
end
if isfield(cfg,'includeSm') && cfg.includeSm && ismember('Sm_ppm', raw.Properties.VariableNames)
    raw.Sm_atoms_g = toAtoms(raw.Sm_ppm, MW_Sm);
    if ismember('Sm_ppm_se', raw.Properties.VariableNames)
        raw.Sm_atoms_g_se = toAtoms(raw.Sm_ppm_se, MW_Sm);
    end
end

% >>>> Paste NEW code block here <<<<
% ---------------- Sm output guard-rails ----------------
if isfield(cfg,'includeSm') && cfg.includeSm && ...
        ismember('Sm_ppm', raw.Properties.VariableNames)
    relSE = NaN(height(raw),1);
    if ismember('Sm_ppm_se', raw.Properties.VariableNames)
        relSE = raw.Sm_ppm_se ./ max(raw.Sm_ppm, eps);
    end

    hasSmCps = ismember('cpssm', raw.Properties.VariableNames);
    weakCps  = false(height(raw),1);
    if hasSmCps
        weakCps = ~isfinite(raw.cpssm) | raw.cpssm < SmMinCps;
    end

    outlier  = isfinite(raw.Sm_ppm)    & raw.Sm_ppm    > SmMaxPlaus;
    unstable = isfinite(relSE)         & relSE         > SmRelSEmax;

    badMask  = outlier | unstable | weakCps;

    if ~ismember('flags', raw.Properties.VariableNames)
        raw.flags = strings(height(raw),1);
    end
    addFlag = strings(height(raw),1);
    addFlag(outlier)  = strtrim(addFlag(outlier)  + " " + "SM_OUTLIER");
    addFlag(unstable) = strtrim(addFlag(unstable) + " " + "SM_UNSTABLE");
    if hasSmCps
        addFlag(weakCps) = strtrim(addFlag(weakCps) + " " + "SM_WEAK");
    end
    raw.flags = strtrim(raw.flags + " " + addFlag);

    raw.Sm_ppm(badMask)    = NaN;
    if ismember('Sm_ppm_se', raw.Properties.VariableNames)
        raw.Sm_ppm_se(badMask) = NaN;
    end
    if ismember('Sm_atoms_g', raw.Properties.VariableNames)
        raw.Sm_atoms_g(badMask) = NaN;
    end
    if ismember('Sm_atoms_g_se', raw.Properties.VariableNames)
        raw.Sm_atoms_g_se(badMask) = NaN;
    end
end
% -------------------------------------------------------
% ---------------------- Output tidy table --------------------------------
want = {'file','type','stdname','t0','t1','seq', ...
        'cpsu','cpsu_se','cpsth','cpsth_se','cpssi','cpssi_se','cpssm','cpssm_se', ... % ← added
        'R_U','R_Th','R_Sm', ...
        'U_ppm','U_ppm_se','Th_ppm','Th_ppm_se','Sm_ppm','Sm_ppm_se', ...
        'U_atoms_g','U_atoms_g_se','Th_atoms_g','Th_atoms_g_se','Sm_atoms_g','Sm_atoms_g_se', ...
        'nist_u_anchor_role','nist_th_anchor_role','nist_sm_anchor_role', ...
        'nist_review_kMAD','nist_review_half_window','nist_auto_exclude_reviews', ...
        'integration_mode','window_override_source', ...
        'flags','notes'};


have = intersect(want, raw.Properties.VariableNames, 'stable');
outtbl = raw(:, have);

% sanity print (post-scaling NIST should remain ~37/37)
if showNistSummary && all(ismember({'U_ppm','Th_ppm'}, outtbl.Properties.VariableNames))
    is612 = strcmpi(outtbl.type,'NIST612');
    if any(is612)
        fprintf('NIST612 check — median U_ppm: %.2f, Th_ppm: %.2f\n', ...
            median(outtbl.U_ppm(is612),'omitnan'), median(outtbl.Th_ppm(is612),'omitnan'));
    end
end

if ~isempty(saveAs)
    writetable(outtbl, saveAs);
    fprintf('Wrote %s (%.0f rows)\n', char(saveAs), double(height(outtbl)));
end

end  % <<< CLOSES main function: function outtbl = reduce_core_pv(...)

% ====================== helpers ==========================================

function [K_clean, seq_clean, n_review, review] = hampel_filter(K, seq, kMAD, half_win, autoExclude)
% HAMPEL_FILTER  Review calibration anchors in a K* time series.
%
%   K        : calibration factor values (e.g. known_ppm / cpsu) per anchor
%   seq      : sequence index of each anchor (monotonically increasing)
%   kMAD     : review threshold in scaled MADs (default 6)
%   half_win : half-width of sliding window in number of anchors (default 12)
%   autoExclude : if true, reviewed anchors are removed; default false
%
% Returns the K/sequence actually used, review count, and a logical review
% mask on the original inputs. Under the default policy all reviewed rows
% remain in K_clean/seq_clean.

if nargin < 3, kMAD     = 6; end
if nargin < 4, half_win = 12; end
if nargin < 5, autoExclude = false; end

n     = numel(K);
review = false(n, 1);
robstd = @(x) 1.4826 * mad(x, 1);

for ii = 1:n
    lo = max(1,   ii - half_win);
    hi = min(n,   ii + half_win);
    neighbours = K([lo:ii-1, ii+1:hi]);   % exclude self
    if numel(neighbours) < 2, continue; end
    m = median(neighbours, 'omitnan');
    s = robstd(neighbours);
    if s > 0 && abs(K(ii) - m) > kMAD * s
        review(ii) = true;
    end
end

n_review = sum(review);
keep = true(n,1);
if autoExclude
    keep = ~review;
end

% Safety: never drop so many that fewer than 3 anchors remain
if sum(keep) < 3
    warning('hampel_filter: would leave <3 anchors; keeping all %d.', n);
    keep(:) = true;
end

K_clean   = K(keep);
seq_clean = seq(keep);
end

function out = ternary_text(tf, ifTrue, ifFalse)
if tf, out = ifTrue; else, out = ifFalse; end
end

function A = aliasSet(main)
a = string(main);
A = unique([a, regexprep(a,'[^0-9a-zA-Z]+',''), upper(a), lower(a), ...
            a+"_cps", regexprep(a+"_cps",'[^0-9a-zA-Z]+','')]);
end

function [t, M, colmap] = read_icpms_csv(f)
L = readlines(f);
assert(numel(L)>=16,'File too short: %s', f);
hdrLine = strtrim(L(14));
cComma = count(hdrLine, ","); cSemi = count(hdrLine, ";"); cTab = count(hdrLine, char(9));
if cTab >= cComma && cTab >= cSemi, delim = "\t";
elseif cSemi >= cComma,              delim = ";";
else,                                delim = ","; end
headers = strtrim(split(hdrLine, delim))'; ncol = numel(headers);
dataLines = L(16:end); dataLines = dataLines(strlength(strtrim(dataLines))>0);
if isempty(dataLines), error('No data rows after header in %s', f); end
S = split(dataLines, delim); maxc = max(width(S));
if maxc < ncol, S(:, end+1:ncol) = ""; elseif maxc > ncol, S = S(:, 1:ncol); end
M = zeros(size(S), 'double'); for j = 1:ncol, M(:, j) = str2double(S(:, j)); end
norm = lower(regexprep(headers, '[^0-9a-zA-Z]+',''));
colmap = containers.Map; for j = 1:ncol, colmap(norm{j}) = j; end
keys = ["time","times","time_s","time_sec","elapsedtime","seconds","sec","t","acquisitiontime","scantime","timepoint"];
ti = []; for k = 1:numel(keys)
    key = lower(regexprep(keys(k),'[^0-9a-zA-Z]+','')); if isKey(colmap,key), ti = colmap(key); break; end
end
if isempty(ti), anyTime = find(contains(norm,'time'), 1, 'first'); assert(~isempty(anyTime),'No time column'); ti = anyTime; end
t = M(:, ti);
end

function col = pickCol(M, colmap, aliases)
col = [];
for a = 1:numel(aliases)
    key = lower(regexprep(aliases(a),'[^0-9a-zA-Z]+',''));
    if isKey(colmap,key), col = M(:, colmap(key)); return; end
end
end

function [i0,i1,notes] = robust_window(x, kMAD)
% Finds the longest continuous segment above threshold.
% Tail guard: trims i1 back where signal drops to <5% of plateau median,
% preventing extension into post-ablation noise (e.g. after ~55s for 65um pits).
x = x(:); n = numel(x); [~,ord] = sort(x);
baseIdx = ord(1:floor(0.3*n)); base = median(x(baseIdx)); madv = mad(x(baseIdx),1);
thr = base + kMAD*madv; above = x > thr; d = diff([false; above; false]);
starts = find(d==1); ends = find(d==-1)-1;
if isempty(starts), i0 = nan; i1 = nan; notes = "NO_SEGMENT"; return; end
[~,k] = max(ends - starts + 1); i0 = starts(k); i1 = ends(k);

% Tail guard: walk i1 back while last 4 rows are <5% of plateau median.
% Stops signal integration from extending into post-ablation noise.
minPlat = 10;
plat_med = median(x(i0:i1), 'omitnan');
while i1 > i0 + minPlat
    tail_seg = x(max(i0, i1-3):i1);
    if median(tail_seg, 'omitnan') < 0.05 * plat_med
        i1 = i1 - 1;
    else
        break
    end
end

notes = "OK";
end

function row = makeRow(mdrow, f, t, i0, i1, ...
                       u, th, IS, RU, RTh, Upp, Thpp, flag, notes, seq, ...
                       cpsU_se, cpsTh_se, cpsIS_se, cpsSm_se, med_sm)


% Build one output row struct with optional sigma fields & sequence index
row = struct();

% metadata (case-insensitive)
fields = lower(fieldnames(mdrow));
row.file = string(mdrow.file);
row.type = string(mdrow.type);

if ismember('stdname', fields)
    row.stdname = string(mdrow.stdname);
else
    row.stdname = "";
end
if ismember('known_u_ppm', fields)
    row.known_u_ppm = double(mdrow.known_u_ppm);
else
    row.known_u_ppm = NaN;
end
if ismember('known_th_ppm', fields)
    row.known_th_ppm = double(mdrow.known_th_ppm);
else
    row.known_th_ppm = NaN;
end
if ismember('known_sm_ppm', fields)
    row.known_sm_ppm = double(mdrow.known_sm_ppm);
else
    row.known_sm_ppm = NaN;
end

% times
if isnan(i0) || isnan(i1) || i1<=i0 || i0<1 || i1>numel(t)
    row.t0 = NaN; row.t1 = NaN;
else
    row.t0 = t(i0); row.t1 = t(i1);
end

% signals & ratios
row.cpsu     = u;
row.cpsth    = th;
row.cpssi    = IS;
row.cpsu_se  = cpsU_se;
row.cpsth_se = cpsTh_se;
row.cpssi_se = cpsIS_se;

% Sm optional
row.cpssm    = med_sm;    % ← keep the measured Sm cps median
row.cpssm_se = cpsSm_se;


row.R_U  = RU;
row.R_Th = RTh;
row.R_Sm = NaN;

% concentrations / atoms placeholders (filled later)
row.U_ppm     = Upp;
row.Th_ppm    = Thpp;
row.U_ppm_se  = NaN;
row.Th_ppm_se = NaN;
row.Sm_ppm    = NaN;
row.Sm_ppm_se = NaN;

row.U_atoms_g = NaN; row.U_atoms_g_se = NaN;
row.Th_atoms_g = NaN; row.Th_atoms_g_se = NaN;
row.Sm_atoms_g = NaN; row.Sm_atoms_g_se = NaN;

row.flags = string(flag);
row.notes = string(notes);
row.seq   = seq;
end
