function cfg = he_default_cfg(mineral)
% HE_DEFAULT_CFG  Build a run-level config struct for He_reduce.
%
% USAGE
%   cfg = he_default_cfg('zircon');    % or 'apatite'
%
% Then pass cfg as the fourth argument to He_reduce along with three
% file paths: He raw file, pit volume file, and air standard file.
%
% FIELDS SET HERE (do not change)
%   cfg.mineral  — 'zircon' or 'apatite' (applies to entire run)
%   cfg.density  — g/cm3 (zircon 4.65, apatite 3.19)
%
% OPTIONAL FIELDS YOU MAY ADD
%   cfg.renameMap — fix grain ID mismatches that survive normalisation
%                   Format: {'raw_id_1','pitvol_id_1'; 'raw_id_2','pitvol_id_2'}
%                   Only needed if he_normID() cannot resolve the mismatch.
%
% NOTE: FourHeAir is no longer a cfg field. It is read per-shot from the
% air standard CSV (airStdFile argument to He_reduce). See he_airstd_template.csv.

switch lower(mineral)
    case 'zircon'
        cfg.mineral = 'zircon';
        cfg.density = 4.65;   % g/cm3
    case 'apatite'
        cfg.mineral = 'apatite';
        cfg.density = 3.19;   % g/cm3
    otherwise
        error('Unknown mineral "%s". Use ''zircon'' or ''apatite''.', mineral);
end

cfg.renameMap = {};  % empty by default; add entries if needed
end
