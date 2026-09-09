function outtbl = ladd_reduce_helium(rawHeFile, hePitVolumeFile, ...
    airCalibrationFile, sampleTypeFile, varargin)
%LADD_REDUCE_HELIUM Reduce ordinary blank-corrected helium measurements.
%
% rawHeFile must contain the blank-corrected signal exported by the
% acquisition software.
%
% Required name-value option:
%   'Mineral'  "apatite" or "zircon"
%
% Optional name-value options:
%   'OutputFile'       output CSV path (default: no file)
%   'RenameMap'        N-by-2 mapping of source to corrected GrainID
%   'BlankReviewRatio' flag rows when median blank / signal exceeds this
%                      ratio (default 1; review only)
%   'StandardReviewMAD' robust standard-review threshold (default 5)
%   'ExcludeStandardIDs' explicitly documented bad standard GrainIDs

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'rawHeFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'hePitVolumeFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'airCalibrationFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'sampleTypeFile', @(x)ischar(x) || isstring(x));
addParameter(p, 'Mineral', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'OutputFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'RenameMap', strings(0,2), ...
    @(x)iscell(x) || isstring(x));
addParameter(p, 'BlankReviewRatio', 1, ...
    @(x)isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'StandardReviewMAD', 5, ...
    @(x)isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ExcludeStandardIDs', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
parse(p, rawHeFile, hePitVolumeFile, airCalibrationFile, sampleTypeFile, ...
    varargin{:});
opt = p.Results;

mineral = lower(strtrim(string(opt.Mineral)));
assert(any(mineral == ["apatite","zircon"]), ...
    'ladd_reduce_helium: Mineral must be explicitly set to apatite or zircon.');

requiredFiles = string({rawHeFile; hePitVolumeFile; ...
    airCalibrationFile; sampleTypeFile});
for i = 1:numel(requiredFiles)
    assert(isfile(requiredFiles(i)), ...
        'ladd_reduce_helium: input file not found: %s', requiredFiles(i));
end

cfg = he_default_cfg(char(mineral));
if ~isempty(opt.RenameMap)
    renameMap = string(opt.RenameMap);
    assert(size(renameMap,2) == 2, ...
        'ladd_reduce_helium: RenameMap must have two columns.');
    cfg.renameMap = cellstr(renameMap);
end

outtbl = He_reduce(char(rawHeFile), char(hePitVolumeFile), ...
    char(airCalibrationFile), char(sampleTypeFile), cfg, ...
    'saveAs', char(string(opt.OutputFile)), ...
    'minBlankFrac', opt.BlankReviewRatio, ...
    'stdOutlierMAD', opt.StandardReviewMAD, ...
    'stdExcludeIDs', opt.ExcludeStandardIDs);
end
