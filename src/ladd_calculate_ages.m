function ages = ladd_calculate_ages(heliumFile, parentFile, varargin)
%LADD_CALCULATE_AGES Match LADD outputs and calculate sub-grain He dates.
%
% Required name-value option:
%   'Mineral'  "apatite" or "zircon"
%
% Optional name-value options:
%   'OutputFile'       output CSV path (default: no file)
%   'ExcludeGrainIDs'  explicitly documented GrainIDs to withhold
%   'Verbose'          print matching and date summary (default true)

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'heliumFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'parentFile', @(x)ischar(x) || isstring(x));
addParameter(p, 'Mineral', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'OutputFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'ExcludeGrainIDs', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'Verbose', true, @(x)islogical(x) && isscalar(x));
parse(p, heliumFile, parentFile, varargin{:});
opt = p.Results;

mineral = lower(strtrim(string(opt.Mineral)));
assert(any(mineral == ["apatite","zircon"]), ...
    'ladd_calculate_ages: Mineral must be explicitly set to apatite or zircon.');
assert(isfile(heliumFile), ...
    'ladd_calculate_ages: helium file not found: %s', heliumFile);
assert(isfile(parentFile), ...
    'ladd_calculate_ages: parent file not found: %s', parentFile);

ages = ladd_age_calc(char(heliumFile), char(parentFile), ...
    'mineral', char(mineral), ...
    'saveAs', char(string(opt.OutputFile)), ...
    'excludeGrainIDs', opt.ExcludeGrainIDs, ...
    'verbose', opt.Verbose);
end
