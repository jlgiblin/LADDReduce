function outtbl = ladd_reduce_apatite(rawFolder, metadataFile, ...
    hePitVolumeFile, uthPitVolumeFile, varargin)
%LADD_REDUCE_APATITE Reduce apatite U-Th-Sm from measured nested pits.
%
% Required name-value options:
%   'BridgeStandardName'  exact stdname used as the mineral bridge
%   'SmReferenceBasis'    "total" or "147isotope"
%
% Optional name-value options:
%   'AnchorMode'          nearest, median, nist_following, or
%                         nist_interpolated (default nearest)
%   'OutputFile'          output CSV path (default: no file)
%   'ExcludeNistFiles'    documented bad NIST measurement filenames
%   'ExcludeBridgeIDs'    documented bad bridge-standard GrainIDs

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'rawFolder', @(x)ischar(x) || isstring(x));
addRequired(p, 'metadataFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'hePitVolumeFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'uthPitVolumeFile', @(x)ischar(x) || isstring(x));
addParameter(p, 'BridgeStandardName', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'SmReferenceBasis', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'AnchorMode', "nearest", @(x)ischar(x) || isstring(x));
addParameter(p, 'OutputFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'ExcludeNistFiles', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'ExcludeBridgeIDs', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
parse(p, rawFolder, metadataFile, hePitVolumeFile, uthPitVolumeFile, ...
    varargin{:});
opt = p.Results;

assert(isfolder(rawFolder), ...
    'ladd_reduce_apatite: rawFolder not found: %s', rawFolder);
for f = string({metadataFile; hePitVolumeFile; uthPitVolumeFile})'
    assert(isfile(f), 'ladd_reduce_apatite: input file not found: %s', f);
end
assert(strlength(strtrim(string(opt.BridgeStandardName))) > 0, ...
    'ladd_reduce_apatite: BridgeStandardName must be explicitly supplied.');
assert(any(lower(strtrim(string(opt.SmReferenceBasis))) == ...
    ["total","147isotope"]), ...
    'ladd_reduce_apatite: SmReferenceBasis must be total or 147isotope.');

outtbl = ladd_reduce_apatite_pvmeas_hybrid(char(rawFolder), ...
    char(metadataFile), char(hePitVolumeFile), char(uthPitVolumeFile), ...
    'bridgeStandardName', char(string(opt.BridgeStandardName)), ...
    'smReferenceBasis', char(string(opt.SmReferenceBasis)), ...
    'anchorMode', char(string(opt.AnchorMode)), ...
    'saveAs', char(string(opt.OutputFile)), ...
    'excludeNistFiles', opt.ExcludeNistFiles, ...
    'excludeBridgeIDs', opt.ExcludeBridgeIDs);
end
