function outtbl = ladd_reduce_apatite(rawFolder, metadataFile, ...
    hePitVolumeFile, uthPitVolumeSource, varargin)
%LADD_REDUCE_APATITE Reduce apatite U-Th-Sm from nested pits.
%
% uthPitVolumeSource can be a per-analysis pit-volume CSV or one positive
% session-average U-Th pit volume in cubic micrometres.
%
% Required name-value options:
%   'BridgeStandardName'  exact stdname used as the mineral bridge
%   'SmReferenceBasis'    "total" or "147isotope"
%
% Optional name-value options:
%   'AnchorMode'          nearest, median, nist_following, or
%                         nist_interpolated (default nearest)
%   'ReferenceLookupFile' user-supplied standard concentrations; optional
%                         when known_* values are already in metadata
%   'RunOrder'            explicit raw-file prefix order for metadata sorting
%   'UthAverage1SD'       required when uthPitVolumeSource is numeric
%   'AllowBridgeOnly'     allow bridge calibration when NIST612 was not run
%   'OutputFile'          output CSV path (default: no file)
%   'ExcludeNistFiles'    documented bad NIST measurement filenames
%   'ExcludeBridgeIDs'    documented bad bridge-standard GrainIDs

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'rawFolder', @(x)ischar(x) || isstring(x));
addRequired(p, 'metadataFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'hePitVolumeFile', @(x)ischar(x) || isstring(x));
addRequired(p, 'uthPitVolumeSource', ...
    @(x)isnumeric(x) || ischar(x) || isstring(x));
addParameter(p, 'BridgeStandardName', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'SmReferenceBasis', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'AnchorMode', "nearest", @(x)ischar(x) || isstring(x));
addParameter(p, 'ReferenceLookupFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'RunOrder', strings(0,1), @(x)iscell(x) || isstring(x));
addParameter(p, 'UthAverage1SD', NaN, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'AllowBridgeOnly', false, @(x)islogical(x) || isnumeric(x));
addParameter(p, 'OutputFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'ExcludeNistFiles', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'ExcludeBridgeIDs', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
parse(p, rawFolder, metadataFile, hePitVolumeFile, uthPitVolumeSource, ...
    varargin{:});
opt = p.Results;

assert(isfolder(rawFolder), ...
    'ladd_reduce_apatite: rawFolder not found: %s', rawFolder);
for f = string({metadataFile; hePitVolumeFile})'
    assert(isfile(f), 'ladd_reduce_apatite: input file not found: %s', f);
end
if isnumeric(uthPitVolumeSource)
    assert(isscalar(uthPitVolumeSource) && isfinite(uthPitVolumeSource) && ...
        uthPitVolumeSource > 0, ...
        'ladd_reduce_apatite: numeric U-Th pit volume must be positive.');
    assert(isfinite(opt.UthAverage1SD) && opt.UthAverage1SD > 0, ...
        ['ladd_reduce_apatite: UthAverage1SD must be supplied when using ', ...
         'an average U-Th pit volume.']);
else
    assert(isfile(uthPitVolumeSource), ...
        'ladd_reduce_apatite: U-Th pit-volume file not found: %s', ...
        uthPitVolumeSource);
end
anchorMode = lower(strtrim(string(opt.AnchorMode)));
assert(any(anchorMode == ["nearest","median","nist_following","nist_interpolated"]), ...
    'ladd_reduce_apatite: unrecognized AnchorMode.');
if any(anchorMode == ["nearest","median"])
    assert(strlength(strtrim(string(opt.BridgeStandardName))) > 0, ...
        'ladd_reduce_apatite: BridgeStandardName is required for bridge modes.');
end
assert(any(lower(strtrim(string(opt.SmReferenceBasis))) == ...
    ["total","147isotope"]), ...
    'ladd_reduce_apatite: SmReferenceBasis must be total or 147isotope.');

metadataForRun = char(metadataFile);
cleanupMetadata = [];
lookupFile = strtrim(string(opt.ReferenceLookupFile));
if strlength(lookupFile) > 0
    assert(isfile(lookupFile), ...
        'ladd_reduce_apatite: reference lookup file not found: %s', lookupFile);
    enriched = ladd_enrich_metadata(metadataFile, ...
        'lookupTable', lookupFile, 'runOrder', opt.RunOrder);
    metadataForRun = [tempname, '.csv'];
    writetable(enriched, metadataForRun);
    cleanupMetadata = onCleanup(@() delete_if_exists(metadataForRun));
end

outtbl = ladd_reduce_apatite_pvmeas_hybrid(char(rawFolder), ...
    metadataForRun, char(hePitVolumeFile), uthPitVolumeSource, ...
    'bridgeStandardName', char(string(opt.BridgeStandardName)), ...
    'smReferenceBasis', char(string(opt.SmReferenceBasis)), ...
    'anchorMode', char(string(opt.AnchorMode)), ...
    'uthAverage1sd', opt.UthAverage1SD, ...
    'allowBridgeOnly', logical(opt.AllowBridgeOnly), ...
    'saveAs', char(string(opt.OutputFile)), ...
    'excludeNistFiles', opt.ExcludeNistFiles, ...
    'excludeBridgeIDs', opt.ExcludeBridgeIDs);
clear cleanupMetadata
end

function delete_if_exists(pathText)
if isfile(pathText), delete(pathText); end
end
