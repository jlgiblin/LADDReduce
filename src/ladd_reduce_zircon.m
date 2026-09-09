function outtbl = ladd_reduce_zircon(rawFolder, metadataFile, ...
    hePitVolumeFile, uthPitVolumeSource, varargin)
%LADD_REDUCE_ZIRCON Reduce zircon U-Th from measured or declared nested pits.
%
% uthPitVolumeSource can be a per-analysis pit-volume CSV or one positive
% session-average U-Th pit volume in cubic micrometres.
%
% Required name-value option:
%   'BridgeStandardName'  exact stdname used as the zircon bridge
%
% Optional name-value options:
%   'AnchorMode'          following, median, or single (default following)
%   'SingleBridgeID'      required when AnchorMode is single
%   'UthAverage1SD'       required when uthPitVolumeSource is numeric
%   'ReferenceLookupFile' user-supplied standard concentrations; optional
%                         when known_* values are already in metadata
%   'RunOrder'            explicit raw-file prefix order for metadata sorting
%   'AllowBridgeOnly'     explicitly allow a run with no NIST glass
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
addParameter(p, 'AnchorMode', "following", @(x)ischar(x) || isstring(x));
addParameter(p, 'SingleBridgeID', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'UthAverage1SD', NaN, ...
    @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'ReferenceLookupFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'RunOrder', strings(0,1), @(x)iscell(x) || isstring(x));
addParameter(p, 'AllowBridgeOnly', false, ...
    @(x)islogical(x) || isnumeric(x));
addParameter(p, 'OutputFile', "", @(x)ischar(x) || isstring(x));
addParameter(p, 'ExcludeNistFiles', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'ExcludeBridgeIDs', strings(0,1), ...
    @(x)ischar(x) || isstring(x) || iscell(x));
parse(p, rawFolder, metadataFile, hePitVolumeFile, uthPitVolumeSource, ...
    varargin{:});
opt = p.Results;

assert(isfolder(rawFolder), ...
    'ladd_reduce_zircon: rawFolder not found: %s', rawFolder);
assert(isfile(metadataFile), ...
    'ladd_reduce_zircon: metadata file not found: %s', metadataFile);
assert(isfile(hePitVolumeFile), ...
    'ladd_reduce_zircon: He pit-volume file not found: %s', hePitVolumeFile);
if ~isnumeric(uthPitVolumeSource)
    assert(isfile(uthPitVolumeSource), ...
        'ladd_reduce_zircon: U-Th pit-volume file not found: %s', ...
        uthPitVolumeSource);
end
assert(strlength(strtrim(string(opt.BridgeStandardName))) > 0, ...
    'ladd_reduce_zircon: BridgeStandardName must be explicitly supplied.');

metadataForRun = char(metadataFile);
cleanupMetadata = [];
lookupFile = strtrim(string(opt.ReferenceLookupFile));
if strlength(lookupFile) > 0
    assert(isfile(lookupFile), ...
        'ladd_reduce_zircon: reference lookup file not found: %s', lookupFile);
    enriched = ladd_enrich_metadata(metadataFile, ...
        'lookupTable', lookupFile, 'runOrder', opt.RunOrder);
    metadataForRun = [tempname, '.csv'];
    writetable(enriched, metadataForRun);
    cleanupMetadata = onCleanup(@() delete_if_exists(metadataForRun));
end

outtbl = ladd_reduce_zircon_pvmeas(char(rawFolder), metadataForRun, ...
    char(hePitVolumeFile), uthPitVolumeSource, ...
    'bridgeStandardName', char(string(opt.BridgeStandardName)), ...
    'anchorMode', char(string(opt.AnchorMode)), ...
    'singleBridgeID', char(string(opt.SingleBridgeID)), ...
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
