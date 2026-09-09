function tests = test_public_parent_options
%TEST_PUBLIC_PARENT_OPTIONS Exercise documented calibration and volume modes.
tests = functiontests(localfunctions);
end

function testLookupPreservesMetadataOrder(testCase)
[fixture, cleanup] = makeFixture(); %#ok<ASGLU>
original = readtable(fixture.apatite.metadata_file, ...
    'VariableNamingRule','preserve','TextType','string');
enriched = ladd_enrich_metadata(fixture.apatite.metadata_file, ...
    'lookupTable',fixture.apatite.reference_lookup_file,'verbose',false);
verifyEqual(testCase,string(enriched.file),string(original.file));
end

function testApatiteBridgeOnly(testCase)
[fixture, cleanup] = makeFixture(); %#ok<ASGLU>
md = readtable(fixture.apatite.metadata_file, ...
    'VariableNamingRule','preserve','TextType','string');
md = md(~strcmpi(string(md.type),'NIST612'),:);
bridgeOnlyMetadata = fullfile(fileparts(fixture.apatite.metadata_file), ...
    'parent_metadata_bridge_only.csv');
writetable(md,bridgeOnlyMetadata);
parentOutput = fullfile(fileparts(fixture.apatite.metadata_file), ...
    'parents_bridge_only.csv');
ageOutput = fullfile(fileparts(fixture.apatite.metadata_file), ...
    'ages_bridge_only.csv');

consoleText = evalc("parents = ladd_reduce_apatite(" + ...
    "fixture.apatite.raw_parent_folder, bridgeOnlyMetadata, " + ...
    "fixture.apatite.he_pit_volume_file, fixture.apatite.uth_pit_volume_file, " + ...
    "'ReferenceLookupFile', fixture.apatite.reference_lookup_file, " + ...
    "'BridgeStandardName', 'ReferenceMaterial', 'SmReferenceBasis', 'total', " + ...
    "'AnchorMode', 'median', 'AllowBridgeOnly', true, " + ...
    "'OutputFile', parentOutput);");
ages = ladd_calculate_ages(fixture.apatite.helium_output,parentOutput, ...
    'Mineral','apatite','OutputFile',ageOutput,'Verbose',false);

verifyTrue(testCase,all(~parents.independent_nist_check_available));
verifyTrue(testCase,all(parents.primary_calibration_path == "BRIDGE_ONLY_NO_NIST"));
verifyEqual(testCase,ages.Age_Ma,fixture.apatite.target_age_Ma,'AbsTol',0.02);
verifyFalse(testCase,contains(consoleText,"NIST612 check — median"));
verifyTrue(testCase,contains(consoleText,"Bridge-only session"));
end

function testZirconBridgeOnlyConsole(testCase)
[fixture, cleanup] = makeFixture(); %#ok<ASGLU>
md = readtable(fixture.zircon.metadata_file, ...
    'VariableNamingRule','preserve','TextType','string');
md = md(~strcmpi(string(md.type),'NIST612'),:);
bridgeOnlyMetadata = fullfile(fileparts(fixture.zircon.metadata_file), ...
    'parent_metadata_bridge_only.csv');
writetable(md,bridgeOnlyMetadata);
parentOutput = fullfile(fileparts(fixture.zircon.metadata_file), ...
    'parents_bridge_only.csv');

consoleText = evalc("parents = ladd_reduce_zircon(" + ...
    "fixture.zircon.raw_parent_folder, bridgeOnlyMetadata, " + ...
    "fixture.zircon.he_pit_volume_file, fixture.zircon.uth_pit_volume_file, " + ...
    "'ReferenceLookupFile', fixture.zircon.reference_lookup_file, " + ...
    "'BridgeStandardName', 'ReferenceMaterial', 'AnchorMode', 'median', " + ...
    "'AllowBridgeOnly', true, 'OutputFile', parentOutput);");

verifyTrue(testCase,all(~parents.independent_nist_check_available));
verifyTrue(testCase,all(parents.parent_calibration_path == "BRIDGE_ONLY_NO_NIST"));
verifyFalse(testCase,contains(consoleText,"NIST612 check — median"));
verifyTrue(testCase,contains(consoleText,"bridge-only signal extraction"));
end

function testSessionAverageVolumeForBothMinerals(testCase)
[fixture, cleanup] = makeFixture(); %#ok<ASGLU>
minerals = ["apatite","zircon"];
for mineral = minerals
    data = fixture.(char(mineral));
    parentOutput = fullfile(fileparts(data.metadata_file), ...
        'parents_session_average.csv');
    ageOutput = fullfile(fileparts(data.metadata_file), ...
        'ages_session_average.csv');
    common = {'ReferenceLookupFile',data.reference_lookup_file, ...
        'BridgeStandardName','ReferenceMaterial','AnchorMode','median', ...
        'UthAverage1SD',20,'OutputFile',parentOutput};
    if mineral == "apatite"
        parents = ladd_reduce_apatite(data.raw_parent_folder,data.metadata_file, ...
            data.he_pit_volume_file,2000,common{:},'SmReferenceBasis','total');
    else
        parents = ladd_reduce_zircon(data.raw_parent_folder,data.metadata_file, ...
            data.he_pit_volume_file,2000,common{:});
    end
    ages = ladd_calculate_ages(data.helium_output,parentOutput, ...
        'Mineral',mineral,'OutputFile',ageOutput,'Verbose',false);
    verifyTrue(testCase,all(parents.uth_pit_volume_mode == "SESSION_AVERAGE"));
    verifyEqual(testCase,ages.Age_Ma,data.target_age_Ma,'AbsTol',0.02);
end
end

function [fixture, cleanup] = makeFixture()
packageRoot = fileparts(fileparts(mfilename('fullpath')));
exampleDir = fullfile(packageRoot,'examples','synthetic');
addpath(exampleDir);
outputRoot = tempname;
cleanup = onCleanup(@() removeIfPresent(outputRoot));
fixture = run_synthetic_examples(outputRoot);
end

function removeIfPresent(folder)
if isfolder(folder), rmdir(folder,'s'); end
end
