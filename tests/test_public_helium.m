function tests = test_public_helium
%TEST_PUBLIC_HELIUM Verify ordinary blank-corrected public He reduction.
tests = functiontests(localfunctions);
end

function testOrdinaryBlankCorrectedInput(testCase)
packageRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(packageRoot, 'src'));

fixtureRoot = tempname;
mkdir(fixtureRoot);
cleanup = onCleanup(@() rmdir(fixtureRoot, 's')); %#ok<NASGU>

raw = table( ...
    ["air-001";"blank-001";"sample-001";"reference-001"], ...
    ["Air";"Blank";"SampleA";"ReferenceMaterial"], ...
    ["air-001";"blank-001";"SampleA-grain-001";"ReferenceMaterial-grain-001"], ...
    [100;10;50;80], [1;1;2;2], ...
    'VariableNames', {'RunID','SampleName','GrainID','He4_cps','He4_1SD'});
rawFile = fullfile(fixtureRoot, 'helium.csv');
writetable(raw, rawFile);

pit = table(["SampleA-grain-001";"ReferenceMaterial-grain-001"], ...
    [1000;1000], [20;20], ...
    'VariableNames', {'GrainID','PitVol_um3','PV1SD_um3'});
pitFile = fullfile(fixtureRoot, 'pit.csv');
writetable(pit, pitFile);

air = table("air-001", 2e8, 1e6, ...
    'VariableNames', {'AirID','FourHeAir','FourHeAir1SD'});
airFile = fullfile(fixtureRoot, 'air.csv');
writetable(air, airFile);

types = table(["Air";"Blank";"SampleA";"ReferenceMaterial"], ...
    [1;2;3;4], 'VariableNames', {'SampleName','RunScript'});
typeFile = fullfile(fixtureRoot, 'types.csv');
writetable(types, typeFile);

out = ladd_reduce_helium(rawFile, pitFile, airFile, typeFile, ...
    'Mineral', 'apatite');
row = find(string(out.GrainID) == "SampleA-grain-001", 1);
verifyNotEmpty(testCase, row);
verifyEqual(testCase, out.He4Unk_atoms(row), 1e10, 'AbsTol', 1e-6);
verifyEqual(testCase, out.BlankMedian_cps(row), 10, 'AbsTol', 1e-12);
end
