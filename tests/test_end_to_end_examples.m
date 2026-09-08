function tests = test_end_to_end_examples
%TEST_END_TO_END_EXAMPLES Exercise all public entry points for both minerals.
tests = functiontests(localfunctions);
end

function testApatiteAndZirconWorkflows(testCase)
packageRoot = fileparts(fileparts(mfilename('fullpath')));
exampleDir = fullfile(packageRoot,'examples','synthetic');
addpath(exampleDir);

outputRoot = tempname;
cleanup = onCleanup(@() removeIfPresent(outputRoot)); %#ok<NASGU>
result = run_synthetic_examples(outputRoot);

verifyEqual(testCase,result.apatite.ages.Age_Ma,[50;65],'AbsTol',0.02);
verifyEqual(testCase,result.zircon.ages.Age_Ma,[80;100],'AbsTol',0.02);
verifyTrue(testCase,all(result.apatite.ages.converged));
verifyTrue(testCase,all(result.zircon.ages.converged));
verifyTrue(testCase,isfile(result.apatite.age_output));
verifyTrue(testCase,isfile(result.zircon.age_output));
end

function removeIfPresent(folder)
if isfolder(folder)
    rmdir(folder,'s');
end
end
