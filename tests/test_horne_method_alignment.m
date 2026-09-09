function tests = test_horne_method_alignment
%TEST_HORNE_METHOD_ALIGNMENT Regression checks against Horne et al. (2016).
tests = functiontests(localfunctions);
end

function testPublishedZirconAge(testCase)
% Horne et al. (2016), Table 4, zircon analysis 1 reports 84.2 nmol/g
% 4He, 1719 nmol/g 238U, 1570 nmol/g 232Th, and an age of 31.3 Ma.
NA = 6.02214076e23;
f238 = 0.992742;
grainID = "Horne2016-z1";

he = table(grainID,"Horne2016",3,84.2e-9*NA,0, ...
    'VariableNames',{'GrainID','SampleName','RunScript', ...
    'He4Unk_atoms_g','He4Unk1SD_atoms_g'});
parents = table(grainID,"Unknown",(1719e-9*NA)/f238,0,1570e-9*NA,0, ...
    'VariableNames',{'GrainID','type','u_atoms_g','u_atoms_g_se', ...
    'th_atoms_g','th_atoms_g_se'});

fixtureRoot = tempname;
mkdir(fixtureRoot);
cleanup = onCleanup(@() removeIfPresent(fixtureRoot)); %#ok<NASGU>
heFile = fullfile(fixtureRoot,'he.csv');
parentFile = fullfile(fixtureRoot,'parents.csv');
writetable(he,heFile);
writetable(parents,parentFile);

ages = ladd_age_calc(heFile,parentFile,'mineral','zircon','verbose',false);
verifyEqual(testCase,ages.Age_Ma,31.3,'AbsTol',0.05);
verifyEqual(testCase,ages.ThU,1570/(1719/f238),'RelTol',1e-12);
end

function removeIfPresent(folder)
if isfolder(folder), rmdir(folder,'s'); end
end
