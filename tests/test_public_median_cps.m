function tests = test_public_median_cps
%TEST_PUBLIC_MEDIAN_CPS Self-contained regression for CPS integration units.
tests = functiontests(localfunctions);
end

function testMedianSignalMinusBackground(testCase)
packageRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(packageRoot, 'src'));

fixtureRoot = tempname;
mkdir(fixtureRoot);
cleanup = onCleanup(@() rmdir(fixtureRoot, 's')); %#ok<NASGU>

names = ["calibration_01.csv"; "calibration_02.csv"; ...
    "calibration_03.csv"; "sample_01.csv"];
uSignal = [1010; 1110; 1210; 710];
thSignal = [510; 560; 610; 310];
for k = 1:numel(names)
    writeSyntheticQtegra(fullfile(fixtureRoot, names(k)), ...
        uSignal(k), thSignal(k));
end

type = ["NIST612"; "NIST612"; "NIST612"; "UNKNOWN"];
stdname = ["NIST612"; "NIST612"; "NIST612"; ""];
known_u_ppm = [40; 40; 40; NaN];
known_th_ppm = [40; 40; 40; NaN];
metadata = table(names, type, stdname, known_u_ppm, known_th_ppm, ...
    'VariableNames', {'file','type','stdname','known_u_ppm','known_th_ppm'});
metadataFile = fullfile(fixtureRoot, 'metadata.csv');
writetable(metadata, metadataFile);

cfg = struct('mineral','zircon','density',4.65,'mode','nois', ...
    'IS_element','Si','IS_ppm_std',3.10e5,'IS_ppm_unk',1.52e5, ...
    'pick',struct('U',"238U",'Th',"232Th",'IS',"29Si"), ...
    'includeSm',false,'matrixScalarName',"NONE");
out = reduce_core_pv(fixtureRoot, metadataFile, cfg, ...
    'autoExcludeNistReviews', false);

sampleRow = find(string(out.file) == "sample_01.csv", 1);
verifyNotEmpty(testCase, sampleRow);
% Synthetic background is 10 cps; signal levels above are absolute cps.
verifyEqual(testCase, out.cpsu(sampleRow), 700, 'AbsTol', 1e-12);
verifyEqual(testCase, out.cpsth(sampleRow), 300, 'AbsTol', 1e-12);
verifyEqual(testCase, string(out.integration_mode(sampleRow)), "median_cps");
end

function writeSyntheticQtegra(fileName, uPlateau, thPlateau)
n = 100;
t = (0:n-1)';
u = 10 * ones(n,1);
th = 10 * ones(n,1);
u(31:80) = uPlateau;
th(31:80) = thPlateau;

lines = strings(n + 15, 1);
for k = 1:13
    lines(k) = "synthetic header";
end
lines(14) = "Time,238U,232Th";
lines(15) = "s,cps,cps";
for k = 1:n
    lines(k + 15) = sprintf('%g,%g,%g', t(k), u(k), th(k));
end
writelines(lines, fileName);
end
