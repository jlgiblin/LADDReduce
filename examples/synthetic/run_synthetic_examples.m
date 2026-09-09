function result = run_synthetic_examples(outputRoot)
%RUN_SYNTHETIC_EXAMPLES Run complete generic apatite and zircon workflows.
%
% result = run_synthetic_examples()
% result = run_synthetic_examples(outputRoot)
%
% The function creates small synthetic raw-style tables and Qtegra-style
% time-series files, runs every public LADDReduce entry point, and checks the
% recovered ages against the ages used to generate the helium signals. No
% laboratory or research data are used.

if nargin < 1 || strlength(string(outputRoot)) == 0
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    outputRoot = fullfile(tempdir, ['LADDReduce_synthetic_', stamp]);
end
outputRoot = string(outputRoot);
assert(~isfolder(outputRoot), ...
    'Example output folder already exists; choose a new folder: %s', outputRoot);
mkdir(outputRoot);

exampleDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(exampleDir));
addpath(fullfile(packageRoot, 'src'));

apatite = runMineralExample(outputRoot, "apatite");
zircon = runMineralExample(outputRoot, "zircon");
result = struct('output_root', outputRoot, 'apatite', apatite, 'zircon', zircon);

fprintf('\nSynthetic LADDReduce examples completed successfully.\n');
fprintf('Inputs and outputs: %s\n', outputRoot);
end

function result = runMineralExample(outputRoot, mineral)
mineralDir = fullfile(outputRoot, mineral);
rawDir = fullfile(mineralDir, 'parent_timeseries');
mkdir(rawDir);

switch mineral
    case "apatite"
        density = 3.19;
        bridgeU = 20;
        bridgeTh = 100;
        bridgeSm = 500;
        bridgeCps = [100000, 200000, 150000];
        samplePpm = [30, 60, 400; 45, 90, 600];
        targetAge = [50; 65];
        nistCps = [200000, 80000, 12000];
    case "zircon"
        density = 4.65;
        bridgeU = 100;
        bridgeTh = 50;
        bridgeSm = NaN;
        bridgeCps = [500000, 100000, 0];
        samplePpm = [40, 80, 0; 60, 100, 0];
        targetAge = [80; 100];
        nistCps = [200000, 80000, 0];
    otherwise
        error('Unknown mineral: %s', mineral);
end

sampleIDs = mineral + "-sample-" + compose('%03d',(1:2)');
bridgeIDs = mineral + "-reference-" + compose('%03d',(1:3)');
nistIDs = "nist-" + compose('%03d',(1:3)');
fileNames = ["nist-001.csv"; "reference-001.csv"; "sample-001.csv"; ...
    "nist-002.csv"; "reference-002.csv"; "sample-002.csv"; ...
    "nist-003.csv"; "reference-003.csv"];
types = ["NIST612";"MineralStd";"Unknown";"NIST612"; ...
    "MineralStd";"Unknown";"NIST612";"MineralStd"];
stdnames = ["NIST612";"ReferenceMaterial";"";"NIST612"; ...
    "ReferenceMaterial";"";"NIST612";"ReferenceMaterial"];
grainids = [nistIDs(1);bridgeIDs(1);sampleIDs(1);nistIDs(2); ...
    bridgeIDs(2);sampleIDs(2);nistIDs(3);bridgeIDs(3)];

knownU = NaN(8,1); knownTh = NaN(8,1); knownSm = NaN(8,1);

residualPV = 500;
bridgePV = 1000;
sampleCps = zeros(2,3);
sampleCps(:,1) = (samplePpm(:,1) ./ bridgeU) .* ...
    (residualPV / bridgePV) .* bridgeCps(1);
sampleCps(:,2) = (samplePpm(:,2) ./ bridgeTh) .* ...
    (residualPV / bridgePV) .* bridgeCps(2);
if mineral == "apatite"
    sampleCps(:,3) = (samplePpm(:,3) ./ bridgeSm) .* ...
        (residualPV / bridgePV) .* bridgeCps(3);
end

sampleIndex = 0;
for k = 1:numel(fileNames)
    if types(k) == "NIST612"
        cps = nistCps;
    elseif types(k) == "MineralStd"
        cps = bridgeCps;
    else
        sampleIndex = sampleIndex + 1;
        cps = sampleCps(sampleIndex,:);
    end
    writeQtegra(fullfile(rawDir,fileNames(k)),cps,mineral == "apatite");
end

metadata = table(fileNames,types,stdnames,grainids,knownU,knownTh,knownSm, ...
    'VariableNames',{'file','type','stdname','grainid','known_u_ppm', ...
    'known_th_ppm','known_sm_ppm'});
metadataFile = fullfile(mineralDir,'parent_metadata.csv');
writetable(metadata,metadataFile);

lookupNames = ["NIST612";"ReferenceMaterial"];
lookupU = [40;bridgeU];
lookupTh = [40;bridgeTh];
if mineral == "apatite"
    lookupSm = [40;bridgeSm];
else
    lookupSm = [NaN;NaN];
end
lookupBasis = ["total";"total"];
lookupReference = ["synthetic fixture";"synthetic fixture"];
lookupU1sd = 0.01 .* lookupU;
lookupTh1sd = 0.01 .* lookupTh;
lookupSm1sd = 0.01 .* lookupSm;
referenceLookup = table(lookupNames,lookupU,lookupTh,lookupSm, ...
    lookupU1sd,lookupTh1sd,lookupSm1sd, ...
    lookupBasis,lookupReference,'VariableNames',{'stdname','known_u_ppm', ...
    'known_th_ppm','known_sm_ppm','known_u_1sd_ppm','known_th_1sd_ppm', ...
    'known_sm_1sd_ppm','sm_reference_basis','reference'});
referenceLookupFile = fullfile(mineralDir,'reference_material_lookup.csv');
writetable(referenceLookup,referenceLookupFile);

hePV = table(sampleIDs,repmat(1000,2,1),repmat(20,2,1), ...
    'VariableNames',{'GrainID','PitVol_um3','PV1SD_um3'});
hePVFile = fullfile(mineralDir,'he_pit_volumes.csv');
writetable(hePV,hePVFile);

uthIDs = [sampleIDs;bridgeIDs];
uthVolumes = [repmat(1500,2,1);repmat(bridgePV,3,1)];
uthPV = table(uthIDs,uthVolumes,repmat(20,5,1), ...
    'VariableNames',{'GrainID','PitVol_um3','PV1SD_um3'});
uthPVFile = fullfile(mineralDir,'uth_pit_volumes.csv');
writetable(uthPV,uthPVFile);

fourHeAir = 1e4;
heAtomsG = zeros(2,1);
for k = 1:2
    heAtomsG(k) = heliumAtAge(targetAge(k),samplePpm(k,:),mineral);
end
grainMass = density * 1000 * 1e-12;
heCps = heAtomsG .* grainMass ./ fourHeAir;

heInput = table( ...
    ["air-001";"blank-001";sampleIDs], ...
    ["Air";"Blank";"SampleA";"SampleA"], ...
    ["air-001";"blank-001";sampleIDs], ...
    [100;10;heCps], [1;1;0.01.*heCps], ...
    'VariableNames',{'RunID','SampleName','GrainID','He4_cps','He4_1SD'});
heInputFile = fullfile(mineralDir,'helium_input.csv');
writetable(heInput,heInputFile);

airCalibration = table("air-001",fourHeAir,100, ...
    'VariableNames',{'AirID','FourHeAir','FourHeAir1SD'});
airFile = fullfile(mineralDir,'air_calibration.csv');
writetable(airCalibration,airFile);

sampleTypes = table(["Air";"Blank";"SampleA"],[1;2;3], ...
    'VariableNames',{'SampleName','RunScript'});
typeFile = fullfile(mineralDir,'helium_metadata.csv');
writetable(sampleTypes,typeFile);

heOutput = fullfile(mineralDir,'helium_reduced.csv');
parentOutput = fullfile(mineralDir,'parents_reduced.csv');
ageOutput = fullfile(mineralDir,'ladd_ages.csv');

ladd_reduce_helium(heInputFile,hePVFile,airFile,typeFile, ...
    'Mineral',mineral,'OutputFile',heOutput);
if mineral == "apatite"
    ladd_reduce_apatite(rawDir,metadataFile,hePVFile,uthPVFile, ...
        'BridgeStandardName','ReferenceMaterial', ...
        'ReferenceLookupFile',referenceLookupFile, ...
        'SmReferenceBasis','total','AnchorMode','median', ...
        'OutputFile',parentOutput);
else
    ladd_reduce_zircon(rawDir,metadataFile,hePVFile,uthPVFile, ...
        'BridgeStandardName','ReferenceMaterial','AnchorMode','median', ...
        'ReferenceLookupFile',referenceLookupFile, ...
        'OutputFile',parentOutput);
end
ages = ladd_calculate_ages(heOutput,parentOutput, ...
    'Mineral',mineral,'OutputFile',ageOutput,'Verbose',false);

assert(height(ages) == 2, '%s example did not return two ages.', mineral);
assert(all(ages.converged), '%s example age solver did not converge.', mineral);
assert(all(abs(ages.Age_Ma-targetAge) < 0.02), ...
    '%s example did not recover its synthetic target ages.', mineral);

result = struct('target_age_Ma',targetAge,'ages',ages, ...
    'helium_output',string(heOutput),'parent_output',string(parentOutput), ...
    'age_output',string(ageOutput),'metadata_file',string(metadataFile), ...
    'reference_lookup_file',string(referenceLookupFile), ...
    'he_pit_volume_file',string(hePVFile), ...
    'uth_pit_volume_file',string(uthPVFile), ...
    'raw_parent_folder',string(rawDir));
end

function he = heliumAtAge(ageMa,ppm,mineral)
NA = 6.02214076e23;
MWU = 238.02891;
MWTh = 232.03806;
MWSm = 150.36;
f238 = 0.992742;
f235 = 0.007204;
f147 = 0.1499;
lambda = [1.55125e-10,9.8485e-10,4.9475e-11,6.54e-12];
t = ageMa * 1e6;
uAtomsG = (ppm(1)/1e6)/MWU*NA;
thAtomsG = (ppm(2)/1e6)/MWTh*NA;
he = 8*f238*uAtomsG*(exp(lambda(1)*t)-1) + ...
    7*f235*uAtomsG*(exp(lambda(2)*t)-1) + ...
    6*thAtomsG*(exp(lambda(3)*t)-1);
if mineral == "apatite"
    smAtomsG = (ppm(3)/1e6)/MWSm*NA;
    he = he + f147*smAtomsG*(exp(lambda(4)*t)-1);
end
end

function writeQtegra(fileName,cps,includeSm)
n = 100;
time = (0:n-1)';
background = 10;
u = background*ones(n,1);
th = background*ones(n,1);
u(31:80) = background+cps(1);
th(31:80) = background+cps(2);
if includeSm
    sm = background*ones(n,1);
    sm(31:80) = background+cps(3);
end

lines = strings(n+15,1);
lines(1:13) = "synthetic Qtegra header";
if includeSm
    lines(14) = "Time,238U,232Th,147Sm";
    lines(15) = "s,cps,cps,cps";
    for k = 1:n
        lines(k+15) = sprintf('%g,%g,%g,%g',time(k),u(k),th(k),sm(k));
    end
else
    lines(14) = "Time,238U,232Th";
    lines(15) = "s,cps,cps";
    for k = 1:n
        lines(k+15) = sprintf('%g,%g,%g',time(k),u(k),th(k));
    end
end
writelines(lines,fileName);
end
