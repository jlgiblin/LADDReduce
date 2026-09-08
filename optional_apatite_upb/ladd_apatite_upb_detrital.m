function result = ladd_apatite_upb_detrital(inputDir, varargin)
%LADD_APATITE_UPB_DETRITAL Correct detrital apatite U-Pb grains individually.
%
% result = ladd_apatite_upb_detrital opens a folder picker, normalizes all
% iolite CSV files, applies an age-dependent Stacey-Kramers common-Pb
% correction to each grain, and models each sample's corrected age
% distribution with one to four error-aware Gaussian components.

if nargin < 1 || strlength(string(inputDir)) == 0
    inputDir = uigetdir(pwd, 'Choose the folder containing detrital apatite U-Pb CSVs');
    if isequal(inputDir, 0), error('No input folder selected.'); end
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'inputDir', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'NormalizationDir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'InputErrors', "2se_abs", @(x) ischar(x) || isstring(x));
addParameter(p, 'ColumnOverrides', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'MaxComponents', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 6 && x == fix(x));
addParameter(p, 'RefreshNormalizedInput', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'RExecutable', "Rscript", @(x) ischar(x) || isstring(x));
parse(p, inputDir, varargin{:});
opt = p.Results;

inputDir = string(inputDir);
assert(isfolder(inputDir), 'Input folder does not exist: %s', inputDir);
if strlength(string(opt.NormalizationDir)) == 0
    normalizationDir = fullfile(inputDir, 'isochron_results');
else
    normalizationDir = string(opt.NormalizationDir);
end
if strlength(string(opt.OutputDir)) == 0
    outputDir = fullfile(inputDir, 'detrital_results');
else
    outputDir = string(opt.OutputDir);
end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

normalizedFile = fullfile(normalizationDir, 'apatite_upb_normalized_input.csv');
if opt.RefreshNormalizedInput || ~isfile(normalizedFile)
    batchArgs = {'OutputDir', normalizationDir, 'InputErrors', string(opt.InputErrors), ...
        'RExecutable', string(opt.RExecutable), 'DisequilibriumMode', 'none'};
    if strlength(string(opt.ColumnOverrides)) > 0
        batchArgs = [batchArgs, {'ColumnOverrides', string(opt.ColumnOverrides)}]; %#ok<AGROW>
    end
    ladd_apatite_upb_batch(inputDir, batchArgs{:});
end
assert(isfile(normalizedFile), 'Normalized ratio file was not created: %s', normalizedFile);

codeRoot = string(fileparts(mfilename('fullpath')));
rScript = fullfile(codeRoot, 'ladd_apatite_upb_detrital.R');
assert(isfile(rScript), 'Detrital R helper is missing: %s', rScript);

parts = [shell_quote(string(opt.RExecutable)), shell_quote(rScript), ...
    "--input", shell_quote(normalizedFile), "--output", shell_quote(outputDir), ...
    "--max_components", string(opt.MaxComponents)];
isochronSummary = fullfile(normalizationDir, 'apatite_upb_isochron_summary.csv');
if isfile(isochronSummary)
    parts = [parts, "--isochron_summary", shell_quote(isochronSummary)]; %#ok<AGROW>
end

[status, message] = system(char(strjoin(parts, " ")));
if status ~= 0, error('Detrital apatite U-Pb reduction failed:\n%s', message); end

summaryFile = fullfile(outputDir, 'apatite_upb_detrital_distribution_summary.csv');
modesFile = fullfile(outputDir, 'apatite_upb_detrital_population_modes.csv');
grainsFile = fullfile(outputDir, 'apatite_upb_detrital_corrected_grains.csv');
assert(isfile(summaryFile) && isfile(modesFile) && isfile(grainsFile), ...
    'Detrital reduction finished without creating all expected outputs.');

result = struct;
result.outputDir = outputDir;
result.summaryFile = summaryFile;
result.populationModesFile = modesFile;
result.modelSelectionFile = fullfile(outputDir, 'apatite_upb_detrital_model_selection.csv');
result.correctedGrainsFile = grainsFile;
result.methodsFile = fullfile(outputDir, 'apatite_upb_detrital_methods.txt');
result.summary = readtable(summaryFile, 'VariableNamingRule', 'preserve');
result.populationModes = readtable(modesFile, 'VariableNamingRule', 'preserve');
result.correctedGrains = readtable(grainsFile, 'VariableNamingRule', 'preserve');
result.plotFiles = dir(fullfile(outputDir, '*_detrital_age_distribution.png'));

fprintf('\nDetrital apatite U-Pb reduction complete.\nOutput folder: %s\n', outputDir);
wanted = {'Sample','N_corrected','Corrected_mean_Ma','Corrected_median_Ma', ...
    'Recommended_K','Candidate_tight_modes','Distribution_interpretation'};
disp(result.summary(:, intersect(wanted, result.summary.Properties.VariableNames, 'stable')));
end

function q = shell_quote(x)
x = string(x);
q = "'" + replace(x, "'", "'\"'\"'") + "'";
end
