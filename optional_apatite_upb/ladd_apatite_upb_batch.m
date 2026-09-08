function result = ladd_apatite_upb_batch(inputDir, varargin)
%LADD_APATITE_UPB_BATCH Reduce every apatite U-Pb CSV in one folder.
%
% result = ladd_apatite_upb_batch opens a folder picker, treats the iolite
% uncertainties as absolute propagated 2SE, and runs one common-Pb isochron
% per Sample value. Raw CSV files are never modified.

if nargin < 1 || strlength(string(inputDir)) == 0
    inputDir = uigetdir(pwd, 'Choose the folder containing apatite U-Pb CSVs');
    if isequal(inputDir, 0), error('No input folder selected.'); end
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'inputDir', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'InputErrors', "2se_abs", @(x) ischar(x) || isstring(x));
addParameter(p, 'ColumnOverrides', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'MinimumN', 4, @(x) isnumeric(x) && isscalar(x) && x >= 3);
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'DisequilibriumMode', "none", @(x) ischar(x) || isstring(x));
addParameter(p, 'InitialTh230U238', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'InitialTh230U238SE', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'RExecutable', "Rscript", @(x) ischar(x) || isstring(x));
parse(p, inputDir, varargin{:});
opt = p.Results;

inputDir = string(inputDir);
assert(isfolder(inputDir), 'Input folder does not exist: %s', inputDir);
if strlength(string(opt.OutputDir)) == 0
    outputDir = fullfile(inputDir, 'isochron_results');
else
    outputDir = string(opt.OutputDir);
end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

columnOverrides = string(opt.ColumnOverrides);
if strlength(columnOverrides) == 0
    savedOverrides = fullfile(outputDir, 'apatite_upb_column_overrides.csv');
    if isfile(savedOverrides), columnOverrides = savedOverrides; end
end

codeRoot = string(fileparts(mfilename('fullpath')));
rScript = fullfile(codeRoot, 'ladd_apatite_upb_batch.R');
assert(isfile(rScript), 'R batch helper is missing: %s', rScript);

parts = [shell_quote(string(opt.RExecutable)), shell_quote(rScript), ...
    "--input_dir", shell_quote(inputDir), "--output_dir", shell_quote(outputDir), ...
    "--input_errors", shell_quote(lower(string(opt.InputErrors))), ...
    "--alpha", string(opt.Alpha), "--minimum_n", string(opt.MinimumN), ...
    "--diseq_mode", shell_quote(lower(string(opt.DisequilibriumMode))), ...
    "--th230_u238", string(opt.InitialTh230U238), ...
    "--th230_u238_se", string(opt.InitialTh230U238SE)];
if strlength(columnOverrides) > 0
    assert(isfile(columnOverrides), 'Column-override file does not exist: %s', columnOverrides);
    parts = [parts, "--column_overrides", shell_quote(columnOverrides)]; %#ok<AGROW>
end

[status, message] = system(char(strjoin(parts, " ")));
if status ~= 0, error('Batch IsoplotR reduction failed:\n%s', message); end

summaryFile = fullfile(outputDir, 'apatite_upb_isochron_summary.csv');
assert(isfile(summaryFile), 'Batch reduction finished without creating its summary.');
result = struct;
result.outputDir = outputDir;
result.summaryFile = summaryFile;
result.correctedAnalysesFile = fullfile(outputDir, 'apatite_upb_commonPb_corrected_analyses.csv');
result.excludedOrInvalidFile = fullfile(outputDir, 'apatite_upb_excluded_or_invalid.csv');
result.inputInventoryFile = fullfile(outputDir, 'apatite_upb_input_inventory.csv');
result.columnMappingFile = fullfile(outputDir, 'apatite_upb_column_mapping.csv');
result.summary = readtable(summaryFile, 'VariableNamingRule', 'preserve');

fprintf('\nApatite U-Pb batch reduction complete.\nOutput folder: %s\n', outputDir);
disp(result.summary(:, intersect({'Sample','N','Status','Age_Ma', ...
    'Age_95CI_preferred_Ma','CommonPb207Pb206','MSWD','P_value'}, ...
    result.summary.Properties.VariableNames, 'stable')));
end

function q = shell_quote(x)
x = string(x);
q = "'" + replace(x, "'", "'\"'\"'") + "'";
end
