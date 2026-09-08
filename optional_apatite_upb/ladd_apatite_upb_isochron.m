function result = ladd_apatite_upb_isochron(inputFile, varargin)
%LADD_APATITE_UPB_ISOCHRON Batch apatite U-Pb common-Pb isochron reduction.
%
% result = ladd_apatite_upb_isochron(inputFile) reads an iolite CSV/XLSX,
% identifies the Wetherill ratios and uncertainties, groups analyses by
% sample, and runs IsoplotR's Ludwig (1998) semitotal-Pb/U regression.
%
% The safest input contains these fields (names need not match exactly):
%   Sample, Analysis, 207Pb/235U, error, 206Pb/238U, error, rho
%
% Name-value options:
%   OutputDir               output folder (default beside input)
%   SampleColumn            explicit sample column name
%   AnalysisColumn          explicit analysis/spot column name
%   Pb207U235Column         explicit ratio column name
%   Pb207U235ErrorColumn    explicit uncertainty column name
%   Pb206U238Column         explicit ratio column name
%   Pb206U238ErrorColumn    explicit uncertainty column name
%   RhoColumn               explicit error-correlation column name
%   IncludeColumn           optional logical include column
%   InputErrors             1se_abs (default), 2se_abs, 1se_pct, 2se_pct
%   MinimumN                minimum valid analyses per sample (default 4)
%   Alpha                   probability cutoff (default 0.05)
%   DisequilibriumMode      none (default) or fixed_th230_u238
%   InitialTh230U238        fixed initial activity ratio (default 1)
%   InitialTh230U238SE      its 1SE uncertainty (default 0)
%   RExecutable             Rscript executable (default Rscript)
%
% With no inputFile, a file picker opens.

if nargin < 1 || strlength(string(inputFile)) == 0
    [name, folder] = uigetfile({'*.csv;*.xlsx;*.xls','Iolite tables (*.csv, *.xlsx, *.xls)'}, ...
        'Choose an apatite U-Pb iolite export');
    if isequal(name, 0)
        error('No input file selected.');
    end
    inputFile = fullfile(folder, name);
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'inputFile', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'SampleColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'AnalysisColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Pb207U235Column', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Pb207U235ErrorColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Pb206U238Column', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Pb206U238ErrorColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'RhoColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'IncludeColumn', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'InputErrors', "1se_abs", @(x) ischar(x) || isstring(x));
addParameter(p, 'MinimumN', 4, @(x) isnumeric(x) && isscalar(x) && x >= 3);
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'DisequilibriumMode', "none", @(x) ischar(x) || isstring(x));
addParameter(p, 'InitialTh230U238', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'InitialTh230U238SE', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'RExecutable', "Rscript", @(x) ischar(x) || isstring(x));
parse(p, inputFile, varargin{:});
opt = p.Results;

inputFile = string(inputFile);
assert(isfile(inputFile), 'Input file does not exist: %s', inputFile);
[inputFolder, inputStem] = fileparts(inputFile);
if strlength(string(opt.OutputDir)) == 0
    runTag = string(datetime('now','Format','yyyyMMdd_HHmmss'));
    outputDir = fullfile(inputFolder, inputStem + "_apatite_upb_isochron_" + runTag);
else
    outputDir = string(opt.OutputDir);
end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

raw = readtable(inputFile, 'VariableNamingRule', 'preserve');
assert(height(raw) > 0, 'The input table contains no data rows.');
names = string(raw.Properties.VariableNames);

sampleCol = resolve_column(names, string(opt.SampleColumn), "sample", false);
analysisCol = resolve_column(names, string(opt.AnalysisColumn), "analysis", false);
ratio75Col = resolve_column(names, string(opt.Pb207U235Column), "ratio75", true);
error75Col = resolve_column(names, string(opt.Pb207U235ErrorColumn), "error75", true);
ratio68Col = resolve_column(names, string(opt.Pb206U238Column), "ratio68", true);
error68Col = resolve_column(names, string(opt.Pb206U238ErrorColumn), "error68", true);
rhoCol = resolve_column(names, string(opt.RhoColumn), "rho", false);
includeCol = resolve_column(names, string(opt.IncludeColumn), "include", false);

if strlength(analysisCol) == 0
    analysis = "row_" + string((1:height(raw))');
else
    analysis = strip(string(raw.(analysisCol)));
end
if strlength(sampleCol) == 0
    sample = repmat(inputStem, height(raw), 1);
    warning('No sample column was found; all analyses will be fit as sample "%s".', inputStem);
else
    sample = strip(string(raw.(sampleCol)));
end
sample(ismissing(sample) | strlength(sample) == 0) = inputStem;
analysis(ismissing(analysis) | strlength(analysis) == 0) = "row_" + string(find(ismissing(analysis) | strlength(analysis) == 0));

ratio75 = numeric_column(raw.(ratio75Col));
error75 = numeric_column(raw.(error75Col));
ratio68 = numeric_column(raw.(ratio68Col));
error68 = numeric_column(raw.(error68Col));
[error75, error68] = convert_errors(ratio75, error75, ratio68, error68, string(opt.InputErrors));

rhoWasAssumed = false(height(raw), 1);
if strlength(rhoCol) == 0
    rho = zeros(height(raw), 1);
    rhoWasAssumed(:) = true;
    warning(['No ratio-error correlation column was found. Rho is set to zero. ' ...
        'The fit will run, but exporting rho from iolite is preferable.']);
else
    rho = numeric_column(raw.(rhoCol));
    badRho = ~isfinite(rho) | abs(rho) > 1;
    rho(badRho) = 0;
    rhoWasAssumed(badRho) = true;
    if any(badRho)
        warning('%d invalid/missing rho values were set to zero.', nnz(badRho));
    end
end

if strlength(includeCol) == 0
    include = true(height(raw), 1);
else
    include = logical_column(raw.(includeCol));
end

normalized = table((1:height(raw))', sample, analysis, ratio75, error75, ...
    ratio68, error68, rho, include, rhoWasAssumed, ...
    'VariableNames', {'InputRow','Sample','Analysis','Pb207U235','SE_Pb207U235', ...
    'Pb206U238','SE_Pb206U238','Rho','Include','RhoWasAssumed'});
normalizedFile = fullfile(outputDir, 'apatite_upb_normalized_input.csv');
writetable(normalized, normalizedFile);

mapping = table(["Sample";"Analysis";"Pb207U235";"SE_Pb207U235"; ...
    "Pb206U238";"SE_Pb206U238";"Rho";"Include"], ...
    [sampleCol;analysisCol;ratio75Col;error75Col;ratio68Col;error68Col;rhoCol;includeCol], ...
    'VariableNames', {'NormalizedField','SourceColumn'});
mapping.SourceColumn(strlength(mapping.SourceColumn) == 0) = "<not found/default used>";
mappingFile = fullfile(outputDir, 'apatite_upb_column_mapping.csv');
writetable(mapping, mappingFile);

diseqMode = lower(string(opt.DisequilibriumMode));
assert(any(diseqMode == ["none","fixed_th230_u238"]), ...
    'DisequilibriumMode must be "none" or "fixed_th230_u238".');

codeRoot = string(fileparts(mfilename('fullpath')));
rScript = fullfile(codeRoot, 'ladd_apatite_upb_isoplotr.R');
assert(isfile(rScript), 'R helper is missing: %s', rScript);
cmd = strjoin([shell_quote(string(opt.RExecutable)), shell_quote(rScript), ...
    "--input", shell_quote(normalizedFile), ...
    "--output", shell_quote(outputDir), ...
    "--alpha", string(opt.Alpha), ...
    "--minimum_n", string(opt.MinimumN), ...
    "--diseq_mode", shell_quote(diseqMode), ...
    "--th230_u238", string(opt.InitialTh230U238), ...
    "--th230_u238_se", string(opt.InitialTh230U238SE)], " ");
[status, message] = system(char(cmd));
if status ~= 0
    error('IsoplotR reduction failed:\n%s', message);
end

summaryFile = fullfile(outputDir, 'apatite_upb_isochron_summary.csv');
correctedFile = fullfile(outputDir, 'apatite_upb_commonPb_corrected_analyses.csv');
excludedFile = fullfile(outputDir, 'apatite_upb_excluded_or_invalid.csv');
methodsFile = fullfile(outputDir, 'apatite_upb_methods.txt');
assert(isfile(summaryFile), 'Reduction finished without creating its summary.');

result = struct;
result.outputDir = outputDir;
result.summaryFile = summaryFile;
result.correctedAnalysesFile = correctedFile;
result.excludedOrInvalidFile = excludedFile;
result.normalizedInputFile = normalizedFile;
result.columnMappingFile = mappingFile;
result.methodsFile = methodsFile;
result.summary = readtable(summaryFile, 'VariableNamingRule', 'preserve');

fprintf('\nApatite U-Pb common-Pb isochron reduction complete.\n');
fprintf('Output folder: %s\n', outputDir);
disp(result.summary(:, intersect({'Sample','N','Status','Age_Ma', ...
    'Age_95CI_preferred_Ma','CommonPb207Pb206','MSWD','P_value'}, ...
    result.summary.Properties.VariableNames, 'stable')));
end

function col = resolve_column(names, explicit, role, required)
if strlength(explicit) > 0
    match = find(strcmpi(names, explicit), 1);
    assert(~isempty(match), 'Requested column "%s" was not found.', explicit);
    col = names(match);
    return
end

n = lower(regexprep(names, '[^A-Za-z0-9]', ''));
isError = contains(n, ["err","error","unc","sigma","se","sd","2s","1s"]);
isAge = contains(n, "age");
switch role
    case "sample"
        exact = ismember(n, ["sample","samplename","sampleid","group","groupname"]);
        score = 100*exact + 10*contains(n,"sample") + 3*contains(n,"group");
    case "analysis"
        exact = ismember(n, ["analysis","analysisname","name","grainid","spot","spotname","label"]);
        score = 100*exact + 10*contains(n,"analysis") + 8*contains(n,"grain") + 6*contains(n,"spot");
    case "ratio75"
        hasRatio = contains(n,"207") & contains(n,"235");
        score = 20*hasRatio - 30*isError - 30*isAge + 5*contains(n,"mean");
    case "error75"
        hasRatio = contains(n,"207") & contains(n,"235");
        score = 20*hasRatio + 20*isError - 30*isAge;
    case "ratio68"
        hasRatio = contains(n,"206") & contains(n,"238");
        score = 20*hasRatio - 30*isError - 30*isAge + 5*contains(n,"mean");
    case "error68"
        hasRatio = contains(n,"206") & contains(n,"238");
        score = 20*hasRatio + 20*isError - 30*isAge;
    case "rho"
        exact = ismember(n, ["rho","rxy","errorcorrelation","correlation"]);
        score = 100*exact + 20*contains(n,"rho") + 10*contains(n,"correlation") + ...
            3*(contains(n,"207") & contains(n,"235") & contains(n,"206") & contains(n,"238"));
    case "include"
        exact = ismember(n, ["include","use","accepted","accept","keep"]);
        score = 100*exact;
    otherwise
        error('Unknown column role: %s', role);
end

[best, idx] = max(score);
threshold = 1;
if any(role == ["ratio75","error75","ratio68","error68"]), threshold = 20; end
if best < threshold
    col = "";
else
    tied = find(score == best);
    if required && numel(tied) > 1
        error('More than one possible %s column was found: %s. Specify it explicitly.', ...
            role, strjoin(names(tied), ', '));
    end
    col = names(idx);
end
if required && strlength(col) == 0
    error('Could not identify the required %s column. Available columns: %s', ...
        role, strjoin(names, ', '));
end
end

function x = numeric_column(v)
if isnumeric(v)
    x = double(v);
else
    x = str2double(strrep(strip(string(v)), ',', ''));
end
x = x(:);
end

function x = logical_column(v)
if islogical(v)
    x = v;
elseif isnumeric(v)
    x = isfinite(v) & v ~= 0;
else
    s = lower(strip(string(v)));
    x = ismember(s, ["1","true","t","yes","y","include","included","keep"]);
end
x = x(:);
end

function [e75, e68] = convert_errors(r75, e75, r68, e68, mode)
mode = lower(mode);
switch mode
    case "1se_abs"
        return
    case "2se_abs"
        e75 = e75 / 2;
        e68 = e68 / 2;
    case "1se_pct"
        e75 = abs(r75) .* e75 / 100;
        e68 = abs(r68) .* e68 / 100;
    case "2se_pct"
        e75 = abs(r75) .* e75 / 200;
        e68 = abs(r68) .* e68 / 200;
    otherwise
        error('InputErrors must be 1se_abs, 2se_abs, 1se_pct, or 2se_pct.');
end
end

function q = shell_quote(x)
x = string(x);
q = "'" + replace(x, "'", "'\"'\"'") + "'";
end
