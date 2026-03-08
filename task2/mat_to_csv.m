function mat_to_csv(matFile, csvFile, varName)
%MAT_TO_CSV Convert a MAT variable to CSV.
%   mat_to_csv() converts NASDAQ.mat -> NASDAQ.csv using the first variable.
%   mat_to_csv(matFile, csvFile, varName) lets you choose file/variable.

if nargin < 1 || isempty(matFile)
    matFile = 'NASDAQ.mat';
end
if nargin < 2 || isempty(csvFile)
    csvFile = 'NASDAQ.csv';
end

S = load(matFile);
fields = fieldnames(S);
if isempty(fields)
    error('No variables found in MAT file: %s', matFile);
end

if nargin < 3 || isempty(varName)
    varName = fields{1};
end

if ~isfield(S, varName)
    error('Variable "%s" not found in %s. Available: %s', ...
        varName, matFile, strjoin(fields, ', '));
end

X = S.(varName);

if istimetable(X)
    T = timetable2table(X, 'ConvertRowTimes', true);
elseif istable(X)
    T = X;
elseif isstruct(X)
    T = struct2table(X(:));
elseif isnumeric(X) || islogical(X)
    T = array2table(X);
elseif iscell(X)
    T = cell2table(X);
else
    error('Unsupported variable type: %s', class(X));
end

writetable(T, csvFile);
fprintf('Converted %s:%s -> %s (%d rows, %d columns).\n', ...
    matFile, varName, csvFile, height(T), width(T));
end
