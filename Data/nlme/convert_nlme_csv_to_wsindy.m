function DataGroups = convert_nlme_csv_to_wsindy(csv_file, varargin)
% CONVERT_NLME_CSV_TO_WSINDY  Convert R nlme-style grouped growth CSVs into
%   WSINDy struct(s), one per experimental group.
%
%   DataGroups = convert_nlme_csv_to_wsindy(csv_file)
%   DataGroups = convert_nlme_csv_to_wsindy(csv_file, 'Name', Value, ...)
%
%   Handles the standard R datasets package growth CSVs:
%     ChickWeight  – weight vs Time, grouped by Diet, individual = Chick
%     Loblolly     – height vs age, individual = Seed (no higher group)
%     Orange       – circumference vs age, individual = Tree (no higher group)
%
%   The function auto-detects column roles from the known R column names.
%   For non-standard CSVs, override with name-value pairs below.
%
%   OUTPUT
%     DataGroups – struct array with one element per group.  Each element is
%       a WSINDy-compatible Data struct (with .tracks, .t, .growth, etc.)
%       plus an extra field .group_label identifying the group.
%
%       If the dataset has no group column, DataGroups is a scalar struct
%       containing all individuals.
%
%   NAME-VALUE PAIRS
%     'id_col'       ""    Column name for individual ID
%     'time_col'     ""    Column name for time / age
%     'obs_col'      ""    Column name for the observed measurement
%     'group_col'    ""    Column name for grouping variable (empty = none)
%     'assumed_model' "vB" Placeholder model for .growth.type
%     'dataset_name' ""    Label (auto-detected from filename if empty)
%
%   EXAMPLES
%     % Auto-detect for each known dataset:
%     D_chick = convert_nlme_csv_to_wsindy('ChickenData.csv');
%     D_lob   = convert_nlme_csv_to_wsindy('LoblollyData.csv');
%     D_oran  = convert_nlme_csv_to_wsindy('OrangeTreeData.csv');
%
%     % Run one group through WSINDy:
%     results = run_wsindy_experiment(D_chick(1), ...
%                   struct('use_data_for_score', true));
%
%     % Override columns for a custom CSV:
%     D = convert_nlme_csv_to_wsindy('mydata.csv', ...
%             'id_col','SubjectID', 'time_col','Day', ...
%             'obs_col','Mass', 'group_col','Treatment');

%% -------------------- Parse inputs --------------------
p = inputParser;
addRequired(p, 'csv_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'id_col',        "",    @(x) ischar(x)||isstring(x));
addParameter(p, 'time_col',      "",    @(x) ischar(x)||isstring(x));
addParameter(p, 'obs_col',       "",    @(x) ischar(x)||isstring(x));
addParameter(p, 'group_col',     "",    @(x) ischar(x)||isstring(x));
addParameter(p, 'assumed_model', "vB",  @(x) ischar(x)||isstring(x));
addParameter(p, 'dataset_name',  "",    @(x) ischar(x)||isstring(x));
parse(p, csv_file, varargin{:});
opts = p.Results;

%% -------------------- Read CSV --------------------
fprintf('Reading %s...\n', opts.csv_file);
T = readtable(opts.csv_file, 'VariableNamingRule', 'preserve', 'TextType', 'string');
fprintf('  Loaded %d rows x %d columns.\n', height(T), width(T));

cols = string(T.Properties.VariableNames);
cols_lower = lower(cols);
fprintf('  Columns: %s\n', strjoin(cols, ', '));

% Drop R row-index column if present (unnamed first column, often "Var1" or "x")
if cols_lower(1) == "var1" || cols_lower(1) == "x" || cols_lower(1) == ""
    % Check if it's just 1:N
    first_col = T{:,1};
    if isnumeric(first_col) && all(diff(first_col) >= 0)
        fprintf('  Dropping R row-index column "%s".\n', cols(1));
        T(:,1) = [];
        cols = string(T.Properties.VariableNames);
        cols_lower = lower(cols);
    end
end

%% -------------------- Auto-detect dataset & columns --------------------
% Try to identify the dataset from filename or columns
[~, fname, ~] = fileparts(char(opts.csv_file));
fname_lower = lower(fname);

dataset_name = string(opts.dataset_name);
id_col_name  = string(opts.id_col);
time_col_name = string(opts.time_col);
obs_col_name  = string(opts.obs_col);
group_col_name = string(opts.group_col);

if dataset_name == ""
    if contains(fname_lower, "chick")
        dataset_name = "ChickWeight";
    elseif contains(fname_lower, "loblolly")
        dataset_name = "Loblolly";
    elseif contains(fname_lower, "orange")
        dataset_name = "Orange";
    else
        dataset_name = fname;
    end
end
fprintf('  Detected dataset: %s\n', dataset_name);

% Apply known column mappings if user didn't override
switch dataset_name
    case "ChickWeight"
        if id_col_name == "",    id_col_name    = detect_col(cols, cols_lower, ["chick"]); end
        if time_col_name == "",  time_col_name  = detect_col(cols, cols_lower, ["time"]); end
        if obs_col_name == "",   obs_col_name   = detect_col(cols, cols_lower, ["weight"]); end
        if group_col_name == "", group_col_name = detect_col(cols, cols_lower, ["diet"]); end
        obs_units = "gm";  time_units = "days";
        species = "Gallus gallus domesticus (chick)";

    case "Loblolly"
        if id_col_name == "",    id_col_name    = detect_col(cols, cols_lower, ["seed"]); end
        if time_col_name == "",  time_col_name  = detect_col(cols, cols_lower, ["age"]); end
        if obs_col_name == "",   obs_col_name   = detect_col(cols, cols_lower, ["height"]); end
        % No higher-level group for Loblolly
        obs_units = "ft";  time_units = "years";
        species = "Pinus taeda (Loblolly pine)";

    case "Orange"
        if id_col_name == "",    id_col_name    = detect_col(cols, cols_lower, ["tree"]); end
        if time_col_name == "",  time_col_name  = detect_col(cols, cols_lower, ["age"]); end
        if obs_col_name == "",   obs_col_name   = detect_col(cols, cols_lower, ["circumference","circumf","circ"]); end
        % No higher-level group for Orange
        obs_units = "mm";  time_units = "days";
        species = "Citrus sinensis (orange tree)";

    otherwise
        obs_units = "unknown";  time_units = "unknown";
        species = "unknown";
end

% Resolve column indices
id_idx    = find(strcmpi(cols, id_col_name), 1);
time_idx  = find(strcmpi(cols, time_col_name), 1);
obs_idx   = find(strcmpi(cols, obs_col_name), 1);
if group_col_name ~= ""
    group_idx = find(strcmpi(cols, group_col_name), 1);
else
    group_idx = [];
end

if isempty(id_idx),   error('Cannot find ID column. Use ''id_col'' to specify. Columns: %s', strjoin(cols,', ')); end
if isempty(time_idx), error('Cannot find time column. Use ''time_col'' to specify. Columns: %s', strjoin(cols,', ')); end
if isempty(obs_idx),  error('Cannot find observation column. Use ''obs_col'' to specify. Columns: %s', strjoin(cols,', ')); end

fprintf('  ID column      : %s\n', cols(id_idx));
fprintf('  Time column    : %s\n', cols(time_idx));
fprintf('  Obs column     : %s\n', cols(obs_idx));
if ~isempty(group_idx)
    fprintf('  Group column   : %s\n', cols(group_idx));
else
    fprintf('  Group column   : (none)\n');
end

%% -------------------- Extract data --------------------
ids  = string(T{:, id_idx});
t_all = T{:, time_idx};
x_all = T{:, obs_idx};

if iscell(t_all) || isstring(t_all), t_all = str2double(t_all); end
if iscell(x_all) || isstring(x_all), x_all = str2double(x_all); end

if ~isempty(group_idx)
    groups = string(T{:, group_idx});
else
    groups = repmat("all", height(T), 1);
end

% Filter valid rows
valid = isfinite(t_all) & isfinite(x_all) & x_all > 0;
ids    = ids(valid);
t_all  = t_all(valid);
x_all  = x_all(valid);
groups = groups(valid);

fprintf('  %d valid rows.\n', sum(valid));

%% -------------------- Build structs per group --------------------
unique_groups = unique(groups, 'stable');
n_groups = numel(unique_groups);
fprintf('  %d group(s): %s\n', n_groups, strjoin(unique_groups, ', '));

% Placeholder model info
assumed_model = string(opts.assumed_model);
switch assumed_model
    case "vB",       nW = 2; lib_terms = ["1", "x"];
    case "vB23",     nW = 2; lib_terms = ["x^(2/3)", "x"];
    case "logistic", nW = 2; lib_terms = ["x", "x^2"];
    case "gompertz", nW = 2; lib_terms = ["x*ln(x)", "x"];
    case "exp",      nW = 1; lib_terms = "x";
    otherwise,       nW = 2; lib_terms = ["1", "x"];
end

for g = 1:n_groups
    gmask = groups == unique_groups(g);
    g_ids  = ids(gmask);
    g_t    = t_all(gmask);
    g_x    = x_all(gmask);

    [u_ids, ~, igrp] = unique(g_ids, 'stable');
    n_indiv = numel(u_ids);

    xcell = cell(n_indiv, 1);
    tcell = cell(n_indiv, 1);
    keep = true(n_indiv, 1);

    for j = 1:n_indiv
        idx = find(igrp == j);
        [t_sorted, ord] = sort(g_t(idx));
        x_sorted = g_x(idx(ord));

        % Remove duplicate times
        [t_uniq, ia] = unique(t_sorted, 'first');
        x_uniq = x_sorted(ia);

        if numel(t_uniq) < 2
            keep(j) = false;
            continue;
        end

        tcell{j} = t_uniq(:);
        xcell{j} = x_uniq(:);
    end

    xcell = xcell(keep);
    tcell = tcell(keep);
    u_ids = u_ids(keep);
    ntraj = numel(xcell);

    D = struct();
    D.tracks   = xcell(:);
    D.t        = tcell(:);
    D.div_type = "none";

    D.growth.type            = char(assumed_model);
    D.growth.weights         = nan(ntraj, nW);
    D.growth.base_param      = nan(ntraj, 1);
    D.growth.base_param_name = "k";
    D.growth.library_terms   = lib_terms;

    all_x0  = cellfun(@(x) x(1), xcell);
    all_xN  = cellfun(@(x) x(end), xcell);
    all_npts = cellfun(@numel, xcell);
    all_span = cellfun(@(t) t(end) - t(1), tcell);

    D.growth.params = struct( ...
        'dataset',       char(dataset_name), ...
        'species',       char(species), ...
        'obs_units',     char(obs_units), ...
        'time_units',    char(time_units), ...
        'group_label',   char(unique_groups(g)), ...
        'group_col',     char(group_col_name), ...
        'n_individuals', ntraj, ...
        'individual_ids', {cellstr(u_ids)}, ...
        'mean_first_obs', mean(all_x0), ...
        'mean_last_obs',  mean(all_xN), ...
        'mean_n_points',  mean(all_npts), ...
        'mean_time_span', mean(all_span) ...
    );

    D.source      = char(dataset_name);
    D.medium      = char(unique_groups(g));
    D.group_label = char(unique_groups(g));

    if g == 1
        DataGroups = D;
    else
        DataGroups(g) = D;
    end

    fprintf('\n  --- Group "%s" ---\n', unique_groups(g));
    fprintf('    N individuals : %d\n', ntraj);
    fprintf('    Mean points   : %.1f\n', mean(all_npts));
    fprintf('    Mean obs range: %.1f – %.1f %s\n', mean(all_x0), mean(all_xN), obs_units);
    fprintf('    Mean time span: %.1f %s\n', mean(all_span), time_units);
end

fprintf('\n  Assumed model : %s (placeholder)\n', assumed_model);
fprintf('  Output: struct array with %d group(s).\n', n_groups);
fprintf('  Usage:  results = run_wsindy_experiment(DataGroups(1), ...\n');
fprintf('              struct(''use_data_for_score'',true));\n\n');

end

%% ========================= HELPER ======================
function name = detect_col(cols, cols_lower, candidates)
% Return the original-case column name matching the first candidate.
name = "";
for c = 1:numel(candidates)
    idx = find(cols_lower == lower(candidates(c)), 1);
    if ~isempty(idx), name = cols(idx); return; end
end
% Partial match fallback
for c = 1:numel(candidates)
    idx = find(contains(cols_lower, lower(candidates(c))), 1);
    if ~isempty(idx), name = cols(idx); return; end
end
end