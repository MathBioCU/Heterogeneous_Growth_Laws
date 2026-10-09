clear all; 
close all;
addpath(genpath('../wsindy_obj_base/'));
addpath(genpath('./Utils/'));

%% artificial Data
wanted_div_type = "ind";
wanted_growth_type = "logistic"; % exp | logistic | gompertz | vB | vB23

opts = struct();

dr = "./Data/";
filename = wanted_div_type + "_"+wanted_growth_type+"_"+"artificial.mat";
try
    load(dr+ filename)        % loads struct "Data"
catch
    filename = wanted_div_type + "_"+wanted_growth_type+"_"+"clean.mat";
    load(dr+filename)
end
opts.noise_ratio = 0.05 ; % additive gaussian noise iid accross trajectories
opts.want_shrinkage = true; % run mixed effect shrinkage
opts.want_linargs = true; % sign constraints on parameters
opts.n_select = 500; % number of individuals

%% Real Data
% Data_choice = 1;
% opts = struct();
% switch Data_choice
% 
%     case 1 % Loblolly pine tree 
%         dr = "./Data/";
%         filename = "nlme/LoblollyData.mat";
%         opts.want_interp = true;
% 
% end
% opts.noise_ratio = 0.0;
% opts.want_shrinkage = true;
% opts.want_linargs = true;
% opts.n_select = 500;

%%=======
try
    load(dr+ filename); % loads struct "Data"
catch
    error("Can't find file"+dr+filename)
end

disp(filename)

% defaults
dflt = struct( ...
    'subsample_k',     1,        ...
    'min_length',      6,        ...
    'n_select',        [],       ...
    'noise_ratio',     0.00,     ...
    'tf_power',        15,       ...
    'want_shrinkage',  false,     ...
    'max_iter',        10,       ...
    'min_valid_frac',  0.5,      ...
    'want_linargs',    true,     ...
    'cell_idx',        [],       ...
    'verbose',         true,     ...
    'want_plots',      true,     ...
    'want_interp',     false,    ...   %  toggle interpolation
    'interp_n_points', 150,      ...   %  target points per trajectory
    'interp_method',   'linear');      % 'pchip' | 'spline' | 'linear' | 'makima'
flds = fieldnames(dflt);
for f = 1:numel(flds)
    if ~isfield(opts, flds{f}), opts.(flds{f}) = dflt.(flds{f}); end
end

%% Prep Data
xcell = Data.tracks;  tcell = Data.t;
k = opts.subsample_k;
xcell = cellfun(@(x) x(1:k:end), xcell, 'UniformOutput', false);
tcell = cellfun(@(t) t(1:k:end), tcell, 'UniformOutput', false);


valid_idx = cellfun(@(x) numel(x) > 5 && all(isfinite(x)), xcell);
xcell = xcell(valid_idx);  tcell = tcell(valid_idx);
track_lens = cellfun(@numel, xcell);
valid_idx  = track_lens >= opts.min_length;
xcell = xcell(valid_idx);  tcell = tcell(valid_idx);
true_weights = Data.growth.weights(valid_idx,:);

[tcell, xcell] = check_uniform_dt(tcell, xcell);
tcell_raw = tcell; xcell_raw = xcell;
% --- Optional uniform-mesh interpolation -----------------------------
if opts.want_interp
    [tcell, xcell] = interpolate_trajectories(tcell, xcell, ...
                        opts.interp_n_points, opts.interp_method);
end
% ---------------------------------------------------------------------

n_total  = numel(xcell);


if isempty(opts.n_select)
    n_select = n_total;
else
    n_select = min(opts.n_select, n_total);
end
rng('shuffle');
if isempty(opts.cell_idx)
    rand_idx = randperm(n_total, n_select);
else
    rand_idx = opts.cell_idx;
end
xcell = xcell(rand_idx);  tcell = tcell(rand_idx);
xcell_raw = xcell_raw(rand_idx);  tcell_raw = tcell_raw(rand_idx);

xcell = cellfun(@(x) x(:), xcell, 'UniformOutput', false);
tcell = cellfun(@(t) t(:), tcell, 'UniformOutput', false);
xcell_raw = cellfun(@(x) x(:), xcell_raw, 'UniformOutput', false);
tcell_raw = cellfun(@(t) t(:), tcell_raw, 'UniformOutput', false);

N = numel(xcell);
Uobj = arrayfun(@(i) wsindy_data(xcell{i}, tcell{i}), (1:N)');
for i = 1:N
    Ui = Uobj(i);
    Ui.addnoise(opts.noise_ratio);
    ob = [Ui.Uobs{:}];
    if any(ob <= 0)
        new_ob = max(ob, eps);
        Ui.noise{1,1} = Ui.noise{1,1} + (new_ob - ob);   % adjust stored noise
        Ui.Uobs = {new_ob};
    end
    Uobj(i) = Ui;
end

%% Lib
t_const = term('fhandle', @(x) ones(size(x)));
t_x     = term('fhandle', @(x) x);
t_x2    = term('fhandle', @(x) x.^2);
t_x23   = term('fhandle', @(x) x.^(2/3));
t_xlogx_grad = term('fhandle', @(x) log(max(x,eps)) + 1, ...
                     'gradterms', term('fhandle', @(x) 1 ./ x));
t_xlogx = term('fhandle', @(x) x .* log(x), ...
               'gradterms', t_xlogx_grad);
lib  = library('terms', {t_const, t_x, t_x2, t_x23, t_xlogx});
tags = ["const", "x", "x^2", "x^{2/3}", "x log(x)"];

%% TF
tf_full = cell(size(Uobj));
for i = 1:numel(Uobj)
    tf_full{i} = testfcn(Uobj(i), 'phifuns', @(v) (1-v.^2).^opts.tf_power);
end

%% ============== Combined WENDy + BIC Model Selection ===============
% For every candidate growth law:
%      (1) build a reduced library containing only that law's support,
%      (2) run the full WENDy estimator on every trajectory (with linargs
%          if opts.want_linargs is true), and
%      (3) score the law by BIC using the WENDy-converged weights.
%  The law with minimum BIC is selected; its WENDy fits become the final
%  per-trajectory parameters consumed by the rest of the script.
%
%  Library term ordering:
%      1 -> 1     2 -> x     3 -> x^2     4 -> x^{2/3}     5 -> x*log(x)
% ====================================================================

candidate_laws = {'exponential','logistic','gompertz','vb_linear','vb_metabolic'};
n_laws = numel(candidate_laws);

term_map = struct();
term_map.exponential   = 2;        % {x}
term_map.logistic      = [2, 3];   % {x, x^2}
term_map.gompertz      = [2, 5];   % {x, x*log(x)}
term_map.vb_linear     = [1, 2];   % {1, x}
term_map.vb_metabolic  = [2, 4];   % {x, x^{2/3}}  (lib order)

wendy_results  = cell(n_laws, 1);
bic_per_law    = inf(n_laws, 1);
rss_per_law    = inf(n_laws, 1);
np_per_law     = zeros(n_laws, 1);
K_per_law      = zeros(n_laws, 1);
n_obs_per_law  = zeros(n_laws, 1);

% --- Additional fixed-effects score diagnostics ---
bic_global_fe_per_law = inf(n_laws, 1);   % current global fixed-effects BIC
bic_sumtraj_per_law   = inf(n_laws, 1);   % sum_i BIC_i
aic_sumtraj_per_law   = inf(n_laws, 1);   % sum_i AIC_i
biceff_sumtraj_per_law= inf(n_laws, 1);   % sum_i BIC with Keff
rss_mean_per_law      = inf(n_laws, 1);   % RSS / weak observation
bic_i_all             = cell(n_laws, 1);  % individual BIC_i values
biceff_i_all          = cell(n_laws, 1);  % individual BIC_i values using Keff
rss_i_all             = cell(n_laws, 1);  % individual RSS_i values
Ki_all                = cell(n_laws, 1);  % weak observations per trajectory
Keff_all              = cell(n_laws, 1);  % effective weak dof per trajectory

% --- Which score drives the final model choice ----------------------------
%  'sumtraj_keff' : sum_i BIC_i with Ki replaced by the effective number of
%                   independent weak equations Keff_i. (Default)
%  'sumtraj'      : sum_i BIC_i with Ki     
%  'global'       : single global fixed-effects BIC over all weak equations
if ~isfield(opts, 'selection_score') || isempty(opts.selection_score)
    opts.selection_score = 'sumtraj_keff';
end

if opts.verbose
    fprintf('\n--- Fitting %d candidate growth laws via WENDy ---\n', n_laws);
end

parfor m = 1:n_laws
    warning('off', 'MATLAB:rankDeficientMatrix');
    law       = candidate_laws{m};
    idx_terms = term_map.(law);
    K_M       = numel(idx_terms);

    % Reduced library and per-law bounds
    lib_M = library('terms', lib.terms(idx_terms));
    if opts.want_linargs
        linargs_M = get_linargs(law);
    else
        linargs_M = {};
    end

    WS_wendy_M = cell(N, 1);
    CovWs_M    = cell(N, 1);
    params_M   = nan(K_M, N);
    RSS_M      = 0;
    n_obs_M    = 0;
    n_failed   = 0;


    BIC_sumtraj_M    = 0;
    AIC_sumtraj_M    = 0;
    BICeff_sumtraj_M = 0;

    bic_i_vec = nan(N,1);
    biceff_i_vec = nan(N,1);
    aic_i_vec = nan(N,1);
    rss_i_vec = nan(N,1);
    Ki_vec    = nan(N,1);
    Keff_vec    = nan(N,1);

    for i = 1:N
            % Build WENDy Model
            optm_i  = WS_opt();
            W_wendy = wendy_model(Uobj(i), lib_M, tf_full{i}, [1,0]);

            % Run WENDy
            if opts.want_linargs
                [W_wendy, ~, ~, ~, CovW, ~] = optm_i.wendy(W_wendy, ...
                                                "linregargs", linargs_M);
            else
                [W_wendy, ~, ~, ~, CovW, ~] = optm_i.wendy(W_wendy);
            end
            WS_wendy_M{i}  = W_wendy;
            CovWs_M{i}     = CovW;
            params_M(:, i) = W_wendy.weights;

            % WENDy outputs
            G_i     = W_wendy.G{1};
            b_i     = W_wendy.b{1};
            w_i     = W_wendy.weights;
            C       = W_wendy.cov;

            resid = b_i - G_i * w_i; % Equation error residual 
            
            % resid = W_wendy.res; % Weighted residual
            
            % Individual weak residual sum of squares
            RSS_i = norm(resid)^2;
            Ki    = numel(b_i);     % Num of weak-form equations for this trajectory
            d_m   = K_M;            % Num of fitted parameters for this law per trajectory

            
            if isequal(size(C), [Ki Ki]) && all(isfinite(C(:)))
                % Effective num of weak-form equations
                Keff = trace(C)^2 / max(sum(C(:).^2), eps); 
            else
                Keff = Ki;
            end
            Keff = min(max(Keff, d_m + 1), Ki);  % keep the logs well posed

            % Aggregate global RSS
            RSS_M   = RSS_M   + RSS_i;
            n_obs_M = n_obs_M + Ki;
            
            % Per-trajectory fixed-effects BIC/AIC
            BIC_i = Ki * log(max(RSS_i, eps) / Ki) + d_m * log(Ki);
            AIC_i = Ki * log(max(RSS_i, eps) / Ki) + 2 * d_m;
            % BICeff_i = Keff * log(max(RSS_i, eps) / Keff) + d_m * log(Keff);
            BICeff_i = Ki * log(max(RSS_i, eps) / Ki) + d_m * log(Keff);
            
            BIC_sumtraj_M = BIC_sumtraj_M + BIC_i;
            AIC_sumtraj_M = AIC_sumtraj_M + AIC_i;
            BICeff_sumtraj_M = BICeff_sumtraj_M + BICeff_i;
            
            % Store trajectory-level diagnostics
            rss_i_vec(i) = RSS_i;
            Ki_vec(i)    = Ki;
            Keff_vec(i)  = Keff;
            bic_i_vec(i) = BIC_i;
            aic_i_vec(i) = AIC_i;
            biceff_i_vec(i) = BICeff_i;

    end

    if n_obs_M == 0
        if opts.verbose
            fprintf('  %-14s : all trajectories failed -- skipped\n', law);
        end
        continue
    end

    n_succ   = N - n_failed;
    n_params = K_M * n_succ;
    
    % Score 1: Global fixed-effects weak-residual BIC
    BIC_global_FE = n_obs_M * log(max(RSS_M, eps) / n_obs_M) ...
                  + n_params * log(n_obs_M);
    
    % Score 2: Sum of per-trajectory fixed-effects BICs
    BIC_sumtraj_FE = BIC_sumtraj_M;
    
    % Optional: summed per-trajectory AIC
    AIC_sumtraj_FE = AIC_sumtraj_M;
    
    % Mean weak RSS per weak equation
    RSS_mean = RSS_M / n_obs_M;
    
    % Store scores
    BIC_M = BIC_global_FE;  
    
    bic_per_law(m)             = BIC_M;
    bic_global_fe_per_law(m)   = BIC_global_FE;
    bic_sumtraj_per_law(m)     = BIC_sumtraj_FE;
    aic_sumtraj_per_law(m)     = AIC_sumtraj_FE;
    rss_mean_per_law(m)        = RSS_mean;
    
    rss_per_law(m)   = RSS_M;
    np_per_law(m)    = n_params;
    K_per_law(m)     = K_M;
    n_obs_per_law(m) = n_obs_M;
    
    biceff_sumtraj_per_law(m)  = BICeff_sumtraj_M;

    bic_i_all{m}    = bic_i_vec;
    biceff_i_all{m} = biceff_i_vec;
    rss_i_all{m}    = rss_i_vec;
    Ki_all{m}       = Ki_vec;
    Keff_all{m}     = Keff_vec;

    bic_per_law(m)   = BIC_M;
    rss_per_law(m)   = RSS_M;
    np_per_law(m)    = n_params;
    K_per_law(m)     = K_M;
    n_obs_per_law(m) = n_obs_M;

    wendy_results{m} = struct( ...
        'WS_wendy', {WS_wendy_M}, ...
        'CovWs',    {CovWs_M},    ...
        'params',   params_M,     ...
        'n_obs',    n_obs_M);

    if opts.verbose
        fprintf('  %-14s : K=%d  n_succ=%d  RSS=%.3e  BIC=%.2f\n', ...
            law, K_M, n_succ, RSS_M, BIC_M);
    end
end

%%  Optional: score every law against a COMMON Keff 
if ~isfield(opts, 'keff_reference'), opts.keff_reference = ''; end

if ~isempty(opts.keff_reference)
    ref = find(strcmp(candidate_laws, opts.keff_reference), 1);
    if isempty(ref) || isempty(Keff_all{ref})
        error('opts.keff_reference = "%s" did not produce usable fits.', ...
              opts.keff_reference);
    end
    Keff_ref = max(Keff_all{ref}, max(K_per_law) + 1);
    Keff_ref(~isfinite(Keff_all{ref})) = NaN; 
    for m = 1:n_laws
        if isempty(rss_i_all{m}), continue; end
        rss_m = rss_i_all{m};
        d_m   = K_per_law(m);
        ok_m  = isfinite(rss_m) & isfinite(Keff_ref);
        bi    = Keff_ref(ok_m) .* log(max(rss_m(ok_m), eps) ./ Keff_ref(ok_m)) ...
              + d_m * log(Keff_ref(ok_m));
        biceff_sumtraj_per_law(m) = sum(bi);
        biceff_i_all{m}(ok_m)     = bi;
        Keff_all{m}(ok_m)         = Keff_ref(ok_m);
    end
    fprintf('Keff taken from the "%s" fits for every law (common-Keff mode).\n', ...
            opts.keff_reference);
end

% ---- Select winning model ------------------------------------------------
switch lower(opts.selection_score)
    case 'sumtraj_keff'
        selection_vec  = biceff_sumtraj_per_law;
        selection_name = 'summed per-trajectory BIC (Keff)';
    case 'sumtraj'
        selection_vec  = bic_sumtraj_per_law;
        selection_name = 'summed per-trajectory BIC (Ki)';
    case 'global'
        selection_vec  = bic_global_fe_per_law;
        selection_name = 'global fixed-effects BIC';
    otherwise
        error(['opts.selection_score must be ''sumtraj_keff'', ''sumtraj'' ' ...
               'or ''global''; got ''%s''.'], opts.selection_score);
end
[~, best_idx] = min(selection_vec);
if ~isfinite(selection_vec(best_idx))
    error('Model selection failed: no candidate law produced a finite BIC.');
end

model_id   = candidate_laws{best_idx};
WS_wendy   = wendy_results{best_idx}.WS_wendy;
CovWs      = wendy_results{best_idx}.CovWs;
params     = wendy_results{best_idx}.params;
p_lib      = K_per_law(best_idx);
n_obs_best = n_obs_per_law(best_idx);

% Mark and drop trajectories that failed under the winning law
valid_traj    = ~any(isnan(params), 1);
n_failed_best = sum(~valid_traj);

% ---- BIC table (sorted) ---------
if strcmpi(opts.selection_score, 'sumtraj_keff') && ~isempty(Keff_all{best_idx})
    n_units_best = sum(Keff_all{best_idx}, 'omitnan');
else
    n_units_best = n_obs_best;
end
delta_bic      = selection_vec - selection_vec(best_idx);
norm_delta_bic = delta_bic / max(n_units_best, eps);
[~, sort_idx]  = sort(selection_vec, 'ascend');

bic_table = table( ...
    candidate_laws(sort_idx)',   ...
    selection_vec(sort_idx),     ...
    norm_delta_bic(sort_idx),    ...
    rss_per_law(sort_idx),       ...
    np_per_law(sort_idx),        ...
    K_per_law(sort_idx),         ...
    'VariableNames', {'Law','BIC','NormDeltaBIC','RSS','n_params','K'});

if opts.verbose
    fprintf('\n');
    disp(bic_table)
    fprintf('Selected model: %s   (failed traj under winner: %d / %d)\n\n', ...
        model_id, n_failed_best, N);
end

% ---- Compare fixed-effects model-selection scores ---------------------

[~, best_global_idx]  = min(bic_global_fe_per_law);
[~, best_sumtraj_idx] = min(bic_sumtraj_per_law);
[~, best_biceff_idx]  = min(biceff_sumtraj_per_law);
[~, best_aic_idx]     = min(aic_sumtraj_per_law);

delta_global  = bic_global_fe_per_law  - min(bic_global_fe_per_law);
delta_sumtraj = bic_sumtraj_per_law    - min(bic_sumtraj_per_law);
delta_biceff  = biceff_sumtraj_per_law - min(biceff_sumtraj_per_law);
delta_aic     = aic_sumtraj_per_law    - min(aic_sumtraj_per_law);

% Median Ki / Keff per law, for the table and for the sanity check below
Ki_med_per_law   = nan(n_laws, 1);
Keff_med_per_law = nan(n_laws, 1);
for m = 1:n_laws
    if isempty(Ki_all{m}), continue; end
    Ki_med_per_law(m)   = median(Ki_all{m},   'omitnan');
    Keff_med_per_law(m) = median(Keff_all{m}, 'omitnan');
end

score_table = table( ...
    candidate_laws(:), ...
    K_per_law(:), ...
    n_obs_per_law(:), ...
    np_per_law(:), ...
    rss_per_law(:), ...
    rss_mean_per_law(:), ...
    Ki_med_per_law(:), ...
    Keff_med_per_law(:), ...
    bic_global_fe_per_law(:), ...
    delta_global(:), ...
    bic_sumtraj_per_law(:), ...
    delta_sumtraj(:), ...
    biceff_sumtraj_per_law(:), ...
    delta_biceff(:), ...
    aic_sumtraj_per_law(:), ...
    delta_aic(:), ...
    'VariableNames', { ...
        'Law', ...
        'K', ...
        'n_obs', ...
        'n_params_global', ...
        'RSS', ...
        'RSS_per_obs', ...
        'Ki_med', ...
        'Keff_med', ...
        'BIC_global_FE', ...
        'Delta_BIC_global_FE', ...
        'BIC_sumtraj_FE', ...
        'Delta_BIC_sumtraj_FE', ...
        'BIC_sumtraj_Keff', ...
        'Delta_BIC_sumtraj_Keff', ...
        'AIC_sumtraj_FE', ...
        'Delta_AIC_sumtraj_FE'});

% Sort by whichever score is actually driving the selection
switch lower(opts.selection_score)
    case 'sumtraj_keff', sort_by = 'BIC_sumtraj_Keff';
    case 'sumtraj',      sort_by = 'BIC_sumtraj_FE';
    otherwise,           sort_by = 'BIC_global_FE';
end
score_table = sortrows(score_table, sort_by, 'ascend');

fprintf('\n=== Fixed-effects WENDy model-selection diagnostics ===\n');
disp(score_table)

fprintf('Winner by global fixed-effects BIC:       %s\n', ...
    candidate_laws{best_global_idx});
fprintf('Winner by summed trajectory BIC (Ki):     %s\n', ...
    candidate_laws{best_sumtraj_idx});
fprintf('Winner by summed trajectory BIC (Keff):   %s\n', ...
    candidate_laws{best_biceff_idx});
fprintf('Winner by summed trajectory AIC:          %s\n', ...
    candidate_laws{best_aic_idx});
fprintf('Selection in force: %s  ->  %s\n', selection_name, model_id);
fprintf('======================================================\n\n');


infl = Ki_med_per_law ./ max(Keff_med_per_law, eps);

% fprintf('=== Ki vs Keff ===\n');
% fprintf('%-14s  %8s  %9s  %11s\n', 'Law', 'med Ki', 'med Keff', 'inflation');
% for m = 1:n_laws
%     if isnan(Ki_med_per_law(m)), continue; end
%     fprintf('%-14s  %8.0f  %9.1f  %10.1fx\n', ...
%         candidate_laws{m}, Ki_med_per_law(m), Keff_med_per_law(m), infl(m));
% end

% if all(abs(infl(isfinite(infl)) - 1) < 1e-6)
%     fprintf(['\n  All inflation factors are 1: W_wendy.cov is not being\n' ...
%              '  returned as a Ki x Ki matrix, so Keff has fallen back to Ki\n' ...
%              '  and the two scores are identical by construction.  Check the\n' ...
%              '  covariance field before concluding the correction has no effect.\n']);
% end

% if best_sumtraj_idx ~= best_biceff_idx
%     a = best_biceff_idx;   
%     b = best_sumtraj_idx;  
%     dA = bic_i_all{b}    - bic_i_all{a};      
%     dB = biceff_i_all{b} - biceff_i_all{a};
%     ok = isfinite(dA) & isfinite(dB);
%     fprintf('\n  The two scores disagree: %s (Ki) vs %s (Keff).\n', ...
%         candidate_laws{b}, candidate_laws{a});
%     fprintf('  Trajectories individually preferring %s: %.1f%% (Ki), %.1f%% (Keff)\n', ...
%         candidate_laws{b}, 100*mean(dA(ok) < 0), 100*mean(dB(ok) < 0));
%     fprintf('  Summed Delta (%s - %s): %+.1f (Ki), %+.1f (Keff)\n', ...
%         candidate_laws{b}, candidate_laws{a}, sum(dA(ok)), sum(dB(ok)));
%     fprintf('  Per-trajectory Delta (Keff) quantiles [5 25 50 75 95]:');
%     fprintf(' %+.2f', quantile(dB(ok), [0.05 0.25 0.50 0.75 0.95]));
%     fprintf('\n');
% end
% fprintf('==================\n\n');

%% Convert to biological parameters
params_converted_struct = weights_to_biological_params(params', model_id);



%% ==================== Empirical Bayes Shrinkage ====================
%  Assumes biological parameters (r, S_max) are normally distributed
%  across the population.  Shrinkage is performed in biological-parameter
%  space using the delta-method covariance propagated from the raw GLS
%  coefficient covariance.
%  ===================================================================

if opts.want_shrinkage

% ---- 1. Raw GLS coefficient covariance (full matrix per indv) ----
CovW_gls_full = cell(N,1);
for i = 1:N
    G_i = WS_wendy{i}.G{1};
    C_i = WS_wendy{i}.cov;
    CovW_gls_full{i} = pinv(G_i' * (C_i \ G_i));
end

% ---- 2. Select parameters to shrink: r always, S_max if finite ----
valid = params_converted_struct.valid;
N_valid = sum(valid);
idx_valid = find(valid);

r_all = params_converted_struct.r;

shrink_names = {'r'};
if isfield(params_converted_struct, 'S_max')
    Smax_vals = params_converted_struct.S_max(valid);
    if all(isfinite(Smax_vals))
        shrink_names{end+1} = 'S_max';
    else
        fprintf('  Skipping S_max shrinkage (non-finite values present; e.g. exponential growth)\n');
    end
end
n_shrink = numel(shrink_names);

% ---- 3. Assemble biological parameter matrix (valid cells only) ----
bio_params = zeros(N_valid, n_shrink);
for f = 1:n_shrink
    bio_params(:,f) = params_converted_struct.(shrink_names{f})(valid);
end

% ---- 4. Propagate covariance via delta method ----
%  Cov(theta_bio) = J * Cov(w) * J'
%  J_{fk} = d(theta_bio_f) / d(w_k)   (central differences)
bio_S_hat = zeros(N_valid, n_shrink);
h_scale = 1e-7;

for ii = 1:N_valid
    i = idx_valid(ii);
    w_i = params(:, i);
    S_full_i = CovW_gls_full{i};

    J = zeros(n_shrink, p_lib);
    for k = 1:p_lib
        h = h_scale * max(abs(w_i(k)), 1e-8);
        w_fwd = w_i; w_fwd(k) = w_fwd(k) + h;
        w_bwd = w_i; w_bwd(k) = w_bwd(k) - h;
        s_fwd = weights_to_biological_params(w_fwd', model_id);
        s_bwd = weights_to_biological_params(w_bwd', model_id);
        for f = 1:n_shrink
            J(f,k) = (s_fwd.(shrink_names{f}) - s_bwd.(shrink_names{f})) / (2*h);
        end
    end

    CovBio = J * S_full_i * J';
    bio_S_hat(ii,:) = max(diag(CovBio)', 0);
end

% Cap non-finite or extreme covariances
for f = 1:n_shrink
    emp_var_f = var(bio_params(:,f));
    cap = 100 * emp_var_f;
    bad = ~isfinite(bio_S_hat(:,f)) | bio_S_hat(:,f) > cap;
    if any(bad)
        fprintf('  Capped %d extreme S_hat entries for "%s" at %.2e\n', ...
            sum(bad), shrink_names{f}, cap);
        bio_S_hat(bad,f) = cap;
    end
end

% ---- 5. REML estimation + BLUP shrinkage (inline, per parameter) ----
mu_hat   = zeros(1, n_shrink);
tau2_hat = zeros(1, n_shrink);
gamma    = zeros(N_valid, n_shrink);
bio_params_shrunk = zeros(N_valid, n_shrink);
V_post   = zeros(N_valid, n_shrink);

for f = 1:n_shrink
    theta = bio_params(:,f);
    s2    = bio_S_hat(:,f);

    % -- Initial tau^2 from method of moments --
    %    Var(theta_hat) = tau^2 + E[s_n^2]
    %    => tau^2 ~ Var(theta_hat) - mean(s_n^2)
    tau2_init = max(var(theta) - mean(s2), var(theta)/4);

    % -- REML profile log-likelihood for tau^2 --
    %    l_R(tau2) = -0.5 * [ sum log(V_n) + log(sum 1/V_n)
    %                        + sum (theta_n - mu(tau2))^2 / V_n ]
    %    where V_n = tau2 + s_n^2,  mu(tau2) = precision-weighted mean
    neg_reml = @(tau2) reml_neg_profile(tau2, theta, s2);

    % Search bounds
    tau2_lo  = 0;
    tau2_hi  = max(10 * var(theta), 10 * max(s2));
    if ~isfinite(tau2_hi) || tau2_hi <= 0
        tau2_hi = 1;
    end

    % Optimize
    [tau2_opt, fval_opt] = fminbnd(neg_reml, tau2_lo, tau2_hi, ...
        optimset('TolX', 1e-12, 'MaxIter', 500, 'Display', 'off'));

    % Check boundary
    if neg_reml(0) <= fval_opt
        tau2_opt = 0;
    end

    % -- Population mean at optimal tau^2 --
    V_f   = max(tau2_opt + s2,eps);
    om_f  = 1 ./ max(V_f, eps);
    mu_f  = (om_f' * theta) / sum(om_f);

    % -- Shrinkage gains and BLUP --
    gam_f = tau2_opt ./ max(V_f, eps);
    theta_shrunk = mu_f + gam_f .* (theta - mu_f);
    vpost_f = gam_f .* s2;

    % Store
    mu_hat(f)   = mu_f;
    tau2_hat(f) = tau2_opt;
    gamma(:,f)  = gam_f;
    bio_params_shrunk(:,f) = theta_shrunk;
    V_post(:,f) = vpost_f;

    % -- Diagnostic --
    if max(gam_f) < 0.05
        warning('shrinkage:collapsed', ...
            ['Parameter "%s": all shrinkage gains < 0.05 (max = %.4f).\n' ...
             'All estimates collapsed to population mean %.4e.\n' ...
             'Likely cause: measurement covariance too large relative to\n' ...
             'true heterogeneity. Check delta-method Jacobian and CovW_gls.'], ...
            shrink_names{f}, max(gam_f), mu_f);
    end
end

% ---- 6. Pack shrunk biological parameters back into struct ----
params_bio_shrunk_struct = params_converted_struct;   % copy metadata + unshrunk fields
for f = 1:n_shrink
    vals = params_bio_shrunk_struct.(shrink_names{f});
    vals(valid) = bio_params_shrunk(:,f);
    params_bio_shrunk_struct.(shrink_names{f}) = vals;
end

% ---- Summary ----
fprintf('\n=== Empirical Bayes Shrinkage Summary ===\n');
fprintf('Model: %s  (%s)\n', params_converted_struct.model_tag, ...
    params_converted_struct.model_rhs_form);
fprintf('N = %d valid indivs\n', N_valid);
for f = 1:n_shrink
    fn = shrink_names{f};
    fprintf('\n  %s:\n', fn);
    fprintf('    mu_hat          = %10.4e\n', mu_hat(f));
    fprintf('    tau_hat (std)    = %10.4e\n', sqrt(tau2_hat(f)));
    fprintf('    naive std       = %10.4e\n', std(bio_params(:,f)));
    fprintf('    mean(s_n^2)     = %10.4e\n', mean(bio_S_hat(:,f)));
    fprintf('    gamma: min=%.3f  med=%.3f  max=%.3f\n', ...
        min(gamma(:,f)), median(gamma(:,f)), max(gamma(:,f)));
end
fprintf('=========================================\n');
params_shrunk_w = bio_to_weights(params_bio_shrunk_struct, model_id, N);

end % want_shrinkage

%% ======== Population Sims ============================================

% ---- Options -----------------------------------------------------------
if ~isfield(opts, 'T_ref_quantile') || isempty(opts.T_ref_quantile)
    opts.T_ref_quantile = 0.25;
end
opts.cluster_space = 'weights';
% Clustering space: 'weights' (raw ODE weights) or 'bio' (r, S_max)
if ~isfield(opts, 'cluster_space')
    opts.cluster_space = 'bio';   % 'bio' | 'weights'
end
if ~isfield(opts, 'M_cohorts')
    opts.M_cohorts = [];          % [] = auto (sqrt heuristic)
end

T_ref_q = opts.T_ref_quantile;

% ---- Canonical parameter matrix: always (p_lib x N), raw ODE weights --
% Both branches produce the same shape so all downstream indexing is
% identical regardless of whether shrinkage is on or off.
if opts.want_shrinkage
    params_full = bio_to_weights(params_bio_shrunk_struct, model_id, N)';
    params_tag  = 'shrunk (empirical Bayes)';
else
    params_full = params;          
    params_tag  = 'raw WENDy';
end

% Tracking old inconsistencies with matrix sizes
assert(size(params_full, 1) == p_lib && size(params_full, 2) == N, ...
    'params_full shape mismatch: expected (%d x %d), got (%d x %d)', ...
    p_lib, N, size(params_full,1), size(params_full,2));

N_full     = N;   % numel(xcell) after subsampling
% Validity: finite params AND the trajectory itself was not a WENDy failure
good_param = all(isfinite(params_full), 1) & valid_traj(:)';  % (1 x N)

% ---- Age-since-birth axis ---------------------------------------------
age_cell   = cellfun(@(t) t - t(1), tcell, 'UniformOutput', false);
final_ages = cellfun(@(a) a(end),   age_cell);
T_ref      = quantile(final_ages, T_ref_q);

if T_ref <= 0 || ~isfinite(T_ref)
    warning('Population-sim: degenerate age axis (T_ref = %.3g); skipping.', T_ref);
    return
end

% ---- Survival mask + empirical sizes at T_ref -------------------------
survives       = false(N_full, 1);
sizes_ref_full = nan(N_full, 1);
for i = 1:N_full
    if age_cell{i}(1) <= T_ref && age_cell{i}(end) >= T_ref
        sizes_ref_full(i) = interp1(age_cell{i}, xcell{i}, T_ref, 'linear');
        survives(i)       = true;
    end
end

% ---- Joint mask -------------------------------------------------------
keep       = good_param(:) & survives(:);   % (N x 1)
N_pop      = sum(keep);
idx_keep   = find(keep);

if N_pop < 10
    warning(['Population-sim: only %d individuals survive to T_ref = %.3g. ' ...
             'Reduce opts.T_ref_quantile.'], N_pop, T_ref);
end

params_pop = params_full(:, keep);          % (p_lib x N_pop)
sizes_emp  = sizes_ref_full(keep);          % (N_pop x 1)

% ---- Fit x0 under the ACTUAL params used for simulation ---------------
% Use the shrunk/raw weights (whichever is in params_pop) so that the
% initial condition is consistent with the RHS we integrate.
x0_vec = nan(N_pop, 1);
for ii = 1:N_pop
    i     = idx_keep(ii);
    w_i   = params_pop(:, ii);
    rhs_i = @(x) growth_rhs(x, w_i', model_id);
    tg_i  = Uobj(i).grid{1}(:);
    y_i   = Uobj(i).Uobs{1}(:);
    x0_vec(ii) = fit_x0([], rhs_i, tg_i, y_i, y_i(1));
    x0_vec(ii) = max(x0_vec(ii),1e-4);
end

% ---- Clustering space -------------------------------------------------
% Build the feature matrix used for k-means.  Two options:
%   'weights' : raw ODE weights  (p_lib x N_pop) -> cluster rows of params_pop'
%   'bio'     : biological params (r, S_max)      -> more interpretable centroids
switch lower(opts.cluster_space)

    case 'weights'
        cluster_features = params_pop';          % (N_pop x p_lib)
        feature_names    = arrayfun(@(k) sprintf('w%d', k), ...
                               (1:p_lib)', 'UniformOutput', false);

    case 'bio'
        % Convert each individual's ODE weights to (r, S_max).
        % weights_to_biological_params expects rows = individuals.
        bio_pop = weights_to_biological_params(params_pop', model_id);

        fn_list = {'r'};
        if isfield(bio_pop, 'S_max') && all(isfinite(bio_pop.S_max))
            fn_list{end+1} = 'S_max';
        end
        cluster_features = zeros(N_pop, numel(fn_list));
        for f = 1:numel(fn_list)
            cluster_features(:, f) = bio_pop.(fn_list{f});
        end
        feature_names = fn_list;

    otherwise
        error('opts.cluster_space must be ''weights'' or ''bio''; got ''%s''.', ...
              opts.cluster_space);
end

% ---- K-means ----------------------------------------------------------
if isempty(opts.M_cohorts)
    M_cohorts = min(20, max(2, round(sqrt(N_pop))));
else
    M_cohorts = min(opts.M_cohorts, N_pop);
end

rng(1);
[cluster_idx, centroids_feat] = kmeans(cluster_features, M_cohorts, ...
    'Replicates', 10, 'Display', 'off');

% Convert centroids back to ODE weights if we clustered in bio space
if strcmpi(opts.cluster_space, 'bio')
    % Rebuild a bio struct from the centroid feature values
    centroid_bio_struct = params_bio_shrunk_struct;  % copy for field metadata
    for f = 1:numel(fn_list)
        centroid_bio_struct.(fn_list{f}) = centroids_feat(:, f);
    end
    % bio_to_weights returns (M_cohorts x p_lib)
    centroids_w = bio_to_weights(centroid_bio_struct, model_id, M_cohorts);
else
    centroids_w = centroids_feat;   % already ODE weights (M_cohorts x p_lib)
end

% Cohort sizes (diagnostic)
cohort_sizes = histcounts(cluster_idx, 1:M_cohorts+1);
fprintf('K-means cohort sizes [%s space]: min=%d  max=%d  mean=%.1f\n', ...
    opts.cluster_space, min(cohort_sizes), max(cohort_sizes), mean(cohort_sizes));

% ---- k-means diagnostic plot -----------------------------------------
n_feat = size(cluster_features, 2);

if n_feat == 1
    Z_plot      = [cluster_features(:,1), 0.02*randn(N_pop,1)];
    Z_cent_plot = [centroids_feat(:,1),   zeros(M_cohorts,1)];
    xlab = feature_names{1};  ylab = 'jitter';
elseif n_feat == 2
    Z_plot      = cluster_features;
    Z_cent_plot = centroids_feat;
    xlab = feature_names{1};  ylab = feature_names{2};
else
    [coeff, score, ~, ~, explained, mu_pca] = pca(cluster_features);
    Z_plot      = score(:,1:2);
    Z_cent_plot = (centroids_feat - mu_pca) * coeff(:,1:2);
    xlab = sprintf('PC1 (%.1f%%)', explained(1));
    ylab = sprintf('PC2 (%.1f%%)', explained(2));
end

figure('Name', sprintf('K-means clusters [%s space]', opts.cluster_space));
hold on; box on; grid on;
scatter(Z_plot(:,1), Z_plot(:,2), 36, cluster_idx, 'filled', ...
    'MarkerFaceAlpha', 0.65);
colormap(lines(M_cohorts));
scatter(Z_cent_plot(:,1), Z_cent_plot(:,2), 144, (1:M_cohorts)', 'x', ...
    'LineWidth', 2.5);
xlabel(xlab); ylabel(ylab);
title(sprintf('K-means  (%s, M=%d, space=%s)', model_id, M_cohorts, opts.cluster_space));

% ---- Size grid for KDE ------------------------------------------------
all_sizes = [x0_vec; sizes_emp];
x_lo      = 0.9 * min(all_sizes(isfinite(all_sizes)));
x_hi      = 1.1 * max(all_sizes(isfinite(all_sizes)));
x_grid    = linspace(x_lo, x_hi, 600)';
bw        = 0.04 * (x_hi - x_lo);

rho_init    = ksdensity(x0_vec,    x_grid, 'Bandwidth', bw);
rho_emp_ref = ksdensity(sizes_emp, x_grid, 'Bandwidth', bw);

% ---- (1) Mixture closure: per-individual under their cohort centroid ---
x_t_mix = zeros(N_pop, 1);
for m = 1:M_cohorts
    members = find(cluster_idx == m);
    if isempty(members), continue; end
    x_t_mix(members) = advect_population( ...
        x0_vec(members), T_ref, centroids_w(m, :), model_id);
end
rho_mix = ksdensity(x_t_mix, x_grid, 'Bandwidth', bw);
rho_mix = rho_mix(:);

% ---- (2a) Mean-weight PDE ---------------------------------------------
w_mean  = mean(params_pop, 2)';           % (1 x p_lib)
x_t_mw  = advect_population(x0_vec, T_ref, w_mean, model_id);
rho_mw  = ksdensity(x_t_mw, x_grid, 'Bandwidth', bw);
rho_mw  = rho_mw(:);

% ---- (2b) Empirical g_eff PDE -----------------------------------------
g_eff   = @(x) eval_g_eff_empirical(x, params_pop, model_id);
x_t_ge  = advect_with_geff(x0_vec, T_ref, g_eff);
rho_ge  = ksdensity(x_t_ge, x_grid, 'Bandwidth', bw);
rho_ge  = rho_ge(:);

% ---- Per-individual oracle prediction ---------------------------------
x_t_indiv = zeros(N_pop, 1);
for j = 1:N_pop
    x_t_indiv(j) = advect_population( ...
        x0_vec(j), T_ref, params_pop(:, j)', model_id);
end
rho_indiv   = ksdensity(x_t_indiv, x_grid, 'Bandwidth', bw);
rho_indiv   = rho_indiv(:);
indiv_resid = x_t_indiv - sizes_emp;
indiv_rrmse = sqrt(mean(indiv_resid.^2)) / max(std(sizes_emp), eps);

% ---- Density plot -----------------------------------------------------
figure('Name', 'Population density evolution', 'Position', [100 100 900 480]);
hold on; box on; grid on;
plot(x_grid, rho_init,    'k-',  'LineWidth', 1.8, ...
    'DisplayName', 'initial (age 0, survivors)');
plot(x_grid, rho_emp_ref, 'r-',  'LineWidth', 2.0, ...
    'DisplayName', sprintf('empirical at age %.2g', T_ref));
plot(x_grid, rho_indiv,   ':',   'LineWidth', 1.8, 'Color', [0.4 0.4 0.4], ...
    'DisplayName', 'per-individual ODE (oracle)');
plot(x_grid, rho_mix,     '--',  'LineWidth', 1.7, 'Color', [0.85 0.65 0.13], ...
    'DisplayName', sprintf('mixture closure (M=%d, %s)', M_cohorts, opts.cluster_space));
plot(x_grid, rho_mw,      ':',   'LineWidth', 2.0, 'Color', [0.50 0.00 0.50], ...
    'DisplayName', 'mean-weight PDE');
plot(x_grid, rho_ge,      '-.',  'LineWidth', 1.7, 'Color', [0.10 0.55 0.30], ...
    'DisplayName', 'empirical g_{eff} PDE');
xlabel('size  x'); ylabel('density  \rho(x)');
title(sprintf('Population density at age T=%.3g  [%s, params: %s, N=%d]', ...
    T_ref, model_id, params_tag, N_pop));
legend('Location', 'best', 'FontSize', 9);

% ---- Quantitative comparison ------------------------------------------
L1 = @(a, b) trapz(x_grid, abs(a(:) - b(:)));
cdf_emp = cumtrapz(x_grid, rho_emp_ref);
cdf_emp = cdf_emp / cdf_emp(end);
W1 = @(rho) trapz(x_grid, abs( ...
    cumtrapz(x_grid, rho(:)) / max(trapz(x_grid, rho(:)), eps) - cdf_emp));

fprintf('\n=== Population-level density at age T=%.3g (q=%.2f) ===\n', T_ref, T_ref_q);
fprintf('  Clustering space : %s\n', opts.cluster_space);
fprintf('  N used           : %d survivors with valid params / %d total\n', N_pop, N_full);
fprintf('  per-individual RRMSE on sizes(T_ref): %.4f\n', indiv_rrmse);
fprintf('                              L1         W1\n');
fprintf('  per-individual ODE  :  %.4f     %.4f\n', ...
    L1(rho_indiv, rho_emp_ref), W1(rho_indiv));
fprintf('  mixture (M=%2d, %s):  %.4f     %.4f\n', ...
    M_cohorts, opts.cluster_space, L1(rho_mix, rho_emp_ref), W1(rho_mix));
fprintf('  mean-weight PDE     :  %.4f     %.4f\n', ...
    L1(rho_mw,  rho_emp_ref), W1(rho_mw));
fprintf('  empirical g_eff PDE :  %.4f     %.4f\n', ...
    L1(rho_ge,  rho_emp_ref), W1(rho_ge));
fprintf('======================================================\n\n'); 


%% ==================== Plotting ====================
if opts.want_plots
    WENDy_color = [0.2 0.4 0.8];
    Shrunk_color = [0.9 0.3 0.2];
    True_color = [0.3 0.8 0.3];

try
    true_weights = true_weights(rand_idx,:);
    growth_type  = Data.growth.type;
catch
    true_weights = [];
end

ode_opts_plot = odeset('RelTol',1e-10,'AbsTol',1e-10);
has_true = ~isempty(true_weights) && ~all(isnan(true_weights(:)));

bio_raw    = params_converted_struct;
if opts.want_shrinkage
bio_shrunk = params_bio_shrunk_struct;
else
    bio_shrunk = [];
end

bio_true   = [];
if has_true    
    bio_true = weights_to_biological_params(true_weights, Data.growth.type);
end

% Fields to plot: only r and S_max (if finite)
plot_fields = {'r'};
if isfield(bio_raw, 'S_max') && any(isfinite(bio_raw.S_max))
    plot_fields{end+1} = 'S_max';
end
n_plot = numel(plot_fields);

% --- RRMSE ---
err_raw    = nan(N,1);
err_shrunk = nan(N,1);
for i = 1:N
    try
        x0 = Uobj(i).get_x0([]);
        tg = Uobj(i).grid{1}(:);
        y  = Uobj(i).Uobs{1}(:);
        ynorm = max(norm(y), eps);

        rhs_raw_i = WS_wendy{i}.get_rhs;
        [~,xr] = ode45(@(t,x) rhs_raw_i(x), tg, x0, ode_opts_plot);
        err_raw(i) = norm(xr(:,1) - y) / ynorm;

        % [~,xs] = ode45(@(t,x) bio_rhs(x, bio_shrunk, model_id, i), ...
        %                 tg, x0, ode_opts_plot);
        % err_shrunk(i) = norm(xs(:,1) - y) / ynorm;
    catch
    end
end

% --- E2 ---
if opts.want_shrinkage
    params_shrunk_w = bio_to_weights(bio_shrunk, model_id, N);
end
if has_true

    e2_raw    = vecnorm(params' - true_weights, 2, 2) ...
              ./ max(vecnorm(true_weights, 2, 2), eps);

    if opts.want_shrinkage
    e2_shrunk = vecnorm(params_shrunk_w - true_weights, 2, 2) ...
              ./ max(vecnorm(true_weights, 2, 2), eps);
    end
end

% ================================================================
%  FIGURE 1: Biological parameter histograms
% ================================================================
figure('Name','Biological parameters');
for f = 1:n_plot
    all_vals =[];
    subplot(1, n_plot, f); hold on; box on;
    fn = plot_fields{f};

    vr = bio_raw.(fn);    vr = vr(isfinite(vr));
    if opts.want_shrinkage
        vs = bio_shrunk.(fn); vs = vs(isfinite(vs));
        all_vals = [vr; vs];

    else
        all_vals = vr;
    end

    if isempty(vr), title(strrep(fn,'_','\_')); continue; end


    vt = [];
    if (~isempty(bio_true) && isfield(bio_true, fn)) && has_true
        vt = bio_true.(fn); vt = vt(isfinite(vt));
        all_vals = [all_vals; vt];
    end
    lo = quantile(all_vals, 0.005);
    hi = quantile(all_vals, 0.995);
    if lo == hi, lo = min(all_vals)-eps; hi = max(all_vals)+eps; end
    edges = linspace(lo, hi, 35);



    histogram(vr, edges', 'Normalization','pdf', ...
        'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.8, ...
        'EdgeColor','none', 'DisplayName','WENDy');

    if opts.want_shrinkage
        collapsed = (std(vs) / max(std(vr),eps)) < 0.01;
    if ~collapsed
        histogram(vs, edges, 'Normalization','pdf', ...
            'FaceColor',[0.9 0.3 0.2], 'FaceAlpha',0.4, ...
            'EdgeColor','none', 'DisplayName','Shrunk');
    else
        xline(mean(vs), 'r--', 'LineWidth',2, ...
            'DisplayName', sprintf('Shrunk (%.3g)', mean(vs)));
    end
    end
    if ~isempty(vt)
        histogram(vt, edges, 'Normalization','pdf', ...
            'FaceColor',[0.3 0.8 0.3], 'FaceAlpha',0.35, ...
            'EdgeColor','none', 'DisplayName','True');
    end

    xlabel(strrep(fn,'_','\_')); ylabel('pdf');
    legend('Location','best','FontSize',8);
    title(strrep(fn,'_','\_'));
end
sgtitle(sprintf('Biological parameters  (%s)', bio_raw.model_tag));

% %% histfit
% figure('Name','Biological parameters');
% for f = 1:n_plot
%     all_vals =[];
%     subplot(1, n_plot, f); hold on; box on;
%     fn = plot_fields{f};
% 
%     vr = bio_raw.(fn);    vr = vr(isfinite(vr));
%     if opts.want_shrinkage
%         vs = bio_shrunk.(fn); vs = vs(isfinite(vs));
%         all_vals = [vr; vs];
%     else
%         all_vals = vr;
%     end
%     if isempty(vr), title(strrep(fn,'_','\_')); continue; end
% 
% 
%     vt = [];
%     if ~isempty(bio_true) && isfield(bio_true, fn)
%         vt = bio_true.(fn); vt = vt(isfinite(vt));
%         all_vals = [all_vals; vt];
%     end
%     lo = quantile(all_vals, 0.005);
%     hi = quantile(all_vals, 0.995);
%     if lo == hi, lo = min(all_vals)-eps; hi = max(all_vals)+eps; end
%     edges = linspace(lo, hi, 35);
% 
% 
%     bins = 50;
%     dist = "normal";
%     h1 = histfit(vr,bins,dist);
%     delete(h1(1))
%     set(h1(2),'Color',WENDy_color);
%     lgnd = "WENDy (var = "+num2str(var(vr))+")";
%     if opts.want_shrinkage
%         collapsed = (std(vs) / max(std(vr),eps)) < 0.01;
%         if ~collapsed
%             h2 = histfit(vs,bins,dist);
%             delete(h2(1))
%             lgnd = [lgnd, "Shrunk (var = "+num2str(var(vs))+")"];
%             set(h2(2),'Color',Shrunk_color);
%         else    
%         end
%     end
%     if ~isempty(vt)
%         h3 = histfit(vt,bins,dist);
%         delete(h3(1))
%         lgnd = [lgnd, "True (var = "+num2str(var(vt))+")"];
%         set(h3(2),'Color',True_color);
%     end
%     xlabel(strrep(fn,'_','\_')); ylabel('pdf');
%     legend(lgnd,'Location','best','FontSize',8);
%     title(strrep(fn,'_','\_'));
% end
% sgtitle(sprintf('Biological parameters  (%s)', bio_raw.model_tag));

%% fitted PDFs
figure('Name','Biological parameters');

for f = 1:n_plot
    subplot(1, n_plot, f); hold on; box on;
    fn = plot_fields{f};

    vr = bio_raw.(fn);
    vr = vr(isfinite(vr));

    all_vals = vr;

    if opts.want_shrinkage
        vs = bio_shrunk.(fn);
        vs = vs(isfinite(vs));
        all_vals = [all_vals; vs];
    else
        vs = [];
    end

    vt = [];
    if ~isempty(bio_true) && isfield(bio_true, fn)
        vt = bio_true.(fn);
        vt = vt(isfinite(vt));
        all_vals = [all_vals; vt];
    end

    if isempty(all_vals)
        title(strrep(fn,'_','\_'));
        continue;
    end

    lo = quantile(all_vals, 0.005);
    hi = quantile(all_vals, 0.995);

    if lo == hi
        lo = min(all_vals) - eps;
        hi = max(all_vals) + eps;
    end

    x = linspace(lo, hi, 500);

    lines = gobjects(0);
    lgnd = strings(0);

    % WENDy
    if ~isempty(vr) && std(vr) > 0
        mu = mean(vr);
        sig = std(vr);
        y = normpdf(x, mu, sig);

        lines(end+1) = plot(x, y, 'LineWidth', 1.5, 'Color', WENDy_color);
        lgnd(end+1) = "WENDy (var = " + num2str(var(vr), '%.3g') + ")";
    end

    % Shrunk
    if opts.want_shrinkage && ~isempty(vs) && std(vs) > 0
        collapsed = (std(vs) / max(std(vr), eps)) < 0.01;

        if ~collapsed
            mu = mean(vs);
            sig = std(vs);
            y = normpdf(x, mu, sig);

            lines(end+1) = plot(x, y, 'LineWidth', 1.5, 'Color', Shrunk_color);
            lgnd(end+1) = "Shrunk (var = " + num2str(var(vs), '%.3g') + ")";
        end
    end

    % True
    if ~isempty(vt) && std(vt) > 0
        mu = mean(vt);
        sig = std(vt);
        y = normpdf(x, mu, sig);

        lines(end+1) = plot(x, y, 'LineWidth', 1.5, 'Color', True_color);
        lgnd(end+1) = "True (var = " + num2str(var(vt), '%.3g') + ")";
    end

    xlabel(strrep(fn,'_','\_'));
    ylabel('pdf');
    title(strrep(fn,'_','\_'));
    xlim([lo hi]);

    if ~isempty(lines)
        legend(lines, lgnd, 'Location', 'best', 'FontSize', 8);
    end
end

sgtitle(sprintf('Biological parameters  (%s)', bio_raw.model_tag));
%% histogram the weights
figure('Name','Weights');
nW = length(WS_wendy{1}.weights);

for i = 1:nW
    subplot(1,nW,i); hold on; box on;
    h1 = histogram(params(i,:),'Normalization','pdf', ...
        'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.8, ...
        'EdgeColor','none', 'DisplayName','WENDy');

    if opts.want_shrinkage
    histogram(params_shrunk_w(:,i),'BinWidth',h1.BinWidth, 'Normalization','pdf', ...
            'FaceColor',[0.9 0.3 0.2], 'FaceAlpha',0.4, ...
            'EdgeColor','none', 'DisplayName','Shrunk');
    end

    if has_true
    histogram(true_weights(:,i), 'BinWidth',h1.BinWidth, 'Normalization','pdf', ...
            'FaceColor',[0.3 0.8 0.3], 'FaceAlpha',0.35, ...
            'EdgeColor','none', 'DisplayName','True');
    end
    xlabel("W_"+num2str(i)); ylabel('pdf')
end
%% ================================================================
%  FIGURE 2: Trajectory comparisons (5th, 50th, 95th pctl RRMSE)
% ================================================================
opts.want_shrinkage = false;
valid_err = find(isfinite(err_raw));
if numel(valid_err) >= 3
    [~, ord] = sort(err_raw(valid_err));
    pct_idx  = [max(1, round(0.05*numel(ord))), ...
                round(0.50*numel(ord)), ...
                min(numel(ord), round(0.95*numel(ord)))];
    comp_idx = valid_err(ord(pct_idx));
    comp_lbl = {'5th pctl (best)','Median','95th pctl (worst)'};
else
    comp_idx = valid_err(:)';
    comp_lbl = arrayfun(@(i) sprintf('traj %d',i), comp_idx, 'Uni',false);
end

figure('Name','Trajectories','Position',[100 50 400*numel(comp_idx) 350]);
% comp_idx = [comp_idx;2];
for j = 1:numel(comp_idx)
    ci = comp_idx(j);
    subplot(1, numel(comp_idx), j); hold on; box on;

    % x0 = Uobj(ci).get_x0([]);
    x0_raw = Uobj(ci).get_x0([]);
    rhs_r  = WS_wendy{ci}.get_rhs;
    tg     = Uobj(ci).grid{1}(:);
    y      = Uobj(ci).Uobs{1}(:);
    x0     = fit_x0([], rhs_r, tg, y, x0_raw);

    tg = Uobj(ci).grid{1}(:);          
    y  = Uobj(ci).Uobs{1}(:);          
    
    % original measurements
    if opts.want_interp
        t_plot = tcell_raw{ci};
        y_plot = xcell_raw{ci};
    else
        t_plot = Uobj(ci).grid{1}(:);
        y_plot = Uobj(ci).Uobs{1}(:);
    end
    
    scatter(t_plot, y_plot, 20, [0.5 0.5 0.5], 'filled', ...
        'MarkerFaceAlpha',0.6, 'DisplayName','Data');

    if has_true
        wt = true_weights(ci,:);
        if all(isfinite(wt))
            try
                [~,xt] = ode45(@(t,x) growth_rhs(x, wt, growth_type), ...
                               tg, x0, ode_opts_plot);
                plot(tg, xt(:,1), 'g-', 'LineWidth',2.5, 'DisplayName','True');
            catch, warning("True ODE failure "),
            end
        end
    end

    try
        rhs_r = WS_wendy{ci}.get_rhs;
        [~,xr] = ode45(@(t,x) rhs_r(x), tg, x0, ode_opts_plot);
        plot(tg, xr(:,1), 'b--', 'LineWidth',2, 'DisplayName','WENDy');
    catch, warning("WENDy ODE failure "),
    end

    % try
    %     [~,xs] = ode45(@(t,x) bio_rhs(x, bio_shrunk, model_id, ci), ...
    %                     tg, x0, ode_opts_plot);
    %     plot(tg, xs(:,1), 'r-.', 'LineWidth',2, 'DisplayName','Shrunk');
    % catch,  warning("Shrunk ODE failure "),
    % end

    legend('Location','best','FontSize',8);
    xlabel('Time'); ylabel('Size');
    ttl = sprintf('%s  (traj %d)\nRRMSE: raw=%.3g', comp_lbl{j}, ci, err_raw(ci));
    if isfinite(err_shrunk(ci))
        ttl = sprintf('%s  sh=%.3g', ttl, err_shrunk(ci));
    end
    title(ttl, 'FontSize',10);
end
sgtitle(sprintf('Trajectories by RRMSE  (%s)', bio_raw.model_tag));

% ================================================================
%  FIGURE 3: Learned vs true scatter
% ================================================================
if has_true
    figure('Name','Learned vs True','Position',[100 400 350*n_plot 320]);
    for f = 1:n_plot
        subplot(1, n_plot, f); hold on; box on;
        fn = plot_fields{f};
        if ~isfield(bio_true, fn), continue; end

        vr = bio_raw.(fn);
        vs = bio_shrunk.(fn);
        vt = bio_true.(fn);
        mask = isfinite(vr) & isfinite(vt);
        if ~any(mask), continue; end

        tmin = min(vt(mask)); tmax = max(vt(mask));
        pad = 0.05 * max(tmax-tmin, eps);
        plot([tmin-pad tmax+pad], [tmin-pad tmax+pad], '-', ...
            'Color',[0.6 0.6 0.6], 'HandleVisibility','off');

        scatter(vt(mask), vr(mask), 16, [0.2 0.4 0.8], 'filled', ...
            'MarkerFaceAlpha',0.4, 'DisplayName','WENDy');
        mask_s = isfinite(vs) & isfinite(vt);
        scatter(vt(mask_s), vs(mask_s), 16, [0.9 0.3 0.2], 'filled', ...
            'MarkerFaceAlpha',0.4, 'DisplayName','Shrunk');

        xlabel([strrep(fn,'_','\_') ' (true)']);
        ylabel([strrep(fn,'_','\_') ' (learned)']);
        legend('Location','best','FontSize',8);
        title(strrep(fn,'_','\_'));

    end
    sgtitle(sprintf('Learned vs True  (%s)', bio_raw.model_tag));
end

% ================================================================
%  FIGURE 4: RRMSE and E2 histograms
% ================================================================
n_panels = 1 + has_true;
figure('Name','Errors','Position',[100 50 450*n_panels 350]);

subplot(1, n_panels, 1); hold on; box on;
vr = err_raw(isfinite(err_raw));
vs = err_shrunk(isfinite(err_shrunk));
if ~isempty(vr)
    all_e = [vr; vs];
    edg = linspace(log10(max(min(all_e),eps)), log10(max(all_e)), 35);
    histogram(log10(vr), edg, 'FaceColor',[0.2 0.4 0.8], ...
        'FaceAlpha',0.5, 'EdgeColor','none', 'DisplayName','WENDy');
    histogram(log10(vs), edg, 'FaceColor',[0.9 0.3 0.2], ...
        'FaceAlpha',0.45, 'EdgeColor','none', 'DisplayName','Shrunk');
    xline(log10(median(vr)),'b--','LineWidth',1.5,'HandleVisibility','off');
    xline(log10(median(vs)),'r--','LineWidth',1.5,'HandleVisibility','off');
end
xlabel('log_{10}(RRMSE)'); ylabel('Count');
title(sprintf('RRMSE  (med: raw=%.2e, sh=%.2e)', median(vr), median(vs)));
legend('Location','best');

if has_true
    subplot(1,2,2); hold on; box on;
    vr2 = e2_raw(isfinite(e2_raw));
    vs2 = e2_shrunk(isfinite(e2_shrunk));
    if ~isempty(vr2)
        all_e2 = [vr2; vs2];
        edg2 = linspace(log10(max(min(all_e2),eps)), log10(max(all_e2)), 35);
        histogram(log10(vr2), edg2, 'FaceColor',[0.2 0.4 0.8], ...
            'FaceAlpha',0.5, 'EdgeColor','none', 'DisplayName','WENDy');
        histogram(log10(vs2), edg2, 'FaceColor',[0.9 0.3 0.2], ...
            'FaceAlpha',0.45, 'EdgeColor','none', 'DisplayName','Shrunk');
        xline(log10(median(vr2)),'b--','LineWidth',1.5,'HandleVisibility','off');
        xline(log10(median(vs2)),'r--','LineWidth',1.5,'HandleVisibility','off');
    end
    xlabel('log_{10}(E_2)'); ylabel('Count');
    title(sprintf('E_2  (med: raw=%.2e, sh=%.2e)', median(vr2), median(vs2)));
    legend('Location','best');
end
sgtitle(sprintf('Error distributions  (%s)', bio_raw.model_tag));

end % want_plots


%% ==================== Helpers ====================


function [tc, xc] = check_uniform_dt(tcell, xcell)
    L = length(tcell); flag = false;
    tc = cell(size(tcell)); xc = cell(size(xcell));
    for i = 1:L
        t_i = tcell{i};  s_i = xcell{i};
        dt_nom = median(diff(t_i));
        if abs(diff(t_i(end-1:end)) - dt_nom) > 1e-10 * dt_nom
            t_i = t_i(1:end-1);  s_i = s_i(1:end-1); flag = true;
        end
        tc{i} = t_i;  xc{i} = s_i;
    end
    if flag, [tc, xc] = check_uniform_dt(tc, xc); end
end

function nll = reml_neg_profile(tau2, theta, s2)
% Negative REML profile log-likelihood for a single parameter.
    V   = max(tau2 + s2,eps);
    om  = 1 ./ max(V, eps);
    mu  = (om' * theta) / sum(om);
    res = theta - mu;
    nll = 0.5 * sum(log(V)) + 0.5 * log(sum(om)) + 0.5 * sum(res.^2 .* om);
end

function dx = growth_rhs(x, w, gtype)
% RHS from raw ODE weights.  Weight ordering matches the library:
%   lib_wendy active terms after model selection.
%
%   exp:         w(1)*x
%   vb_linear:   w(1)*1 + w(2)*x          (const + x terms)
%   logistic:    w(1)*x + w(2)*x^2
%   gompertz:    w(1)*x + w(2)*x*log(x)
%   vb23:        w(1)*x + w(2)*x^(2/3)    (x + x^{2/3} terms)
    switch lower(string(gtype))
        case {"exp","exponential"}
            dx = w(1) .* x;
        case {"vb","vb_linear"}
            dx = w(1) + w(2) .* x;
        case "logistic"
            dx = w(1).*x + w(2).*x.^2;
        case {"gompertz","gomp"}
            dx = w(1).*x + w(2).*(x.*log(max(x,eps)));
        case {"vb23","vb_metabolic"}
            dx = w(1).*x + w(2).*x.^(2/3);
        otherwise
            error('growth_rhs: unknown type "%s"', gtype);
    end
end

function dx = bio_rhs(x, bio_struct, model_tag, idx)
% RHS from biological parameters (r, S_max) for individual idx.
%
%   All models use r = intrinsic growth rate, S_max = asymptotic size.
%
%   model           ODE in bio params
%   -----           -----------------
%   exp             dx/dt = r * x
%   vb_linear       dx/dt = r * (S_max - x)
%   logistic        dx/dt = r * x * (1 - x/S_max)
%   gompertz        dx/dt = r * x * log(S_max / x)
%   vb23            dx/dt = r * S_max^{1/3} * x^{2/3}  -  r * x
%                         = r * (S_max^{1/3} * x^{2/3} - x)
    r = bio_struct.r(idx);
    switch lower(string(model_tag))
        case {"exp","exponential"}
            dx = r .* x;
        case {"vb","vb_linear"}
            Sm = bio_struct.S_max(idx);
            dx = r .* (Sm - x);
        case "logistic"
            Sm = bio_struct.S_max(idx);
            dx = r .* x .* (1 - x./Sm);
        case {"gompertz","gomp"}
            Sm = bio_struct.S_max(idx);
            dx = r .* x .* log(Sm ./ max(x, eps));
        case {"vb23","vb_metabolic"}
            Sm = bio_struct.S_max(idx);
            dx = r .* (Sm.^(1/3) .* x.^(2/3) - x);
        otherwise
            error('bio_rhs: unknown model "%s"', model_tag);
    end
end

function W = bio_to_weights(bio_struct, model_tag, N)
% Convert (r, S_max) back to ODE weights for E2 computation.
%
%   model           w1                  w2
%   -----           --                  --
%   exp             r                   (none)
%   vb_linear       r*S_max             -r
%   logistic        r                   -r/S_max
%   gompertz        r*log(S_max)        -r
%   vb23            -r                  r*S_max^{1/3}
    r = bio_struct.r(:);
    switch lower(string(model_tag))
        case {"exp","exponential"}
            W = r;
        case {"vb","vb_linear"}
            Sm = bio_struct.S_max(:);
            W = [r.*Sm,  -r];
        case "logistic"
            Sm = bio_struct.S_max(:);
            W = [r,  -r./Sm];
        case {"gompertz","gomp"}
            Sm = bio_struct.S_max(:);
            W = [r.*log(Sm),  -r];
        case {"vb23","vb_metabolic"}
            Sm = bio_struct.S_max(:);
            W = [-r,  r.*Sm.^(1/3)];
        otherwise
            error('bio_to_weights: unknown model "%s"', model_tag);
    end
end


function la = get_linargs(model_tag)
switch lower(string(model_tag))
    case {"exponential"}
        la = {"LB",[0],"UB",[inf]};
    case {"logistic"}
        la = {'LB', [0   -inf], 'UB', [inf   eps]};
    case {"gomp", "gompertz"}
        % la = {'LB', [ -inf 0 ], 'UB', [ 0  inf]};
        la = {'LB', [0 -inf], 'UB', [inf 0]};
    case {"vb","vb_linear"}
        la = {'LB', [0   -inf], 'UB', [inf   0]};
    case {"vb23","vb_metabolic"}
        la = {'LB', [ -inf 0], 'UB', [ 0  inf]};
end
end

function [tcell_out, xcell_out] = interpolate_trajectories(tcell, xcell, ...
                                                           n_target, method)
% INTERPOLATE_TRAJECTORIES  Resample trajectories onto a uniform mesh.
%
%   Defaults: N_TARGET = 150, METHOD = 'pchip'.


    if nargin < 4 || isempty(method),   method   = 'pchip'; end
    if nargin < 3 || isempty(n_target), n_target = 150;     end

    L = numel(tcell);
    tcell_out = cell(size(tcell));
    xcell_out = cell(size(xcell));

    for i = 1:L
        t_i = tcell{i}(:);
        x_i = xcell{i}(:);

        % Guard: need at least 2 points to interpolate
        if numel(t_i) < 2
            tcell_out{i} = t_i;
            xcell_out{i} = x_i;
            continue
        end

        % Skip if already as dense as (or denser than) requested
        if numel(t_i) >= n_target
            tcell_out{i} = t_i;
            xcell_out{i} = x_i;
            continue
        end

        t_new = linspace(t_i(1), t_i(end), n_target).';
        x_new = interp1(t_i, x_i, t_new, method);

        tcell_out{i} = t_new;
        xcell_out{i} = x_new;
    end
end

%% ==================== Population-Sim Helpers ====================
 
function x_T = advect_population(x0, T, w, model_id)
% ADVECT_POPULATION  Advance all initial sizes x0 to age T under
%   dx/ds = growth_rhs(x, w, model_id), as a single vectorized ode45 call.
    if T == 0, x_T = x0; return; end
    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1:numel(x0));
    try
        [~, X] = ode45(@(s, x) growth_rhs(x, w, model_id), [0 T], x0(:), ode_opts);
        x_T    = X(end, :)';
    catch
        % Fallback: per-particle loop in case a model misbehaves at a few x0
        x_T = zeros(size(x0));
        for j = 1:numel(x0)
            try
                [~, xs] = ode45(@(s, x) growth_rhs(x, w, model_id), [0 T], x0(j), ode_opts);
                x_T(j) = xs(end);
            catch
                x_T(j) = x0(j);
            end
        end
    end
end
 
function x_T = advect_with_geff(x0, T, g_eff_handle)
% ADVECT_WITH_GEFF  Advance all initial sizes by an arbitrary g_eff(x)
%   handle (vectorized in x).
    if T == 0, x_T = x0; return; end
    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1:numel(x0));
    try
        [~, X] = ode45(@(s, x) g_eff_handle(x), [0 T], x0(:), ode_opts);
        x_T    = X(end, :)';
    catch
        x_T = zeros(size(x0));
        for j = 1:numel(x0)
            try
                [~, xs] = ode45(@(s, x) g_eff_handle(x), [0 T], x0(j), ode_opts);
                x_T(j) = xs(end);
            catch
                x_T(j) = x0(j);
            end
        end
    end
end
 
function g = eval_g_eff_empirical(x, params_mat, model_id)
% EVAL_G_EFF_EMPIRICAL  Empirical conditional-mean velocity field
%   g_eff(x) = (1/N) sum_i growth_rhs(x, p^(i), model_id),  vectorized in x.
    g = zeros(size(x));
    Nloc = size(params_mat, 2);
    for ii = 1:Nloc
        g = g + growth_rhs(x, params_mat(:, ii)', model_id);
    end
    g = g / Nloc;
end

%--- initial condition-----------
function x0_opt = fit_x0(w_i, rhs_handle, tgrid, xdata, x0_guess)
% FIT_X0  Find the initial condition that best matches data given fixed ODE weights.
%
%   x0_opt = fit_x0(w_i, rhs_handle, tgrid, xdata, x0_guess)
%
%   Given a fixed RHS (from WENDy), solve:
%     x0* = argmin_{x0} sum_j (x(t_j; x0) - xdata_j)^2
%
%   This is a 1D bounded optimization (x0 > 0 for growth data).

    ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    % Search bounds: 50% to 200% of the guess
    lo = 0.5 * x0_guess;
    hi = 2.0 * x0_guess;

    cost = @(x0) trajectory_cost(x0, rhs_handle, tgrid, xdata, ode_opts);

    x0_opt = fminbnd(cost, lo, hi, ...
        optimset('TolX', 1e-10, 'MaxIter', 200, 'Display', 'off'));
end

function J = trajectory_cost(x0, rhs_handle, tgrid, xdata, ode_opts)
    try
        [~, X] = ode45(@(t,x) rhs_handle(x), tgrid, x0, ode_opts);
        J = sum((X(:,1) - xdata(:)).^2);
    catch
        J = inf;
    end
end