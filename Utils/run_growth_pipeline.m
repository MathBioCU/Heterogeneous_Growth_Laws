function res = run_growth_pipeline(filename, opts)
%RUN_GROWTH_PIPELINE  WENDy model selection + empirical Bayes shrinkage
%   for one individual-growth dataset.  Adapted directly from the user's
%   main pipeline script with the population-sims and per-dataset plots
%   stripped out; results are returned in a struct so that the same
%   pipeline can be looped across multiple datasets.
%
%   res = run_growth_pipeline(filename, opts)
%
%   filename  : path to a .mat file containing a struct Data with fields
%               .tracks (cell array of size trajectories) and .t (cell
%               array of corresponding time vectors).
%   opts      : same options struct as the main script (see defaults).
%
%   res       : struct with all post-shrinkage outputs needed for
%               downstream plotting/analysis.

% ---- Default options (same as main script) ----
dflt = struct( ...
    'subsample_k',     1,        ...
    'min_length',      6,        ...
    'n_select',        [],       ...
    'noise_ratio',     0.00,     ...
    'tf_power',        15,       ...
    'want_shrinkage',  true,     ...
    'max_iter',        10,       ...
    'min_valid_frac',  0.5,      ...
    'want_linargs',    true,     ...
    'cell_idx',        [],       ...
    'verbose',         true,     ...
    'want_interp',     false,    ...
    'interp_n_points', 150,      ...
    'interp_method',   'linear', ...
    'rng_seed',        [],       ...   % [] = shuffle; otherwise seed value
    'selection_score', 'sumtraj_keff', ... % only supported value; see below
    'force_model',     '');            % '' = use BIC winner; else override
                                       %       (e.g. 'exponential')
flds = fieldnames(dflt);
for f = 1:numel(flds)
    if ~isfield(opts, flds{f}), opts.(flds{f}) = dflt.(flds{f}); end
end
if ~strcmpi(opts.selection_score, 'sumtraj_keff')
    error('run_growth_pipeline:Score', ...
          ['This version scores laws only with the effective-dof BIC; ' ...
           'opts.selection_score must be ''sumtraj_keff'', got ''%s''.'], ...
          opts.selection_score);
end

% ---- Load data ----
try
    S = load(filename);
    Data = S.Data;
catch ME
    error('run_growth_pipeline:LoadFailed', ...
        'Could not load %s: %s', filename, ME.message);
end
if opts.verbose
    [~, fname_only, ext] = fileparts(filename);
    fprintf('\nLoaded %s%s\n', fname_only, ext);
end

% ---- Prep Data ----
xcell = Data.tracks;  tcell = Data.t;
k = opts.subsample_k;
xcell = cellfun(@(x) x(1:k:end), xcell, 'UniformOutput', false);
tcell = cellfun(@(t) t(1:k:end), tcell, 'UniformOutput', false);

valid_idx = cellfun(@(x) numel(x) > 5 && all(isfinite(x)), xcell);
xcell = xcell(valid_idx);  tcell = tcell(valid_idx);
track_lens = cellfun(@numel, xcell);
valid_idx  = track_lens >= opts.min_length;
xcell = xcell(valid_idx);  tcell = tcell(valid_idx);

[tcell, xcell] = check_uniform_dt(tcell, xcell);
tcell_raw = tcell; xcell_raw = xcell;

if opts.want_interp
    [tcell, xcell] = interpolate_trajectories(tcell, xcell, ...
                        opts.interp_n_points, opts.interp_method);
end

n_total  = numel(xcell);
if isempty(opts.n_select)
    n_select = n_total;
else
    n_select = min(opts.n_select, n_total);
end

if isempty(opts.rng_seed)
    rng('shuffle');
else
    rng(opts.rng_seed, 'twister');
end
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
        Ui.noise{1,1} = Ui.noise{1,1} + (new_ob - ob);
        Ui.Uobs = {new_ob};
    end
    Uobj(i) = Ui;
end

% ---- Library and test functions ----
t_const = term('fhandle', @(x) ones(size(x)));
t_x     = term('fhandle', @(x) x);
t_x2    = term('fhandle', @(x) x.^2);
t_x23   = term('fhandle', @(x) x.^(2/3));
t_xlogx_grad = term('fhandle', @(x) log(max(x,eps)) + 1, ...
                     'gradterms', term('fhandle', @(x) 1 ./ x));
t_xlogx = term('fhandle', @(x) x .* log(x), ...
               'gradterms', t_xlogx_grad);
lib  = library('terms', {t_const, t_x, t_x2, t_x23, t_xlogx});

tf_full = cell(size(Uobj));
for i = 1:numel(Uobj)
    tf_full{i} = testfcn(Uobj(i), 'phifuns', @(v) (1-v.^2).^opts.tf_power);
end

% ---- WENDy + BIC model selection ----
candidate_laws = {'exponential','logistic','gompertz','vb_linear','vb_metabolic'};
n_laws = numel(candidate_laws);

term_map = struct();
term_map.exponential   = 2;
term_map.logistic      = [2, 3];
term_map.gompertz      = [2, 5];
term_map.vb_linear     = [1, 2];
term_map.vb_metabolic  = [2, 4];

wendy_results         = cell(n_laws, 1);
bic_per_law           = inf(n_laws, 1);
rss_per_law           = inf(n_laws, 1);
np_per_law            = zeros(n_laws, 1);
K_per_law             = zeros(n_laws, 1);
n_obs_per_law         = zeros(n_laws, 1);
n_eff_per_law         = zeros(n_laws, 1);
bic_global_fe_per_law = inf(n_laws, 1);
bic_sumtraj_per_law   = inf(n_laws, 1);   % effective-dof BIC, summed
aic_sumtraj_per_law   = inf(n_laws, 1);
rss_mean_per_law      = inf(n_laws, 1);
Ki_all                = cell(n_laws, 1);
Keff_all              = cell(n_laws, 1);
bic_i_all             = cell(n_laws, 1);

if opts.verbose
    fprintf('Fitting %d candidate growth laws via WENDy on N = %d trajectories...\n', ...
        n_laws, N);
end

parfor m = 1:n_laws
    warning('off', 'MATLAB:rankDeficientMatrix');
    law       = candidate_laws{m};
    idx_terms = term_map.(law);
    K_M       = numel(idx_terms);

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
    n_eff_M    = 0;
    BIC_sumtraj_M = 0;
    AIC_sumtraj_M = 0;
    Ki_vec    = nan(N, 1);
    Keff_vec  = nan(N, 1);
    bic_i_vec = nan(N, 1);

    for i = 1:N
        optm_i  = WS_opt();
        W_wendy = wendy_model(Uobj(i), lib_M, tf_full{i}, [1,0]);
        if opts.want_linargs
            [W_wendy, ~, ~, ~, CovW, ~] = optm_i.wendy(W_wendy, ...
                                            "linregargs", linargs_M);
        else
            [W_wendy, ~, ~, ~, CovW, ~] = optm_i.wendy(W_wendy);
        end
        WS_wendy_M{i}  = W_wendy;
        CovWs_M{i}     = CovW;
        params_M(:, i) = W_wendy.weights;

        G_i = W_wendy.G{1};
        b_i = W_wendy.b{1};
        w_i = W_wendy.weights;
        resid = b_i - G_i * w_i; % Equation error residual
        % resid = W_wendy.res; % WENDy residual

        RSS_i = norm(resid)^2;
        Ki    = numel(b_i);
        d_m   = K_M;

        % ---------------------------------------------------------
        % Effective number of independent weak equations.
        % ---------------------------------------------------------
        C_i = W_wendy.cov;
        if isequal(size(C_i), [Ki Ki]) && all(isfinite(C_i(:)))
            Keff = trace(C_i)^2 / max(sum(C_i(:).^2), eps);
        else
            Keff = Ki;
        end
        Keff = min(max(Keff, d_m + 1), Ki);   % keep the logs well posed

        RSS_M   = RSS_M   + RSS_i;
        n_obs_M = n_obs_M + Ki;
        n_eff_M = n_eff_M + Keff;

        BIC_i = Keff * log(max(RSS_i, eps) / Keff) + d_m * log(Keff);
        AIC_i = Keff * log(max(RSS_i, eps) / Keff) + 2 * d_m;

        BIC_sumtraj_M = BIC_sumtraj_M + BIC_i;
        AIC_sumtraj_M = AIC_sumtraj_M + AIC_i;

        Ki_vec(i)    = Ki;
        Keff_vec(i)  = Keff;
        bic_i_vec(i) = BIC_i;
    end

    if n_obs_M == 0
        if opts.verbose
            fprintf('  %-14s : all trajectories failed -- skipped\n', law);
        end
        continue
    end

    n_params      = K_M * N;
    BIC_global_FE = n_eff_M * log(max(RSS_M, eps) / n_eff_M) ...
                  + n_params * log(n_eff_M);

    bic_per_law(m)             = BIC_global_FE;
    bic_global_fe_per_law(m)   = BIC_global_FE;
    bic_sumtraj_per_law(m)     = BIC_sumtraj_M;
    aic_sumtraj_per_law(m)     = AIC_sumtraj_M;
    rss_mean_per_law(m)        = RSS_M / n_obs_M;
    rss_per_law(m)             = RSS_M;
    np_per_law(m)              = n_params;
    K_per_law(m)               = K_M;
    n_obs_per_law(m)           = n_obs_M;
    n_eff_per_law(m)           = n_eff_M;
    Ki_all{m}                  = Ki_vec;
    Keff_all{m}                = Keff_vec;
    bic_i_all{m}               = bic_i_vec;

    wendy_results{m} = struct( ...
        'WS_wendy', {WS_wendy_M}, ...
        'CovWs',    {CovWs_M},    ...
        'params',   params_M,     ...
        'n_obs',    n_obs_M);

    if opts.verbose
        fprintf('  %-14s : K=%d  RSS=%.3e  med Keff=%.1f  BIC_sumtraj(Keff)=%.2f\n', ...
            law, K_M, RSS_M, median(Keff_vec, 'omitnan'), BIC_sumtraj_M);
    end
end

% ---- Select winner: BIC by default, or override with opts.force_model ----
if ~isempty(opts.force_model)
    forced_idx = find(strcmp(candidate_laws, opts.force_model), 1);
    if isempty(forced_idx)
        error(['opts.force_model = "%s" is not in the candidate list. ' ...
               'Valid choices: %s'], opts.force_model, ...
              strjoin(candidate_laws, ', '));
    end
    if ~isfinite(bic_sumtraj_per_law(forced_idx))
        error('Forced model "%s" failed to fit any trajectories.', ...
              opts.force_model);
    end
    best_idx = forced_idx;
    if opts.verbose
        [~, bic_winner_idx] = min(bic_sumtraj_per_law);
        fprintf(['\nForcing model = %s  (BIC winner would have been %s, ' ...
                 'DeltaBIC_sumtraj = %.2f)\n'], ...
                opts.force_model, candidate_laws{bic_winner_idx}, ...
                bic_sumtraj_per_law(forced_idx) - bic_sumtraj_per_law(bic_winner_idx));
    end
else
    [~, best_idx] = min(bic_sumtraj_per_law);
end
if ~isfinite(bic_sumtraj_per_law(best_idx))
    error('Model selection failed: no candidate produced a finite BIC.');
end

model_id    = candidate_laws{best_idx};
WS_wendy    = wendy_results{best_idx}.WS_wendy;
CovWs       = wendy_results{best_idx}.CovWs;
params      = wendy_results{best_idx}.params;
p_lib       = K_per_law(best_idx);
n_obs_best  = n_obs_per_law(best_idx);
n_eff_best  = n_eff_per_law(best_idx);

valid_traj    = ~any(isnan(params), 1);
n_failed_best = sum(~valid_traj);

% ---- Score tables ----
%  Everything below is in effective-dof units: the summed per-trajectory
%  BIC is the criterion, and Delta is normalised by sum_i Keff_i (the
%  number of independent weak equations the score actually counts), not by
%  the number of test functions.
delta_bic      = bic_sumtraj_per_law - bic_sumtraj_per_law(best_idx);
norm_delta_bic = delta_bic / max(n_eff_best, eps);
[~, sort_idx]  = sort(bic_sumtraj_per_law, 'ascend');

Ki_med_per_law   = nan(n_laws, 1);
Keff_med_per_law = nan(n_laws, 1);
for m = 1:n_laws
    if isempty(Ki_all{m}), continue; end
    Ki_med_per_law(m)   = median(Ki_all{m},   'omitnan');
    Keff_med_per_law(m) = median(Keff_all{m}, 'omitnan');
end

bic_table = table( ...
    candidate_laws(sort_idx)',   ...
    bic_sumtraj_per_law(sort_idx), ...
    norm_delta_bic(sort_idx),    ...
    rss_per_law(sort_idx),       ...
    Keff_med_per_law(sort_idx),  ...
    np_per_law(sort_idx),        ...
    K_per_law(sort_idx),         ...
    'VariableNames', ...
    {'Law','BIC_Keff','NormDeltaBIC','RSS','Keff_med','n_params','K'});

delta_global  = bic_global_fe_per_law - min(bic_global_fe_per_law);
delta_sumtraj = bic_sumtraj_per_law   - min(bic_sumtraj_per_law);
delta_aic     = aic_sumtraj_per_law   - min(aic_sumtraj_per_law);

score_table = table( ...
    candidate_laws(:), ...
    K_per_law(:), ...
    n_obs_per_law(:), ...
    n_eff_per_law(:), ...
    Ki_med_per_law(:), ...
    Keff_med_per_law(:), ...
    np_per_law(:), ...
    rss_per_law(:), ...
    rss_mean_per_law(:), ...
    bic_global_fe_per_law(:), ...
    delta_global(:), ...
    bic_sumtraj_per_law(:), ...
    delta_sumtraj(:), ...
    aic_sumtraj_per_law(:), ...
    delta_aic(:), ...
    'VariableNames', { ...
        'Law','K','n_obs','n_eff','Ki_med','Keff_med', ...
        'n_params_global','RSS','RSS_per_obs', ...
        'BIC_global_Keff','Delta_BIC_global_Keff', ...
        'BIC_sumtraj_Keff','Delta_BIC_sumtraj_Keff', ...
        'AIC_sumtraj_Keff','Delta_AIC_sumtraj_Keff'});
score_table = sortrows(score_table, 'BIC_sumtraj_Keff', 'ascend');

if opts.verbose
    infl = median(Ki_med_per_law ./ max(Keff_med_per_law, eps), 'omitnan');
    fprintf('\nSelected model: %s  (failed under winner: %d / %d)\n', ...
        model_id, n_failed_best, N);
    fprintf('Weak dof: median Ki = %.0f, median Keff = %.1f  (inflation %.1fx)\n', ...
        median(Ki_med_per_law, 'omitnan'), ...
        median(Keff_med_per_law, 'omitnan'), infl);
    if abs(infl - 1) < 1e-6
        fprintf(['  NOTE: Keff fell back to Ki everywhere, so WS_wendy.cov is\n' ...
                 '  not a Ki x Ki matrix in this build and no correction is\n' ...
                 '  being applied.  Check that field.\n']);
    end
    fprintf('\n');
    disp(bic_table)

    % Per-trajectory split between the top two laws.  A summed preference
    % carried by a minority of cells is a different finding from one shared
    % by most of them, and only this split separates the two.
    if numel(sort_idx) >= 2
        a = sort_idx(1);  b = sort_idx(2);
        if ~isempty(bic_i_all{a}) && ~isempty(bic_i_all{b})
            d  = bic_i_all{b} - bic_i_all{a};    % >0 : cell prefers the winner
            ok = isfinite(d);
            fprintf(['Trajectories individually preferring %s over %s: ' ...
                     '%.1f%%  (summed Delta = %+.1f)\n\n'], ...
                candidate_laws{a}, candidate_laws{b}, ...
                100*mean(d(ok) > 0), sum(d(ok)));
        end
    end
end

% ---- Convert to biological parameters ----
params_converted_struct = weights_to_biological_params(params', model_id);

% ---- Empirical Bayes shrinkage ----
params_shrunk_w          = [];
params_bio_shrunk_struct = [];
mu_hat                   = [];
tau2_hat                 = [];
tau_hat                  = [];
gamma                    = [];
shrink_names             = {};
bio_S_hat                = [];

if opts.want_shrinkage
    CovW_gls_full = cell(N,1);
    for i = 1:N
        G_i = WS_wendy{i}.G{1};
        C_i = WS_wendy{i}.cov;
        CovW_gls_full{i} = pinv(G_i' * (C_i \ G_i));
    end

    valid     = params_converted_struct.valid;
    N_valid   = sum(valid);
    idx_valid = find(valid);

    shrink_names = {'r'};
    if isfield(params_converted_struct, 'S_max')
        Smax_vals = params_converted_struct.S_max(valid);
        if all(isfinite(Smax_vals))
            shrink_names{end+1} = 'S_max';
        elseif opts.verbose
            fprintf('  Skipping S_max shrinkage (non-finite values present)\n');
        end
    end
    n_shrink = numel(shrink_names);

    bio_params = zeros(N_valid, n_shrink);
    for f = 1:n_shrink
        bio_params(:,f) = params_converted_struct.(shrink_names{f})(valid);
    end

    bio_S_hat = zeros(N_valid, n_shrink);
    h_scale = 1e-7;
    for ii = 1:N_valid
        i = idx_valid(ii);
        w_i = params(:, i);
        S_full_i = CovW_gls_full{i};

        J = zeros(n_shrink, p_lib);
        for k_idx = 1:p_lib
            h = h_scale * max(abs(w_i(k_idx)), 1e-8);
            w_fwd = w_i; w_fwd(k_idx) = w_fwd(k_idx) + h;
            w_bwd = w_i; w_bwd(k_idx) = w_bwd(k_idx) - h;
            s_fwd = weights_to_biological_params(w_fwd', model_id);
            s_bwd = weights_to_biological_params(w_bwd', model_id);
            for f = 1:n_shrink
                J(f,k_idx) = (s_fwd.(shrink_names{f}) - s_bwd.(shrink_names{f})) / (2*h);
            end
        end
        CovBio = J * S_full_i * J';
        bio_S_hat(ii,:) = max(diag(CovBio)', 0);
    end

    for f = 1:n_shrink
        emp_var_f = var(bio_params(:,f));
        cap = 100 * emp_var_f;
        bad = ~isfinite(bio_S_hat(:,f)) | bio_S_hat(:,f) > cap;
        if any(bad)
            if opts.verbose
                fprintf('  Capped %d extreme S_hat entries for "%s" at %.2e\n', ...
                    sum(bad), shrink_names{f}, cap);
            end
            bio_S_hat(bad,f) = cap;
        end
    end

    mu_hat            = zeros(1, n_shrink);
    tau2_hat          = zeros(1, n_shrink);
    gamma             = zeros(N_valid, n_shrink);
    bio_params_shrunk = zeros(N_valid, n_shrink);

    for f = 1:n_shrink
        theta = bio_params(:,f);
        s2    = bio_S_hat(:,f);

        neg_reml = @(tau2) reml_neg_profile(tau2, theta, s2);
        tau2_lo  = 0;
        tau2_hi  = max(10 * var(theta), 10 * max(s2));
        if ~isfinite(tau2_hi) || tau2_hi <= 0, tau2_hi = 1; end

        [tau2_opt, fval_opt] = fminbnd(neg_reml, tau2_lo, tau2_hi, ...
            optimset('TolX', 1e-12, 'MaxIter', 500, 'Display', 'off'));
        if neg_reml(0) <= fval_opt, tau2_opt = 0; end

        V_f = tau2_opt + s2;
        om_f = 1 ./ max(V_f, eps);
        mu_f = (om_f' * theta) / sum(om_f);
        gam_f = tau2_opt ./ max(V_f, eps);
        theta_shrunk = mu_f + gam_f .* (theta - mu_f);

        mu_hat(f)              = mu_f;
        tau2_hat(f)            = tau2_opt;
        gamma(:,f)             = gam_f;
        bio_params_shrunk(:,f) = theta_shrunk;
    end
    tau_hat = sqrt(tau2_hat);

    params_bio_shrunk_struct = params_converted_struct;
    for f = 1:n_shrink
        vals = params_bio_shrunk_struct.(shrink_names{f});
        vals(valid) = bio_params_shrunk(:,f);
        params_bio_shrunk_struct.(shrink_names{f}) = vals;
    end
    params_shrunk_w = bio_to_weights(params_bio_shrunk_struct, model_id, N);

    if opts.verbose
        fprintf('\n--- Empirical Bayes Summary ---\n');
        for f = 1:n_shrink
            fprintf('  %s: mu_hat = %.4e, tau_hat = %.4e (naive std = %.4e), CV = %.3f\n', ...
                shrink_names{f}, mu_hat(f), tau_hat(f), ...
                std(bio_params(:,f)), tau_hat(f)/mu_hat(f));
        end
    end
end

% ---- Package output ----
res = struct();
res.filename             = filename;
res.opts                 = opts;
res.model_id             = model_id;
res.candidate_laws       = {candidate_laws};
res.bic_table            = bic_table;
res.score_table          = score_table;
res.bic_per_law          = bic_per_law;
res.bic_sumtraj_per_law  = bic_sumtraj_per_law;   % effective-dof BIC
res.bic_global_per_law   = bic_global_fe_per_law;
res.aic_sumtraj_per_law  = aic_sumtraj_per_law;
res.bic_i_per_law        = bic_i_all;             % per-trajectory BIC_i
res.Ki_per_law           = Ki_all;
res.Keff_per_law         = Keff_all;
res.n_eff_per_law        = n_eff_per_law;
res.params               = params;
res.params_shrunk_w      = params_shrunk_w;
res.params_converted     = params_converted_struct;
res.params_bio_shrunk    = params_bio_shrunk_struct;
res.mu_hat               = mu_hat;
res.tau_hat              = tau_hat;
res.tau2_hat             = tau2_hat;
res.gamma                = gamma;
res.shrink_names         = {shrink_names};
res.tracks               = xcell;
res.t                    = tcell;
res.tracks_raw           = xcell_raw;
res.t_raw                = tcell_raw;
res.N                    = N;
res.rand_idx             = rand_idx;
res.valid_traj           = valid_traj;
res.p_lib                = p_lib;
res.WS_wendy             = {WS_wendy};
res.CovWs                = {CovWs};

end % main function


% ======================== Local Helpers ========================

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
    V   = tau2 + s2;
    om  = 1 ./ max(V, eps);
    mu  = (om' * theta) / sum(om);
    res = theta - mu;
    nll = 0.5 * sum(log(V)) + 0.5 * log(sum(om)) + 0.5 * sum(res.^2 .* om);
end

function W = bio_to_weights(bio_struct, model_tag, N)
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
            la = {'LB', [0 -inf], 'UB', [inf 0]};
        case {"vb","vb_linear"}
            la = {'LB', [0   -inf], 'UB', [inf   0]};
        case {"vb23","vb_metabolic"}
            la = {'LB', [ -inf 0], 'UB', [ 0  inf]};
    end
end

function [tcell_out, xcell_out] = interpolate_trajectories(tcell, xcell, ...
                                                           n_target, method)
    if nargin < 4 || isempty(method),   method   = 'pchip'; end
    if nargin < 3 || isempty(n_target), n_target = 150;     end

    L = numel(tcell);
    tcell_out = cell(size(tcell));
    xcell_out = cell(size(xcell));

    for i = 1:L
        t_i = tcell{i}(:);
        x_i = xcell{i}(:);
        if numel(t_i) < 2
            tcell_out{i} = t_i;  xcell_out{i} = x_i;  continue
        end
        if numel(t_i) >= n_target
            tcell_out{i} = t_i;  xcell_out{i} = x_i;  continue
        end
        t_new = linspace(t_i(1), t_i(end), n_target).';
        x_new = interp1(t_i, x_i, t_new, method);
        tcell_out{i} = t_new;
        xcell_out{i} = x_new;
    end
end