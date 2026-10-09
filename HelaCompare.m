%% run_hela_phase_comparison.m
%  Run WENDy + BIC model selection on HeLa G1 and S-G2 phases separately,
%  then perform a paired t-test on growth rates for cells observed in both.
%
%  Prerequisites:
%    - All4Exp_df.csv, All4Exp_table.csv from Dryad doi:10.6078/D12M6C
%    - convert_hela_fxm_to_wsindy.m on the path
%    - wsindy_obj_base/ on the path

clear all; close all;
addpath(genpath('./wsindy_obj_base/'));

%% ==================== 1. Convert both phases ============================
fprintf('=== Converting G1 phase ===\n');
D = load('./Data/HeLa/HeLa_G1.mat');
Data_G1 = D.Data;

fprintf('\n=== Converting S-G2 phase ===\n');
D = load('./Data/HeLa/HeLa_SG2.mat');
Data_SG2 = D.Data;

%% ==================== 2. Shared options ========================
opts = struct();
opts.noise_ratio     = 0.0;
opts.want_linargs    = true;
opts.min_length      = 6;
opts.subsample_k     = 1;
opts.tf_power        = 15;
opts.want_interp     = false;
opts.interp_n_points = 150;
opts.interp_method   = 'linear';
opts.verbose         = true;

%% ==================== 3. Library & candidate laws =======================
t_const = term('fhandle', @(x) ones(size(x)));
t_x     = term('fhandle', @(x) x);
t_x2    = term('fhandle', @(x) x.^2);
t_x23   = term('fhandle', @(x) x.^(2/3));
t_xlogx_grad = term('fhandle', @(x) log(max(x,eps)) + 1, ...
                     'gradterms', term('fhandle', @(x) 1 ./ x));
t_xlogx = term('fhandle', @(x) x .* log(x), ...
               'gradterms', t_xlogx_grad);
lib  = library('terms', {t_const, t_x, t_x2, t_x23, t_xlogx});

candidate_laws = {'exponential','logistic','gompertz','vb_linear','vb_metabolic'};
n_laws = numel(candidate_laws);

term_map = struct();
term_map.exponential   = 2;
term_map.logistic      = [2, 3];
term_map.gompertz      = [2, 5];
term_map.vb_linear     = [1, 2];
term_map.vb_metabolic  = [2, 4];

%% ==================== 4. Run WENDy on both phases =======================
phase_names = {'G1', 'SG2'};
phase_data  = {Data_G1, Data_SG2};
phase_results = struct();

for ph = 1:2
    Data = phase_data{ph};
    pname = phase_names{ph};
    fprintf('\n############################################################\n');
    fprintf('  Phase: %s  (%d trajectories)\n', pname, numel(Data.tracks));
    fprintf('############################################################\n');

    % ---- Prep data ----
    xcell = Data.tracks;  tcell = Data.t;
    k = opts.subsample_k;
    xcell = cellfun(@(x) x(1:k:end), xcell, 'UniformOutput', false);
    tcell = cellfun(@(t) t(1:k:end), tcell, 'UniformOutput', false);

    valid_idx = cellfun(@(x) numel(x) > 5 && all(isfinite(x)), xcell);
    xcell = xcell(valid_idx);  tcell = tcell(valid_idx);

    % Also filter cell_ids to match
    cell_ids_all = Data.growth.params.cell_ids;
    cell_ids_all = cell_ids_all(valid_idx);

    track_lens = cellfun(@numel, xcell);
    valid_idx2 = track_lens >= opts.min_length;
    xcell = xcell(valid_idx2);  tcell = tcell(valid_idx2);
    cell_ids_all = cell_ids_all(valid_idx2);

    [tcell, xcell] = check_uniform_dt(tcell, xcell);

    if opts.want_interp
        [tcell, xcell] = interpolate_trajectories(tcell, xcell, ...
                            opts.interp_n_points, opts.interp_method);
    end

    xcell = cellfun(@(x) x(:), xcell, 'UniformOutput', false);
    tcell = cellfun(@(t) t(:), tcell, 'UniformOutput', false);

    N = numel(xcell);
    if N == 0
        warning('Phase %s: no valid trajectories.', pname);
        continue;
    end
    fprintf('  %d valid trajectories after filtering.\n', N);

    % ---- Build wsindy_data objects ----
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

    % ---- Test functions ----
    tf_full = cell(size(Uobj));
    for i = 1:numel(Uobj)
        tf_full{i} = testfcn(Uobj(i), 'phifuns', @(v) (1-v.^2).^opts.tf_power);
    end

    % ---- WENDy + BIC model selection ----
    wendy_results  = cell(n_laws, 1);
    bic_sumtraj_per_law = inf(n_laws, 1);
    rss_per_law    = inf(n_laws, 1);
    K_per_law      = zeros(n_laws, 1);
    n_obs_per_law  = zeros(n_laws, 1);

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
        BIC_sumtraj_M = 0;

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

            G_i   = W_wendy.G{1};
            b_i   = W_wendy.b{1};
            w_i   = W_wendy.weights;
            C     = W_wendy.cov;
            resid = b_i - G_i * w_i;

            RSS_i = norm(resid)^2;
            Ki    = numel(b_i);
            d_m   = K_M;

            if isequal(size(C), [Ki Ki]) && all(isfinite(C(:)))
                % Effective num of weak-form equations
                Keff = trace(C)^2 / max(sum(C(:).^2), eps); 
            else
                Keff = Ki;
            end
            Keff = min(max(Keff, d_m + 1), Ki);  % keep the logs well posed

            RSS_M   = RSS_M   + RSS_i;
            n_obs_M = n_obs_M + Ki;

            BIC_i = Ki * log(max(RSS_i, eps) / Ki) + d_m * log(Keff);
            BIC_sumtraj_M = BIC_sumtraj_M + BIC_i;
        end

        bic_sumtraj_per_law(m) = BIC_sumtraj_M;
        rss_per_law(m)         = RSS_M;
        K_per_law(m)           = K_M;
        n_obs_per_law(m)       = n_obs_M;

        wendy_results{m} = struct( ...
            'WS_wendy', {WS_wendy_M}, ...
            'CovWs',    {CovWs_M},    ...
            'params',   params_M,     ...
            'n_obs',    n_obs_M);

        if opts.verbose
            fprintf('  %s | %-14s : K=%d  RSS=%.3e  BIC_sum=%.2f\n', ...
                pname, law, K_M, RSS_M, BIC_sumtraj_M);
        end
    end

    % ---- Select winner ----
    [~, best_idx] = min(bic_sumtraj_per_law);
    model_id = candidate_laws{best_idx};
    WS_wendy = wendy_results{best_idx}.WS_wendy;
    params   = wendy_results{best_idx}.params;
    p_lib    = K_per_law(best_idx);

    % BIC table
    [~, sort_idx] = sort(bic_sumtraj_per_law, 'ascend');
    bic_table = table( ...
        candidate_laws(sort_idx)', ...
        bic_sumtraj_per_law(sort_idx), ...
        rss_per_law(sort_idx), ...
        K_per_law(sort_idx), ...
        'VariableNames', {'Law','BIC_sumtraj','RSS','K'});
    fprintf('\n  %s BIC table:\n', pname);
    disp(bic_table);
    fprintf('  %s => Selected model: %s\n', pname, model_id);

    % ---- Convert to biological parameters ----
    bio = weights_to_biological_params(params', model_id);

    % ---- GLS covariances for precision-weighted test ----
    CovW_gls = cell(N, 1);
    for i = 1:N
        G_i = WS_wendy{i}.G{1};
        C_i = WS_wendy{i}.cov;
        CovW_gls{i} = pinv(G_i' * (C_i \ G_i));
    end

    % ---- Propagate to get per-cell variance of r ----
    r_var = nan(N, 1);
    h_scale = 1e-7;
    for i = 1:N
        w_i = params(:, i);
        S_i = CovW_gls{i};
        J = zeros(1, p_lib);
        for kk = 1:p_lib
            h = h_scale * max(abs(w_i(kk)), 1e-8);
            w_fwd = w_i; w_fwd(kk) = w_fwd(kk) + h;
            w_bwd = w_i; w_bwd(kk) = w_bwd(kk) - h;
            s_fwd = weights_to_biological_params(w_fwd', model_id);
            s_bwd = weights_to_biological_params(w_bwd', model_id);
            J(1,kk) = (s_fwd.r - s_bwd.r) / (2*h);
        end
        r_var(i) = max(J * S_i * J', 0);
    end

    % ---- Store ----
    % Convert cell_ids to numeric for matching
    if iscell(cell_ids_all)
        cids_num = cellfun(@(x) double(x), cell_ids_all);
    else
        cids_num = double(cell_ids_all);
    end

    phase_results(ph).name      = pname;
    phase_results(ph).model_id  = model_id;
    phase_results(ph).N         = N;
    phase_results(ph).bio       = bio;
    phase_results(ph).params    = params;
    phase_results(ph).WS_wendy  = WS_wendy;
    phase_results(ph).cell_ids  = cids_num(:);
    phase_results(ph).r         = bio.r(:);
    phase_results(ph).r_var     = r_var(:);
    phase_results(ph).bic_table = bic_table;
end

%% ==================== 5. Match cells across phases ======================
ids_G1  = phase_results(1).cell_ids;
ids_SG2 = phase_results(2).cell_ids;
r_G1    = phase_results(1).r;
r_SG2   = phase_results(2).r;
rv_G1   = phase_results(1).r_var;
rv_SG2  = phase_results(2).r_var;

[shared_ids, idx_G1, idx_SG2] = intersect(ids_G1, ids_SG2);
n_paired = numel(shared_ids);

fprintf('\n=== Cell matching ===\n');
fprintf('  G1 cells  : %d\n', numel(ids_G1));
fprintf('  SG2 cells : %d\n', numel(ids_SG2));
fprintf('  Paired    : %d\n', n_paired);
fprintf('  G1-only   : %d\n', numel(ids_G1) - n_paired);
fprintf('  SG2-only  : %d\n', numel(ids_SG2) - n_paired);

if n_paired < 3
    error('Too few paired cells (%d) for a meaningful test.', n_paired);
end

r_G1_paired  = r_G1(idx_G1);
r_SG2_paired = r_SG2(idx_SG2);
delta_r      = r_SG2_paired - r_G1_paired;

% Filter out any NaN/Inf pairs
good = isfinite(delta_r);
r_G1_paired  = r_G1_paired(good);
r_SG2_paired = r_SG2_paired(good);
delta_r      = delta_r(good);
idx_G1_good  = idx_G1(good);
idx_SG2_good = idx_SG2(good);
n_paired     = sum(good);

fprintf('  Valid pairs after NaN removal: %d\n', n_paired);

%% ==================== 6. Paired t-test ==================================
[h, p_ttest, ci, stats] = ttest(delta_r);

fprintf('\n=== Paired t-test: lambda_SG2 - lambda_G1 ===\n');
fprintf('  n_paired               = %d\n', n_paired);
fprintf('  mean(lambda_G1)        = %.6f / hour\n', mean(r_G1_paired));
fprintf('  mean(lambda_SG2)       = %.6f / hour\n', mean(r_SG2_paired));
fprintf('  mean(delta)            = %.6f / hour\n', mean(delta_r));
fprintf('  std(delta)             = %.6f\n', std(delta_r));
fprintf('  t(%d)                  = %.4f\n', stats.df, stats.tstat);
fprintf('  p-value                = %.4e\n', p_ttest);
fprintf('  95%% CI for delta      = [%.6f, %.6f]\n', ci(1), ci(2));
fprintf('  Significant at 0.05?   %s\n', string(h == 1));

%% ==================== 7. Nonparametric confirmation =====================
[p_wilcox, ~, stats_w] = signrank(delta_r);
fprintf('\n=== Wilcoxon signed-rank (nonparametric) ===\n');
fprintf('  p-value                = %.4e\n', p_wilcox);

%% ==================== 8. Precision-weighted test ========================
s2_G1_paired  = rv_G1(idx_G1_good);
s2_SG2_paired = rv_SG2(idx_SG2_good);
s2_delta      = s2_G1_paired + s2_SG2_paired;
omega         = 1 ./ max(s2_delta, eps);
mu_weighted   = (omega' * delta_r) / sum(omega);
se_weighted   = sqrt(1 / sum(omega));
z_stat        = mu_weighted / se_weighted;
p_weighted    = 2 * (1 - normcdf(abs(z_stat)));

fprintf('\n=== Precision-weighted z-test ===\n');
fprintf('  weighted mean(delta)   = %.6f / hour\n', mu_weighted);
fprintf('  SE                     = %.6f\n', se_weighted);
fprintf('  z                      = %.4f\n', z_stat);
fprintf('  p-value                = %.4e\n', p_weighted);

%% ==================== 9. Effect size ====================================
d_cohen = mean(delta_r) / std(delta_r);
pct_increase = 100 * mean(delta_r) / max(mean(r_G1_paired), eps);

fprintf('\n=== Effect size ===\n');
fprintf('  Cohen''s d (paired)    = %.4f\n', d_cohen);
fprintf('  Relative increase     = %.1f%%\n', pct_increase);

%% ==================== 10. Plots =========================================

% ---- Figure 1: Paired growth rates scatter ----
figure('Name','Paired growth rates','Position',[50 300 500 450]);
hold on; box on; grid on;
scatter(r_G1_paired, r_SG2_paired, 30, [0.2 0.5 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.5);
mn = min([r_G1_paired; r_SG2_paired]);
mx = max([r_G1_paired; r_SG2_paired]);
pad = 0.05 * (mx - mn);
plot([mn-pad mx+pad], [mn-pad mx+pad], 'k--', 'LineWidth', 1);
xlabel('\lambda_{G1} (hour^{-1})');
ylabel('\lambda_{S-G2} (hour^{-1})');
title(sprintf('Paired growth rates (n = %d cells)', n_paired));

% ---- Figure 2: Histogram of delta_r ----
figure('Name','Delta lambda','Position',[570 300 500 400]);
hold on; box on;
histogram(delta_r, 30, 'Normalization','pdf', ...
    'FaceColor', [0.4 0.6 0.9], 'EdgeColor','none', 'FaceAlpha', 0.7);
xline(0, 'k--', 'LineWidth', 1.5);
xline(mean(delta_r), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('mean = %.4f', mean(delta_r)));
xlabel('\Delta\lambda = \lambda_{S-G2} - \lambda_{G1}  (hour^{-1})');
ylabel('pdf');
title(sprintf('\\Delta\\lambda distribution   (p_{paired-t} = %.2e,  d = %.2f)', ...
    p_ttest, d_cohen));
legend('Location','best');

% ---- Figure 3: Overlaid growth rate distributions ----
figure('Name','Growth rate distributions','Position',[50 50 550 400]);
hold on; box on;

% Use ALL cells (not just paired) for the marginal distributions
r_G1_all  = phase_results(1).r;
r_SG2_all = phase_results(2).r;
r_G1_all  = r_G1_all(isfinite(r_G1_all));
r_SG2_all = r_SG2_all(isfinite(r_SG2_all));

all_r = [r_G1_all; r_SG2_all];
lo = quantile(all_r, 0.005);
hi = quantile(all_r, 0.995);
edges = linspace(lo, hi, 40);

histogram(r_G1_all, edges, 'Normalization','pdf', ...
    'FaceColor', [0.3 0.5 0.9], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName', sprintf('G1 (n=%d, \\mu=%.4f)', numel(r_G1_all), mean(r_G1_all)));
histogram(r_SG2_all, edges, 'Normalization','pdf', ...
    'FaceColor', [0.9 0.4 0.3], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName', sprintf('S-G2 (n=%d, \\mu=%.4f)', numel(r_SG2_all), mean(r_SG2_all)));
xlabel('\lambda  (hour^{-1})');
ylabel('pdf');
legend('Location','best','FontSize',10);
title('Growth rate distributions: G1 vs S-G2');

% ---- Figure 4: Fitted normal PDFs ----
figure('Name','Fitted PDFs','Position',[620 50 550 400]);
hold on; box on;
x_pdf = linspace(lo, hi, 500);
plot(x_pdf, normpdf(x_pdf, mean(r_G1_all), std(r_G1_all)), '-', ...
    'Color', [0.3 0.5 0.9], 'LineWidth', 2, ...
    'DisplayName', sprintf('G1 (\\mu=%.4f, \\sigma=%.4f)', mean(r_G1_all), std(r_G1_all)));
plot(x_pdf, normpdf(x_pdf, mean(r_SG2_all), std(r_SG2_all)), '-', ...
    'Color', [0.9 0.4 0.3], 'LineWidth', 2, ...
    'DisplayName', sprintf('S-G2 (\\mu=%.4f, \\sigma=%.4f)', mean(r_SG2_all), std(r_SG2_all)));
xlabel('\lambda  (hour^{-1})');
ylabel('pdf');
legend('Location','best','FontSize',10);
title('Fitted normal PDFs: G1 vs S-G2');

% ---- Figure 5: Per-cell paired comparison (slopegraph) ----
figure('Name','Slopegraph','Position',[100 100 350 500]);
hold on; box on;
for j = 1:n_paired
    if delta_r(j) > 0
        clr = [0.8 0.2 0.2 0.15];
    else
        clr = [0.2 0.2 0.8 0.15];
    end
    plot([1 2], [r_G1_paired(j) r_SG2_paired(j)], '-', ...
        'Color', clr, 'LineWidth', 0.8);
end
% Overlay means
plot([1 2], [mean(r_G1_paired) mean(r_SG2_paired)], 'ko-', ...
    'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
set(gca, 'XTick', [1 2], 'XTickLabel', {'G1','S-G2'}, 'XLim', [0.7 2.3]);
ylabel('\lambda  (hour^{-1})');
title(sprintf('Per-cell growth rate shift  (n = %d)', n_paired));

%% ==================== Summary table =====================================
fprintf('\n\n================== SUMMARY ==================\n');
fprintf('%-20s %-12s %-12s\n', '', 'G1', 'S-G2');
fprintf('%-20s %-12s %-12s\n', '--------', '--------', '--------');
fprintf('%-20s %-12s %-12s\n', 'Selected model', ...
    phase_results(1).model_id, phase_results(2).model_id);
fprintf('%-20s %-12d %-12d\n', 'N trajectories', ...
    phase_results(1).N, phase_results(2).N);
fprintf('%-20s %-12.5f %-12.5f\n', 'mean(lambda)', ...
    mean(r_G1_all), mean(r_SG2_all));
fprintf('%-20s %-12.5f %-12.5f\n', 'std(lambda)', ...
    std(r_G1_all), std(r_SG2_all));
fprintf('\n%-20s %d\n', 'Paired cells', n_paired);
fprintf('%-20s %.5f\n', 'mean(delta)', mean(delta_r));
fprintf('%-20s %.4f\n', 'Cohen''s d', d_cohen);
fprintf('%-20s %.1f%%\n', 'Relative increase', pct_increase);
fprintf('%-20s %.4e\n', 'p (paired t)', p_ttest);
fprintf('%-20s %.4e\n', 'p (Wilcoxon)', p_wilcox);
fprintf('%-20s %.4e\n', 'p (weighted z)', p_weighted);
fprintf('=============================================\n\n');


%% ==================== Helper functions ==================================

function [tc, xc] = check_uniform_dt(tcell, xcell)
    L = length(tcell); flag = false;
    tc = cell(size(tcell)); xc = cell(size(xcell));
    for i = 1:L
        t_i = tcell{i};  s_i = xcell{i};
        if numel(t_i) < 2, tc{i}=t_i; xc{i}=s_i; continue; end
        dt_nom = median(diff(t_i));
        if abs(diff(t_i(end-1:end)) - dt_nom) > 1e-10 * dt_nom
            t_i = t_i(1:end-1);  s_i = s_i(1:end-1); flag = true;
        end
        tc{i} = t_i;  xc{i} = s_i;
    end
    if flag, [tc, xc] = check_uniform_dt(tc, xc); end
end

function la = get_linargs(model_tag)
    switch lower(string(model_tag))
        case {"exponential"}
            la = {"LB",[0],"UB",[inf]};
        case {"logistic"}
            la = {'LB', [0   -inf], 'UB', [inf   0]};
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
        t_i = tcell{i}(:);  x_i = xcell{i}(:);
        if numel(t_i) < 2 || numel(t_i) >= n_target
            tcell_out{i} = t_i;  xcell_out{i} = x_i;  continue;
        end
        t_new = linspace(t_i(1), t_i(end), n_target).';
        x_new = interp1(t_i, x_i, t_new, method);
        tcell_out{i} = t_new;  xcell_out{i} = x_new;
    end
end