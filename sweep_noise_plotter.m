%% plot_sweep_results.m
% Visualize WENDy growth-law-ID sweeps across multiple datasets.
%
% Figures:
%   1. Selection success rate vs noise (one solid line per model).
%   2. E2 mean OR median vs noise with +/- 1 std sleeve.
%   3. Relative error in learned vs true heterogeneity for r and S_max
%      (one panel per parameter), raw and shrunk overlaid.  Values above
%      the cap (default 1) are clipped to the cap and the y-axis is
%      relabeled ">1" at the cap to flag the clipping.  Exponential is
%      omitted from the S_max panel.

clear; close all;

%% ----------------- Configuration -----------------
results_dir = './sweep_results/';
file_list   = {};            % leave empty to glob sweep_*.mat in results_dir

% Figure 2 toggles
e2_metric              = 'median';   % 'mean' or 'median'
e2_yscale              = 'log';      % 'linear' or 'log'
e2_sleeve_clamp_factor = 0.1;        % log-y only: lower sleeve >= ctr*this
fit_powerlaw           = true;       % fit & overlay a global E2 vs noise power law

% Figure 3 toggles
het_agg = 'mean';   % 'mean' or 'median' (which across-realization aggregate)
het_cap = 3.0;        % cap on |rel diff|; values above shown at the cap

% Styling
LW           = 1.8;
MS           = 7;
sleeve_alpha = 0.12;

%% ----------------- Discover & load -----------------
if isempty(file_list)
    L = dir(fullfile(results_dir, 'sweep_*.mat'));
    [~, ord] = sort({L.name});
    L = L(ord);
    file_list = {L.name};
end
if isempty(file_list)
    error('No sweep_*.mat files found in %s.', results_dir);
end
n_files = numel(file_list);
fprintf('Loading %d sweep file(s):\n', n_files);

styles      = build_style_map();
datasets    = cell(n_files, 1);
seen_models = {};

for k = 1:n_files
    fp = file_list{k};
    if ~exist(fp, 'file')
        fp = fullfile(results_dir, file_list{k});
    end
    fprintf('  [%d/%d] %s\n', k, n_files, fp);
    S = load(fp);

    tm = normalize_tag(S.agg.true_model);
    if any(strcmp(seen_models, tm))
        warning('plot:dup', ...
            'Multiple files for true model "%s"; lines will overlap.', tm);
    end
    seen_models{end+1} = tm; %#ok<SAGROW>

    [rs, ss] = pool_e2_stats(S.results_flat, S.noise_levels, S.n_realizations);

    datasets{k} = struct( ...
        'file',          file_list{k}, ...
        'true_model',    tm, ...
        'display_name',  display_name(tm), ...
        'noise_levels',  S.noise_levels(:), ...
        'agg',           S.agg, ...
        'e2_raw',        rs, ...
        'e2_shrunk',     ss);
end

%% ----------------- Figure 1: Success rate -----------------
figure('Name','Selection success rate');
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax, 'XScale', 'log');

for k = 1:n_files
    D  = datasets{k};
    st = styles.(D.true_model);
    nl = D.noise_levels;
    yv = D.agg.success_rate;
    ok = nl > 0 & isfinite(yv);
    if ~any(ok), continue; end
    plot(nl(ok), yv(ok), ...
        'Color', st.color, 'LineStyle', '-', ...
        'Marker', st.marker, 'LineWidth', LW, 'MarkerSize', MS, ...
        'MarkerFaceColor', st.color, 'DisplayName', D.display_name);
end
xlabel('noise ratio  (log scale; noise=0 omitted)');
ylabel('selection success rate');
title('Model selection accuracy vs noise level');
ylim([-0.02 1.05]);
legend('Location', 'best', 'FontSize', 10);


%% ----------------- Figure 2: E2 raw ---------------
figure('Name','E2 vs noise');
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax, 'XScale', 'log', 'YScale', e2_yscale);

% Style key anchors (plotted first so they sit first in the legend).
plot(NaN, NaN, 'k-',  'LineWidth', LW, 'DisplayName', 'raw WENDy');
plot(NaN, NaN, 'k--', 'LineWidth', LW, 'DisplayName', 'shrunk');

for k = 1:n_files
    D  = datasets{k};
    st = styles.(D.true_model);
    nl = D.noise_levels;

    for tag = {'raw'}
        which = tag{1};
        S2  = D.(['e2_' which]);
        ctr = S2.(e2_metric);
        sd  = S2.std;
        ok  = nl > 0 & isfinite(ctr) & isfinite(sd);
        if ~any(ok), continue; end

        nl_p  = nl(ok);
        ctr_p = ctr(ok);
        sd_p  = sd(ok);

        up = ctr_p + sd_p;
        lo = ctr_p - sd_p;
        if strcmp(e2_yscale, 'log')
            lo = max(lo, ctr_p * e2_sleeve_clamp_factor);
        end

        % Sleeve (no legend entry).
        fill([nl_p; flipud(nl_p)], [up; flipud(lo)], st.color, ...
            'FaceAlpha', sleeve_alpha, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Center line.  Only the raw line gets a legend entry; the
        % shrunk line reuses the same color/marker so the legend stays
        % at n_models + 2 entries.
        if strcmp(which, 'raw')
            ls = '-';  dn = D.display_name;  hvis = 'on';
        else
            ls = '--'; dn = '';              hvis = 'off';
        end
        plot(nl_p, ctr_p, ...
            'Color', st.color, 'LineStyle', ls, ...
            'Marker', st.marker, 'LineWidth', LW, 'MarkerSize', MS, ...
            'MarkerFaceColor', st.color, ...
            'DisplayName', dn, 'HandleVisibility', hvis);
    end
end
xlabel('noise ratio  (log scale; noise=0 omitted)');
ylabel(sprintf('%s E_2  (\\pm 1 std)', e2_metric));
title(sprintf('Parameter error E_2 (pooled across cells \\times realizations, %s)', ...
    e2_metric));
% --- Global power-law fit -------------------------------------------------
% Fit log10(E2_center) = m * log10(noise) + b using every (noise, center)
% point from every (model, raw/shrunk) curve.  This is what the eye is
% averaging when it reports "all the lines have the same slope".
if fit_powerlaw
    [m_fit, b_fit, slope_ci, r2_fit, n_pts, per_group] = ...
        fit_global_powerlaw(datasets, e2_metric);

    if ~isnan(m_fit)
        % Overlay reference line across the visible x range
        xl = xlim;
        x_ref = logspace(log10(xl(1)), log10(xl(2)), 100);
        man_power = 1;
        y_ref = 1.8*10^(b_fit) * x_ref.^man_power;
        plot(x_ref, y_ref, 'k--', 'LineWidth', 2.5, ...
            'DisplayName', ...
            sprintf('E_2 \\propto \\sigma^{%.2f}', man_power));

        % Console diagnostics
        fprintf('\n=== Global E_2 power-law fit (using %s center curves) ===\n', e2_metric);
        fprintf('  log10(E_2) = %.4f * log10(noise) + %.4f\n', m_fit, b_fit);
        fprintf('  slope  = %.4f   (95%% CI: [%.4f, %.4f])\n', ...
            m_fit, slope_ci(1), slope_ci(2));
        fprintf('  E_2    = %.4g * noise^{%.4f}\n', 10^b_fit, m_fit);
        fprintf('  R^2    = %.4f   (n = %d points)\n', r2_fit, n_pts);
        fprintf('  Per-group slopes (sanity check that they agree):\n');
        for jj = 1:numel(per_group)
            fprintf('    %-32s %-6s : %.3f   (n=%d)\n', ...
                per_group(jj).name, per_group(jj).which, ...
                per_group(jj).slope, per_group(jj).n);
        end
        fprintf('========================================================\n');
    end
end

legend('Location', 'best', 'FontSize', 9);


%% ----------------- Figure 3: Heterogeneity rel error -------------------
bio_params = {'r', 'S_max'};
n_bio      = numel(bio_params);

figure('Name','Heterogeneity rel error');
for bi = 1:n_bio
    pname    = bio_params{bi};
    skip_exp = strcmp(pname, 'S_max');

    ax = subplot(1, n_bio, bi);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    set(ax, 'XScale', 'log');

    % Style key anchors.
    plot(NaN, NaN, 'k-',  'LineWidth', LW, 'DisplayName', 'raw WENDy');
    plot(NaN, NaN, 'k--', 'LineWidth', LW, 'DisplayName', 'shrunk');

    any_data = false;
    for k = 1:n_files
        D = datasets{k};
        if skip_exp && strcmp(D.true_model, 'exponential')
            continue;   % no S_max for exponential -> no legend entry
        end
        st = styles.(D.true_model);
        nl = D.noise_levels;

        for tag = {'raw','shrunk'}
            which = tag{1};
            fn = sprintf('rel_diff_hetero_%s_%s_%s', which, het_agg, pname);
            if ~isfield(D.agg, fn), continue; end
            y = D.agg.(fn);
            ok = nl > 0 & isfinite(y);
            if ~any(ok), continue; end

            % Clip large positive blow-ups (WENDy collapsing toward an
            % exponential under heavy noise) so they don't blow out the
            % y-axis.  The y-tick at the cap is relabeled ">cap" below.
            y_p  = min(y(ok), het_cap);
            nl_p = nl(ok);
            any_data = true;

            if strcmp(which, 'raw')
                ls = '-';  dn = D.display_name;  hvis = 'on';
            else
                ls = '--'; dn = '';              hvis = 'off';
            end
            plot(nl_p, y_p, ...
                'Color', st.color, 'LineStyle', ls, ...
                'Marker', st.marker, 'LineWidth', LW, 'MarkerSize', MS, ...
                'MarkerFaceColor', st.color, ...
                'DisplayName', dn, 'HandleVisibility', hvis);
        end
    end

    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    xlabel('noise ratio');
    ylabel('E_\tau'); yscale("log");
    title(sprintf('%s', tex_param(pname)));
    if any_data
        apply_cap_yaxis(ax, het_cap);
    end
    legend('Location', 'best', 'FontSize', 8);
end
sgtitle(sprintf(['Heterogeneity recovery (%s across realizations)' newline ...
    '  values >%g clipped to the cap'], het_agg, het_cap));


%% =====================================================================
%  Helpers
%  =====================================================================
function S = build_style_map()
% Each model: unique color + unique marker.  Line style is reserved for
% raw vs shrunk (solid / dashed) and is NOT set here.
S = struct();
S.exponential  = struct('color', [0.00 0.45 0.74], 'marker', 'o');
S.logistic     = struct('color', [0.85 0.33 0.10], 'marker', 's');
S.gompertz     = struct('color', [0.93 0.69 0.13], 'marker', '^');
S.vb_linear    = struct('color', [0.49 0.18 0.56], 'marker', 'd');
S.vb_metabolic = struct('color', [0.47 0.67 0.19], 'marker', 'v');
end

function tag = normalize_tag(s)
s = lower(string(s));
switch s
    case {"exp", "exponential"},        tag = 'exponential';
    case {"vb", "vb_linear"},           tag = 'vb_linear';
    case {"logistic"},                  tag = 'logistic';
    case {"gomp", "gompertz"},          tag = 'gompertz';
    case {"vb23", "vb_metabolic"},      tag = 'vb_metabolic';
    otherwise,                          tag = char(s);
end
end

function n = display_name(tag)
switch tag
    case 'exponential',  n = 'Exponential';
    case 'logistic',     n = 'Logistic';
    case 'gompertz',     n = 'Gompertz';
    case 'vb_linear',    n = 'von Bertalanffy (linear)';
    case 'vb_metabolic', n = 'von Bertalanffy (metabolic)';
    otherwise,           n = tag;
end
end

function s = tex_param(p)
switch p
    case 'r',     s = 'r';
    case 'S_max', s = 'S_{max}';
    otherwise,    s = strrep(p, '_', '\_');
end
end

function apply_cap_yaxis(ax, cap)
% Cap the upper y-limit at `cap`, force a tick at `cap`, and relabel that
% tick as ">cap" so the reader knows clipping has occurred.  Also draws a
% light dotted reference line at y = cap.
yl = ylim(ax);
ylim(ax, [yl(1), cap * 1.05]);

yt = get(ax, 'YTick');
yt = yt(yt <= cap + 1e-9);
if ~any(abs(yt - cap) < 1e-9)
    yt = sort([yt, cap]);
end
ytl = arrayfun(@(v) num2str(v, '%g'), yt, 'UniformOutput', false);
idx = find(abs(yt - cap) < 1e-9, 1);
if ~isempty(idx)
    ytl{idx} = sprintf('>%g', cap);
end
set(ax, 'YTick', yt, 'YTickLabel', ytl);

yline(ax, cap, ':', 'Color', [0.45 0.45 0.45], ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
end

function [rs, ss] = pool_e2_stats(results_flat, noise_levels, n_realizations)
% Pool e2_raw_vec / e2_shrunk_vec across realizations at each noise level
% and compute mean, median, std on the pooled set.
n_noise = numel(noise_levels);
rs = struct('mean', nan(n_noise, 1), 'median', nan(n_noise, 1), 'std', nan(n_noise, 1));
ss = rs;
results = reshape(results_flat, [n_noise, n_realizations]);
for ni = 1:n_noise
    pool_r = []; pool_s = [];
    for ri = 1:n_realizations
        r = results{ni, ri};
        if ~isstruct(r) || isfield(r, 'failed'), continue; end
        if isfield(r, 'e2_raw_vec')
            pool_r = [pool_r, r.e2_raw_vec(:).']; %#ok<AGROW>
        end
        if isfield(r, 'e2_shrunk_vec')
            pool_s = [pool_s, r.e2_shrunk_vec(:).']; %#ok<AGROW>
        end
    end
    pr = pool_r(isfinite(pool_r));
    ps = pool_s(isfinite(pool_s));
    if ~isempty(pr)
        rs.mean(ni)   = mean(pr);
        rs.median(ni) = median(pr);
        rs.std(ni)    = std(pr);
    end
    if ~isempty(ps)
        ss.mean(ni)   = mean(ps);
        ss.median(ni) = median(ps);
        ss.std(ni)    = std(ps);
    end
end
end

function [m, b, slope_ci, r2, n, per_group] = fit_global_powerlaw(datasets, metric)
% Fit log10(E2_center) = m * log10(noise) + b using every (noise, center)
% point from every (model, raw/shrunk) curve in the dataset list.  Also
% fits each (model, raw/shrunk) curve individually and returns the slopes
% so the caller can sanity-check that the global slope is representative.
%
%   m         : global slope (power-law exponent)
%   b         : global intercept (so E2 = 10^b * noise^m)
%   slope_ci  : 1x2 95% CI for m from OLS
%   r2        : coefficient of determination of the global fit
%   n         : number of points used in the global fit
%   per_group : struct array with .name, .which, .slope, .n

    all_x = [];
    all_y = [];
    per_group = struct('name', {}, 'which', {}, 'slope', {}, 'n', {});

    for k = 1:numel(datasets)
        D  = datasets{k};
        nl = D.noise_levels;
        for tag = {'raw', 'shrunk'}
            which = tag{1};
            S2  = D.(['e2_' which]);
            ctr = S2.(metric);
            ok  = nl > 0 & isfinite(ctr) & ctr > 0;
            if ~any(ok), continue; end

            xg = log10(nl(ok));
            yg = log10(ctr(ok));
            all_x = [all_x; xg]; %#ok<AGROW>
            all_y = [all_y; yg]; %#ok<AGROW>

            % Per-group slope (only if enough points to fit a line)
            slope_g = NaN;
            if numel(xg) >= 2
                Ag = [xg, ones(size(xg))];
                bg = Ag \ yg;
                slope_g = bg(1);
            end
            per_group(end+1) = struct( ...
                'name',  D.display_name, ...
                'which', which, ...
                'slope', slope_g, ...
                'n',     numel(xg)); %#ok<AGROW>
        end
    end

    n = numel(all_x);
    if n < 3
        m = NaN; b = NaN; slope_ci = [NaN NaN]; r2 = NaN;
        return
    end

    A    = [all_x, ones(n, 1)];
    beta = A \ all_y;
    m = beta(1);
    b = beta(2);

    resid  = all_y - A * beta;
    dof    = n - 2;
    sigma2 = sum(resid.^2) / max(dof, 1);
    covB   = sigma2 * pinv(A' * A);
    se_m   = sqrt(max(covB(1, 1), 0));
    % 95% CI (z=1.96; close enough to t for moderate dof)
    slope_ci = [m - 1.96 * se_m, m + 1.96 * se_m];

    ss_tot = sum((all_y - mean(all_y)).^2);
    ss_res = sum(resid.^2);
    r2 = 1 - ss_res / max(ss_tot, eps);
end