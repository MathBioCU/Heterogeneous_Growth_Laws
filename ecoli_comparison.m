
clear; close all;
addpath(genpath('./wsindy_obj_base/'));
addpath(genpath('./Utils/'));

% =================== Load  ===================
manifest_file = "./Data/elife/ecoli_qc_manifest.mat";
if ~isfile(manifest_file)
    error('run_ecoli_model_selection:NoManifest', ...
          'Manifest %s not found. Run prune_ecoli_data.m first.', manifest_file);
end
load(char(manifest_file));
media      = MF.media;        % tag | display name | canonical <T_d> (min)
qc         = MF.qc;
qc_info    = MF.qc_info;
fname_used = MF.fname_used;
n_media    = size(media, 1);

fprintf('Loaded QC manifest %s (created %s)\n', manifest_file, string(MF.created));
if qc.enable
    fprintf('  QC: %s, threshold %.3g um, window %d, growth_mult %g, mode %s, trim [%g %g]\n', ...
            qc.direction, qc.abs_thresh, qc.window, qc.growth_mult, qc.mode, ...
            qc.trim_frac(1), qc.trim_frac(2));
else
    fprintf('  QC was disabled in STEP 1: running on the raw data.\n');
end

% =================== WENDy+Shrinkage options ===================
opts = struct();
opts.noise_ratio    = 0.0;
opts.want_shrinkage = true;
opts.want_linargs   = true;
opts.n_select       = 500;
opts.want_interp    = false;
opts.verbose        = true;
opts.rng_seed       = 42;     % reproducible subsample across runs
opts.min_length     = qc.min_length;   % same cut-off as the QC step
opts.force_model    = '';               % '' = let the criterion choose
opts.selection_score = 'sumtraj_keff';  % effective-dof BIC (see header)

% =================== Run WENDy per medium ===================
results = cell(n_media, 1);

for m = 1:n_media
    nice = media{m, 2};
    Td   = media{m, 3};

    fname_m = fname_used(m);
    if strlength(fname_m) == 0 || ~isfile(fname_m)
        warning('Skipping %s (no data file in manifest or file missing: %s)', ...
                nice, fname_m);
        continue
    end

    fprintf('\n============== %s  (<T_d> = %d min) ==============\n', ...
            nice, Td);
    fprintf('  data: %s\n', fname_m);

    R = run_growth_pipeline(char(fname_m), opts);
    R.tag          = media{m, 1};
    R.nice         = nice;
    R.Td_canonical = Td;
    R.mu_ref       = log(2) / Td;
    R.qc           = qc_info{m};
    results{m}     = R;
end

% =================== comparison figure ===================
palette = [ ...
    0.00 0.45 0.70;   % blue   - alanine     (slow)
    0.90 0.60 0.00;   % orange - glycerol    (medium)
    0.00 0.62 0.45  ];  % green  - glucose-cas (fast)

fig = figure('Name', 'E. coli cross-media comparison', ...
             'Position', [100 100 1200 760], 'Color', 'w');
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Row 1: example trajectories per medium ---
n_show = 25;
for m = 1:n_media
    nexttile(tl, m); hold on; box on;
    if isempty(results{m})
        title(media{m,2});  continue
    end
    R = results{m};
    npick = min(n_show, R.N);
    rng(0);
    idx = randperm(R.N, npick);
    for j = 1:npick
        i   = idx(j);
        t   = R.t{i};
        x   = R.tracks{i};
        age = t - t(1);
        plot(age, x, 'Color', [palette(m,:) 0.35], 'LineWidth', 0.8);
    end
    xlabel('age (min)');
    ylabel('length (\mum)');
    title(sprintf('%s  (\\langle T_d\\rangle = %d min, n = %d, %s)', ...
                  R.nice, R.Td_canonical, R.N, R.model_id), ...
          'FontWeight', 'normal', 'FontSize', 11);
    set(gca, 'FontSize', 10);
end

% --- Row 2, panel A: KDE overlay of shrunken r distributions ---
nexttile(tl, 4); hold on; box on; grid on;

r_all = [];
for m = 1:n_media
    if isempty(results{m}), continue; end
    rv = results{m}.params_bio_shrunk.r;
    r_all = [r_all; rv(isfinite(rv))];
end
r_lo   = max(min(r_all) * 0.5, 1e-5);
r_hi   = max(r_all) * 1.5;
r_grid = logspace(log10(r_lo), log10(r_hi), 400)';

for m = 1:n_media
    if isempty(results{m}), continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    rv = R.params_bio_shrunk.r;
    rv = rv(isfinite(rv));
    bw = 0.06 * R.mu_hat(fi);
    if bw <= 0, bw = 0.05 * std(rv); end
    pdf = ksdensity(rv, r_grid, 'Bandwidth', bw);
    plot(r_grid, pdf, '-', 'Color', palette(m,:), 'LineWidth', 2.0, ...
         'DisplayName', R.nice);
    xline(R.mu_ref, '--', 'Color', palette(m,:), 'Alpha', 0.6, ...
          'LineWidth', 1.0, 'HandleVisibility', 'off');
end
xlabel('growth rate r  (min^{-1})');
ylabel('density (shrunk)');
title('Growth-rate distributions', 'FontWeight', 'normal');
legend('Location', 'northwest', 'FontSize', 9, 'Box', 'off');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);
text(0.98, 0.04, 'dashed: ln 2 / \langle T_d\rangle', ...
     'Units', 'normalized', 'HorizontalAlignment', 'right', ...
     'FontSize', 8, 'Color', [0.4 0.4 0.4]);

% --- Row 2, panel B: forest plot of (mu_hat +/- tau_hat) vs canonical ---
nexttile(tl, 5); hold on; box on; grid on;

y_pos = n_media:-1:1;
labels_y = cell(n_media, 1);
for m = 1:n_media
    if isempty(results{m})
        labels_y{m} = media{m,2};  continue
    end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    mu = R.mu_hat(fi);
    ta = R.tau_hat(fi);

    % Filled marker = mu_hat, error bars = +/- tau_hat (population spread)
    errorbar(mu, y_pos(m), [], [], ta, ta, ...
             'o', 'Color', palette(m,:), 'MarkerFaceColor', palette(m,:), ...
             'LineWidth', 1.5, 'MarkerSize', 9, 'CapSize', 10);

    % Cross = canonical ln(2) / <T_d> reference
    plot(R.mu_ref, y_pos(m), 'x', 'Color', [0.2 0.2 0.2], ...
         'MarkerSize', 11, 'LineWidth', 1.5);

    labels_y{m} = R.nice;
end
yticks(1:n_media);
yticklabels(flip(labels_y));
ylim([0.3, n_media + 0.7]);
xlabel('growth rate  (min^{-1})');
title('$\hat\mu_r \pm \hat\tau_r$ vs.\ canonical $\ln 2/\langle T_d\rangle$', ...
      'Interpreter', 'latex', 'FontWeight', 'normal');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);
text(0.02, 0.07, '\bullet  \mu_r \pm \tau_r     \times  ln 2 / \langle T_d\rangle', ...
     'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% --- Row 2, panel C: biological CV vs mean growth rate ---
nexttile(tl, 6); hold on; box on; grid on;

for m = 1:n_media
    if isempty(results{m}), continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    mu = R.mu_hat(fi);
    ta = R.tau_hat(fi);
    cv = ta / mu;
    plot(mu, cv, 'o', ...
         'Color', palette(m,:), 'MarkerFaceColor', palette(m,:), ...
         'MarkerSize', 12, 'LineWidth', 1.0);
    text(mu * 1.05, cv, R.nice, ...
         'FontSize', 9, 'VerticalAlignment', 'middle');
end
xlabel('population mean  \mu_r  (min^{-1})');
ylabel('biological CV  \tau_r / \mu_r');
title('Heterogeneity vs. mean growth rate', 'FontWeight', 'normal');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);

% --- Overall title ---
sgtitle('Single-cell {\it E. coli} growth across carbon sources  (Kar et al. 2021 data)', ...
        'FontSize', 13, 'FontWeight', 'bold');

% =================== Summary table to stdout ===================
fprintf('\n\n=== Cross-media summary ===\n');
fprintf('%-12s  %5s  %14s  %12s  %12s  %12s  %6s\n', ...
        'Medium', 'N', 'selected model', 'mu_hat', 'ln2/Td', 'tau_hat', 'CV');
fprintf('%-12s  %5s  %14s  %12s  %12s  %12s  %6s\n', ...
        repmat('-',1,12), '-----', repmat('-',1,14), repmat('-',1,12), ...
        repmat('-',1,12), repmat('-',1,12), '------');
for m = 1:n_media
    if isempty(results{m}), continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    mu = R.mu_hat(fi);
    ta = R.tau_hat(fi);
    fprintf('%-12s  %5d  %14s  %12.4e  %12.4e  %12.4e  %6.3f\n', ...
            R.nice, R.N, R.model_id, mu, R.mu_ref, ta, ta/mu);
end
fprintf('==========================\n\n');


%% plot results only
palette = [0.00 0.45 0.70;     % alanine
           0.90 0.60 0.00;     % glycerol
           0.00 0.62 0.45];    % glucose-cas
n_media = numel(results);
 
fig = figure('Position', [100 100 1300 380], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
 
% ---- Panel A: KDE overlay of shrunken r distributions ----
ax_A = nexttile(tl, 1); hold on; box on; grid on;
r_all = [];
for m = 1:n_media
    if isempty(results{m}), continue; end
    rv = results{m}.params_bio_shrunk.r;
    r_all = [r_all; rv(isfinite(rv))];
end
r_grid = logspace(log10(min(r_all)*0.5), log10(max(r_all)*1.5), 400)';
for m = 1:n_media
    if isempty(results{m}), continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    rv = R.params_bio_shrunk.r;  rv = rv(isfinite(rv));
    bw = max(0.06 * R.mu_hat(fi), 0.05 * std(rv));
    pdf = ksdensity(rv, r_grid, 'Bandwidth', bw);
    plot(r_grid, pdf, '-', 'Color', palette(m,:), 'LineWidth', 2.0, ...
         'DisplayName', R.nice);
    xline(R.mu_ref, '--', 'Color', palette(m,:), 'Alpha', 0.6, ...
          'LineWidth', 1.0, 'HandleVisibility', 'off');
end
xlabel('growth rate r  (min^{-1})');
ylabel('density');
legend('Location', 'northeast', 'FontSize', 9, 'Box', 'off');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);
 
% ---- Panel B: forest plot of mu_hat +/- tau_hat ----
ax_B = nexttile(tl, 2); hold on; box on; grid on;
y_pos = n_media:-1:1;  labels_y = cell(n_media, 1);
for m = 1:n_media
    if isempty(results{m}), labels_y{m} = '';  continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    errorbar(R.mu_hat(fi), y_pos(m), [], [], R.tau_hat(fi), R.tau_hat(fi), ...
             'o', 'Color', palette(m,:), 'MarkerFaceColor', palette(m,:), ...
             'LineWidth', 1.5, 'MarkerSize', 9, 'CapSize', 10);
    plot(R.mu_ref, y_pos(m), 'x', 'Color', [0.2 0.2 0.2], ...
         'MarkerSize', 11, 'LineWidth', 1.5);
    labels_y{m} = R.nice;
end
yticks(1:n_media);  yticklabels(flip(labels_y));
ylim([0.3, n_media + 0.7]);
xlabel('growth rate  (min^{-1})');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);
text(0.02, 0.07, '\bullet  \mu_r \pm \tau_r     \times  ln 2 / \langle T_d\rangle', ...
     'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);
 
% ---- Panel C: CV vs mean rate ----
ax_C = nexttile(tl, 3); hold on; box on; grid on;
for m = 1:n_media
    if isempty(results{m}), continue; end
    R  = results{m};
    fi = find(strcmp(R.shrink_names{1}, 'r'), 1);
    mu = R.mu_hat(fi);  ta = R.tau_hat(fi);
    plot(mu, ta/mu, 'o', 'Color', palette(m,:), ...
         'MarkerFaceColor', palette(m,:), 'MarkerSize', 12);
    text(mu*1.05, ta/mu, R.nice, 'FontSize', 9, 'VerticalAlignment', 'middle');
end
xlabel('population mean  \mu_r  (min^{-1})');
ylabel('biological CV  \tau_r / \mu_r');
set(gca, 'XScale', 'log', 'FontSize', 10, 'GridAlpha', 0.15);
 
% ---- A / B / C panel labels above each axis ----
for h = [ax_A, ax_B, ax_C]
    set(h, 'Units', 'normalized');
end
text(ax_A, -0.12, 1.05, 'A', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'FontSize', 16, 'VerticalAlignment', 'bottom');
text(ax_B, -0.12, 1.05, 'B', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'FontSize', 16, 'VerticalAlignment', 'bottom');
text(ax_C, -0.12, 1.05, 'C', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'FontSize', 16, 'VerticalAlignment', 'bottom');
 

% exportgraphics(fig, 'ecoli_cross_media_bottom.pdf', 'ContentType', 'vector');

