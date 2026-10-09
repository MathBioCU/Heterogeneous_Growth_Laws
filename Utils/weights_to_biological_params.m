function params = weights_to_biological_params(weights, model_tag, varargin)
% WEIGHTS_TO_BIOLOGICAL_PARAMS  Convert WSINDy/WENDy weights to
%   interpretable biological quantities for standard growth laws.
%
%   params = weights_to_biological_params(weights, model_tag)
%   params = weights_to_biological_params(weights, model_tag, 'Name', Value, ...)
%   params = weights_to_biological_params(R, [])   % struct shorthand;
%                                                  % uses R.model_tag
%
%   INPUTS
%     weights   : [N x p] matrix of per-cell weight vectors, OR a struct
%                 produced by run_wsindy_experiment_v2 / _bagged_stage1
%                 (in which case model_tag may be omitted and is pulled
%                 from the struct).  Also accepts a 1 x p row.
%     model_tag : string / char -- one of "exp", "logistic", "vB",
%                 "vB23", "gomp" / "gompertz".  Case-insensitive.
%                 For anything else, returns an empty struct with a
%                 warning.
%
%   NAME-VALUE
%     'BasisTags'  [1 x p] string array of basis-term tags matching the
%                  columns of `weights`.  Required when `weights` comes
%                  from the decomp library (columns are {const, x, x^2,
%                  x^(2/3), x*log(x)}, only a subset of which may be
%                  active).  Ignored when the model is given in its
%                  canonical two-weight form (e.g. v1's WENDy output).
%     'Verbose'    logical, default true.  Print a summary table.
%
%   OUTPUT
%     params : struct with model-dependent fields.  All fields are
%              numeric column vectors of length N (one entry per cell),
%              with NaN for cells where the parameter is unphysical
%              (e.g. negative growth rate, or no finite carrying
%              capacity).  Structure:
%
%       model_tag        : echoed tag
%       model_rhs_form   : human-readable RHS in w_1, w_2 convention
%       valid            : [N x 1] logical, true when derived params
%                          are physically meaningful
%       <model-specific fields below>
%
%   MODEL CONVENTIONS (match v1 true_growth_rhs)
%     exp       dx/dt = w_1 * x
%         fields:  mu, doubling_time
%     logistic  dx/dt = w_1*x + w_2*x^2
%         fields:  r, K  (r = w_1 intrinsic rate; K = -w_1/w_2)
%     vB        dx/dt = w_1 + w_2*x
%         fields:  k, S_max (k = -w_2; S_max = -w_1/w_2)
%     vB23      dx/dt = w_1*x^(2/3) + w_2*x
%         fields:  a, b, S_max (a = w_1; b = -w_2; S_max = (a/b)^3)
%     gomp      dx/dt = w_1*x*log(x) + w_2*x
%         fields:  b, K, r0  (b = -w_1; K = exp(-w_2/w_1);
%                             r0 = w_2 = effective initial rate)
%
%   EXAMPLES
%     % Typical use with run_wsindy_experiment output
%     R = run_wsindy_experiment(Data);
%     params = weights_to_biological_params(R.params, R.model_tag);
%
%     % With decomp-library output from v2
%     R = run_wsindy_experiment_v2(Data);
%     params = weights_to_biological_params(R.params, R.model_tag, ...
%                 'BasisTags', R.lib_tags(R.sel_support));
%
%   See also RUN_WSINDY_EXPERIMENT, RUN_WSINDY_EXPERIMENT_V2,
%            RUN_WSINDY_BAGGED_STAGE1.

%% ---- Parse ----
p = inputParser;
p.addParameter('BasisTags', strings(0));
p.addParameter('Verbose', false);
p.parse(varargin{:});
basis_tags = p.Results.BasisTags;
verbose    = p.Results.Verbose;


% Allow a struct-form input as shorthand
if isstruct(weights)
    R = weights;
    if nargin < 2 || isempty(model_tag) || (isstring(model_tag) && strlength(model_tag)==0)
        if isfield(R, 'model_tag'), model_tag = R.model_tag;
        else, error('model_tag missing from struct input.'); end
    end
    if isfield(R, 'params')
        W = R.params;
    elseif isfield(R, 'w_shrunk')
        W = R.w_shrunk;
    elseif isfield(R, 'wendy_pooled') && isfield(R.wendy_pooled, 'W')
        W = R.wendy_pooled.W;
    else
        error('Cannot find a per-cell weights field in the input struct.');
    end
    if isempty(basis_tags) && isfield(R, 'lib_tags') && isfield(R, 'sel_support')
        basis_tags = R.lib_tags(R.sel_support);
    end
else
    W = weights;
end


[N, p_w] = size(W);
tag = lower(string(model_tag));

%% ---- Map decomp-library columns to canonical (w_1, w_2) for the tag ----
[w1, w2, canonical_form] = map_to_canonical(W, tag, basis_tags);



%% ---- Derive biological parameters ----
params = struct();
params.model_tag      = string(model_tag);
params.model_rhs_form = canonical_form;

switch tag
    case {"exp","exponential"}
        mu = w1;
        params.r            = mu;
        params.valid         = mu > 0;
        params.S_max = inf .* ones(size(mu));

    case "logistic"
        r = w1;
        K = -w1 ./ w2;
        K(w2 >= 0 | w1 <= 0) = NaN;
        params.r     = r;
        params.S_max    = K;
        params.valid = (w1 > 0) & (w2 < 0);

    case {"vb","vb_linear","vB"}
        k     = -w2;
        S_max = -w1 ./ w2;
        S_max(w2 >= 0 | w1 <= 0) = NaN;
        params.r     = k;
        params.S_max = S_max;
        params.valid = (w1 > 0) & (w2 < 0);

    case {"vb23","vB23","vb_metabolic"}
        a     = w2;
        b     = -w1;
        S_max = (-w2 ./ w1).^3;
        % S_max(b <= 0 | a <= 0) = NaN;
        params.a     = a;
        params.r     = b;
        params.S_max = S_max;
        params.valid = (w1 < 0) & (w2 > 0);

    case {"gomp", "gompertz"}
        b = -w1;
        % K = exp(-w_2/w_1) = exp(w_2/b)
        K = exp(-w1 ./ w2);
        % K(-w2 >= 0) = NaN;   % growth rate must be positive
        params.r     = -w2;
        params.S_max = K;
        params.valid = (w2 < 0);

    otherwise
        if verbose
            warning(['Non-standard model tag "%s".  Returning empty ', ...
                     'parameter struct.  Supported tags: exp, logistic, ', ...
                     'vB, vB23, gomp/gompertz.'], model_tag);
        end
        params = struct('model_tag', string(model_tag), 'note', ...
                        'non-standard model, no interpretable map');
        return
end

%% ---- Summary ----
if verbose
    print_summary(params, tag, N);
end
end


%% ====================================================================
%% =====================   HELPER: COEFFICIENT MAP   ==================
%% ====================================================================
function [w1, w2, form] = map_to_canonical(W, tag, basis_tags)
% Return the canonical (w_1, w_2) vectors for the given model tag from
% the raw WSINDy/WENDy output `W`.
%
% There are two layouts:
%   A) Canonical two-weight layout from named-model library (v1 / v2
%      named mode).  Columns of W directly correspond to (w_1, w_2) in
%      the growth_type-specific ordering used by v1's true_growth_rhs.
%   B) Decomp-library layout, where columns of W correspond to a
%      subset of {const, x, x^2, x^(2/3), x*log(x)} identified by
%      basis_tags.  Requires basis_tags to match size(W, 2).
N = size(W, 1);
p = size(W, 2);
w1 = nan(N, 1);
w2 = nan(N, 1);
form = "";

% Decide layout.  If basis_tags has exactly as many entries as W has
% columns, use the decomp layout; otherwise assume canonical.
use_decomp = ~isempty(basis_tags) && numel(basis_tags) == p;

if use_decomp
    bt = string(basis_tags(:)');
    idx_const = find(bt == "const",    1);
    idx_x     = find(bt == "x",        1);
    idx_x2    = find(bt == "x^2",      1);
    idx_x23   = find(bt == "x^(2/3)",  1);
    idx_xln   = find(bt == "x*log(x)", 1);

    switch lower(string(tag))
        case {"exp","exponential"}
            if ~isempty(idx_x), w1 = W(:); end
            form = "dx/dt = w_1 * x";
        case "logistic"
            if ~isempty(idx_x),  w1 = W(:, idx_x);  end
            if ~isempty(idx_x2), w2 = W(:, idx_x2); end
            form = "dx/dt = w_1 * x + w_2 * x^2";
        case {"vb","vb_linear","vB"}
            if ~isempty(idx_const), w1 = W(:, idx_const); end
            if ~isempty(idx_x),     w2 = W(:, idx_x);     end
            form = "dx/dt = w_1 + w_2 * x";
        case {"vb23","vB23","vb_metabolic"}
            if ~isempty(idx_x23), w1 = W(:, idx_x23); end
            if ~isempty(idx_x),   w2 = W(:, idx_x);   end
            form = "dx/dt = w_1 * x + w_2 * x^(2/3)";
        case {"gomp", "gompertz"}
            if ~isempty(idx_xln), w2 = W(:, idx_xln); end
            if ~isempty(idx_x),   w1 = W(:, idx_x);   end
            form = "dx/dt = w_1 * x + w_2 * x*log(x)";
    end
else
    % Canonical two-weight layout (matches v1 true_growth_rhs).
    switch lower(string(tag))
        case {"exp","exponential"}
            if p >= 1, w1 = W(:, 1); end
            form = "dx/dt = w_1 * x";
        case "logistic"
            if p >= 2
                w1 = W(:, 1); w2 = W(:, 2);
            end
            form = "dx/dt = w_1 * x + w_2 * x^2";
        case {"vb","vb_linear","vB"}
            if p >= 2
                w1 = W(:, 1); w2 = W(:, 2);
            end
            form = "dx/dt = w_1 + w_2 * x";
        case {"vb23","vB23","vb_metabolic"}
            if p >= 2
                w1 = W(:, 1); w2 = W(:, 2);
            end
            form = "dx/dt = w_1 * x + w_2 * x^(2/3)";
        case {"gomp", "gompertz"}
            if p >= 2
                w1 = W(:, 1); w2 = W(:, 2);
            end
            form = "dx/dt = w_1 * x + w_2 * x*log(x)";
    end
end

if all(isnan(w1)) && all(isnan(w2))
    warning(['Could not map weights to canonical form for model "%s".  ', ...
             'Check that BasisTags matches size(weights, 2) or that ', ...
             'weights are in the canonical two-weight layout.'], tag);
end
end


%% ====================================================================
%% =========================   HELPER: SUMMARY   ======================
%% ====================================================================
function print_summary(params, tag, N)
    fprintf('\n  Biological parameters  (model = %s, N = %d cells):\n', tag, N);
    fprintf('  %s\n', params.model_rhs_form);
    nv = sum(params.valid);
    fprintf('  Physically valid: %d / %d (%.1f%%)\n', nv, N, 100*nv/N);

    flds = setdiff(fieldnames(params), ...
                   {'model_tag','model_rhs_form','valid'});
    fprintf('  %-14s %10s %10s %10s %10s\n', ...
            'field', 'median', 'mean', '5%', '95%');
    for f = 1:numel(flds)
        v = params.(flds{f});
        if ~isnumeric(v), continue, end
        v = v(params.valid);
        if isempty(v)
            fprintf('  %-14s       --         --         --         --\n', flds{f});
            continue
        end
        fprintf('  %-14s %10.4g %10.4g %10.4g %10.4g\n', flds{f}, ...
                median(v,'omitnan'), mean(v,'omitnan'), ...
                quantile(v, 0.05), quantile(v, 0.95));
    end
    fprintf('\n');
end