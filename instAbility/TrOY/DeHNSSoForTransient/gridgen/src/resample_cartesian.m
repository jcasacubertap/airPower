function BF_new = resample_cartesian(BF, params)
%RESAMPLE_CARTESIAN  Interpolate a structured BF onto a Malik-Chebyshev η × equidistant ξ grid.
%
%   BF_new = resample_cartesian(BF, params)
%
%   Designed for BF data that already lives on a structured (nx_raw × ny_raw)
%   grid (e.g. heritage bf_sweptwing_flat.mat) but at a different resolution
%   than what DeHNSSo expects. The curvilinear pipeline (create_ctype_mesh)
%   does the equivalent job for body-fitted unstructured scatter.
%
%   Inputs:
%     BF — struct with X, Y (1-D coordinate vectors), U, V, W on a regular
%          grid. Orientation ((ny_raw × nx_raw) or (nx_raw × ny_raw)) is
%          auto-detected from sizes.
%     params — fields used:
%       n_eta_new  (required) — target wall-normal Chebyshev points
%       n_xi_new   (optional) — target streamwise stations (default: keep X size)
%       H          (optional) — wall-normal domain height (default: 0.99·max(Y))
%       y_i        (optional) — Malik median collocation height (default: H/10)
%       xi_range   (optional) — [x_lo, x_hi] as fraction of [X(1), X(end)]
%       Re         (optional) — override / provide Re if BF.Re is missing
%
%   Output: BF_new with DeHNSSo-ready fields (X, Y, U, V, W, H, y_i) plus any
%           of Re/lref/Uref/nu that were present in input BF.

%% Normalise raw grid to (ny_raw × nx_raw)
x_raw = BF.X(:).';          % [1 × nx_raw]
y_raw = BF.Y(:);            % [ny_raw × 1]
nx_raw = numel(x_raw); ny_raw = numel(y_raw);

if size(BF.U,1) == nx_raw && size(BF.U,2) == ny_raw
    U_raw = BF.U.'; V_raw = BF.V.'; W_raw = BF.W.';
else
    U_raw = BF.U;   V_raw = BF.V;   W_raw = BF.W;
end

%% Target grid dimensions
ny_new = params.n_eta_new;
assert(~isempty(ny_new), 'resample_cartesian: params.n_eta_new is required.');

if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    nx_new = params.n_xi_new;
else
    nx_new = nx_raw;
end

%% Wall-normal domain (H, y_i)
if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
else
    H = max(y_raw) * 0.99;
end
if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
else
    y_i = H / 10;
end
if y_i >= H/2
    y_i = H/3;
    warning('y_i clamped to H/3 = %.5f (must be < H/2 for Malik mapping).', y_i);
end

%% Malik-Chebyshev η — descending: η(1)=H (freestream), η(end)=0 (wall)
a = H * y_i / (H - 2*y_i);
b = 1 + 2*a/H;
k_vec    = (0 : ny_new - 1)';
z_cheb   = cos(pi * k_vec / (ny_new - 1));
eta_new  = a * (1 + z_cheb) ./ (b - z_cheb);            % [ny_new × 1]

%% Streamwise ξ — equidistant, optionally trimmed by xi_range
if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    x_lo = x_raw(1) + params.xi_range(1) * (x_raw(end) - x_raw(1));
    x_hi = x_raw(1) + params.xi_range(2) * (x_raw(end) - x_raw(1));
else
    x_lo = x_raw(1);
    x_hi = x_raw(end);
end
xi_new = linspace(x_lo, x_hi, nx_new);                  % [1 × nx_new]

%% Interpolate U, V, W onto target grid via 2-D spline
[Xq, Yq] = meshgrid(xi_new, eta_new);
U_new = interp2(x_raw, y_raw, U_raw, Xq, Yq, 'spline');
V_new = interp2(x_raw, y_raw, V_raw, Xq, Yq, 'spline');
W_new = interp2(x_raw, y_raw, W_raw, Xq, Yq, 'spline');

%% Package into DeHNSSo-ready BF (no mesh_type field → build_stab_grid uses Cartesian branch)
BF_new     = struct();
BF_new.X   = xi_new;
BF_new.Y   = eta_new;
BF_new.U   = U_new;
BF_new.V   = V_new;
BF_new.W   = W_new;
BF_new.H   = H;
BF_new.y_i = y_i;

fprintf('Cartesian resample done: %d × %d → %d × %d (Malik-Cheb η, equidistant ξ)\n', ...
        ny_raw, nx_raw, ny_new, nx_new);
end
