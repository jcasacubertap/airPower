function BF_new = resample_curvilinear(BF, params)
%RESAMPLE_CURVILINEAR  Interpolate a structured body-fitted BF onto a Malik-Chebyshev grid.
%
%   BF_new = resample_curvilinear(BF, params)
%
%   For BF data on a curvilinear *structured* mesh where X and Y are 2-D
%   arrays of the same shape, but Y is not constant along a row (e.g. a
%   hump-following grid such as BF_SweptWing_hump.mat). The wall is
%   trivially the last row of the grid — no polygon detection is needed.
%
%   Required inputs (params):
%     n_eta_new  — target wall-normal Chebyshev points
%
%   Optional (params):
%     n_xi_new   — target streamwise stations           (default: keep nx_raw)
%     H          — wall-normal domain height            (default: 0.99·max column arc-length)
%     y_i        — Malik median collocation height      (default: H/10)
%     xi_range   — [lo, hi] as fraction of wall arc-length (default: full)
%     xi_trim_inflow / xi_trim_outflow — stations to drop from each end
%                 (integer count ≥1, or fraction 0<f<1)
%     rescale    — if true, divide X,Y,U,V,W by lref/Uref (default false).
%                  Only set true for dimensional input. Heritage BFs are
%                  already non-dimensional — keep rescale=false for those.
%     Re, lref, Uref, nu — propagated into BF_new via ensure_metadata.
%
%   Output: BF_new compatible with build_stab_grid (mesh_type='curvilinear').

%% --- 1. Normalise orientation to [ny × nx], wall at row end --------------

[X, Y, U, V, W, ny_raw, nx_raw] = orient_structured_grid(BF);

if isfield(params,'rescale') && ~isempty(params.rescale) && params.rescale
    assert(isfield(params,'lref') && ~isempty(params.lref) && ...
           isfield(params,'Uref') && ~isempty(params.Uref), ...
        'resample_curvilinear: rescale=true requires params.lref and params.Uref.');
    X = X / params.lref;  Y = Y / params.lref;
    U = U / params.Uref;  V = V / params.Uref;  W = W / params.Uref;
    fprintf('Rescaled input by Uref = %.4f, lref = %.6e\n', params.Uref, params.lref);
end

fprintf('Structured curvilinear input: %d × %d (ny × nx)\n', ny_raw, nx_raw);
fprintf('  X ∈ [%.4f, %.4f],  Y ∈ [%.4f, %.4f]\n', ...
    min(X(:)), max(X(:)), min(Y(:)), max(Y(:)));

%% --- 2. Wall = last row (DeHNSSo convention: eta(end)=0) -----------------

X_wall = X(end, :);                % [1 × nx_raw]
Y_wall = Y(end, :);

%% --- 3. Column arc-length from wall -------------------------------------
% For each column, compute distance from the wall point (last entry) to each
% grid point along the column. Works for vertical columns and wall-normal
% columns alike.

d_raw = zeros(ny_raw, nx_raw);          % d_raw(1,:)=top, d_raw(end,:)=0 (wall)
for j = 1:nx_raw
    dx   = diff(X(:, j));
    dy   = diff(Y(:, j));
    ds   = sqrt(dx.^2 + dy.^2);
    scol = [0; cumsum(ds)];              % ascending from row 1 downward
    d_raw(:, j) = scol(end) - scol;      % flip: now d=0 at wall, d=max at top
end

%% --- 4. Wall-normal domain height H and median-collocation y_i ----------

if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
    fprintf('Using prescribed H = %.5f\n', H);
else
    H = 0.99 * median(max(d_raw, [], 1));
    fprintf('Estimated H = %.5f  (0.99 × median column length)\n', H);
end

if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
else
    y_i = H / 10;
end
if y_i >= H / 2
    y_i = H / 3;
    warning('y_i clamped to H/3 = %.5f (must be < H/2 for Malik mapping).', y_i);
end
fprintf('Using y_i = %.5f\n', y_i);

%% --- 5. Malik-Chebyshev η target distribution ---------------------------
% Descending convention: η(1)=H (top), η(end)=0 (wall).

a = H * y_i / (H - 2*y_i);
b = 1 + 2 * a / H;
k_vec    = (0 : params.n_eta_new - 1)';
z_cheb   = cos(pi * k_vec / (params.n_eta_new - 1));
eta_cheb = a * (1 + z_cheb) ./ (b - z_cheb);          % [n_eta_new × 1]
n_eta_new = numel(eta_cheb);

%% --- 6. ξ stations along the wall ---------------------------------------
% Default: uniform wall arc-length (linspace).
% Optional: Gaussian-weighted clustering around one or more hump centres,
% generalising heritage's "xrefined" to multiple humps.
%
%     params.mug : scalar or vector of hump centres (physical x)
%     params.sig : Gaussian width per hump (fraction of half-domain);
%                  scalar broadcasts to all humps
%     params.ag  : step size at each peak (<1 → denser, >1 → sparser);
%                  scalar broadcasts to all humps
%
% Single hump with (mug = x_m, sig = 0.2, ag = 0.5) reproduces heritage.
% Empty mug or ag = 1 for all humps → uniform spacing (same as before).

ds_wall = sqrt(diff(X_wall).^2 + diff(Y_wall).^2);
s_wall  = [0, cumsum(ds_wall)];                       % [1 × nx_raw]

if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    n_xi_new = params.n_xi_new;
else
    n_xi_new = nx_raw;
end

s_lo = 0;
s_hi = s_wall(end);
if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    s_lo = params.xi_range(1) * s_wall(end);
    s_hi = params.xi_range(2) * s_wall(end);
end

mug = get_opt(params, 'mug', []);
sig = get_opt(params, 'sig', 0.2);
ag  = get_opt(params, 'ag',  1);

if isempty(mug) || all(ag(:) == 1)
    % Uniform wall-arc-length distribution
    xi_new = linspace(s_lo, s_hi, n_xi_new);
    fprintf('ξ distribution: uniform wall arc-length (no refinement)\n');
else
    % Gaussian clustering around one or more hump centres
    mug = mug(:).';
    if isscalar(sig);  sig = sig * ones(size(mug));  end
    if isscalar(ag);   ag  = ag  * ones(size(mug));  end
    assert(numel(sig) == numel(mug) && numel(ag) == numel(mug), ...
        'params.mug/sig/ag must be scalars or same-length vectors.');

    % Map each hump centre (physical x) into the normalised [-1, 1] axis,
    % consistent with heritage grid_gen.m (line 95).
    x_lo_p = interp1(s_wall, X_wall, s_lo, 'pchip');
    x_hi_p = interp1(s_wall, X_wall, s_hi, 'pchip');
    mug_n  = 2 * (mug - x_lo_p) / (x_hi_p - x_lo_p) - 1;

    % Build the local step-size modulation fg(ξ) by taking, at each point,
    % the strongest reduction contributed by any single hump. Using max
    % (rather than sum) keeps fg ∈ [min(ag), 1] even for overlapping humps.
    XC        = linspace(-1, 1, n_xi_new);
    fg_reduce = zeros(size(XC));
    for k = 1:numel(mug)
        G_k       = exp(-0.5 * ((XC - mug_n(k)) / sig(k)).^2);   % peak = 1
        fg_reduce = max(fg_reduce, (1 - ag(k)) * G_k);
    end
    fg = 1 - fg_reduce;                                          % fg ∈ [min(ag), 1]

    % Accumulate variable step sizes into non-uniform stations, then
    % normalise back to the arc-length interval [s_lo, s_hi].
    xi_cumsum = [0, cumsum(fg(2:end))];
    xi_cumsum = xi_cumsum / xi_cumsum(end);
    xi_new    = s_lo + xi_cumsum * (s_hi - s_lo);

    % Summary line, one per hump
    fprintf('ξ distribution: Gaussian clustering, %d hump%s\n', ...
            numel(mug), repmat('s', 1, numel(mug) > 1));
    for k = 1:numel(mug)
        fprintf('    hump %d: x_c = %-8.3f  σ = %-5.3f  a_g = %-5.3f  (Δx ratio = %.2f×)\n', ...
                k, mug(k), sig(k), ag(k), 1 / ag(k));
    end
end

%% --- 7. Resample η per raw column, then resample ξ ----------------------

% Pass 1: η direction (per raw column), raw grid → Malik-Cheb η.
U_mid = zeros(n_eta_new, nx_raw);
V_mid = zeros(n_eta_new, nx_raw);
W_mid = zeros(n_eta_new, nx_raw);
X_mid = zeros(n_eta_new, nx_raw);
Y_mid = zeros(n_eta_new, nx_raw);

for j = 1:nx_raw
    eta_col       = d_raw(:, j);                      % 0 at wall → max at top
    [eta_asc, o]  = sort(eta_col);                    % ascending for interp1
    eta_q         = max(min(eta_cheb, eta_asc(end)), eta_asc(1));

    U_mid(:, j) = interp1(eta_asc, U(o, j), eta_q, 'pchip');
    V_mid(:, j) = interp1(eta_asc, V(o, j), eta_q, 'pchip');
    W_mid(:, j) = interp1(eta_asc, W(o, j), eta_q, 'pchip');
    X_mid(:, j) = interp1(eta_asc, X(o, j), eta_q, 'pchip');
    Y_mid(:, j) = interp1(eta_asc, Y(o, j), eta_q, 'pchip');
end

% Pass 2: ξ direction (per η level), raw → target wall arc-length.
U_new = zeros(n_eta_new, n_xi_new);
V_new = zeros(n_eta_new, n_xi_new);
W_new = zeros(n_eta_new, n_xi_new);
X_new = zeros(n_eta_new, n_xi_new);
Y_new = zeros(n_eta_new, n_xi_new);

for i = 1:n_eta_new
    U_new(i, :) = interp1(s_wall, U_mid(i, :), xi_new, 'pchip');
    V_new(i, :) = interp1(s_wall, V_mid(i, :), xi_new, 'pchip');
    W_new(i, :) = interp1(s_wall, W_mid(i, :), xi_new, 'pchip');
    X_new(i, :) = interp1(s_wall, X_mid(i, :), xi_new, 'pchip');
    Y_new(i, :) = interp1(s_wall, Y_mid(i, :), xi_new, 'pchip');
end

%% --- 8. Optional inflow/outflow station trim ----------------------------

n_in  = resolve_trim(params, 'xi_trim_inflow',  n_xi_new);
n_out = resolve_trim(params, 'xi_trim_outflow', n_xi_new);
if n_in + n_out > 0
    if n_in + n_out >= n_xi_new
        error('xi_trim_inflow + xi_trim_outflow (%d+%d) leaves no stations.', n_in, n_out);
    end
    keep     = (1 + n_in) : (n_xi_new - n_out);
    xi_new   = xi_new(keep);
    U_new    = U_new(:, keep);
    V_new    = V_new(:, keep);
    W_new    = W_new(:, keep);
    X_new    = X_new(:, keep);
    Y_new    = Y_new(:, keep);
    n_xi_new = numel(keep);
    fprintf('ξ trimmed (inflow %d, outflow %d): %d stations retained\n', ...
            n_in, n_out, n_xi_new);
end

%% --- 9. Enforce exact no-slip at the wall row ---------------------------

U_new(end, :) = 0;
V_new(end, :) = 0;
W_new(end, :) = 0;

fprintf('max |U| input:        %.6f\n', max(abs(U(:))));
fprintf('max |U| interpolated: %.6f\n', max(abs(U_new(:))));

%% --- 10. Assemble output struct -----------------------------------------

BF_new = struct();
BF_new.X          = X_new;
BF_new.Y          = Y_new;
BF_new.U          = U_new;
BF_new.V          = V_new;
BF_new.W          = W_new;
BF_new.xi1D       = xi_new;
BF_new.eta1D      = eta_cheb;
BF_new.nx         = n_xi_new;
BF_new.ny         = n_eta_new;
BF_new.H          = H;
BF_new.y_i        = y_i;
BF_new.mesh_type  = 'curvilinear';

fprintf('Curvilinear resample done: %d η × %d ξ  (H=%.4f, y_i=%.4f)\n', ...
        n_eta_new, n_xi_new, H, y_i);
end


function [X, Y, U, V, W, ny, nx] = orient_structured_grid(BF)
% Ensure the grid is [ny × nx] with the wall at the last row.

X = BF.X; Y = BF.Y; U = BF.U; V = BF.V; W = BF.W;

if range(X(1, :)) < range(X(:, 1))
    X = X.'; Y = Y.'; U = U.'; V = V.'; W = W.';
    fprintf('Input transposed so streamwise is dim 2.\n');
end

[ny, nx] = size(X);

% Wall = row with lower mean velocity magnitude (no-slip).
vel_top = mean(abs(U(1,   :)) + abs(V(1,   :)) + abs(W(1,   :)));
vel_bot = mean(abs(U(end, :)) + abs(V(end, :)) + abs(W(end, :)));
if vel_top < vel_bot
    X = flipud(X); Y = flipud(Y); U = flipud(U); V = flipud(V); W = flipud(W);
    fprintf('Grid flipped so wall is at row end (DeHNSSo convention).\n');
end
end
