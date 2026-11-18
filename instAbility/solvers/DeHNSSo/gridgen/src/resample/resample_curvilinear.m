function BF_new = resample_curvilinear(BF, params)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Interpolate a structured body-fitted BF onto a Malik-Chebyshev grid.
%
%   BF_new = resample_curvilinear(BF, params)
%
%   For a curvilinear structured mesh: X and Y both 2-D, Y varying along a row
%   (a hump-following grid). The wall is the last row, so no polygon detection.
%
%   params: n_eta_new (required) wall-normal Chebyshev points; n_xi_new (default
%   nx_raw); H domain height (default 0.99 x max column arc length); y_i Malik
%   median height (default H/10); xi_range as a fraction of wall arc length;
%   xi_trim_inflow / xi_trim_outflow, a count or fraction dropped from each end;
%   rescale (default false, only for dimensional input, since heritage BFs are
%   already non-dimensional); Re, lref, Uref, nu passed on via ensure_metadata.
%
%   Output: BF_new for build_stab_grid (mesh_type='curvilinear').

%% Inputs
% Name   size         unit   explanation
% BF     (1,1) struct [-]    base flow on its own grid
% params (1,1) struct [-]    grid-generation parameters

%% Outputs
% Name   size         unit   explanation
% BF_new (1,1) struct [-]    base flow resampled onto the stability grid

%% --- 1. Normalise orientation to [ny × nx], wall at row end --------------

[X, Y, U, V, W, ny_raw, nx_raw, Draw] = orient_structured_grid(BF);

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

%% --- 6. xi stations along the wall --------------------------------------
% Uniform wall arc-length, or Gaussian clustering around params.mug (see
% cluster_xi). spread (refinement half-width) sizes sig, as on the flat plate.
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

mug   = get_opt(params, 'mug', []);
ag    = get_opt(params, 'ag',  1);
spread = get_opt(params, 'spread', []);

% xi is arc-length; centres are physical x, so map the range ends to wall x.
x_lo_p = interp1(s_wall, X_wall, s_lo, 'pchip');
x_hi_p = interp1(s_wall, X_wall, s_hi, 'pchip');
if ~isempty(spread)
    sig = 3 * spread ./ (0.5 * (x_hi_p - x_lo_p));
else
    sig = get_opt(params, 'sig', 0.2);
end
xi_new = cluster_xi(s_lo, s_hi, n_xi_new, mug, sig, ag, x_lo_p, x_hi_p);
report_clustering(mug, sig, ag, xi_new);

%% --- 7. Resample η per raw column, then resample ξ ----------------------

% Pass 1: η direction (per raw column), raw grid → Malik-Cheb η.
U_mid = zeros(n_eta_new, nx_raw);
V_mid = zeros(n_eta_new, nx_raw);
W_mid = zeros(n_eta_new, nx_raw);
X_mid = zeros(n_eta_new, nx_raw);
Y_mid = zeros(n_eta_new, nx_raw);
dnf   = fieldnames(Draw)';
D_mid = struct();
for q = dnf; D_mid.(q{1}) = zeros(n_eta_new, nx_raw); end

bim = get_opt(params, 'bf_interp', 'pchip');   % base-flow/geometry interp; 'spline' -> smooth C2 derivatives
for j = 1:nx_raw
    eta_col       = d_raw(:, j);                      % 0 at wall → max at top
    [eta_asc, o]  = sort(eta_col);                    % ascending for interp1
    eta_q         = max(min(eta_cheb, eta_asc(end)), eta_asc(1));

    U_mid(:, j) = interp1(eta_asc, U(o, j), eta_q, bim);
    V_mid(:, j) = interp1(eta_asc, V(o, j), eta_q, bim);
    W_mid(:, j) = interp1(eta_asc, W(o, j), eta_q, bim);
    X_mid(:, j) = interp1(eta_asc, X(o, j), eta_q, bim);
    Y_mid(:, j) = interp1(eta_asc, Y(o, j), eta_q, bim);
    for q = dnf; D_mid.(q{1})(:, j) = interp1(eta_asc, Draw.(q{1})(o, j), eta_q, bim); end
end

% Pass 2: ξ direction (per η level), raw → target wall arc-length.
U_new = zeros(n_eta_new, n_xi_new);
V_new = zeros(n_eta_new, n_xi_new);
W_new = zeros(n_eta_new, n_xi_new);
X_new = zeros(n_eta_new, n_xi_new);
Y_new = zeros(n_eta_new, n_xi_new);
D_new = struct();
for q = dnf; D_new.(q{1}) = zeros(n_eta_new, n_xi_new); end

for i = 1:n_eta_new
    U_new(i, :) = interp1(s_wall, U_mid(i, :), xi_new, bim);
    V_new(i, :) = interp1(s_wall, V_mid(i, :), xi_new, bim);
    W_new(i, :) = interp1(s_wall, W_mid(i, :), xi_new, bim);
    X_new(i, :) = interp1(s_wall, X_mid(i, :), xi_new, bim);
    Y_new(i, :) = interp1(s_wall, Y_mid(i, :), xi_new, bim);
    for q = dnf; D_new.(q{1})(i, :) = interp1(s_wall, D_mid.(q{1})(i, :), xi_new, bim); end
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
    for q = dnf; D_new.(q{1}) = D_new.(q{1})(:, keep); end
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
for q = dnf; BF_new.(q{1}) = D_new.(q{1}); end
BF_new.xi1D       = xi_new;
BF_new.eta1D      = eta_cheb;
BF_new.nx         = n_xi_new;
BF_new.ny         = n_eta_new;
BF_new.H          = H;
BF_new.y_i        = y_i;
BF_new.mesh_type  = 'curvilinear';

% Clustering metadata in the xi (arc-length) frame; build_stab_grid maps xc_xi
% to the grid with the same offset it applies to xi1D. See add_wall_forcing.
xc_xi = [];
if ~isempty(mug)
    xc_xi = interp1(X_wall, s_wall, mug(:).', 'pchip');   % physical x -> arc-length
end
BF_new.cluster = struct('mug', mug, 'sig', sig, 'ag', ag, 'xc_xi', xc_xi, 'w', spread);

fprintf('Curvilinear resample done: %d η × %d ξ  (H=%.4f, y_i=%.4f)\n', ...
        n_eta_new, n_xi_new, H, y_i);
end


function [X, Y, U, V, W, ny, nx, Draw] = orient_structured_grid(BF)
% Ensure the grid is [ny × nx] with the wall at the last row.

X = BF.X; Y = BF.Y; U = BF.U; V = BF.V; W = BF.W;
Draw = struct();
for dn = bf_deriv_fields()
    if isfield(BF, dn{1}) && ~isempty(BF.(dn{1})); Draw.(dn{1}) = BF.(dn{1}); end
end

if (max(X(1, :)) - min(X(1, :))) < (max(X(:, 1)) - min(X(:, 1)))
    X = X.'; Y = Y.'; U = U.'; V = V.'; W = W.';
    for dn = fieldnames(Draw)'; Draw.(dn{1}) = Draw.(dn{1}).'; end
    fprintf('Input transposed so streamwise is dim 2.\n');
end

[ny, nx] = size(X);

% Wall = row with lower mean velocity magnitude (no-slip).
vel_top = mean(abs(U(1,   :)) + abs(V(1,   :)) + abs(W(1,   :)));
vel_bot = mean(abs(U(end, :)) + abs(V(end, :)) + abs(W(end, :)));
if vel_top < vel_bot
    X = flipud(X); Y = flipud(Y); U = flipud(U); V = flipud(V); W = flipud(W);
    for dn = fieldnames(Draw)'; Draw.(dn{1}) = flipud(Draw.(dn{1})); end
    fprintf('Grid flipped so wall is at row end (DeHNSSo convention).\n');
end
end
