function [BF_new] = resample_unstructured(BF, params)
%RESAMPLE_UNSTRUCTURED  Build a HNS-ready body-fitted grid from unstructured scatter.
%
%   For BF data that comes in as a flat point cloud (e.g. OpenFOAM CSV
%   exports, M3J midPlane.csv). No grid structure is assumed — the wall
%   is detected from the boundary polygon using the velocity field, then
%   smoothed, parameterised, and used to launch wall-normal columns.
%
% Pipeline:
%   1. Optional rescaling (when params.rescale=true, with lref/Uref)
%   2. Detect wall via boundary polygon + velocity criterion
%   3. Smooth and parameterise wall curve
%   4. Compute wall normals
%   5. Estimate H and delta_99
%   6. Build Malik-Chebyshev target grid
%   7. Column-by-column 1D wall-normal interpolation
%
% Inputs
%   BF      - struct with flat fields X, Y, U, V, W (column vectors)
%   params  - struct (see caller examples in gridgen/benchmark/)
%
% Output
%   BF_new  - struct with [n_eta x n_xi] arrays on Malik-Chebyshev grid,
%             plus xi1D, eta1D, H, y_i, mesh_type='curvilinear'.

%% =====================================================================
%  0. Parse parameters
%  =====================================================================
n_eta_new = 101;
n_xi_new  = [];

if isfield(params, 'n_eta_new') && ~isempty(params.n_eta_new)
    n_eta_new = params.n_eta_new;
end
if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    n_xi_new = params.n_xi_new;
end

X_flat = BF.X(:);
Y_flat = BF.Y(:);
U_flat = BF.U(:);
V_flat = BF.V(:);
W_flat = BF.W(:);

%% =====================================================================
%  1. Rescale (opt-in, for dimensional input)
%  =====================================================================
if isfield(params, 'rescale') && ~isempty(params.rescale) && params.rescale
    assert(isfield(params,'lref') && ~isempty(params.lref) && ...
           isfield(params,'Uref') && ~isempty(params.Uref), ...
        'resample_unstructured: rescale=true requires params.lref and params.Uref.');
    fprintf('Rescaling input: Uref = %.4f, lref = %.6e\n', params.Uref, params.lref);
    X_flat = X_flat / params.lref;
    Y_flat = Y_flat / params.lref;
    U_flat = U_flat / params.Uref;
    V_flat = V_flat / params.Uref;
    W_flat = W_flat / params.Uref;
end

fprintf('Input data: %d points, x=[%.4f, %.4f], y=[%.4f, %.4f]\n', ...
    numel(X_flat), min(X_flat), max(X_flat), min(Y_flat), max(Y_flat));

%% =====================================================================
%  2. Wall detection via boundary polygon + velocity
%  =====================================================================
% boundary() gives the convex/concave hull of the point cloud.
% The wall is the longest contiguous segment of low-velocity boundary points.
bnd     = boundary(X_flat, Y_flat, 1);
bnd     = bnd(1:end-1);               % remove closure duplicate, keep order
n_bnd   = numel(bnd);

uvw_bnd = abs(U_flat(bnd)) + abs(V_flat(bnd)) + abs(W_flat(bnd));
thr     = median(uvw_bnd);
is_low  = uvw_bnd < thr;

% Find longest contiguous run of low-velocity points (with wraparound)
is_low2    = [is_low; is_low];
d_low      = diff([0; double(is_low2); 0]);
run_starts = find(d_low ==  1);
run_ends   = find(d_low == -1) - 1;
run_lens   = min(run_ends - run_starts + 1, n_bnd);
[~, best]  = max(run_lens);
seg        = run_starts(best) : run_starts(best) + run_lens(best) - 1;
seg        = mod(seg - 1, n_bnd) + 1;

wall_idx = bnd(seg);

% Peel endpoints where the boundary polygon detours along an inflow/outflow
% face: either tightly clustered points (small Euclidean spacing) or
% perpendicular excursions (small projection onto the local wall tangent).
X_wi  = X_flat(wall_idx);  Y_wi = Y_flat(wall_idx);
dx_wi = diff(X_wi);
dy_wi = diff(Y_wi);
d_wi  = hypot(dx_wi, dy_wi);
n_seg = numel(d_wi);

robust_med = @(v) median(v(max(1,round(0.10*numel(v))) : ...
                          max(round(0.10*numel(v))+1, round(0.90*numel(v)))));
d_ref = robust_med(sort(d_wi));

% Reference wall tangent: median segment vector over the interior 60%
i_lo_t  = max(1, round(0.20 * n_seg));
i_hi_t  = min(n_seg, round(0.80 * n_seg));
t_ref_x = median(dx_wi(i_lo_t:i_hi_t));
t_ref_y = median(dy_wi(i_lo_t:i_hi_t));
t_ref_n = max(eps, hypot(t_ref_x, t_ref_y));

proj = @(i) (dx_wi(i)*t_ref_x + dy_wi(i)*t_ref_y) / t_ref_n;
ok_seg = @(i) d_wi(i)    > 0.3*d_ref && d_wi(i) < 3*d_ref ...
           && abs(proj(i)) > 0.3*d_ref;

lo = 1;  hi = numel(wall_idx);
while lo < hi - 1 && ~ok_seg(lo)
    lo = lo + 1;
end
while hi > lo + 1 && ~ok_seg(hi - 1)
    hi = hi - 1;
end
wall_idx = wall_idx(lo:hi);
n_wall   = numel(wall_idx);
if isempty(n_xi_new)
    n_xi_new = n_wall;
end
fprintf('Wall detection: %d boundary points (after inflow/outflow trim)\n', n_wall);

X_wall = X_flat(wall_idx);
Y_wall = Y_flat(wall_idx);

% --- Plot wall detection and wait for confirmation ---
figure;
subplot(2,1,1);
scatter(X_flat, Y_flat, 2, abs(U_flat)+abs(V_flat)+abs(W_flat), 'filled');
colormap(jet); colorbar; daspect([1 1 1]);
title('|U|+|V|+|W|');

subplot(2,1,2);
plot(X_flat, Y_flat, 'k.', 'MarkerSize', 1); hold on;
plot(X_wall, Y_wall, 'r.', 'MarkerSize', 6);
daspect([1 1 1]);
legend('input', 'detected wall', 'Location', 'best');
title(sprintf('Wall detection — %d points', n_wall));

if usejava('desktop')
    input('Wall detection OK? Press Enter to continue or Ctrl+C to abort: ', 's');
end

%% =====================================================================
%  3. Sort and smooth wall curve
%  =====================================================================
% Sort wall points along the curve via nearest-neighbour traversal
[~, i0] = min(X_wall);
order      = zeros(n_wall, 1);
visited    = false(n_wall, 1);
order(1)   = i0;
visited(i0) = true;

for k = 2:n_wall
    prev = order(k-1);
    dx   = X_wall - X_wall(prev);
    dy   = Y_wall - Y_wall(prev);
    d    = dx.^2 + dy.^2;
    d(visited) = inf;

    if k > 2
        prev2 = order(k-2);
        dir_x = X_wall(prev) - X_wall(prev2);
        dir_y = Y_wall(prev) - Y_wall(prev2);
        backward = (dx * dir_x + dy * dir_y) < 0;
        d_fwd = d;  d_fwd(backward) = inf;
        if any(~isinf(d_fwd))
            [~, next] = min(d_fwd);
        else
            [~, next] = min(d);
        end
    else
        [~, next] = min(d);
    end
    order(k)     = next;
    visited(next) = true;
end

X_wall = X_wall(order);
Y_wall = Y_wall(order);

% Smooth cell-centre zigzag with a smoothing spline
s_raw  = [0; cumsum(sqrt(diff(X_wall).^2 + diff(Y_wall).^2))];
ppX_sm = csaps(s_raw, X_wall, 1 - 1e-4);
ppY_sm = csaps(s_raw, Y_wall, 1 - 1e-4);
X_wall = fnval(ppX_sm, s_raw);
Y_wall = fnval(ppY_sm, s_raw);

% Recompute arc-length on smoothed curve
dxi_wall = sqrt(diff(X_wall).^2 + diff(Y_wall).^2);
xi_wall  = [0; cumsum(dxi_wall)];   % [n_wall x 1]
fprintf('Wall arc-length: %.4f (%d points)\n', xi_wall(end), n_wall);

%% =====================================================================
%  4. Wall tangent and outward normal
%  =====================================================================
% Resample to uniform arc-length for clean normals
n_xi_wall = max(n_xi_new, n_wall);
xi_uniform = linspace(0, xi_wall(end), n_xi_wall)';

X_wu = interp1(xi_wall, X_wall, xi_uniform, 'pchip');
Y_wu = interp1(xi_wall, Y_wall, xi_uniform, 'pchip');

% Central FD tangent (2nd order)
tx = zeros(n_xi_wall, 1);
ty = zeros(n_xi_wall, 1);
ds = xi_uniform(2) - xi_uniform(1);
tx(1) = (X_wu(2) - X_wu(1)) / ds;
ty(1) = (Y_wu(2) - Y_wu(1)) / ds;
tx(end) = (X_wu(end) - X_wu(end-1)) / ds;
ty(end) = (Y_wu(end) - Y_wu(end-1)) / ds;
tx(2:end-1) = (X_wu(3:end) - X_wu(1:end-2)) / (2*ds);
ty(2:end-1) = (Y_wu(3:end) - Y_wu(1:end-2)) / (2*ds);

t_len = sqrt(tx.^2 + ty.^2);
tx = tx ./ t_len;
ty = ty ./ t_len;

% Outward normal: rotate tangent 90 deg (assumes wall below freestream)
nx_w = -ty;
ny_w =  tx;

% Verify outward direction: normal should point away from wall into fluid.
% Check at midpoint: the mean data centroid should be on the normal side.
mid   = round(n_xi_wall / 2);
x_cen = mean(X_flat);
y_cen = mean(Y_flat);
dot_check = (x_cen - X_wu(mid)) * nx_w(mid) + (y_cen - Y_wu(mid)) * ny_w(mid);
if dot_check < 0
    nx_w = -nx_w;
    ny_w = -ny_w;
    fprintf('Normal direction flipped (pointing into fluid domain)\n');
end

%% =====================================================================
%  5. Estimate H and delta_99
%  =====================================================================
if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
    fprintf('Using prescribed H = %.5f\n', H);
else
    % Project all data onto wall-normal at several stations, take max
    test_j = round(linspace(1, n_xi_wall, min(50, n_xi_wall)));
    H_col  = zeros(numel(test_j), 1);
    for k = 1:numel(test_j)
        j = test_j(k);
        proj = (X_flat - X_wu(j)) .* nx_w(j) + (Y_flat - Y_wu(j)) .* ny_w(j);
        H_col(k) = max(proj);
    end
    H = median(H_col);
    fprintf('Estimated H = %.5f\n', H);
end

% Estimate delta_99 at a mid-domain station
j_est    = round(n_xi_wall / 2);
dx_all   = X_flat - X_wu(j_est);
dy_all   = Y_flat - Y_wu(j_est);
eta_all  = dx_all * nx_w(j_est) + dy_all * ny_w(j_est);
tau_all  = abs(-dx_all * ny_w(j_est) + dy_all * nx_w(j_est));
band_est = (xi_uniform(end) - xi_uniform(1)) / n_xi_wall * 2;
near_est = tau_all < band_est & eta_all >= 0 & eta_all <= H;
[eta_s, si] = sort(eta_all(near_est));
U_s      = U_flat(near_est);  U_s = U_s(si);
U_edge   = U_s(end);
idx99    = find(abs(U_s) >= 0.99 * abs(U_edge), 1, 'first');
if ~isempty(idx99)
    delta_99 = eta_s(idx99);
else
    delta_99 = H / 3;
end
fprintf('Estimated delta_99 = %.5f  (at mid-domain station)\n', delta_99);
fprintf('  Suggestion: H   ~ %.5f  (4 * delta_99)\n', 4 * delta_99);
fprintf('  Suggestion: y_i ~ %.5f  (delta_99)\n',     delta_99);

if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
    fprintf('Using prescribed y_i = %.5f\n', y_i);
else
    y_i = min(delta_99, H / 3);
    fprintf('Using y_i = %.5f\n', y_i);
end

% Malik mapping requires y_i < H/2
if y_i >= H / 2
    y_i = H / 3;
    warning('y_i clamped to H/3 = %.5f (must be < H/2 for valid Malik mapping).', y_i);
end

%% =====================================================================
%  6. Build Malik-Chebyshev target grid
%  =====================================================================
a = H * y_i / (H - 2.0 * y_i);
b = 1.0 + 2.0 * a / H;

% Descending: eta(1)=H (freestream), eta(end)=0 (wall) — DeHNSSo convention
k_vec    = (0 : n_eta_new - 1)';
z_cheb   = cos(pi * k_vec / (n_eta_new - 1));
eta_cheb = a * (1 + z_cheb) ./ (b - z_cheb);   % [n_eta_new x 1]

% Streamwise: uniform arc-length stations
xi1D_new = linspace(0, xi_uniform(end), n_xi_new);   % [1 x n_xi_new]

% Interpolate wall geometry to target xi stations
X_wall_t = interp1(xi_uniform, X_wu, xi1D_new, 'pchip');
Y_wall_t = interp1(xi_uniform, Y_wu, xi1D_new, 'pchip');
nx_t     = interp1(xi_uniform, nx_w, xi1D_new, 'pchip');
ny_t     = interp1(xi_uniform, ny_w, xi1D_new, 'pchip');

% Renormalise normals
n_len = sqrt(nx_t.^2 + ny_t.^2);
nx_t  = nx_t ./ n_len;
ny_t  = ny_t ./ n_len;

% Physical target coordinates: [n_eta_new x n_xi_new]
X_target = X_wall_t + eta_cheb .* nx_t;
Y_target = Y_wall_t + eta_cheb .* ny_t;

%% Step 6b: Optional xi trimming
if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    xi_lo = params.xi_range(1) * xi1D_new(end);
    xi_hi = params.xi_range(2) * xi1D_new(end);
    keep  = xi1D_new >= xi_lo & xi1D_new <= xi_hi;
    xi1D_new = xi1D_new(keep);
    X_wall_t = X_wall_t(keep);
    Y_wall_t = Y_wall_t(keep);
    nx_t     = nx_t(keep);
    ny_t     = ny_t(keep);
    X_target = X_target(:, keep);
    Y_target = Y_target(:, keep);
    n_xi_new = sum(keep);
    fprintf('xi trimmed: %d stations retained\n', n_xi_new);
end

% Additional inflow/outflow trim by count or fraction
n_in  = resolve_trim(params, 'xi_trim_inflow',  numel(xi1D_new));
n_out = resolve_trim(params, 'xi_trim_outflow', numel(xi1D_new));
if n_in + n_out > 0
    if n_in + n_out >= numel(xi1D_new)
        error('xi_trim_inflow + xi_trim_outflow (%d+%d) leaves no stations (have %d).', ...
              n_in, n_out, numel(xi1D_new));
    end
    keep_idx = (1 + n_in) : (numel(xi1D_new) - n_out);
    xi1D_new = xi1D_new(keep_idx);
    X_wall_t = X_wall_t(keep_idx);
    Y_wall_t = Y_wall_t(keep_idx);
    nx_t     = nx_t(keep_idx);
    ny_t     = ny_t(keep_idx);
    X_target = X_target(:, keep_idx);
    Y_target = Y_target(:, keep_idx);
    n_xi_new = numel(keep_idx);
    fprintf('xi trimmed (inflow %d, outflow %d): %d stations retained\n', ...
            n_in, n_out, n_xi_new);
end

%% Step 6c: Elliptic mesh smoothing (optional)
if isfield(params, 'elliptic') && params.elliptic
    ell_opts = struct('max_iter', 5000, 'tol', 1e-6, 'omega', 1.7, 'n_uniform', 201);
    if isfield(params, 'elliptic_opts')
        flds = fieldnames(params.elliptic_opts);
        for f = 1:numel(flds), ell_opts.(flds{f}) = params.elliptic_opts.(flds{f}); end
    end
    [X_target, Y_target] = solve_elliptic_mesh(X_target, Y_target, ...
        X_wall_t, Y_wall_t, nx_t, ny_t, eta_cheb, H, ell_opts);
end

%% =====================================================================
%  7. 2-D scattered interpolation onto the body-fitted target grid
%  =====================================================================
% Cell-centred CFD exports (e.g. OpenFOAM CSV) do not include the wall
% itself — the lowest data row sits half a cell above it. We pick that row
% as our "wall" target and let the natural cell-centre values stand. Do
% NOT force U=0: the resulting step would be amplified by the near-wall
% stencil into a spurious O(10) wall shear. u'=0 is enforced by the
% stability operator independently.

F_U = scatteredInterpolant(X_flat, Y_flat, U_flat, 'linear', 'nearest');
F_V = F_U;  F_V.Values = V_flat;    % reuse triangulation, swap values
F_W = F_U;  F_W.Values = W_flat;

U_new = reshape(F_U(X_target(:), Y_target(:)), n_eta_new, n_xi_new);
V_new = reshape(F_V(X_target(:), Y_target(:)), n_eta_new, n_xi_new);
W_new = reshape(F_W(X_target(:), Y_target(:)), n_eta_new, n_xi_new);

fprintf('max U input:        %.6f\n', max(abs(U_flat)));
fprintf('max U interpolated: %.6f\n', max(abs(U_new(:))));

%% =====================================================================
%  8. Diagnostic plots
%  =====================================================================

% --- Fig: wall-normal profiles with raw-data overlay (sanity check) ---
if isfield(params, 'plot') && params.plot
if isfield(params, 'plot_opts') && isfield(params.plot_opts, 'xi_slices')
    xi_fracs = params.plot_opts.xi_slices;
else
    xi_fracs = linspace(0, 1, 8);
end
j_slices = max(1, min(n_xi_new, round(xi_fracs * n_xi_new)));

figure('Name', 'Wall-normal profiles with raw overlay');
colors_p = lines(numel(j_slices));
h_lines  = gobjects(numel(j_slices), 1);
leg_str  = cell(numel(j_slices), 1);
for k = 1:numel(j_slices)
    j = j_slices(k);

    % Raw data near this station (for comparison).
    % Use a tangential band sized to the local ξ spacing so the scatter
    % overlay captures the nearby raw points.
    xi_spacing_j = max(1e-12, ...
        0.5 * (norm([X_wall_t(min(j+1,end)); Y_wall_t(min(j+1,end))] - ...
                    [X_wall_t(max(j-1,1));   Y_wall_t(max(j-1,1))])));
    dx_r  = X_flat - X_wall_t(j);
    dy_r  = Y_flat - Y_wall_t(j);
    eta_r = dx_r * nx_t(j) + dy_r * ny_t(j);
    tau_r = abs(-dx_r * ny_t(j) + dy_r * nx_t(j));
    near_r = tau_r < 0.6 * xi_spacing_j & eta_r >= 0 & eta_r <= H;

    subplot(1,3,1);
    plot(eta_r(near_r), U_flat(near_r), '.', 'Color', colors_p(k,:)*0.5+0.5, ...
         'MarkerSize', 4, 'HandleVisibility', 'off'); hold on;
    h_lines(k) = plot(eta_cheb, U_new(:,j), '-', 'Color', colors_p(k,:), 'LineWidth', 1.5);

    subplot(1,3,2);
    plot(eta_r(near_r), V_flat(near_r), '.', 'Color', colors_p(k,:)*0.5+0.5, ...
         'MarkerSize', 4, 'HandleVisibility', 'off'); hold on;
    plot(eta_cheb, V_new(:,j), '-', 'Color', colors_p(k,:), 'LineWidth', 1.5);

    subplot(1,3,3);
    plot(eta_r(near_r), W_flat(near_r), '.', 'Color', colors_p(k,:)*0.5+0.5, ...
         'MarkerSize', 4, 'HandleVisibility', 'off'); hold on;
    plot(eta_cheb, W_new(:,j), '-', 'Color', colors_p(k,:), 'LineWidth', 1.5);

    leg_str{k} = sprintf('\\xi = %.3f', xi1D_new(j));
end
subplot(1,3,1); xlabel('\eta'); ylabel('U'); title('U profiles');
legend(h_lines, leg_str, 'Location', 'best');
subplot(1,3,2); xlabel('\eta'); ylabel('V'); title('V profiles');
subplot(1,3,3); xlabel('\eta'); ylabel('W'); title('W profiles');
end

% --- Fig: interpolated mesh and velocity ---
if isfield(params, 'plot') && params.plot
    vel_names = {'U', 'V', 'W'};
    vel_new   = {U_new, V_new, W_new};

    figure('Name', 'Interpolated mesh');
    for ip = 1:3
        subplot(4,1,ip);
        scatter(X_target(:), Y_target(:), 4, vel_new{ip}(:), 'filled');
        clim([min(vel_new{ip}(:)) max(vel_new{ip}(:))]);
        colormap(jet); colorbar; daspect([1 1 1]);
        title(sprintf('Interpolated %s', vel_names{ip}));
    end

    subplot(4,1,4);
    plot(X_flat, Y_flat, 'k.', 'MarkerSize', 1); hold on;
    plot(X_wall_t, Y_wall_t, 'r-', 'LineWidth', 2);
    plot(X_target,  Y_target,  'b-', 'LineWidth', 0.3);
    plot(X_target', Y_target', 'b-', 'LineWidth', 0.3);
    daspect([1 1 1]);
    legend('input', 'wall', 'target mesh', 'Location', 'best');
    title(sprintf('Target mesh: %d \\eta \\times %d \\xi', n_eta_new, n_xi_new));
end

%% =====================================================================
%  9. Assemble output
%  =====================================================================
BF_new.X         = X_target;
BF_new.Y         = Y_target;
BF_new.U         = U_new;
BF_new.V         = V_new;
BF_new.W         = W_new;
BF_new.xi1D      = xi1D_new;
BF_new.eta1D     = eta_cheb;
BF_new.nx        = n_xi_new;
BF_new.ny        = n_eta_new;
BF_new.H         = H;
BF_new.y_i       = y_i;
BF_new.mesh_type = 'curvilinear';

fprintf('Unstructured resample done: %d η × %d ξ  (H=%.4f, y_i=%.4f)\n', ...
    n_eta_new, n_xi_new, H, y_i);
end
