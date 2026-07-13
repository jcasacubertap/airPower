function StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes)
%POSTPROCESS  Normalize shape functions, compute alpha, delete conjugates.

iu = Init.iu;

%% Physical ∂u/∂x via chain rule on the (possibly curvilinear, non-uniform ξ) grid
%   ∂u/∂x = ξx · ∂u/∂ξ + ηx · ∂u/∂η
%   ∂u/∂ξ : sparse 4th-order FD matrix built with Fornberg weights (handles
%           non-uniform ξ correctly; reduces to classical coefficients when ξ is uniform)
%   ∂u/∂η : StabGrid.D1 (Chebyshev or FD, consistent with the stability operator)
D_xi = fd1_matrix_nonuniform(StabGrid.x(:));
for j = Mode.RunJ
    u_j     = squeeze(StabRes.u(j,:,:));          % [ny × nx]
    du_dxi  = u_j * D_xi.';                       % derivative along rows (ξ)
    du_deta = StabGrid.D1 * u_j;                  % spectral η-derivative
    dux(j,:,:) = StabGrid.xix .* du_dxi + StabGrid.etax .* du_deta;
end

%% Normalize shape functions and compute amplitudes
% Interpolate all stations at once using matrix interp1 — one spline per mode
% instead of 2*nx scalar spline calls per mode.
yinteg = linspace(0, max(StabGrid.eta1D), 4000);

for j = Mode.RunJ
    % Interpolate u and dux over all stations simultaneously [4000 x nx]
    U_fine   = interp1(StabGrid.etaun, squeeze(StabRes.u(j,:,:)), yinteg, 'spline');
    DUX_fine = interp1(StabGrid.eta1D, squeeze(dux(j,:,:)),       yinteg, 'spline');

    % Find wall-normal location of max |u| per station (integer index)
    absU_fine   = abs(U_fine);
    [~, k_all]  = max(absU_fine, [], 1);              % 1 × nx
    n_fine      = size(U_fine, 1);
    k_all       = max(2, min(n_fine-1, k_all));       % keep k±1 in range

    % Quadratic sub-index refinement of the peak location.
    %   Fit |u| ≈ a·δ² + b·δ + c on (k-1, k, k+1);  δ_peak = -b/(2a)
    % This converts the discrete row-index jumps into a smooth continuous
    % η-location, eliminating the stair-step wiggle in α_i.
    jj   = 1:StabGrid.nx;
    idxm = sub2ind(size(absU_fine), k_all-1, jj);
    idx0 = sub2ind(size(absU_fine), k_all,   jj);
    idxp = sub2ind(size(absU_fine), k_all+1, jj);
    um = absU_fine(idxm);  u0 = absU_fine(idx0);  up = absU_fine(idxp);
    denom = um - 2*u0 + up;
    denom(abs(denom) < eps) = eps;                    % avoid div-by-zero on flat peaks
    delta = 0.5 * (um - up) ./ denom;                 % ∈ [-1, 1] for well-behaved peaks
    delta = max(-1, min(1, delta));
    k_cont = k_all + delta;                           % continuous peak row

    % Evaluate U_fine and DUX_fine at the continuous peak via 1-D linear interp
    k_floor = floor(k_cont);  k_ceil = min(k_floor+1, n_fine);  frac = k_cont - k_floor;
    idxf = sub2ind(size(U_fine), k_floor, jj);
    idxc = sub2ind(size(U_fine), k_ceil,  jj);
    umax_vec  = (1-frac) .* U_fine(idxf)  + frac .* U_fine(idxc);
    dumax_vec = (1-frac) .* DUX_fine(idxf) + frac .* DUX_fine(idxc);

    % Alpha a posteriori (vectorised)
    StabRes.alpha(j,:) = (dumax_vec ./ umax_vec) / iu;

    % Normalise shape functions and compute amplitudes (vectorised over stations)
    absU = abs(umax_vec);                           % 1 x nx
    StabRes.u(j,:,:) = squeeze(StabRes.u(j,:,:)) ./ absU;
    StabRes.v(j,:,:) = squeeze(StabRes.v(j,:,:)) ./ absU;
    StabRes.w(j,:,:) = squeeze(StabRes.w(j,:,:)) ./ absU;
    StabRes.p(j,:,:) = squeeze(StabRes.p(j,:,:)) ./ absU;
    if isfield(StabRes,'rho')
        StabRes.rho(j,:,:) = squeeze(StabRes.rho(j,:,:)) ./ absU;
        StabRes.T(j,:,:)   = squeeze(StabRes.T(j,:,:))   ./ absU;
    end

    if j == round(Mode.nf/2)
        StabRes.A(j,:) = absU;
    else
        StabRes.A(j,:) = 2 * absU;
    end
end

%% Disturbance energy per mode (kinetic, or Chu/Mack for compressible)
StabRes.E = zeros(size(StabRes.A));
half_nf   = round(Mode.nf/2);
is_compr  = isfield(StabRes,'rho');
if is_compr
    rho_bar = StabGrid.rho;
    T_bar   = StabGrid.T;
    Ma2     = StabGrid.Ma^2;
    ga      = StabGrid.gamma;
end
for j = Mode.RunJ
    if j == half_nf
        scale = StabRes.A(j, :);          % MFD: A = |u_max|
    else
        scale = StabRes.A(j, :) / 2;      % non-MFD: A = 2|u_max|
    end
    u_full = squeeze(StabRes.u(j, :, :)) .* scale;
    v_full = squeeze(StabRes.v(j, :, :)) .* scale;
    w_full = squeeze(StabRes.w(j, :, :)) .* scale;
    if is_compr
        rho_full  = squeeze(StabRes.rho(j, :, :)) .* scale;
        T_full    = squeeze(StabRes.T  (j, :, :)) .* scale;
        integrand = rho_bar .* (abs(u_full).^2 + abs(v_full).^2 + abs(w_full).^2) ...
                  + T_bar   ./ (ga*rho_bar*Ma2)          .* abs(rho_full).^2 ...
                  + rho_bar ./ (ga*(ga-1)*Ma2*T_bar)     .* abs(T_full).^2;
    else
        integrand = abs(u_full).^2 + abs(v_full).^2 + abs(w_full).^2;
    end
    StabRes.E(j, :) = 0.5 * trapz(StabGrid.y, integrand, 1);
end

%% Delete conjugates
half = round(Mode.nf/2);
StabRes.omegavec(:, 1:half-1) = [];
StabRes.betavec(:, 1:half-1)  = [];
StabRes.A(1:half-1, :)        = [];
StabRes.E(1:half-1, :)        = [];
StabRes.alpha(1:half-1, :)    = [];
StabRes.u(1:half-1, :, :)     = [];
StabRes.v(1:half-1, :, :)     = [];
StabRes.w(1:half-1, :, :)     = [];
StabRes.p(1:half-1, :, :)     = [];
if isfield(StabRes,'rho')
    StabRes.rho(1:half-1, :, :) = [];
    StabRes.T(1:half-1, :, :)   = [];
end

if Opt.Sweep
    StabRes.Asweep(1:half-1, :, :)       = [];
    StabRes.alphasweep(1:half-1, :, :)   = [];
    StabRes.usweep(1:half-1, :, :, :)    = [];
    StabRes.vsweep(1:half-1, :, :, :)    = [];
    StabRes.wsweep(1:half-1, :, :, :)    = [];
    StabRes.psweep(1:half-1, :, :, :)    = [];
end

%% Remove working arrays
StabRes = rmfield(StabRes, {'q', 'phi', 'R'});

end
