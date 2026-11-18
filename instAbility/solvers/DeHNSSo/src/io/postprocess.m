function StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes)
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
% Normalize shape functions, recover alpha, drop conjugate modes.
%
%   StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes)
%
%   Recovers the physical streamwise derivative of each mode via the chain
%   rule, normalizes the shape functions by max(|u|), extracts the streamwise
%   wavenumber alpha a posteriori (i.e., after the solve, from the computed
%   shape), evaluates the disturbance energy per mode, and removes the stored
%   complex conjugates.
%

%% Inputs
% Name     size unit   explanation
% StabGrid      [-]    stability grid (metrics, base flow, D1, coordinates)
% Init          [-]    initialization data (Init.iu is the imaginary unit)
% Mode          [-]    mode indexing (RunJ, nf)
% Opt           [-]    run options (Opt.Sweep)
% StabRes       [-]    normalized shape functions, alpha, energy, conjugates removed

%% Outputs
% Name    size unit   explanation
% StabRes      [-]    normalized shape functions, alpha, energy, conjugates removed

iu = Init.iu;

%% Physical streamwise derivative via the chain rule: du/dx = xix*du_dxi + etax*du_deta
% D_xi is the 4th-order Fornberg finite-difference matrix in xi; D1 is the spectral derivative matrix in eta.
D_xi = fd1_matrix_nonuniform(StabGrid.x(:));
for j = Mode.RunJ
    u_j     = squeeze(StabRes.u(j,:,:));          % [ny x nx]
    du_dxi  = u_j * D_xi.';                       % streamwise (xi) derivative
    du_deta = StabGrid.D1 * u_j;                  % spectral eta-derivative
    dux(j,:,:) = StabGrid.xix .* du_dxi + StabGrid.etax .* du_deta;
end

%% Normalize shape functions and compute amplitudes
% makima (non-overshooting) only when a masked step makes u(y) jump at the
% masked region; a spline overshoots there. Gated on hasMask so other cases are unchanged.
eta_interp = 'spline';
if isfield(StabGrid,'hasMask') && StabGrid.hasMask;  eta_interp = 'makima';  end
% Interpolate all stations in one matrix call to interp1 (one spline per mode)
% instead of looping over 2*nx individual scalar spline calls per mode.
yinteg = linspace(0, max(StabGrid.eta1D), 4000);

for j = Mode.RunJ
    % Interpolate u and dux over all stations simultaneously [4000 x nx]
    U_fine   = interp1(StabGrid.etaun, squeeze(StabRes.u(j,:,:)), yinteg, eta_interp);
    DUX_fine = interp1(StabGrid.eta1D, squeeze(dux(j,:,:)),       yinteg, 'spline');

    % Find wall-normal location of max |u| per station (integer index)
    absU_fine   = abs(U_fine);
    [~, k_all]  = max(absU_fine, [], 1);              % 1 x nx
    n_fine      = size(U_fine, 1);
    k_all       = max(2, min(n_fine-1, k_all));       % keep k-1 and k+1 inside the array bounds

    % Refine the peak position with a quadratic fit through k-1, k, k+1, so
    % alpha_i varies smoothly between grid points instead of jumping in integer steps.
    jj   = 1:StabGrid.nx;
    idxm = sub2ind(size(absU_fine), k_all-1, jj);
    idx0 = sub2ind(size(absU_fine), k_all,   jj);
    idxp = sub2ind(size(absU_fine), k_all+1, jj);
    um = absU_fine(idxm);  u0 = absU_fine(idx0);  up = absU_fine(idxp);
    denom = um - 2*u0 + up;
    denom(abs(denom) < eps) = eps;                    % avoid div-by-zero on flat peaks
    delta = 0.5 * (um - up) ./ denom;                 % stays in [-1, 1] for well-behaved peaks
    delta = max(-1, min(1, delta));
    k_cont = k_all + delta;                           % continuous peak row

    % Evaluate U_fine and DUX_fine at the continuous peak via 1-D linear interp
    k_floor = floor(k_cont);  k_ceil = min(k_floor+1, n_fine);  frac = k_cont - k_floor;
    idxf = sub2ind(size(U_fine), k_floor, jj);
    idxc = sub2ind(size(U_fine), k_ceil,  jj);
    umax_vec  = (1-frac) .* U_fine(idxf)  + frac .* U_fine(idxc);
    dumax_vec = (1-frac) .* DUX_fine(idxf) + frac .* DUX_fine(idxc);

    % Compute alpha a posteriori for all stations at once; guard against 0/0
    % producing NaN at stations with no disturbance (umax = 0).
    umax_safe = umax_vec;  umax_safe(~(abs(umax_vec) > 0)) = 1;   % 0 or NaN -> 1
    StabRes.alpha(j,:) = (dumax_vec ./ umax_safe) / iu;

    % Normalize shape functions and compute amplitudes (vectorized over stations)
    absU      = abs(umax_vec);                       % 1 x nx (amplitude; 0 where no disturbance)
    absU_safe = absU;  absU_safe(~(absU > 0)) = 1;   % avoid 0/0 in the shape normalization
    StabRes.u(j,:,:) = squeeze(StabRes.u(j,:,:)) ./ absU_safe;
    StabRes.v(j,:,:) = squeeze(StabRes.v(j,:,:)) ./ absU_safe;
    StabRes.w(j,:,:) = squeeze(StabRes.w(j,:,:)) ./ absU_safe;
    StabRes.p(j,:,:) = squeeze(StabRes.p(j,:,:)) ./ absU_safe;
    if isfield(StabRes,'rho')
        StabRes.rho(j,:,:) = squeeze(StabRes.rho(j,:,:)) ./ absU_safe;
        StabRes.T(j,:,:)   = squeeze(StabRes.T(j,:,:))   ./ absU_safe;
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
        scale = StabRes.A(j, :);          % mean-flow distortion (MFD): A = |u_max|
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

%% Mark the amplitude convention
% StabRes.A is now the physical peak 2|u_hat| (the (0,0) mode keeps |u_hat|:
% it is its own conjugate). A checkpoint written by Opt.save_checkpoint is
% saved before this function runs, so its StabRes.A is |u_hat| and it carries
% no A_physical field. Test isfield(StabRes,'A_physical') to tell them apart.
StabRes.A_physical = true;

%% Remove working arrays
StabRes = rmfield(StabRes, {'q', 'phi', 'R'});

end
