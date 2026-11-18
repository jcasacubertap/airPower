function Grid = build_stab_grid(BF, Grid)
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
% Build the DeHNSSo StabGrid struct from a resampled BF.
%
%   Cartesian (default): BF.X [1 x nx], BF.Y [ny x 1], BF.U/V/W [ny x nx], with
%   optional Malik parameters BF.H and BF.y_i estimated from the data if absent.
%   The metrics are the analytic identity transform.
%
%   Curvilinear (BF.mesh_type = 'curvilinear'): BF is the output of
%   resample_curvilinear or resample_unstructured, and the metrics come from
%   finite differences.
%
%   A compressible BF is detected from BF.rho and BF.T: both are copied to Grid
%   and Sutherland's law fills Grid.mu, Grid.dmuT, Grid.d2muT2 for the energy
%   equation and LHS_compressible. BF.Ma, BF.Pr, BF.gamma and BF.Sstar pass through.

%% Inputs
% Name size         unit   explanation
% BF   (1,1) struct [-]    base flow on its own grid
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

%% Outputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

curvilinear = isfield(BF, 'mesh_type') && strcmp(BF.mesh_type, 'curvilinear');

%% Reynolds number
if isfield(BF, 'Re') && ~isempty(BF.Re)
    Grid.Re = BF.Re;
elseif ~isfield(Grid, 'Re') || isempty(Grid.Re)
    error('build_stab_grid: Reynolds number not provided. Set ctype_params.Re or BF.Re.');
end

% Pass through reference scales if present (needed by downstream callers
% to convert physical λ, f into non-dim β_0, ω_0).
for f_name = {'lref','Uref','nu'}
    if isfield(BF, f_name{1}) && ~isempty(BF.(f_name{1}))
        Grid.(f_name{1}) = BF.(f_name{1});
    end
end

%% Dimensions and Malik parameters
if curvilinear
    Grid.nx  = BF.nx;
    Grid.ny  = BF.ny;
    Grid.H   = BF.H;
    Grid.y_i = BF.y_i;
    Grid.x   = BF.xi1D(:)' + Grid.Re;   % offset so x(1) = Re (ILST convention)
    Grid.y   = BF.eta1D(:);
else
    Grid.nx = length(BF.X);
    Grid.ny = length(BF.Y);
    if isfield(BF, 'H') && ~isempty(BF.H)
        Grid.H = BF.H;
    else
        Grid.H = max(BF.Y);
    end
    if isfield(BF, 'y_i') && ~isempty(BF.y_i)
        Grid.y_i = BF.y_i;
    else
        Grid.y_i = Grid.H / 20;
    end
    Grid.x = BF.X(:)';
    % Use the computational η (wall-distance, starts at 0) if BF carries it
    % — this is what the Malik mapping in compute_cheb_operators expects.
    % Fall back to BF.Y for legacy cartesian inputs where Y already starts at 0.
    if isfield(BF, 'eta1D') && ~isempty(BF.eta1D)
        Grid.y = BF.eta1D(:);
    else
        Grid.y = BF.Y(:);
    end
end

Grid.xi1D    = Grid.x;
Grid.eta1D   = Grid.y;
Grid.dxi     = Grid.x(2) - Grid.x(1);
Grid.S       = Grid.x(1);
Grid.L       = Grid.x(end) - Grid.x(1);
Grid.xun     = Grid.x;          % unmapped streamwise coordinate
Grid.etaun   = Grid.y;          % unmapped wall-normal coordinate
% Solid on the wall. Grid.step is the only thing written
if isfield(BF, 'step') && ~isempty(BF.step)
    Grid.step = BF.step;
end
Grid = step_geometry(Grid);

%% Clustering metadata: centre(s) into the xun frame (a pure shift of xi1D)
if isfield(BF, 'cluster') && ~isempty(BF.cluster) ...
        && isfield(BF.cluster, 'mug') && ~isempty(BF.cluster.mug)
    Grid.cluster    = BF.cluster;
    Grid.cluster.xc = BF.cluster.xc_xi + (Grid.x(1) - BF.xi1D(1));
end

%% Build wall-normal differentiation operators (D1, D2)
% Used consistently for both metric tensors and the stability operator.
use_fd = isfield(Grid, 'FD_eta_method') && strcmpi(Grid.FD_eta_method, 'fd');
% ('cheb' or empty → Chebyshev; 'fd' → sparse FD)
if use_fd
    if ~isfield(Grid, 'FD_eta_order') || isempty(Grid.FD_eta_order)
        Grid.FD_eta_order = 4;
    end
    fprintf('Wall-normal method: FD (order %d, sparse)\n', Grid.FD_eta_order);
    Grid = compute_fd_operators(Grid);
    fprintf('D1/D2 sparsity: %.1f%% nonzeros\n', 100 * nnz(Grid.D1) / numel(Grid.D1));
else
    fprintf('Wall-normal method: Chebyshev (spectral, dense)\n');
    Grid = compute_cheb_operators(Grid);
end

%% Metric tensors
if curvilinear
    xi_2d = repmat(Grid.x, Grid.ny, 1);
    X = BF.X;
    Y = BF.Y;

    % Persist the 2-D physical coordinates (useful for post-processing /
    % debug plots). Cartesian case keeps the 1-D Grid.x, Grid.y only.
    Grid.X = X;
    Grid.Y = Y;

    % First derivatives of forward mapping (xi,eta) -> (x,y)
    % xi-direction: FD along rows (streamwise)
    xxi = zeros(Grid.ny, Grid.nx);   % dx/dxi
    yxi = zeros(Grid.ny, Grid.nx);   % dy/dxi
    for j = 1:Grid.ny
        xxi(j,:) = deriv1_uneven(xi_2d(j,:)', X(j,:)', Grid.FD_xi_order_1)';
        yxi(j,:) = deriv1_uneven(xi_2d(j,:)', Y(j,:)', Grid.FD_xi_order_1)';
    end
    % eta-direction: D1 (Chebyshev or FD, consistent with stability operator)
    xeta = Grid.D1 * X;   % dx/deta  [ny x nx]
    yeta = Grid.D1 * Y;   % dy/deta  [ny x nx]

    J = xxi .* yeta - xeta .* yxi;
    if any(abs(J(:)) < eps)
        warning('build_stab_grid: near-zero Jacobian detected');
    end

    % First-order inverse metrics
    Grid.xix  =  yeta ./ J;
    Grid.etax = -yxi  ./ J;
    Grid.xiy  = -xeta ./ J;
    Grid.etay =  xxi  ./ J;

    % Second derivatives of forward mapping
    % xi-direction: FD along rows
    xxi_xi  = zeros(Grid.ny, Grid.nx);   % d2x/dxi2
    yxi_xi  = zeros(Grid.ny, Grid.nx);   % d2y/dxi2
    xxi_eta = zeros(Grid.ny, Grid.nx);   % d2x/dxi deta  (= d/dxi of xeta)
    yxi_eta = zeros(Grid.ny, Grid.nx);   % d2y/dxi deta  (= d/dxi of yeta)
    for j = 1:Grid.ny
        xxi_xi(j,:)  = deriv2_uneven(xi_2d(j,:)', X(j,:)',    Grid.FD_xi_order_2)';
        yxi_xi(j,:)  = deriv2_uneven(xi_2d(j,:)', Y(j,:)',    Grid.FD_xi_order_2)';
        xxi_eta(j,:) = deriv1_uneven(xi_2d(j,:)', xeta(j,:)', Grid.FD_xi_order_1)';
        yxi_eta(j,:) = deriv1_uneven(xi_2d(j,:)', yeta(j,:)', Grid.FD_xi_order_1)';
    end
    % eta-direction: spectral Chebyshev D2
    xeta_eta = Grid.D2 * X;   % d2x/deta2  [ny x nx]
    yeta_eta = Grid.D2 * Y;   % d2y/deta2  [ny x nx]

    % Second-order inverse metrics — exact formula via implicit differentiation
    % Differentiating xi(x(xi,eta), y(xi,eta)) = xi twice yields a 3x3 system
    % [xixx, xixy, xiyy] with coefficient matrix det = -J^3.  Likewise for eta.
    % Let a=xxi, b=xeta, c=yxi, d=yeta.
    a = xxi; b = xeta; c = yxi; d = yeta;
    J3 = J.^3;

    Grid.xixx  = (-d.^3.*xxi_xi  + b.*d.^2.*yxi_xi  - c.^2.*d.*xeta_eta + b.*c.^2.*yeta_eta + 2*c.*d.^2.*xxi_eta - 2*b.*c.*d.*yxi_eta) ./ J3;
    Grid.xiyy  = (-b.^2.*d.*xxi_xi + b.^3.*yxi_xi   - a.^2.*d.*xeta_eta + a.^2.*b.*yeta_eta + 2*a.*b.*d.*xxi_eta - 2*a.*b.^2.*yxi_eta) ./ J3;
    Grid.etaxx = ( c.*d.^2.*xxi_xi - a.*d.^2.*yxi_xi + c.^3.*xeta_eta   - a.*c.^2.*yeta_eta - 2*c.^2.*d.*xxi_eta + 2*a.*c.*d.*yxi_eta) ./ J3;
    Grid.etayy = ( b.^2.*c.*xxi_xi - a.*b.^2.*yxi_xi + a.^2.*c.*xeta_eta - a.^3.*yeta_eta   - 2*a.*b.*c.*xxi_eta + 2*a.^2.*b.*yxi_eta) ./ J3;

    % Mixed inverse metrics, needed only by the compressible operator's explicit
    % d_x d_y cross-derivatives; the incompressible one uses just the Laplacian.
    % Obtained by chain-ruling the closed-form first-order metrics:
    %   xixy = xiy.*d_dxi(xix) + etay.*d_deta(xix),  likewise etaxy from etax.
    dxix_dxi  = zeros(Grid.ny, Grid.nx);   dxix_deta  = Grid.D1 * Grid.xix;
    detax_dxi = zeros(Grid.ny, Grid.nx);   detax_deta = Grid.D1 * Grid.etax;
    for j = 1:Grid.ny
        dxix_dxi(j,:)  = deriv1_uneven(xi_2d(j,:)', Grid.xix(j,:)',  Grid.FD_xi_order_1)';
        detax_dxi(j,:) = deriv1_uneven(xi_2d(j,:)', Grid.etax(j,:)', Grid.FD_xi_order_1)';
    end
    Grid.xixy  = Grid.xiy .* dxix_dxi  + Grid.etay .* dxix_deta;
    Grid.etaxy = Grid.xiy .* detax_dxi + Grid.etay .* detax_deta;
else
    Grid.xix   = ones(Grid.ny,  Grid.nx);
    Grid.etay  = ones(Grid.ny,  Grid.nx);
    Grid.etax  = zeros(Grid.ny, Grid.nx);
    Grid.xiy   = zeros(Grid.ny, Grid.nx);
    Grid.xixx  = zeros(Grid.ny, Grid.nx);
    Grid.etaxx = zeros(Grid.ny, Grid.nx);
    Grid.xiyy  = zeros(Grid.ny, Grid.nx);
    Grid.etayy = zeros(Grid.ny, Grid.nx);
    Grid.xixy  = zeros(Grid.ny, Grid.nx);
    Grid.etaxy = zeros(Grid.ny, Grid.nx);
end

%% Base flow and derivatives
Grid.U = BF.U;
Grid.V = BF.V;
Grid.W = BF.W;

% Compressible base flow: copy rho, T and thermo metadata, evaluate Sutherland.
is_compr = isfield(BF, 'rho') && isfield(BF, 'T');
if is_compr
    Grid.rho = BF.rho;
    Grid.T   = BF.T;
    for f_name = {'Ma','Pr','gamma','Sstar'}
        if isfield(BF, f_name{1}) && ~isempty(BF.(f_name{1}))
            Grid.(f_name{1}) = BF.(f_name{1});
        end
    end
    if isfield(Grid, 'Sstar') && ~isempty(Grid.Sstar)
        [Grid.mu, Grid.dmuT, Grid.d2muT2] = sutherland(Grid.T, Grid.Sstar);
    end
end

for dn = bf_deriv_fields()
    if isfield(BF, dn{1}) && ~isempty(BF.(dn{1}))
        Grid.(dn{1}) = BF.(dn{1});
    end
end

Grid = compute_bf_derivatives(Grid);

if curvilinear
    fprintf('Curvilinear grid ready: %d streamwise x %d wall-normal\n', Grid.nx, Grid.ny);
else
    fprintf('Cartesian grid ready: %d streamwise x %d wall-normal\n', Grid.nx, Grid.ny);
end
if is_compr
    fprintf('Compressible BF: ρ, T, μ, ∂μ/∂T, ∂²μ/∂T² populated (Ma=%.3f, Pr=%.3f)\n', ...
        getfielddef(Grid,'Ma',NaN), getfielddef(Grid,'Pr',NaN));
end
end

function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
