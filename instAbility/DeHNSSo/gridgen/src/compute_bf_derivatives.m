function Grid = compute_bf_derivatives(Grid)
%COMPUTE_BF_DERIVATIVES  First (and second) physical derivatives of base-flow fields.
%
%   Always differentiates the velocity components U, V, W and stores
%   ∂U/∂x, ∂U/∂y, … as Grid.dxU, Grid.dyU, etc.
%   If Grid.rho / Grid.T are present (compressible base flow), it also
%   differentiates them and stores Grid.dxrho/dyrho, Grid.dxT/dyT.
%   These are the physical-space derivatives the HNS LHS assembly expects,
%   not the raw computational-space derivatives ∂·/∂ξ, ∂·/∂η.
%
%   Second streamwise derivatives ∂²U/∂x², ∂²V/∂x², ∂²W/∂x², ∂²T/∂x² are
%   also computed (Grid.dxxU, …) when the corresponding field is present
%   — required by LHS_compressible for the viscous terms.
%   The chain rule is currently Cartesian-only (xixx = etaxx = 0); on
%   curvilinear grids a warning is emitted and ∂²/∂ξ² is stored verbatim.
%
%   On curvilinear grids the first-derivative chain rule is applied using
%   the inverse metric tensors set earlier in build_stab_grid:
%       ∂·/∂x = ξ_x · ∂·/∂ξ + η_x · ∂·/∂η
%       ∂·/∂y = ξ_y · ∂·/∂ξ + η_y · ∂·/∂η
%   On Cartesian grids (ξ_x = η_y = 1, ξ_y = η_x = 0) this reduces to the
%   identity: ∂·/∂x = ∂·/∂ξ, ∂·/∂y = ∂·/∂η.
%
%   Discretisation:
%     ∂/∂η : D1_bf — Grid.D1 in the interior, local one-sided Fornberg
%            stencils at the wall and freestream rows. The dense Chebyshev
%            boundary rows in Grid.D1 are necessary for the stability
%            operator (spectral accuracy at the (0,0) MFD mode), but they
%            ring on tiny non-smoothness in interpolated CFD base flows
%            (linear scattered interp has piecewise-constant slope) and
%            blow up the wall shear by O(N²). A local FD stencil avoids
%            this without affecting the stability operator.
%     ∂/∂ξ : FD1_uneven on the (possibly non-uniform) Grid.xi1D
%     ∂²/∂ξ²: FD2_uneven (order Grid.FD_xi_order_2, default 4)

% ---- Boundary-safe D1 for base-flow differentiation -----------------------
N     = Grid.ny;
y     = Grid.y(:);
p_bd  = 4;                              % 5-point one-sided stencil at boundaries
D1_bf = Grid.D1;
% Wall row (i=N): use one-sided stencil from the wall up
idx_w = (N - p_bd) : N;
w_w   = fornberg_weights_loc(y(N), y(idx_w), 1);
D1_bf(N, :)     = 0;
D1_bf(N, idx_w) = w_w;
% Freestream row (i=1): one-sided stencil from the top down
idx_t = 1 : (1 + p_bd);
w_t   = fornberg_weights_loc(y(1), y(idx_t), 1);
D1_bf(1, :)     = 0;
D1_bf(1, idx_t) = w_t;

% ---- Differentiate each base-flow field present on Grid -------------------
fields = {'U','V','W'};
if isfield(Grid, 'rho'); fields{end+1} = 'rho'; end
if isfield(Grid, 'T');   fields{end+1} = 'T';   end

xi             = Grid.xi1D(:);
FD_xi_order_2  = getfielddef(Grid, 'FD_xi_order_2', 4);   % 4th-order matches heritage FD1d4o∘FD1d4o
% Streamwise BF second derivatives are only needed by LHS_compressible;
% skip them entirely for incompressible BF (no ρ, no T on Grid).
is_compr       = isfield(Grid, 'rho') && isfield(Grid, 'T');
if is_compr
    need_dxx = {'U','V','W','T'};
else
    need_dxx = {};
end
is_curvilinear = isfield(Grid, 'xixx') && any(Grid.xixx(:) ~= 0);

for k = 1:numel(fields)
    fname = fields{k};
    F     = Grid.(fname);

    % Computational-space first derivatives
    dF_deta = D1_bf * F;
    dF_dxi  = zeros(Grid.ny, Grid.nx);
    for j = 1:Grid.ny
        dF_dxi(j,:) = FD1_uneven(xi, F(j,:)', Grid.FD_xi_order_1)';
    end

    % First derivatives in physical coordinates (Cartesian: identity)
    Grid.(['dx' fname]) = Grid.xix .* dF_dxi + Grid.etax .* dF_deta;
    Grid.(['dy' fname]) = Grid.xiy .* dF_dxi + Grid.etay .* dF_deta;

    % Streamwise second derivative ∂²F/∂x² (Cartesian only; consumed by LHS_compressible)
    if any(strcmp(fname, need_dxx))
        d2F_dxi2 = zeros(Grid.ny, Grid.nx);
        for j = 1:Grid.ny
            d2F_dxi2(j,:) = FD2_uneven(xi, F(j,:)', FD_xi_order_2)';
        end
        Grid.(['dxx' fname]) = d2F_dxi2;
    end
end

if is_curvilinear && is_compr
    warning('compute_bf_derivatives:CurvilinearDxx', ...
        'Streamwise BF second derivatives (dxxU/...) stored as ∂²/∂ξ² — full curvilinear chain rule not yet implemented.');
end
end


function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end


function w = fornberg_weights_loc(z, x, m)
% Local copy of the Fornberg weights helper from compute_fd_operators.m.
% Returns a row vector of weights for the m-th derivative at z on stencil x.
n  = numel(x) - 1;
c  = zeros(n+1, m+1);
c1 = 1;
c4 = x(1) - z;
c(1,1) = 1;

for i = 1:n
    mn  = min(i, m);
    c2  = 1;
    c5  = c4;
    c4  = x(i+1) - z;

    for j = 0:i-1
        c3 = x(i+1) - x(j+1);
        c2 = c2 * c3;

        for k = mn:-1:1
            c(i+1, k+1) = c1 * (k * c(i, k) - c5 * c(i, k+1)) / c2;
        end
        c(i+1, 1) = -c1 * c5 * c(i, 1) / c2;

        for k = mn:-1:1
            c(j+1, k+1) = (c4 * c(j+1, k+1) - k * c(j+1, k)) / c3;
        end
        c(j+1, 1) = c4 * c(j+1, 1) / c3;
    end

    c1 = c2;
end

w = c(:, m+1)';
end
