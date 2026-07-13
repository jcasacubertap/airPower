function Grid = compute_bf_derivatives(Grid)
%COMPUTE_BF_DERIVATIVES  First **physical** derivatives of base-flow velocities U, V, W.
%
%   Returns ∂U/∂x and ∂U/∂y (and likewise for V, W) stored as Grid.dxU, Grid.dyU, etc.
%   These are the physical-space derivatives the HNS LHS assembly expects, not
%   the raw computational-space derivatives ∂U/∂ξ, ∂U/∂η.
%
%   On curvilinear grids the chain rule is applied using the inverse metric
%   tensors set earlier in build_stab_grid:
%       ∂U/∂x = ξ_x · ∂U/∂ξ + η_x · ∂U/∂η
%       ∂U/∂y = ξ_y · ∂U/∂ξ + η_y · ∂U/∂η
%   On Cartesian grids (ξ_x = η_y = 1, ξ_y = η_x = 0) this reduces to the
%   identity: ∂U/∂x = ∂U/∂ξ, ∂U/∂y = ∂U/∂η.
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

% ---- Computational-space derivatives --------------------------------------
dU_deta = D1_bf * Grid.U;
dV_deta = D1_bf * Grid.V;
dW_deta = D1_bf * Grid.W;

xi = Grid.xi1D(:);
dU_dxi = zeros(Grid.ny, Grid.nx);
dV_dxi = zeros(Grid.ny, Grid.nx);
dW_dxi = zeros(Grid.ny, Grid.nx);
for j = 1:Grid.ny
    dU_dxi(j,:) = FD1_uneven(xi, Grid.U(j,:)', Grid.FD1_order)';
    dV_dxi(j,:) = FD1_uneven(xi, Grid.V(j,:)', Grid.FD1_order)';
    dW_dxi(j,:) = FD1_uneven(xi, Grid.W(j,:)', Grid.FD1_order)';
end

% ---- Chain rule to physical coordinates -----------------------------------
% Metrics live on the same (ny × nx) grid as the BF arrays.
Grid.dxU = Grid.xix .* dU_dxi + Grid.etax .* dU_deta;
Grid.dyU = Grid.xiy .* dU_dxi + Grid.etay .* dU_deta;

Grid.dxV = Grid.xix .* dV_dxi + Grid.etax .* dV_deta;
Grid.dyV = Grid.xiy .* dV_dxi + Grid.etay .* dV_deta;

Grid.dxW = Grid.xix .* dW_dxi + Grid.etax .* dW_deta;
Grid.dyW = Grid.xiy .* dW_dxi + Grid.etay .* dW_deta;
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
