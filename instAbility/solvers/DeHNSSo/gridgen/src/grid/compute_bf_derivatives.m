function Grid = compute_bf_derivatives(Grid)
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
% First (and second) physical derivatives of base-flow fields.
%
%   Differentiates U, V, W (and rho, T when present) into Grid.dxU, Grid.dyU,
%   ..., the physical-space derivatives the HNS assembly expects rather than the
%   raw ∂·/∂ξ, ∂·/∂η. Second derivatives of U,V,W,T follow for the viscous terms
%   of LHS_compressible.
%
%   Chain rule with the inverse metrics set in build_stab_grid:
%       ∂·/∂x   = ξ_x ∂_ξ + η_x ∂_η
%       ∂²·/∂x² = ξ_x² ∂_ξξ + 2ξ_xη_x ∂_ξη + η_x² ∂_ηη + ξ_xx ∂_ξ + η_xx ∂_η
%   and analogously in y; both reduce to the identity on a Cartesian grid.
%
%   ∂/∂η uses D1_bf: Grid.D1 in the interior, local one-sided Fornberg stencils
%   at the wall and freestream. The dense Chebyshev boundary rows are needed by
%   the stability operator, but they ring on interpolated CFD base flows and
%   inflate the wall shear by O(N²). ∂/∂ξ uses deriv1_uneven on Grid.xi1D.

%% Inputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

%% Outputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

% ---- Boundary-safe D1 for base-flow differentiation -----------------------
% Replace only the two boundary rows of Grid.D1 by one-sided Fornberg stencils.
% The interior keeps the operator the solver uses, so this affects the base flow
% alone and never the stability operator.
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

% rho and T are differentiated only when present, so an incompressible Grid
% carries no thermodynamic derivatives at all.
xi             = Grid.xi1D(:);
FD_xi_order_2  = getfielddef(Grid, 'FD_xi_order_2', 4);   % 4th-order matches heritage FD1d4o∘FD1d4o
% Second base-flow derivatives are only needed by LHS_compressible;
% skip them entirely for incompressible BF (no ρ, no T on Grid).
is_compr       = isfield(Grid, 'rho') && isfield(Grid, 'T');

n_provided = 0;
for k = 1:numel(fields)
    fname = fields{k};
    F     = Grid.(fname);

    has_dx  = isfield(Grid, ['dx' fname]) && ~isempty(Grid.(['dx' fname]));
    has_dy  = isfield(Grid, ['dy' fname]) && ~isempty(Grid.(['dy' fname]));
    have_d1 = has_dx && has_dy;
    if xor(has_dx, has_dy)
        given = 'dx'; if has_dy; given = 'dy'; end
        warning('compute_bf_derivatives:PartialDeriv', ...
            'Base-flow field %s: %s%s provided without its pair; both dx%s and dy%s recomputed by FD.', ...
            fname, given, fname, fname, fname);
    end
    if have_d1
        n_provided = n_provided + 2;
    else
        dF_deta = D1_bf * F;
        dF_dxi  = zeros(Grid.ny, Grid.nx);
        for j = 1:Grid.ny
            dF_dxi(j,:) = deriv1_uneven(xi, F(j,:)', Grid.FD_xi_order_1)';
        end
        Grid.(['dx' fname]) = Grid.xix .* dF_dxi + Grid.etax .* dF_deta;
        Grid.(['dy' fname]) = Grid.xiy .* dF_dxi + Grid.etay .* dF_deta;
    end
end

% Second physical derivatives for the compressible viscous terms, by the chain
% rule in the header. Precomputed here so the solver never differentiates the
% base flow; read as Grid.dxxU/ddyU/dxyU/dxyV by LHS_compressible and CLST.
if is_compr
    for nm = {'U','V','W','T'}
        fn = nm{1};
        F  = Grid.(fn);
        % computational derivatives of the base-flow field
        F_eta    = D1_bf * F;                          % d/deta
        F_etaeta = Grid.D2 * F;                        % d2/deta2
        F_xi     = zeros(Grid.ny, Grid.nx);
        F_xixi   = zeros(Grid.ny, Grid.nx);
        for j = 1:Grid.ny
            F_xi(j,:)   = deriv1_uneven(xi, F(j,:)', Grid.FD_xi_order_1)';
            F_xixi(j,:) = deriv2_uneven(xi, F(j,:)', FD_xi_order_2)';
        end
        F_xieta = D1_bf * F_xi;                         % d2/dxi deta (= d/deta of d/dxi)
        % physical second derivatives (chain rule)
        provided = @(f) isfield(Grid,f) && ~isempty(Grid.(f));
        if provided(['dxx' fn]); n_provided=n_provided+1; else
            Grid.(['dxx' fn]) = Grid.xix.^2.*F_xixi + 2*Grid.xix.*Grid.etax.*F_xieta ...
                + Grid.etax.^2.*F_etaeta + Grid.xixx.*F_xi + Grid.etaxx.*F_eta;
        end
        if provided(['ddy' fn]); n_provided=n_provided+1; else
            Grid.(['ddy' fn]) = Grid.xiy.^2.*F_xixi + 2*Grid.xiy.*Grid.etay.*F_xieta ...
                + Grid.etay.^2.*F_etaeta + Grid.xiyy.*F_xi + Grid.etayy.*F_eta;
        end
        if any(strcmp(fn,{'U','V'}))
            if provided(['dxy' fn]); n_provided=n_provided+1; else
                Grid.(['dxy' fn]) = Grid.xix.*Grid.xiy.*F_xixi ...
                    + (Grid.xix.*Grid.etay + Grid.xiy.*Grid.etax).*F_xieta ...
                    + Grid.etax.*Grid.etay.*F_etaeta + Grid.xixy.*F_xi + Grid.etaxy.*F_eta;
            end
        end
    end
end

if n_provided > 0
    fprintf('  compute_bf_derivatives: using %d provided base-flow derivative field(s)\n', n_provided);
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
