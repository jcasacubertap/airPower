function [X_out, Y_out] = solve_elliptic_mesh(X_init, Y_init, ...
    X_wall, Y_wall, nx_w, ny_w, eta_cheb, H, opts)
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
% Smooth a structured mesh via elliptic (Poisson) PDE.
%
%   Takes an algebraic (ray-shot) mesh as the initial guess, iterates the Poisson
%   equations on a uniform computational grid with point-by-point SOR, then
%   redistributes the points to Malik-Chebyshev spacing.
%
%   In: X_init, Y_init [n_eta x n_xi]; the wall X_wall, Y_wall and its outward
%   normals nx_w, ny_w [1 x n_xi]; eta_cheb [n_eta x 1], descending; H the domain
%   height. Out: X_out, Y_out on the Malik-Chebyshev nodes.
%
%   opts: max_iter (5000), tol (1e-6), omega SOR factor (1.7), n_uniform eta
%   points for the solve (201), P_strength source strength (0 = Laplace),
%   P_decay source decay (5).

%% Inputs
% Name     size         unit   explanation
% X_init   (n_eta,n_xi) [-]    algebraic ray-shot mesh, x coordinate, used as the initial guess
% Y_init   (n_eta,n_xi) [-]    algebraic ray-shot mesh, y coordinate, used as the initial guess
% X_wall   (1,n_xi)     [-]    wall x coordinate
% Y_wall   (1,n_xi)     [-]    wall y coordinate
% nx_w     (1,n_xi)     [-]    outward wall normal, x component
% ny_w     (1,n_xi)     [-]    outward wall normal, y component
% eta_cheb (n_eta,1)    [-]    Malik-Chebyshev wall-normal nodes, descending
% H        (1,1)        [-]    domain height
% opts     (1,1) struct [-]    solver settings: max_iter, tol, SOR factor omega, n_uniform, P_strength and P_decay

%% Outputs
% Name  size         unit   explanation
% X_out (n_eta,n_xi) [-]    smoothed mesh, x coordinate, on the Malik-Chebyshev nodes
% Y_out (n_eta,n_xi) [-]    smoothed mesh, y coordinate, on the Malik-Chebyshev nodes

%% Parse options
max_iter   = 5000;
tol        = 1e-6;
omega      = 1.7;
n_uniform  = 201;
P_strength = 0;
P_decay    = 5;

if nargin >= 9 && isstruct(opts)
    if isfield(opts, 'max_iter'),   max_iter   = opts.max_iter;   end
    if isfield(opts, 'tol'),        tol        = opts.tol;        end
    if isfield(opts, 'omega'),      omega      = opts.omega;      end
    if isfield(opts, 'n_uniform'),  n_uniform  = opts.n_uniform;  end
    if isfield(opts, 'P_strength'), P_strength = opts.P_strength; end
    if isfield(opts, 'P_decay'),    P_decay    = opts.P_decay;    end
end

[n_eta, n_xi] = size(X_init);

%% Step 1: Interpolate algebraic mesh onto uniform eta grid

eta_asc  = flip(eta_cheb);
eta_uni  = linspace(0, H, n_uniform)';
M = n_uniform;
N = n_xi;

X_uni = zeros(M, N);
Y_uni = zeros(M, N);

for j = 1:N
    X_uni(:,j) = interp1(eta_asc, flip(X_init(:,j)), eta_uni, 'pchip');
    Y_uni(:,j) = interp1(eta_asc, flip(Y_init(:,j)), eta_uni, 'pchip');
end

% Convention: row 1 = wall (eta=0), row M = outer (eta=H)

%% Step 2: Fix boundary conditions
% Wall (eta=0) is pinned to the wall locus. The outer boundary (eta=H) is the
% wall offset along the MEAN outward normal, not the local one: the local-normal
% offset self-intersects on concave walls, leaving a non-simple boundary that
% forces the harmonic map to keep folds. A mean-normal offset stays simple for
% graph-like walls, so the Laplace solution is a diffeomorphism and the folds go.
% Near-wall orthogonality comes from the wall pinning and the Malik-Chebyshev
% redistribution, not from the outer-boundary shape.
X_uni(1,:) = X_wall;
Y_uni(1,:) = Y_wall;
nm = [mean(nx_w), mean(ny_w)];  nm = nm / norm(nm);
X_uni(M,:) = X_wall + H * nm(1);
Y_uni(M,:) = Y_wall + H * nm(2);

%% Precompute Poisson source term Q(eta)
eta_norm = eta_uni / H;
Q_src = -P_strength * exp(-P_decay * eta_norm);

if P_strength > 0
    fprintf('Poisson source: strength=%.2f, decay=%.2f\n', P_strength, P_decay);
else
    fprintf('Laplace mode (no source terms)\n');
end

%% Step 3: SOR iteration on Poisson equations
%
%  alpha * x_ee - 2*beta * x_xe + gamma * x_xx + J^2 * Q * x_e = 0
%
%  Point-wise update (solving for x(i,j)):
%    x_new = [ alpha*(x(i,j+1)+x(i,j-1))
%            + gamma*(x(i+1,j)+x(i-1,j))
%            - beta/2*(x(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1))
%            + J^2*Q/2*(x(i+1,j)-x(i-1,j)) ]
%            / [ 2*(alpha+gamma) ]
%
%  SOR: x(i,j) = (1-omega)*x(i,j) + omega*x_new

fprintf('Elliptic mesh smoothing (SOR): M=%d, N=%d, omega=%.2f\n', M, N, omega);
fprintf('  Iterating');

res_history = zeros(max_iter, 1);

for iter = 1:max_iter
    X_old = X_uni;
    Y_old = Y_uni;

    for i = 2:M-1
        for j = 2:N-1
            % Metric coefficients
            X_x = 0.5 * (X_uni(i,j+1) - X_uni(i,j-1));
            Y_x = 0.5 * (Y_uni(i,j+1) - Y_uni(i,j-1));
            X_e = 0.5 * (X_uni(i+1,j) - X_uni(i-1,j));
            Y_e = 0.5 * (Y_uni(i+1,j) - Y_uni(i-1,j));

            al = X_e^2 + Y_e^2;       % alpha
            be = X_x*X_e + Y_x*Y_e;   % beta
            ga = X_x^2 + Y_x^2;       % gamma
            J2 = (X_x*Y_e - X_e*Y_x)^2;

            denom = 2 * (al + ga);
            if denom < eps; continue; end

            % Cross-derivative term
            cross_X = 0.5 * be * ( ...
                X_uni(i+1,j+1) - X_uni(i+1,j-1) ...
              - X_uni(i-1,j+1) + X_uni(i-1,j-1));
            cross_Y = 0.5 * be * ( ...
                Y_uni(i+1,j+1) - Y_uni(i+1,j-1) ...
              - Y_uni(i-1,j+1) + Y_uni(i-1,j-1));

            % Source term: J^2 * Q * x_e  (central difference for x_e)
            src_X = 0.5 * J2 * Q_src(i) * (X_uni(i+1,j) - X_uni(i-1,j));
            src_Y = 0.5 * J2 * Q_src(i) * (Y_uni(i+1,j) - Y_uni(i-1,j));

            % Point-wise Gauss-Seidel update
            X_new = ( al * (X_uni(i,j+1) + X_uni(i,j-1)) ...
                    + ga * (X_uni(i+1,j) + X_uni(i-1,j)) ...
                    - cross_X + src_X ) / denom;
            Y_new = ( al * (Y_uni(i,j+1) + Y_uni(i,j-1)) ...
                    + ga * (Y_uni(i+1,j) + Y_uni(i-1,j)) ...
                    - cross_Y + src_Y ) / denom;

            % SOR relaxation
            X_uni(i,j) = (1 - omega) * X_uni(i,j) + omega * X_new;
            Y_uni(i,j) = (1 - omega) * Y_uni(i,j) + omega * Y_new;
        end
    end

    % Check convergence
    res = max(abs(X_uni(:) - X_old(:))) + max(abs(Y_uni(:) - Y_old(:)));
    res_history(iter) = res;

    if mod(iter, 500) == 0
        fprintf('.');
    end

    if res < tol
        fprintf('\n  Converged at iteration %d (residual = %.2e)\n', iter, res);
        break
    end
end

n_iter = min(iter, max_iter);
res_history = res_history(1:n_iter);

if iter == max_iter
    fprintf('\n  Warning: did not converge after %d iterations (residual = %.2e)\n', max_iter, res);
end

%% Step 4: Jacobian quality check
Jac = zeros(M-2, N-2);
for i = 2:M-1
    for j = 2:N-1
        X_x = 0.5 * (X_uni(i,j+1) - X_uni(i,j-1));
        Y_x = 0.5 * (Y_uni(i,j+1) - Y_uni(i,j-1));
        X_e = 0.5 * (X_uni(i+1,j) - X_uni(i-1,j));
        Y_e = 0.5 * (Y_uni(i+1,j) - Y_uni(i-1,j));
        Jac(i-1,j-1) = X_x * Y_e - X_e * Y_x;
    end
end

fprintf('  Jacobian: min = %.4e, mean = %.4e\n', min(Jac(:)), mean(Jac(:)));
if any(Jac(:) < 0)
    warning('solve_elliptic_mesh: %d points with negative Jacobian (mesh crossing).', sum(Jac(:) < 0));
end

%% Step 5: Redistribute to Malik-Chebyshev spacing
%  The elliptic solve lives on eta_uni (uniform 0..H, ascending).
%  We need values at eta_cheb (Malik-Chebyshev, descending H..0).
%  Direct interpolation in eta ensures the computational coordinate
%  matches exactly what the Chebyshev D1/D2 operators expect.

eta_cheb_asc = flip(eta_cheb);   % ascending for interp1

X_out = zeros(n_eta, n_xi);
Y_out = zeros(n_eta, n_xi);

for j = 1:N
    X_out(:,j) = flip(interp1(eta_uni, X_uni(:,j), eta_cheb_asc, 'pchip'));
    Y_out(:,j) = flip(interp1(eta_uni, Y_uni(:,j), eta_cheb_asc, 'pchip'));
end

fprintf('  Redistribution to Malik-Chebyshev nodes done.\n');

end

%% =========================================================================
function x = tdma(a, b, c, d)
%TDMA  Thomas tridiagonal matrix algorithm.
n = length(b);
c_star = zeros(n, 1);
d_star = zeros(n, 1);

c_star(1) = c(1) / b(1);
d_star(1) = d(1) / b(1);

for i = 2:n
    m = b(i) - a(i) * c_star(i-1);
    c_star(i) = c(i) / m;
    d_star(i) = (d(i) - a(i) * d_star(i-1)) / m;
end

x = zeros(n, 1);
x(n) = d_star(n);
for i = n-1:-1:1
    x(i) = d_star(i) - c_star(i) * x(i+1);
end
end
