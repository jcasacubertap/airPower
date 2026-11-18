function c = fornberg_weights(z, x, m)
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
% Finite-difference weights on non-uniform grids.
%
%   c = fornberg_weights(z, x, m)
%
%   Returns an (n × (m+1)) weight matrix c such that, for any function f
%   sampled at stencil points x(1..n),
%       f^(k)(z) ≈ sum_j c(j, k+1) * f(x(j))
%   for k = 0, 1, ..., m.
%
%   Implements B. Fornberg, "Generation of finite difference formulas on
%   arbitrarily spaced grids", Math. Comp. 51 (1988) 699-706.

%% Inputs
% Name size  unit   explanation
% z    (1,1) [-]    point at which the derivatives are approximated
% x    (n,1) [-]    stencil points, need not be uniformly spaced
% m    (1,1) [-]    highest derivative order for which weights are returned

%% Outputs
% Name size    unit   explanation
% c    (n,m+1) [-]    weights: f^(k)(z) is approximated by sum_j c(j,k+1)*f(x(j))

%% Initialise
x = x(:);
n = length(x) - 1;
c = zeros(n+1, m+1);

c1 = 1;
c4 = x(1) - z;
c(1, 1) = 1;

%% Fornberg recurrence over the stencil points
for i = 1:n
    mn = min(i, m);
    c2 = 1;
    c5 = c4;
    c4 = x(i+1) - z;
    for j = 0:i-1
        c3 = x(i+1) - x(j+1);
        c2 = c2 * c3;
        if j == i-1
            for k = mn:-1:1
                c(i+1, k+1) = c1 * (k * c(i, k) - c5 * c(i, k+1)) / c2;
            end
            c(i+1, 1) = -c1 * c5 * c(i, 1) / c2;
        end
        for k = mn:-1:1
            c(j+1, k+1) = (c4 * c(j+1, k+1) - k * c(j+1, k)) / c3;
        end
        c(j+1, 1) = c4 * c(j+1, 1) / c3;
    end
    c1 = c2;
end
end
