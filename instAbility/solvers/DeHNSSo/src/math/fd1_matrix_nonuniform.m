function D = fd1_matrix_nonuniform(x)
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
% Sparse 4th-order 1st-derivative matrix for a non-uniform 1-D grid.
%
%   D = fd1_matrix_nonuniform(x)   with x a column vector of n stations.
%
%   D is an (n × n) sparse matrix such that y = D*f gives a 4th-order-accurate
%   estimate of df/dx at each station, including near the boundaries where
%   the stencil becomes biased. On a uniform grid D reduces to the classical
%   4th-order finite-difference coefficients. Stencils mirror those used by
%   the HNS LHS assembly:
%       i = 1        : 5-point forward, x(1:5)
%       i = 2        : 5-point forward-biased, x(1:5)
%       3 ≤ i ≤ n-2  : 5-point centred
%       i = n-1      : 4-point backward-biased, x(n-3:n)
%       i = n        : 3-point backward, x(n-2:n)

%% Inputs
% Name size  unit   explanation
% x    (n,1) [-]    station coordinates, need not be uniformly spaced

%% Outputs
% Name size         unit   explanation
% D    (n,n) sparse [-]    fourth-order first-derivative matrix, biased near the boundaries

%% Validate input
x = x(:);
n = length(x);
if n < 5
    error('fd1_matrix_nonuniform: need at least 5 stations (got %d).', n);
end

%% Assemble stencil weights row by row
% Triplet lists (row, col, val); at most 5 non-zeros per row.
nnz_est = 5 * n;
row = zeros(nnz_est, 1);
col = zeros(nnz_est, 1);
val = zeros(nnz_est, 1);
k   = 0;

for i = 1:n
    if i == 1
        cols  = 1:5;
    elseif i == 2
        cols  = 1:5;
    elseif i <= n-2
        cols  = (i-2) : (i+2);
    elseif i == n-1
        cols  = (n-3) : n;
    else
        cols  = (n-2) : n;
    end
    % Fornberg weights for the 1st derivative (column 2 of c) at station i
    c  = fornberg_weights(x(i), x(cols), 1);
    w1 = c(:, 2).';
    idx = k + (1:numel(cols));
    row(idx) = i;
    col(idx) = cols;
    val(idx) = w1;
    k = k + numel(cols);
end

%% Build sparse derivative matrix
row = row(1:k); col = col(1:k); val = val(1:k);
D = sparse(row, col, val, n, n);
end
