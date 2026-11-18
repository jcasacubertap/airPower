function [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
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
% Map a Chebyshev grid [-1,1] to a clustered physical grid [0,y_max].
%
%   [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
%
%   Malik's algebraic mapping, with half the nodes below the median height y_i.
%   Also returns the physical differentiation matrices for the mapped grid.
%
%     D1eta, D2eta  — Chebyshev differentiation matrices
%     D1y, D2y      — first / second physical differentiation matrices

%% Inputs
% Name  size  unit   explanation
% y_max       [-]    domain height
% y_i         [-]    median height (half the nodes lie in [0, y_i])
% eta         [-]    Chebyshev nodes on [-1, 1]
% D1eta (n,n) [-]    Chebyshev first-derivative matrix on the computational grid
% D2eta (n,n) [-]    Chebyshev second-derivative matrix on the computational grid

%% Outputs
% Name size  unit   explanation
% ynd        [-]    physical wall-normal grid [0, y_max]
% D1y  (n,n) [-]    first-derivative matrix on the mapped physical grid
% D2y  (n,n) [-]    second-derivative matrix on the mapped physical grid

ynd = (y_i*y_max*(1+eta) ./ (y_max - eta*(y_max - 2*y_i)))';

if nargout > 1
    detady   =  (y_max - eta*(y_max - 2*y_i)).^2 / (2*y_i*y_max*(y_max - y_i));
    d2etady2 = -2*detady.^2*(y_max - 2*y_i) ./ (y_max - eta*(y_max - 2*y_i));
    D1y = diag(detady)*D1eta;
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end
end
