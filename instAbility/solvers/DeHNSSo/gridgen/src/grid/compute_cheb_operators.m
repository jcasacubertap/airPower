function Grid = compute_cheb_operators(Grid)
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
% Chebyshev D1, D2 with Malik wall-normal mapping.
%
%   Malik mapping:  z = -1 + (1+b)*eta / (a + eta)
%   where  a = H*y_i / (H - 2*y_i),  b = 1 + 2*a/H
%
%   D1 = diag(dz/deta) * Dn        ->  d/d(eta)
%   D2 = diag(d2z/deta2)*Dn + diag((dz/deta)^2)*D2n

%% Inputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

%% Outputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms

a = Grid.H * Grid.y_i / (Grid.H - 2.0 * Grid.y_i);
b = 1.0 + 2.0 * a / Grid.H;

dzdy   = zeros(Grid.ny, 1);
d2zdy2 = zeros(Grid.ny, 1);
for k = 1:Grid.ny
    dzdy(k)   =  a * (1 + b) / (a + Grid.y(k))^2;
    d2zdy2(k) = -2 * a * (1 + b) / (a + Grid.y(k))^3;
end
Grid.dzdy   = dzdy;
Grid.d2zdy2 = d2zdy2;

Dn  = -chebDiff(Grid.ny);
D2n = Dn * Dn;
Grid.D1 = diag(dzdy)   * Dn;
Grid.D2 = diag(d2zdy2) * Dn + diag(dzdy.^2) * D2n;
end
