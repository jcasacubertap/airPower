function [ ynd,D1y,D2y ] = MappingMalik( y_max,y_i,eta,D1eta,D2eta )
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
% Map the Chebyshev grid onto the physical wall-normal domain.
%
%   [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
%
%   Maps the Chebyshev grid on [-1, 0, 1] to the physical domain
%   [0, y_i, y_max] with the algebraic (Malik) stretching, and transforms the
%   pseudo-spectral differentiation matrices onto that grid.
%

%% Inputs
% Name  size unit   explanation
% y_max      [-]    outer domain height
% y_i        [-]    half of the nodes are clustered below this height
% eta        [-]    Chebyshev grid on [-1, 1]
% D1eta      [-]    first-derivative Chebyshev matrix
% D2eta      [-]    second-derivative Chebyshev matrix

%% Outputs
% Name size unit   explanation
% ynd       [-]    physical (non-dimensional) grid
% D1y       [-]    first-derivative matrix on the physical grid
% D2y       [-]    second-derivative matrix on the physical grid

% Transform Chebyshev to physical grid
ynd = (y_i*y_max*(1+eta)./(y_max - eta*(y_max - 2*y_i)))';

if nargout > 1
    % Chain-rule scale factors for the coordinate transform
    detady   = (y_max - eta*(y_max - 2*y_i)).^2/(2*y_i*y_max*(y_max - y_i));
    d2etady2 = -2*detady.^2*(y_max - 2*y_i)./(y_max - eta*(y_max - 2*y_i));

    D1y = diag(detady)*D1eta;    
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end

end
