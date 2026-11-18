function names = bf_deriv_fields()
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
% Physical base-flow derivative fields that may be supplied on
%   the raw BF (analytic base flow) instead of being computed in compute_bf_derivatives.

%% Inputs
% None.

%% Outputs
% Name  size unit   explanation
% names cell [-]    names of the base-flow derivative fields the solver expects, dxU through ddyT

names = {'dxU','dyU','dxxU','ddyU','dxyU', 'dxV','dyV','dxxV','ddyV','dxyV', ...
         'dxW','dyW','dxxW','ddyW', 'dxrho','dyrho', 'dxT','dyT','dxxT','ddyT'};
end
