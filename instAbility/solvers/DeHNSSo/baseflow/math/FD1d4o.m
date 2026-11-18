function DD = FD1d4o(D, d)
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
% Fourth-order first derivative on a uniform grid.
%
%   DD = FD1d4o(D, d)
%
%   Differentiates D along its rows (uniform spacing d), taking the coordinate
%   to ascend as the column index decreases. Central 5-point differences in the
%   interior, fourth-order one-sided differences at the two end columns.
%

%% Inputs
% Name size unit   explanation
% D         [-]    data, differentiated along its rows
% d         [-]    grid spacing

%% Outputs
% Name size unit   explanation
% DD        [-]    first derivative, same size as D

if size(D, 2) < 6
    error('FD1d4o: need at least 6 columns.')
end

DD = 0*D;

% Central differences (interior columns)
DD(:, 3:end-2) = ( 1/12*D(:, 5:end)   - 2/3*D(:, 4:end-1) ...
                 + 2/3*D(:, 2:end-3)  - 1/12*D(:, 1:end-4) ) / d;

% Fourth-order one-sided differences at the two end columns
DD(:, end-1:end) = (-25/12*D(:, end-1:end) + 4*D(:, end-2:end-1) ...
                   -    3*D(:, end-3:end-2) + 4/3*D(:, end-4:end-3) ...
                   -  1/4*D(:, end-5:end-4)) / d;
DD(:, 1:2)       = (-25/12*D(:, 1:2)       + 4*D(:, 2:3) ...
                   -    3*D(:, 3:4)         + 4/3*D(:, 4:5) ...
                   -  1/4*D(:, 5:6)) / (-d);
end
