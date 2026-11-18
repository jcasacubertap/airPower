function [ DD ] = FD1d4o( D,d )
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
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
%   Differentiates the field D along its rows on a grid with uniform spacing
%   d, assuming the coordinate ascends as the column index decreases. Central
%   differences (finite differences, FD) are used in the interior; biased
%   five-point stencils handle the two boundary columns on each side. At
%   least six columns are needed for these stencils to fit.
%

%% Inputs
% Name size unit   explanation
% D         [-]    field sampled on a uniform grid (derivative taken along rows)
% d         [-]    grid spacing

%% Outputs
% Name size unit   explanation
% DD        [-]    first derivative, same size as D

%% Check enough points for the stencil
if length(D(1,:)) < 6
    error('--- Function length too small; less than 6')
end

%% Differentiate
DD = 0*D;

% Central differences on the interior columns
DD(:,  (1+2):(end-2)) =  (1/12 * D(:,  (1+4):(end-0))...
                      -    2/3 * D(:,  (1+3):(end-1))...
                      +      0 * D(:,  (1+2):(end-2))...
                      +    2/3 * D(:,  (1+1):(end-3))...
                      -   1/12 * D(:,  (1+0):(end-4)))/ d;

% Forward differences on the two smallest-coordinate columns
DD(:,(end-1):(end-0)) =(-25/12 * D(:,(end-1):(end-0))...
                      +      4 * D(:,(end-2):(end-1))...
                      -      3 * D(:,(end-3):(end-2))...
                      +    4/3 * D(:,(end-4):(end-3))...
                      -    1/4 * D(:,(end-5):(end-4)))/ d;

% Backward differences on the two largest-coordinate columns
DD(:,  (1+0):  (1+1)) =(-25/12 * D(:,  (1+0):  (1+1))...
                      +      4 * D(:,  (1+1):  (1+2))...
                      -      3 * D(:,  (1+2):  (1+3))...
                      +    4/3 * D(:,  (1+3):  (1+4))...
                      -    1/4 * D(:,  (1+4):  (1+5)))/(-d);

end
