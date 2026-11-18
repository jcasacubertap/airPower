function [mu, dmuT, d2muT2] = sutherland(T, Sstar)
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
% Sutherland's viscosity law (non-dimensional).
%
%   [mu, dmuT, d2muT2] = sutherland(T, Sstar)
%
%   Non-dimensional form: mu = T^(3/2) (1 + S*) / (T + S*)
%   where S* = S / T_inf  (S = 110.4 K for air).
%

%% Inputs
% Name  size unit   explanation
% T          [-]    non-dimensional temperature T/T_inf  (array)
% Sstar      [-]    S/T_inf  (scalar, e.g. 110.4/288 for air at T_inf = 288 K)

%% Outputs
% Name   size unit   explanation
% mu          [-]    dynamic viscosity mu(T)
% dmuT        [-]    first derivative d(mu)/dT
% d2muT2      [-]    second derivative d^2(mu)/dT^2

%% Viscosity
mu = T.^(3/2) .* (1 + Sstar) ./ (T + Sstar);

%% First derivative
if nargout >= 2
    % dmu/dT = mu/T * (3/2 - T/(T + S*))
    dmuT = mu ./ T .* (3/2 - T ./ (T + Sstar));
end

%% Second derivative
if nargout >= 3
    % obtained by differentiating dmu/dT
    TpS = T + Sstar;
    d2muT2 = mu ./ T.^2 .* (3/4 + Sstar * (3*Sstar - T) ./ (2 * TpS.^2));
end

end
