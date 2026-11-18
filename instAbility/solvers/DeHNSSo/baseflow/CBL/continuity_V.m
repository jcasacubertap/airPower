function v = continuity_V(y_phys, x_phys, u, rho)
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
% Wall-normal base-flow velocity from compressible continuity.
%   v = continuity_V(y_phys, x_phys, u, rho) returns the wall-normal velocity
%   V(y,x) enforcing  d(rho u)/dx + d(rho v)/dy = 0  (2D or infinite-swept,
%   d/dz = 0), i.e.
%
%       rho*v(y) = - int_wall^y  d(rho u)/dx  dy'     [ v(wall) = 0 ].
%
%   This is the EXACT compressible wall-normal velocity for any base flow and
%   replaces the low-Mach Lees-Dorodnitsyn closed form
%       v = sqrt(nu_e Ue/2x) (eta_LD f' - f)/rho ,
%   which is ~correct only as Ma->0 and is ~6x too small at Ma=4.5 (the strong
%   temperature/density variation is not captured). The old form left the base
%   flow violating discrete continuity by ~10x and made the HNS/PSE march
%   over-predict M4.5 2nd-mode growth by ~17%.
%
%   Inputs
%     y_phys  (ny x 1) physical wall-normal grid (any ordering)
%     x_phys  (1 x nx) streamwise stations (physical)
%     u       (ny x nx) streamwise velocity
%     rho     (ny x nx) density (same normalisation used consistently)
%   Output
%     v       (ny x nx) wall-normal velocity, same units/order as u
%
%   d(rho u)/dx is a centred FD across stations (gradient); the wall-normal
%   integral is cumulative-trapezoid from the wall up.

%% Inputs
% Name   size    unit   explanation
% y_phys (ny,nx) [-]    physical wall-normal coordinate
% x_phys (ny,nx) [-]    physical streamwise coordinate
% u      (ny,nx) [-]    streamwise base-flow velocity
% rho    (ny,nx) [-]    base-flow density

%% Outputs
% Name size    unit   explanation
% v    (ny,nx) [-]    wall-normal velocity from compressible continuity, with v = 0 at the wall

y = y_phys(:);
[ny, nx] = size(u);
rhou = rho .* u;
drhoudx = zeros(ny, nx);
for j = 1:ny
    drhoudx(j,:) = gradient(rhou(j,:), x_phys(:).');
end
[ys, is] = sort(y, 'ascend');          % integrate from wall upward
v = zeros(ny, nx);
for i = 1:nx
    rhov = -cumtrapz(ys, drhoudx(is,i));
    v(is,i) = rhov ./ rho(is,i);
end
end
