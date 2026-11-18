function [E_tot, E_loc] = mode_energy(StabGrid, Bound, Mode, StabRes, x_loc)
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
% Mode-summed perturbation kinetic energy.
%
%   [E_tot, E_loc] = MODE_ENERGY(StabGrid, Bound, Mode, StabRes, x_loc)
%
%   Evaluates the measure DeHNSSo-v1 uses for grid convergence (Eqs. 24 and 25
%   there), summed over the active modes:
%
%       E_tot = int_{x0}^{xb} int_0^H 1/2 (|u'|^2+|v'|^2+|w'|^2) dy dx
%       E_loc = int_0^H       1/2 (|u'|^2+|v'|^2+|w'|^2) dy  at x = x_loc
%
%   The streamwise integral stops at the buffer onset, so the attenuated region
%   never enters the norm. x_loc is optional: supply it to evaluate the local
%   measure at that station, omit it and E_loc is returned as NaN. Sec. 5.1 and
%   DeHNSSo-v1 both evaluate the local measure at, so the two sections report
%   the same quantity at the same place.
%
%   Unlike the peak amplitude max_y |u'|, this weighs the whole field rather
%   than one wall-normal maximum of one variable, and it weighs the modes by
%   their energy rather than reporting whichever mode differs most.

%% Inputs
% Name     size         unit   explanation
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Bound    (1,1) struct [-]    boundary conditions and outflow-buffer coefficients
% Mode     (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% StabRes  (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% x_loc    (1,1)        [-]    optional station for the local measure; omit it and E_loc is NaN

%% Outputs
% Name  size  unit   explanation
% E_tot (1,1) [-]    mode-summed perturbation kinetic energy over the buffer-free domain
% E_loc (1,1) [-]    mode-summed perturbation kinetic energy at x_loc, interpolated to that station; NaN when x_loc is omitted

ny = StabGrid.ny;  nx = StabGrid.nx;
y  = StabGrid.y;  if ~isvector(y); y = y(:,1); end;  y = y(:);
x  = StabGrid.x;  if ~isvector(x); x = x(1,:); end;  x = x(:).';

% The wall-normal grid is stored from the top boundary down to the wall, so it
% has to be flipped before integrating: trapz over a descending coordinate
% returns the integral with the wrong sign.
iy = 1:ny;
if y(end) < y(1); iy = ny:-1:1; y = y(iy); end

ib = min(Bound.ibuf, nx);            % buffer onset (ibuf = nx+1 when off)
xr = 1:ib;

phys = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));
e_xy = zeros(ny, numel(xr));
for j = phys(:).'
    u = reshape(StabRes.u(j,:,:), ny, nx);
    v = reshape(StabRes.v(j,:,:), ny, nx);
    w = reshape(StabRes.w(j,:,:), ny, nx);
    e_xy = e_xy + 0.5*(abs(u(iy,xr)).^2 + abs(v(iy,xr)).^2 + abs(w(iy,xr)).^2);
end

e_x   = trapz(y, e_xy, 1);           % wall-normal integral, one value per station
E_tot = trapz(x(xr), e_x);

% E_loc is optional: no station, no local measure. A fixed default cannot suit
% every domain, and silently moving it would change what E_loc means.
if nargin < 5 || isempty(x_loc)
    E_loc = NaN;
else
    assert(x_loc >= x(1) && x_loc <= x(ib), ...
        'x_loc = %g is outside the buffer-free range [%g, %g]', x_loc, x(1), x(ib));
    % Interpolated, not snapped: the nearest node sits at a different fraction of
    % dx on every grid, which would add a sampling error to the discretization one.
    E_loc = interp1(x(xr), e_x, x_loc, 'spline');
end
end
