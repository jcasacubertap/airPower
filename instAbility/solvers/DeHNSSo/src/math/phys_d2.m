function [dxx,dyy,dxy] = phys_d2(q, qxi, qeta, D_xi2, D1, D2, m)
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
% Curvilinear second physical derivatives of a field q(ny,nx).
%
%   [dxx, dyy, dxy] = phys_d2(q, qxi, qeta, D_xi2, D1, D2, m)
%
%   Applies the full curvilinear chain rule (first AND second inverse metrics),
%   matching gridgen/compute_bf_derivatives:
%     d_xx q = xix^2 q_xixi + 2 xix etax q_xieta + etax^2 q_etaeta + xixx q_xi + etaxx q_eta
%     d_yy q = xiy^2 q_xixi + 2 xiy etay q_xieta + etay^2 q_etaeta + xiyy q_xi + etayy q_eta
%     d_xy q = xix xiy q_xixi + (xix etay + xiy etax) q_xieta + etax etay q_etaeta + xixy q_xi + etaxy q_eta
%
%     D1, D2    wall-normal first / second derivative operators
%     dxx, dyy, dxy   second physical derivatives (dxy only if requested)

%% Inputs
% Name  size           unit   explanation
% q     (ny,nx)        [-]    field to differentiate
% qxi   (ny,nx)        [-]    first derivative of q with respect to xi
% qeta  (ny,nx)        [-]    first derivative of q with respect to eta
% D_xi2 (nx,nx) sparse [-]    second-derivative matrix in xi
% D1    (ny,ny)        [-]    first-derivative matrix in eta
% D2    (ny,ny)        [-]    second-derivative matrix in eta
% m     (1,1) struct   [-]    first and second inverse metric terms of the mapping

%% Outputs
% Name size    unit   explanation
% dxx  (ny,nx) [-]    second physical derivative of q with respect to x
% dyy  (ny,nx) [-]    second physical derivative of q with respect to y
% dxy  (ny,nx) [-]    mixed second physical derivative of q with respect to x and y

% Computational second derivatives
qxixi   = q * D_xi2.';        % d^2/dxi^2  (4th-order streamwise operator)
qetaeta = D2 * q;            % d^2/deta^2 (wall-normal spectral/FD operator)
qxieta  = D1 * qxi;         % d/deta of d/dxi

% Assemble physical derivatives via the curvilinear chain rule
dxx = m.xix.^2.*qxixi + 2*m.xix.*m.etax.*qxieta + m.etax.^2.*qetaeta + m.xixx.*qxi + m.etaxx.*qeta;
dyy = m.xiy.^2.*qxixi + 2*m.xiy.*m.etay.*qxieta + m.etay.^2.*qetaeta + m.xiyy.*qxi + m.etayy.*qeta;
if nargout > 2
    dxy = m.xix.*m.xiy.*qxixi + (m.xix.*m.etay + m.xiy.*m.etax).*qxieta ...
        + m.etax.*m.etay.*qetaeta + m.xixy.*qxi + m.etaxy.*qeta;
end
end
