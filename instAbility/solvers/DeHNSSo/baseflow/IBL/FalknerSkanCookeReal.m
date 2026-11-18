function [y,u,v] = FalknerSkanCookeReal(eta_inf,n,beta,nu,x,Ue,halfsine,tolerance)
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
% Falkner-Skan-Cooke solution in physical coordinates.
%
%   [y,u,v] = FalknerSkanCookeReal(eta_inf, n, beta, nu, x, Ue, halfsine, tolerance)
%
%   Solves the similarity problem (FalknerSkanCooke) and maps it to physical
%   coordinates at station x, with m = beta/(2-beta):
%     y = eta sqrt(2 nu x / ((m+1) Ue)),   u = f' Ue,
%     v = -sqrt(nu Ue/(2 x (m+1))) ((m+1) f + (m-1) f' eta).
%
%     eta_inf, n, halfsine, tolerance  — passed to FalknerSkanCooke
%     y, u, v  — physical wall-normal grid and velocity components

%% Inputs
% Name      size          unit   explanation
% eta_inf   (1,1)         [-]    similarity coordinate at the freestream edge, passed to FalknerSkanCooke
% n         (1,1)         [-]    number of similarity nodes, passed to FalknerSkanCooke
% beta                    [-]    Hartree parameter,  nu — kinematic viscosity
% nu        (1,1)         [m^2/s] kinematic viscosity
% x                       [-]    station,            Ue — local streamwise edge velocity
% Ue        (1,1)         [m/s]  edge velocity at station x
% halfsine  (1,1) logical [-]    use the half-sine node distribution, passed to FalknerSkanCooke
% tolerance (1,1)         [-]    bvp5c solver tolerance, passed to FalknerSkanCooke

%% Outputs
% Name size  unit   explanation
% y    (n,1) [m]    physical wall-normal grid at station x
% u    (n,1) [m/s]  streamwise velocity profile
% v    (n,1) [m/s]  wall-normal velocity profile

if nargin < 8
    halfsine = false;
end
if nargin < 9
    tolerance = 1e-8;
end
if beta==2
    error('Beta cant be equal to 2!')
end

% Get similar solution
[sol] = FalknerSkanCooke(eta_inf,n,beta,halfsine,tolerance);

% Transform Hartree to Falkner-Skan parameter
m = beta/(2-beta); % beta = 2*m/(1+m)

% Scale to real coordinates
y = sol.x * sqrt(nu*x/Ue * 2/(m+1));
u = sol.y(2,:) * Ue;
v = -sqrt(nu*Ue/(2*x*(m+1))) * ((m+1)*sol.y(1,:)+(m-1)*sol.y(2,:).*sol.x);
end
