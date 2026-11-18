function [sol] = FalknerSkanCooke(eta_inf,n,beta,halfsine,tol)
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
% Solve the Falkner-Skan-Cooke similarity equations (bvp5c).
%
%   sol = FalknerSkanCooke(eta_inf, n, beta, halfsine, tol)
%
%   Falkner-Skan: f''' + f f'' + beta(1 - f'^2) = 0;  Cooke: g'' + f' g = 0,
%   with eta = y sqrt((m+1)/2 Ue/(nu x)), f' = u/Ue, g = w/We.
%

%% Inputs
% Name     size unit   explanation
% eta_inf       [-]    freestream location where f' = g = 1 is imposed
% n             [-]    number of points
% beta          [-]    Hartree parameter
% halfsine      [-]    (optional) half-sine node distribution, clustered at the wall
% tol           [-]    (optional) absolute tolerance (default 1e-8)

%% Outputs
% Name size unit   explanation
% sol       [-]    bvp5c solution struct

% Half-sine distribution (bias at the wall, at low eta)
if nargin==4 && halfsine
    etas = ( 1-cos(linspace(0,pi/2,n)) ) * eta_inf;
else
    etas = linspace(0,eta_inf,n);
end
if nargin<5
    tol = 1e-8;
end
if beta>50000
    warning('Beta reduced to 50000!')
    beta = 50000;
end

% Tolerances
% Only Absolute tolerance is relevant, Relative wont work if a boundary
% condition is zero, such as with the FSC equations.
options = bvpset('AbsTol',tol,'RelTol',tol*100);  

% Initial solution guess
solinit = bvpinit(etas,@FSCinit,[]);

% Solution Falkner-Skan
warning('off','MATLAB:bvp5c:RelTolNotMet')
sol = bvp5c(@FSCfunction,@FSCbc,solinit,options);
warning('on','MATLAB:bvp5c:RelTolNotMet')

    function dydt = FSCfunction(~,y)
    % Falkner-Skan-Cooke equation in matrix form. 
    % y = [f f' f'']
    dydt = [y(2) 
            y(3) 
            -y(1)*y(3)];
    end

end

function yinit = FSCinit(t)
% Initial guess
yinit = [ t - log(t+1)
          1 - (1/(t+1))
          (t+1)^(-2)];
end

function res = FSCbc(ya,yb) 
% ya is at eta=0, yb is at eta=infinity
res = [ya(1)        % f(0) = 0
       ya(2)        % f'(0) = 0
       yb(2)-1];      % f'(inf) = 1 
end

