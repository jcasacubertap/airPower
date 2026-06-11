function [sol] = FalknerSkanCooke(eta_inf,n,beta,halfsine,tol)
% Solve the Falkner-Skan-Cooke equations. 
% - eta_inf is where the boundary condition f'=g=1 is imposed. 
% - n is the number of points
% - beta is the Hartree parameter
% - halfsine toggles a half-sine distribution of points along eta
%
% Definitions: 
%   Falkner-Skan: f''' + f*f'' + beta(1-f'^2) = 0
%   Cooke: g'' + f'g = 0
% Transformations: 
%   eta = y*sqrt((m+1)/2 * Ue/(nu*x))
%   f' = u/Ue
%   g = w/We

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
% 
% The FSC equation is nested so that it shares the workspace of the main
% function, which includes the variable beta.

    function dydt = FSCfunction(~,y)
    % Falkner-Skan-Cooke equation in matrix form. 
    % y = [f f' f'' g g']
    dydt = [y(2) 
            y(3) 
            -y(1)*y(3) - beta*(1-y(2)^2)
            y(5)
            -y(1)*y(5)];
    end

end

function yinit = FSCinit(t)
% Initial guess
yinit = [ t - log(t+1)
          1 - (1/(t+1))
          (t+1)^(-2)
          1 - (1/(t+1))
          (t+1)^(-2)];
end

function res = FSCbc(ya,yb) 
% ya is at eta=0, yb is at eta=infinity
res = [ya(1)        % f(0) = 0
       ya(2)        % f'(0) = 0
       ya(4)        % g(0) = 0
       yb(2)-1      % f'(inf) = 1
       yb(4)-1];    % g(inf) = 1
end
