function [y,u,v,w] = FalknerSkanCookeReal(eta_inf,n,beta,nu,x,Ue,We,halfsine,tolerance)
% Produces the Falkner-Skan-Cooke solution in real coordinates.
% - eta_inf: input to Falkner-Skan-Cooke function
% - n: input to Falkner-Skan-Cooke function
% - beta: Hartree parameter
% - nu: kinematic viscosity
% - x: station being evaluated
% - Ue: local streamwise reference velocity
% - We: spanwise reference velocity
% - halfsine: input to Falkner-Skan-Cooke function
%
% Using transform:
%  eta = y*sqrt((m+1)/2 * Ue/(nu*x))
%
% inversed relations:
%  y = eta*sqrt(nu*x/Ue * 2/(m+1))
%  u = f' * Ue
%  v = -sqrt(nu*Ue/(2*x*(m+1))) * ((m+1)*f+(m-1)*f'*eta)
%  w = g  * We


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
w = sol.y(4,:) * We;

% Gradient guess using manually derived relations (by JYB)
ddf = 1.1*(beta+0.2).^0.53;
dg = (beta+0.2).^0.05-0.45;

end
