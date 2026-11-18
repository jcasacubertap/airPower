function [prof, alpha0, out] = solver_IPSE(omega, beta, StabGrid, Opt)
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
% Incompressible linear PSE (parabolized stability equations) inflow condition (de Groot & Kotsonis).
%
%   [prof, alpha0, out] = solver_IPSE(omega, beta, StabGrid, Opt)
%
%   Incompressible analogue of solver_CPSE, with state (u,v,w,p). A self-
%   similar precursor first stretches the inflow plane onto [Rstart, R0];
%   a local ILST (incompressible LST, Orr-Sommerfeld/Squire) eigenmode
%   seeds an incompressible PSE march (local function solver_PSE, below)
%   up to R0, normalized by max|u|.
%

%% Inputs
% Name     size unit   explanation
% omega         [-]    angular frequency of the mode
% beta          [-]    spanwise wavenumber of the mode
% StabGrid      [-]    grid/base-flow struct (U, W, y, dzdy, Re, ny; H, y_i opt.)
% Opt           [-]    options (pse_Rstart, pse_npre, pse_ymax, pse_yi)

%% Outputs
% Name   size unit   explanation
% prof        [-]    inflow shape at R0, normalised by max|u| (u,v,w,p)
% alpha0      [-]    complex streamwise wavenumber at R0
% out         [-]    march history (Rv, alpha)

%% Options and defaults
if nargin < 4; Opt = struct(); end
R0 = StabGrid.Re;  ny = StabGrid.ny;  y = StabGrid.y(:);
Rstart = getopt(Opt,'pse_Rstart', max(0.6*R0, R0-200));
Np     = getopt(Opt,'pse_npre',   60);
ymax   = getopt(Opt,'pse_ymax',   getf(StabGrid,'H',   max(y)));
yi     = getopt(Opt,'pse_yi',     getf(StabGrid,'y_i', ymax/20));

%% Self-similar precursor: stretch the inflow profile onto [Rstart, R0]
U0 = StabGrid.U(:,1);  W0 = StabGrid.W(:,1);
Rv = linspace(Rstart, R0, Np);
xv = Rv.^2 / R0;
[~,kw] = min(U0);  ed = 1;  if kw == 1; ed = ny; end
Ur = zeros(ny,Np);  Wr = Ur;
for i = 1:Np
    ys = y*(R0/Rv(i));
    Ur(:,i) = interp1(y,U0,ys,'pchip',U0(ed));
    Wr(:,i) = interp1(y,W0,ys,'pchip',W0(ed));
end
Vr = contV(y, xv, Ur);

%% Base-flow derivatives on the precursor grid
[~,DM] = chebdif(ny,2);
D1y = diag(StabGrid.dzdy(:))*DM(:,:,1);
dyUr = D1y*Ur;  dyVr = D1y*Vr;  dyWr = D1y*Wr;
dxUr = ddx(Ur,xv);  dxVr = ddx(Vr,xv);  dxWr = ddx(Wr,xv);

%% March, then map the developed shape at R0 onto the StabGrid (unit max|u|)
[alpha,u,v,w,p,ye] = solver_PSE(R0, Ur,Vr,Wr, dxUr,dxVr,dxWr, dyUr,dyVr,dyWr, ...
                                xv, y, omega, beta, ny, ymax, yi);
out    = struct('Rv',Rv, 'alpha',alpha);
alpha0 = alpha(end);

uS = interp1(ye(:), u(:,end), y, 'pchip', 0);
vS = interp1(ye(:), v(:,end), y, 'pchip', 0);
wS = interp1(ye(:), w(:,end), y, 'pchip', 0);
pS = interp1(ye(:), p(:,end), y, 'pchip', 0);
un = max(abs(uS));  if un == 0; un = 1; end
prof = struct('u',(uS/un).', 'v',(vS/un).', 'w',(wS/un).', 'p',(pS/un).');

end

%% Local helpers

function V = contV(y, x, U)
%CONTV  Wall-normal velocity from continuity: v = -int_wall^y dU/dx dy', v(wall) = 0.
[ny,nx] = size(U);
dUdx = zeros(ny,nx);
for k = 1:ny; dUdx(k,:) = gradient(U(k,:), x); end

[ys,is] = sort(y(:));
V = zeros(ny,nx);
for i = 1:nx
    I = zeros(ny,1);
    I(is) = cumtrapz(ys, dUdx(is,i));
    V(:,i) = -I;
end
end

function dF = ddx(F, x)
%DDX  Streamwise gradient of each wall-normal row of F w.r.t. x.
dF = zeros(size(F));
for k = 1:size(F,1); dF(k,:) = gradient(F(k,:), x); end
end

function v = getf(S, f, dv)
%GETF  Return S.f if present and non-empty, else default dv.
if isfield(S,f) && ~isempty(S.(f)); v = S.(f); else; v = dv; end
end

function v = getopt(O, f, dv)
%GETOPT  Return O.f if present and non-empty, else default dv.
if isfield(O,f) && ~isempty(O.(f)); v = O.(f); else; v = dv; end
end

%% ------------------------------------------------------------------------
%  PSE marcher (de Groot & Kotsonis): seeded by solver_ILST, marches the
%  linear PSE from the precursor inlet to R0. Local to solver_IPSE.
function [alpha,u,v,w,p,y] = solver_PSE(Re,Ur,Vr,Wr,dxUr,dxVr,dxWr,dyUr,dyVr,dyWr,x_md,y_md,omega,beta,N,y_max,y_i)
%SOLVER_PSE  March the linear PSE from the inflow station to the outflow.
%
%   [alpha, u, v, w, p, y] = solver_PSE(Re, Ur, Vr, Wr, dxUr, dxVr, dxWr, dyUr, dyVr, dyWr, x_md, y_md, omega, beta, N, y_max, y_i)
%
%   Executes a PSE march from min(x_md) to x_md(:,end) over the mean flow
%   (Ur, Vr, Wr) given on a uniform grid, with the fields ordered from the
%   freestream to the wall as the row index increases and downstream stations
%   as the column index increases. The inflow is seeded by an ILST
%   (Orr-Sommerfeld) eigenmode, and each station is iterated on alpha through
%   the auxiliary normalisation condition.
%
%   Inputs:
%     Re                Reynolds number
%     Ur,Vr,Wr          base-flow velocity components (ny x nx)
%     dxUr,dxVr,dxWr    their streamwise derivatives
%     dyUr,dyVr,dyWr    their wall-normal derivatives
%     x_md,y_md         base-flow streamwise / wall-normal grids
%     omega             angular frequency
%     beta              spanwise wavenumber
%     N                 number of Chebyshev nodes
%     y_max             outer domain height
%     y_i               node-clustering height
%   Outputs:
%     alpha             complex streamwise wavenumber alpha(x), based on max|u|
%     u,v,w,p           complex shape functions, normalised by local max|u|
%     y                 physical wall-normal grid
%
%   Written by Koen de Groot and Marios Kotsonis

[~,nx] = size(Ur);

%% Spectral basis functions

% Set up Chebyshev grid and associated Pseudo-Spectral Differentiation
% Matrices (PSDMs) in the computational domain: i.e. eta in [-1,1]
[eta,DM] = chebdif(N,2);

% Transform the PSDMs to the physical domain: eta |-> y
[y,D1,D2] = MappingMalik(y_max,y_i,eta',DM(:,:,1),DM(:,:,2));
clear DM eta

% The PSDMs generally have size (N-2)x(N-2); define the matching identity
% and zero column/row accordingly
I      = eye(N-2);
Zc     = zeros(N-2,1);
Zr     = zeros(1,N-2);

%% Initial Condition Generation with ILST (expanded Orr-Sommerfeld)

% Execute ILST analysis for initialization
tic
fprintf(' Running ILST for inflow conditions...\n')
[EVOS,uOS,vOS,wOS,pOS,yOS] = solver_ILST(...
    Re,Ur(:,1),Wr(:,1),dyUr(:,1),dyWr(:,1),y_md,...
    omega,beta,N,y_max,y_i);
fprintf('\b\b\b\b, finished in %.2f seconds.\n',toc)

% Find boundary layer thickness for filtering
i_edge = find(Ur(:,1)<0.999*max(Ur(:,1)),1); % index where U starts decreasing
d99 = interp1(Ur(i_edge:end,1),y_md(i_edge:end),0.99*max(Ur(i_edge:end,1)),'spline');

% Filter eigenfunctions for correct eigenvalue
[OSeigval,index] = EVfilter(EVOS,vOS,yOS,omega,beta,d99);

% Initialize solution vector phi = [u v w p]'(x,y) and alpha = alpha(x)
phi    = zeros(4*(N-2)+2,nx);
alpha  = zeros(1,nx);

% Set ILST eigenfunction in first phi station
phi(:,1) = ([uOS(2:end-1,index);
            vOS(2:end-1,index);
            wOS(2:end-1,index);
            pOS(:,index)]);

% Set ILST eigenvalue in first alpha station
alpha(1) = OSeigval;

%% Execute marching algorithm
iu=sqrt(-1);
stabon = false;
percentage = round(100*2/nx);
fprintf(' Starting marching method.')
fprintf(' Progress:  %.0f percent',0)
for i=2:nx    
    % Interpolate base flow components and create diagonal matrix form
      U = diag(interp1(y_md,  Ur(:,i),y(2:N-1),'spline'));
    dxU = diag(interp1(y_md,dxUr(:,i),y(2:N-1),'spline'));
    dyU = diag(interp1(y_md,dyUr(:,i),y(2:N-1),'spline'));
    
      V = diag(interp1(y_md,  Vr(:,i),y(2:N-1),'spline'));
    dxV = diag(interp1(y_md,dxVr(:,i),y(2:N-1),'spline')); % note that this term is order 1/Re^2
    dyV = diag(interp1(y_md,dyVr(:,i),y(2:N-1),'spline'));
    
      W = diag(interp1(y_md,  Wr(:,i),y(2:N-1),'spline'));
    dxW = diag(interp1(y_md,dxWr(:,i),y(2:N-1),'spline'));
    dyW = diag(interp1(y_md,dyWr(:,i),y(2:N-1),'spline'));
    
    % Initialize alpha at next station
    alpha(i) = alpha(i-1);
    % Reset convergence measure
    dal = 1;
    while abs(dal)*Re >= 1e-7
        
        % Compile common convection-diffusion terms
        Del = -iu*omega*I + iu*alpha(i)*U + iu*beta*W + I*(alpha(i)^2 + beta^2)/Re+V*D1(2:end-1,2:end-1)-D2(2:end-1,2:end-1)/Re;
        
        % Define A in A*q + B*(q-qp)/dx + C*q*dalpha/dx = 0
        A = [Del+dxU               dyU          0*I [Zc iu*alpha(i)*I Zc]
                  Zr       -D2(1,  2:end-1)/Re   Zr       D1(1,  :)
                 dxV           Del+dyV          0*I     D1(2:end-1,:)
                  Zr       -D2(end,2:end-1)/Re   Zr       D1(end,:)
                 dxW               dyW          Del   [Zc iu*beta*I Zc]
             I*alpha(i)*iu D1(2:end-1,2:end-1) I*beta*iu [Zc 0*I Zc]];
                 
        % Define B in A*q + B*(q-qp)/dx + C*q*dalpha/dx = 0
        B = [U-I*(2*alpha(i)*iu)/Re    I*0  I*0  [Zc  I  Zc]
                  Zr                    Zr   Zr  [0   Zr  0]
                 I*0 U-I*(2*alpha(i)*iu)/Re I*0  [Zc I*0 Zc] 
                  Zr  Zr                     Zr  [0   Zr  0]
                 I*0 I*0  U-I*(2*alpha(i)*iu)/Re [Zc I*0 Zc]
                  I  I*0  I*0                    [Zc I*0 Zc]];
              
        % Define C in A*q + B*(q-qp)/dx + C*q*dalpha/dx = 0
        C = [-I*(iu/Re)                I*0  I*0  [Zc I*0 Zc] 
                  Zr                    Zr   Zr  [0   Zr  0]
                 I*0 -I*(iu/Re)             I*0  [Zc I*0 Zc] 
                  Zr  Zr                     Zr  [0   Zr  0]
                 I*0 I*0 -I*(iu/Re)              [Zc I*0 Zc] 
                 I*0 I*0  I*0                    [Zc I*0 Zc]];
        
        % Determine dalpha/dx
	    dx  =   x_md(1,i) - x_md(1,i-1);
        dax = (alpha(i) - alpha(i-1))/dx;
        
        % Stabilization amount s: nonzero only when the step is too small
        % relative to the streamwise wavelength
        s = (1/abs(real(alpha(i))) - dx)*(real(alpha(i)) ~= 0 && dx <= 1.01/abs(real(alpha(i))));
        if s && ~stabon
            stabon = i; 
        end 
        
        % Apply stabilization if required
        if s B = B + s*(A + C*dax); end
        
        % Form LHS matrix and RHS vector
        Ml  = A + B/dx + C*dax;
        Mr  =     B/dx*phi(:,i-1);

        % Solve problem
        phi(:,i) = Ml\Mr;

        % Use the whole eigenfunction for the auxiliary condition; this moves
        % the residual growth of the shape function into the exponential
        % growth rate, which is optimal
        uold = phi(:,i-1);
        unew = phi(:,i);
        yint = [repmat(y(2:end-1),3,1);y];

        % Enforce auxiliary condition to maximize growth contribution in
        % the exponential term in the ansatz partition
        dal  = -iu*trapz(yint,conj(unew).*(unew-uold)/dx)/trapz(yint,abs(unew).^2);
        
        % Update alpha
        alpha(i)=alpha(i)+dal;

    % Notify progress
    if round(100*i/nx) == percentage
        for j = 1:ceil(log10(percentage))+9; fprintf('\b'); end
        fprintf('%.0f percent ',percentage)
        percentage = percentage + 1;
    end
    
    end
end

% Augment Dirichlet conditions to perturbation velocity components
phi = [zeros(1,nx); phi(0*(N-2)+1:1*(N-2)  ,:); zeros(1,nx);
       zeros(1,nx); phi(1*(N-2)+1:2*(N-2)  ,:); zeros(1,nx);
       zeros(1,nx); phi(2*(N-2)+1:3*(N-2)  ,:); zeros(1,nx);
                      phi(3*(N-2)+1:4*(N-2)+2,:)];
                 
% Extract the individual eigenvectors
u   = phi(0*N+1:1*N,:);
v   = phi(1*N+1:2*N,:);
w   = phi(2*N+1:3*N,:);
p   = phi(3*N+1:4*N,:);
         
%% Normalise shape functions and correct the growth rate

% Calculate the streamwise derivative of u
dux = zeros(size(u));           % march grid is non-uniform in x
for k = 1:size(u,1); dux(k,:) = gradient(u(k,:), x_md(1,:)); end

% Residual growth of the eigenfunction evaluated at umax (Bertolotti thesis,
% p.28). 
alphaphi = zeros(1,nx);

yinteg=linspace(0,max(y),4000);

for i = 1:nx

    [~,k]=max(abs(interp1(y,u(:,i),yinteg,'spline')));
    umax=interp1(y,u(:,i),yinteg(k),'spline');
    dumax=interp1(y,dux(:,i),yinteg(k),'spline');

    alphaphi(1,i) = (dumax/umax)*sqrt(-1);

    % Normalise all shape functions with the local u maximum
    u(:,i)=u(:,i)/abs(umax);
    v(:,i)=v(:,i)/abs(umax);
    w(:,i)=w(:,i)/abs(umax);
    p(:,i)=p(:,i)/abs(umax);

end

% Correct the growth rate with the residual growth rate based on umax
alpha=alpha-alphaphi;

% Report whether stabilization was activated during any of the steps
if stabon
    fprintf('\n Stabilization was activated.\n')
end

end
