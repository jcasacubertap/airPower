function [EV,u,v,w,p,y,fn]=solver_ILST(Re,Ur,Wr,dUr,dWr,y_md,omega,beta,N,y_max,y_i)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
%     Modified (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
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
% Solve the Orr-Sommerfeld eigenvalue problem for the inflow mode.
%
%   [EV, u, v, w, p, y, fn] = solver_ILST(Re, Ur, Wr, dUr, dWr, y_md, omega, beta, N, y_max, y_i)
%
%   Solves the Orr-Sommerfeld problem on the domain y in [0, y_i, y_max],
%   clustering half the nodes between y = 0 and y_i with the same scaling as
%   y_md. Clamped conditions are imposed as Dirichlet conditions at the
%   boundary nodes. The base flow is interpolated onto the Chebyshev grid,
%   extrapolating where y_max exceeds max(y_md). The quadratic eigenvalue
%   problem is recast into companion (generalised linear) form and solved for
%   the full spectrum; the calling routine filters it for the physical mode.
%
%     Ur,Wr      streamwise / spanwise base-flow velocity profiles
%     dUr,dWr    their wall-normal derivatives
%     u,v,w,p    eigenfunctions on the N nodes stored in y

%% Inputs
% Name  size            unit   explanation
% Re                    [-]    Reynolds number
% Ur    (numel(y_md),1) [-]    streamwise base-flow velocity profile on the base-flow grid y_md
% Wr    (numel(y_md),1) [-]    spanwise base-flow velocity profile on the base-flow grid y_md
% dUr   (numel(y_md),1) [-]    wall-normal derivative of Ur
% dWr   (numel(y_md),1) [-]    wall-normal derivative of Wr
% y_md                  [-]    base-flow wall-normal grid
% omega                 [-]    angular frequency
% beta                  [-]    spanwise wavenumber
% N                     [-]    number of Chebyshev nodes
% y_max                 [-]    outer domain height
% y_i                   [-]    node-clustering height

%% Outputs
% Name size    unit   explanation
% EV           [-]    eigenvalue spectrum (2*(4*N-6) eigenvalues)
% u    (N,4*N) [-]    streamwise eigenfunctions, one column per eigenvalue
% v    (N,4*N) [-]    wall-normal eigenfunctions, one column per eigenvalue
% w    (N,4*N) [-]    spanwise eigenfunctions, one column per eigenvalue
% p    (N,4*N) [-]    pressure eigenfunctions, one column per eigenvalue
% y            [-]    physical wall-normal grid
% fn           [-]    Frobenius norm of the operator (eig error ~ eps*fn)

%% Spectral basis functions

% Set up Chebyshev grid and associated Pseudo-Spectral Differentiation
% Matrices (PSDMs) in the computational domain: i.e. eta in [-1,1]
[eta,DM] = chebdif(N,2);

% Transform the PSDMs to the physical domain: eta |-> y
[y,D1,D2] = MappingMalik(y_max,y_i,eta',DM(:,:,1),DM(:,:,2));
clear DM eta

% The transformed matrices represent the full operators. Removing the outer
% elements yields the Dirichlet conditions
D1e = D1(1,  :);
D1w = D1(end,:);
D1x = D1(2:N-1,:);
D1  = D1(2:N-1,2:N-1);
D2e = D2(1,  :);
D2w = D2(end,:);
D2  = D2(2:N-1,2:N-1);
% These matrices now have size (N-2)x(N-2); define the matching identity,
% zero matrix, and zero column/row
I   = eye(N-2);
Z   = zeros(N-2);
Zc  = zeros(N-2,1);
Zr  = zeros(1,N-2);

%% Base flow on the Chebyshev grid

% Interpolate the base flow coefficients U and W and their derivatives onto
% the (physical non-dimensional) Chebyshev grid y
[~,i_top] = max(y_md);
U  = diag( interp1(y_md, Ur,y(2:N-1) ,'spline',  Ur(i_top)) );
DU = diag( interp1(y_md,dUr,y(2:N-1) ,'spline', dUr(i_top)) );
W  = diag( interp1(y_md, Wr,y(2:N-1) ,'spline',  Wr(i_top)) );
DW = diag( interp1(y_md,dWr,y(2:N-1) ,'spline', dWr(i_top)) );
% Where y_max exceeds max(y_md), interp1 extrapolates with the same spline;
% safest is to supply the base flow to a very large max(y_md) and keep
% y_max < max(y_md).

%% Operator assembly

% Imaginary unit
iu=sqrt(-1);

% Diagonal convection-diffusion terms
CD = (beta^2*I - D2)/Re + iu*beta*W - iu*omega*I;

%% Define A in A*Xi = alpha*B*Xi + alpha^2*C*Xi

A11 = CD;
A12 = DU;

A22 = CD;
PC2u2 = -D2e(2:end-1)/Re;
PC2l2 = -D2w(2:end-1)/Re;
A24 = D1x;
PC2u4 = D1e(1,:);
PC2l4 = D1w(end,:);

A32 = DW;
A33 = CD;
A34 = iu*beta*I;
clear CD

A42 = D1;
A43 = iu*beta*I;

%% Define B in A*Xi = alpha*B*Xi + alpha^2*C*Xi

B11 = -iu*U;
B14 = -iu*I;

B22 = -iu*U;

B33 = -iu*U;

B41 = -iu*I;

%% Define C in A*Xi = alpha*B*Xi + alpha^2*C*Xi

C11 = -I/Re;

C22 = -I/Re;

C33 = -I/Re;

%% Generate matrices

% LHS matrix
%     u    v    w   pe  p  pw
A = [A11  A12   Z  [Zc  Z  Zc]   % x-momentum
      Zr PC2u2  Zr    PC2u4      % y-compatibility (edge boundary)
      Z   A22   Z      A24       % y-momentum
      Zr PC2l2  Zr    PC2l4      % y-compatibility (wall boundary)
      Z   A32  A33 [Zc A34 Zc]   % z-momentum
      Z   A42  A43 [Zc  Z  Zc]]; % continuity

% RHS matrix (alpha)
B = [B11   Z    Z  [Zc B14 Zc]
      Zr   Zr   Zr [0   Zr  0]
      Z   B22   Z  [Zc  Z  Zc]
      Zr   Zr   Zr [0   Zr  0]
      Z    Z   B33 [Zc  Z  Zc]
     B41   Z    Z  [Zc  Z  Zc] ];

% RHS matrix (alpha^2)
C = [C11   Z    Z  [Zc  Z  Zc]
      Zr   Zr   Zr [0   Zr  0]
      Z   C22   Z  [Zc  Z  Zc]
      Zr   Zr   Zr [0   Zr  0]
      Z    Z   C33 [Zc  Z  Zc]
      Z    Z    Z  [Zc  Z  Zc] ];

%% Construct companion matrices

% A-sized zero and identity matrices for the companion linearisation
Z = zeros(size(A));
I =   eye(size(A));

% Build the companion matrices (quadratic to generalised linear problem)
A = [A -B; Z I];
B = [Z  C; I Z];

% Frobenius norm of A; error of eig = machine precision * ||A||_F
fn = norm(A,'fro');

%% Solve the eigenvalue problem

[V,D] = eig(A,B);

% Extract the eigenvalues
EV = diag(D);

% Extract the unique eigenvector part
V = V(1:end/2,:);

% Augment the Dirichlet boundary conditions to the functions
V = [zeros(1,2*(4*(N-2)+2)); V(0*(N-2)+1:1*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
     zeros(1,2*(4*(N-2)+2)); V(1*(N-2)+1:2*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
     zeros(1,2*(4*(N-2)+2)); V(2*(N-2)+1:3*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
                             V(3*(N-2)+1:4*(N-2)+2,:)];

% Extract the eigenvectors
u   = V(0*N+1:1*N,:);
v   = V(1*N+1:2*N,:);
w   = V(2*N+1:3*N,:);
p   = V(3*N+1:4*N,:);

end
