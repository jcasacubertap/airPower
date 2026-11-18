function [Eval, Evec] = solver_CLST(StabGrid, i, omega, beta, bc_top, drop_alpha2)
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
% Assemble and solve the compressible LST (linear stability theory) eigenvalue problem.
%
%   [Eval, Evec] = solver_CLST(StabGrid, i, omega, beta, bc_top, drop_alpha2)
%
%   Builds the spatial compressible stability operator at station i in
%   variables (rho, u, v, w, T) and solves the quadratic eigenvalue problem
%   (AA - alpha*BB - alpha^2*CCxx) q = 0 via polyeig. Matrix assembly is
%   delegated to CLST_Matrix (one call per grid point), using the ideal-gas
%   equation of state (EOS) and Sutherland's viscosity law. The caller
%   filters Eval to find the physical mode (see EVfilter).
%
%                   Re, Ma, Pr, gamma, Sstar)
%     omega, beta   angular frequency, spanwise wavenumber

%% Inputs
% Name        size  unit   explanation
% StabGrid          [-]    grid structure (U, W, rho, T and y-derivatives, D-metrics,
% i                 [-]    streamwise station
% omega       (1,1) [-]    angular frequency of the mode
% beta        (1,1) [-]    spanwise wavenumber of the mode
% bc_top            [-]    freestream BC: 'H_DR' (default) or 'H_NM'
% drop_alpha2       [-]    true zeroes CCxx (old linear GEP, paper comparison)

%% Outputs
% Name size unit   explanation
% Eval      [-]    eigenvalue spectrum (finite alpha only)
% Evec      [-]    eigenvectors [5*ny x n_eig], ordered [rho,u,v,w,T] per node

if ~exist('bc_top','var') || isempty(bc_top); bc_top = 'H_DR'; end
if ~exist('drop_alpha2','var') || isempty(drop_alpha2); drop_alpha2 = false; end

%% Parameters
ny    = StabGrid.ny;
Re    = StabGrid.Re;
Ma    = StabGrid.Ma;
Pr    = StabGrid.Pr;
ga    = StabGrid.gamma;
Sstar = StabGrid.Sstar;
Ec    = (ga - 1) * Ma^2;

%% Base flow at station i
U    = StabGrid.U(:,i);
W    = StabGrid.W(:,i);
Rb   = StabGrid.rho(:,i);
Tb   = StabGrid.T(:,i);
dyU  = StabGrid.dyU(:,i);
dyW  = StabGrid.dyW(:,i);
dyRb = StabGrid.dyrho(:,i);
dyTb = StabGrid.dyT(:,i);

d2yU  = StabGrid.ddyU(:,i);
d2yTb = StabGrid.ddyT(:,i);

[mu_b, dmuT, d2muT] = sutherland(Tb, Sstar);
mu_y = dmuT .* dyTb;

%% Chebyshev operators
[~, DM] = chebdif(ny, 2);
Dn  = DM(:,:,1);
Dn2 = Dn * Dn;
I5  = eye(5);
Total_Dn  = kron(Dn,  I5);
Total_DDn = kron(Dn2, I5);

dzdy   = StabGrid.dzdy;
d2zdy2 = StabGrid.d2zdy2;

%% Assemble coefficient matrices
zz = 1i;
N5 = 5 * ny;
AAL1 = zeros(N5);  AAL2 = zeros(N5);  AAL3 = zeros(N5);
BBR1 = zeros(N5);  BBR2 = zeros(N5);
CCxx = zeros(N5);   % alpha^2 (Vxx, streamwise 2nd-deriv viscous) coefficient,
                    % block-diagonal (no y-derivative). Kept via polyeig below.

gM2 = ga * Ma^2;
Rg  = 1 / gM2;
cv  = Rg / (ga - 1);  % de/dT = cv (ideal gas, nondimensional)
eT  = cv;

for k = 1:ny
    pk = Rb(k) * Tb(k) / gM2;  % base pressure (nondimensional)

    [Lt, Lx, Ly, Lz, Lq, Vxx, Vxy, Vyy, Vxz, Vyz, Vzz] = CLST_Matrix( ...
        Re, Pr, Ec, ...
        Rb(k), dyRb(k), Tb(k), dyTb(k), d2yTb(k), ...
        U(k), dyU(k), d2yU(k), W(k), dyW(k), ...
        mu_b(k), dmuT(k), d2muT(k), mu_y(k), pk, eT, gM2);

    % Row scaling: multiply energy row (5) by Ec for conditioning
    Lt(5,:) = Lt(5,:) * Ec;   Lx(5,:) = Lx(5,:) * Ec;
    Ly(5,:) = Ly(5,:) * Ec;   Lz(5,:) = Lz(5,:) * Ec;
    Lq(5,:) = Lq(5,:) * Ec;
    Vxx(5,:) = Vxx(5,:) * Ec; Vyy(5,:) = Vyy(5,:) * Ec;
    Vzz(5,:) = Vzz(5,:) * Ec;
    Vxy(5,:) = Vxy(5,:) * Ec; Vxz(5,:) = Vxz(5,:) * Ec;
    Vyz(5,:) = Vyz(5,:) * Ec;

    % Spatial eigenvalue problem, quadratic in alpha: (AA - alpha*BB - alpha^2*Vxx) q = 0.
    % A1..B2 are the alpha^0 and alpha^1 coefficients; Vxx is the alpha^2 (streamwise
    % 2nd-deriv viscous) coefficient, assembled into CCxx and solved via polyeig below.
    A1 = -zz*omega*Lt + zz*beta*Lz + Lq - beta^2*Vzz;
    A2 = Ly + zz*beta*Vyz;
    A3 = Vyy;
    B1 = -zz*Lx + beta*Vxz;
    B2 = -zz*Vxy;

    % Apply Malik mapping
    AA1 = A1;
    AA2 = A2 * dzdy(k) + A3 * d2zdy2(k);
    AA3 = A3 * dzdy(k)^2;
    BB1 = B1;
    BB2 = B2 * dzdy(k);

    icc = (k - 1) * 5;
    AAL1(icc+1:icc+5, icc+1:icc+5) = AA1;
    AAL2(icc+1:icc+5, icc+1:icc+5) = AA2;
    AAL3(icc+1:icc+5, icc+1:icc+5) = AA3;
    BBR1(icc+1:icc+5, icc+1:icc+5) = BB1;
    BBR2(icc+1:icc+5, icc+1:icc+5) = BB2;
    CCxx(icc+1:icc+5, icc+1:icc+5) = Vxx;   % d2/dx2 -> -alpha^2*Vxx (Vxx already Ec-scaled on row 5)
end

% Global matrices
AA = AAL1 + AAL2 * Total_Dn + AAL3 * Total_DDn;
BB = BBR1 + BBR2 * Total_Dn;

%% Boundary conditions
% Index 1 = freestream (y = H), index ny = wall (y = 0)
% Ordering within each point: [rho, u, v, w, T]

% Freestream: u = v = w = T = 0  ('H_DR') OR d/dy = 0  ('H_NM')
% The Neumann condition acts on the wall-normal derivative operator Total_Dn.
use_NM_top = strcmpi(bc_top, 'H_NM');
for m = 2:5
    if use_NM_top
        AA(m,:) = Total_Dn(m,:);   % Neumann: dq^/dy = 0 at top
    else
        AA(m,:) = 0;  AA(m,m) = 1; % Dirichlet (default)
    end
    BB(m,:) = 0;   CCxx(m,:) = 0;
end

% Wall: u = v = w = T = 0  (last 4 rows of last point)
for m = 0:3
    row = N5 - 3 + m;
    AA(row,:) = 0;  AA(row,row) = 1;
    BB(row,:) = 0;   CCxx(row,:) = 0;
end

%% Solve QUADRATIC eigenvalue problem  (AA - alpha*BB - alpha^2*CCxx) q = 0
%  via polyeig(AA,-BB,-CCxx). The alpha^2 (streamwise 2nd-deriv viscous) term is
%  essential at finite Mach: dropping it under-damps the 2nd Mack mode (Im(alpha)
%  ~11% low at Ma=4.5). drop_alpha2=true zeroes CCxx to reproduce the old
%  linear generalised eigenvalue problem (GEP) used for paper comparison.
if drop_alpha2
    CCxx = 0*CCxx;
end
[Evec, Eval] = polyeig(AA, -BB, -CCxx);
fin = isfinite(Eval);  Eval = Eval(fin);  Evec = Evec(:, fin);

end

%% ------------------------------------------------------------------------
function [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz] = CLST_Matrix(Re, Pr, Ec, ...
    rk, ryk, tk, tyk, tyyk, uk, uyk, uyyk, wk, wyk, ...
    muk, mutk, muttk, muyk, pk, eT, gM2)
%CLST_MATRIX  5x5 coefficient matrices for one grid point, variables (rho,u,v,w,T).
%   Ideal-gas equation of state (EOS): p = rho*T/(gamma*Ma^2), e = cv*T.
%     muyk = mutk*tyk (y-derivative of viscosity); pk = rk*tk/gM2; eT = de/dT = cv.

% EOS derivatives (ideal gas: p = rho*T/gM2)
p_rho = tk / gM2;
p_T   = rk / gM2;
p_rhoT = 1 / gM2;

%% Lt (temporal: -iw coefficient)
Lt = diag([1, rk, rk, rk, rk*eT]);

%% Lx (streamwise: ia coefficient)
Lx = zeros(5);
Lx(1,1) = uk;              Lx(1,2) = rk;
Lx(2,1) = p_rho;           Lx(2,2) = rk*uk;
Lx(2,3) = -muyk/Re;        Lx(2,5) = p_T;
Lx(3,2) = 2/3*muyk/Re;     Lx(3,3) = rk*uk;
Lx(3,5) = -mutk*uyk/Re;
Lx(4,4) = rk*uk;
Lx(5,2) = pk;              Lx(5,3) = -2*muk*uyk/Re;
Lx(5,5) = rk*uk*eT;

%% Ly (wall-normal: d/dy coefficient)
Ly = zeros(5);
Ly(1,3) = rk;
Ly(2,2) = -muyk/Re;        Ly(2,5) = -mutk*uyk/Re;
Ly(3,1) = p_rho;            Ly(3,3) = -4/3*muyk/Re;
Ly(3,5) = p_T;
Ly(4,4) = -muyk/Re;
Ly(5,2) = -2*muk*uyk/Re;   Ly(5,3) = pk;
Ly(5,5) = -(muyk + mutk*tyk) / (Re*Pr*Ec);

%% Lz (spanwise: ib coefficient)
Lz = zeros(5);
Lz(1,4) = rk;
Lz(3,4) = 2/3*muyk/Re;
Lz(4,1) = p_rho;           Lz(4,3) = -muyk/Re;
Lz(4,5) = p_T;
Lz(5,4) = pk;

%% Lq (zero-derivative terms)
Lq = zeros(5);
Lq(1,3) = ryk;
Lq(2,3) = rk*uyk;
Lq(2,5) = -(mutk*uyyk + muttk*uyk*tyk) / Re;
Lq(3,1) = p_rhoT*tyk;      Lq(3,5) = p_rhoT*ryk;
Lq(4,3) = rk*wyk;
Lq(5,3) = rk*eT*tyk;
Lq(5,5) = -(tyyk*mutk + tyk^2*muttk)/(Re*Pr*Ec) ...
          - mutk*(uyk^2 + wyk^2)/Re;

%% Viscous second-derivative matrices
Vxx = diag([0, -4/3*muk/Re, -muk/Re,     -muk/Re,     -muk/(Re*Pr*Ec)]);
Vyy = diag([0, -muk/Re,     -4/3*muk/Re,  -muk/Re,     -muk/(Re*Pr*Ec)]);
Vzz = diag([0, -muk/Re,     -muk/Re,      -4/3*muk/Re, -muk/(Re*Pr*Ec)]);

Vxy = zeros(5);  Vxy(2,3) = -1/3*muk/Re;  Vxy(3,2) = -1/3*muk/Re;
Vxz = zeros(5);  Vxz(2,4) = -1/3*muk/Re;  Vxz(4,2) = -1/3*muk/Re;
Vyz = zeros(5);  Vyz(3,4) = -1/3*muk/Re;  Vyz(4,3) = -1/3*muk/Re;

end
