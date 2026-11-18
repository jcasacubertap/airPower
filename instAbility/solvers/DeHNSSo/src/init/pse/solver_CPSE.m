function [prof, alpha0, out] = solver_CPSE(omega, beta, StabGrid, Opt)
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
% Compressible linear PSE (parabolized stability equations) inflow condition (Chang & Malik 1993).
%
%   [prof, alpha0, out] = solver_CPSE(omega, beta, StabGrid, Opt)
%
%   Builds a non-parallel inflow shape by marching the linear parabolized
%   stability equations up to the inflow station R0. A self-similar
%   precursor first stretches the inflow plane onto [Rstart, R0]; a local
%   compressible LST (linear stability theory) eigenmode seeds the march,
%   and a first-order backward scheme (BDF1) integrates the density-form
%   operator (state [rho;u;v;w;T], with pressure eliminated through the
%   ideal-gas equation of state). The Chang & Malik normalization supplies
%   the auxiliary closure for alpha.
%
%                 Re, Ma, Pr, gamma, Sstar)
%                 pse_imag_max, pse_alpha0, mode_select, bc_top)

%% Inputs
% Name     size unit   explanation
% omega         [-]    angular frequency of the mode
% beta          [-]    spanwise wavenumber of the mode
% StabGrid      [-]    grid/base-flow struct (U, W, rho, T, y, dzdy, d2zdy2,
% Opt           [-]    options (pse_Rstart, pse_npre, pse_nx, pse_maxit, pse_tol,

%% Outputs
% Name   size unit   explanation
% prof        [-]    inflow shape at R0, normalised by max|u| (u,v,w,rho,T,p)
% alpha0      [-]    complex streamwise wavenumber at R0
% out         [-]    march history (Rv, alpha, phi)

%% Options and defaults
if nargin < 4; Opt = struct(); end
Opt.pse_Rstart  = getopt(Opt,'pse_Rstart', max(0.6*StabGrid.Re, StabGrid.Re-200));
Opt.pse_npre    = getopt(Opt,'pse_npre',    60);
Opt.pse_maxit   = getopt(Opt,'pse_maxit',   40);
Opt.pse_tol     = getopt(Opt,'pse_tol',     1e-9);
Opt.mode_select = getopt(Opt,'mode_select', 'auto');

% Top boundary: reduce a cell bc_top to a scalar tag (Dirichlet unless Neumann)
bc = getopt(Opt,'bc_top','H_DR');
if iscell(bc)
    bc = 'H_DR';
    if any(strcmpi(getopt(Opt,'bc_top','H_DR'),'H_NM')); bc = 'H_NM'; end
end
Opt.bc_top = bc;

%% Precursor and march
SG = precursor(StabGrid, Opt.pse_Rstart, Opt.pse_npre);
[prof, alpha0, out] = march(SG, omega, beta, Opt);

end

%% Local helpers

function SG = precursor(G, Rstart, Np)
%PRECURSOR  Self-similar upstream extension of the inflow plane.
%   Stretches station 1 of G onto Np stations over [Rstart, R0] using the
%   scaling y -> y*(R0/R), and recovers V from continuity.

R0 = G.Re;  ny = G.ny;  y = G.y(:);

% Inflow-plane profiles (station 1)
U0 = G.U(:,1);  W0 = G.W(:,1);  r0 = G.rho(:,1);  T0 = G.T(:,1);

% Streamwise stations and edge node (opposite the wall)
Rv = linspace(Rstart, R0, Np);
xn = Rv.^2 / R0;
[~,kw] = min(U0);  ed = 1;  if kw == 1; ed = ny; end

% Self-similar stretch: thinner BL upstream, so sample the profile at larger y
SG.U = zeros(ny,Np);  SG.W = SG.U;  SG.rho = SG.U;  SG.T = SG.U;
for i = 1:Np
    ys = y*(R0/Rv(i));
    SG.U(:,i)   = interp1(y,U0,ys,'pchip',U0(ed));
    SG.W(:,i)   = interp1(y,W0,ys,'pchip',W0(ed));
    SG.rho(:,i) = interp1(y,r0,ys,'pchip',r0(ed));
    SG.T(:,i)   = interp1(y,T0,ys,'pchip',T0(ed));
end
SG.V = contV(y, xn, SG.U, SG.rho);

% Carry grid metrics and fluid properties
[SG.ny, SG.y, SG.Re, SG.x]         = deal(ny, y, R0, xn);
[SG.Ma, SG.Pr, SG.gamma, SG.Sstar] = deal(G.Ma, G.Pr, G.gamma, G.Sstar);
[SG.dzdy, SG.d2zdy2]               = deal(G.dzdy, G.d2zdy2);

end

function V = contV(y, x, U, R)
%CONTV  Wall-normal velocity from continuity: rho*v = -int_wall^y d(rho u)/dx dy'.
%   v(wall) = 0; integrated in ascending physical y.

[ny,nx] = size(U);
dRUdx = zeros(ny,nx);
for k = 1:ny; dRUdx(k,:) = gradient(R(k,:).*U(k,:), x); end

[ys,is] = sort(y(:));
V = zeros(ny,nx);
for i = 1:nx
    I = zeros(ny,1);
    I(is) = cumtrapz(ys, dRUdx(is,i));
    V(:,i) = -I ./ max(R(:,i), eps);
end

end

function [prof, alpha0, out] = march(SG, omega, beta, Opt)
%MARCH  BDF1 PSE march from Rstart to R0 with the alpha closure.

%% Setup
N  = SG.ny;  Re = SG.Re;  Ma = SG.Ma;  y = SG.y(:);  ga = SG.gamma;
prm = struct('Re',Re, 'kM2',ga*Ma^2, 'Ec',(ga-1)*Ma^2, 'sig',SG.Pr, ...
             'N',N, 'I',eye(N), 'Z',zeros(N), 'iu',1i, 'Omega',ones(N,1));

% Chebyshev wall-normal derivatives, mapped to physical y
[~,DM] = chebdif(N,2);  Dn = DM(:,:,1);
D1 = diag(SG.dzdy(:))*Dn;
D2 = diag(SG.d2zdy2(:))*Dn + diag(SG.dzdy(:).^2)*(Dn*Dn);

%% PSE stations and base flow
R_fine = sqrt(SG.x(:).'*Re);  Rstart = R_fine(1);  R0 = R_fine(end);
nxp = getopt(Opt,'pse_nx', max(6, 1+floor((R0^2-Rstart^2)/Re/(1.5/0.04))));
Rv  = linspace(Rstart, R0, nxp);
xv  = Rv.^2 / Re;
bf  = interp_bf(SG, R_fine, Rv, D1, D2);

%% Inflow seed: local LST eigenmode of L0 + a*L1 + a^2*L2  (a = i*alpha)
E0 = egcol(bf,1, 0,omega,beta,D1,D2,prm);
Ep = egcol(bf,1,+1,omega,beta,D1,D2,prm);
Em = egcol(bf,1,-1,omega,beta,D1,D2,prm);
[L0,L1,L2] = bcevp(E0, (Ep-Em)/2, (Ep+Em)/2-E0, N, Opt.bc_top, D1);
[V,av] = polyeig(L0,L1,L2);
a = -1i*av;  k = isfinite(a) & abs(a) < 1e3;  a = a(k);  V = V(:,k);

% Physical-band filter (pse_imag_max caps damping; 0.05 keeps near-neutral modes)
imx = getopt(Opt,'pse_imag_max',0.05);
b = real(a) > 0.005 & real(a) < 0.5 & imag(a) > -0.1 & imag(a) < imx;
if any(b); a = a(b);  V = V(:,b); end

% Boundary-layer edge (99% U) for EVfilter, then select the physical mode
Uc = bf.U(:,1);
kd = find(Uc <= 0.99*max(Uc), 1);  if isempty(kd); kd = 1; end
if ~isempty(getopt(Opt,'pse_alpha0',[]))
    [~,i] = min(abs(a - Opt.pse_alpha0));  alpha = a(i);
else
    [alpha,i] = EVfilter(a, V(2*N+1:3*N,:), y, omega, beta, y(kd), Ma, Opt.mode_select);
end
phi = V(:,i);

%% Promote LST -> PSE eigenmode (d/dx -> 0 self-consistent)
dx = xv(2) - xv(1);
for it = 1:20
    [E,G] = egcol(bf,1,1i*alpha,omega,beta,D1,D2,prm);
    [A,r] = bcstep(E+G/dx, (G/dx)*phi, N);  pn = A\r;
    da = closure(pn, (pn-phi)/dx, bf.rho(:,1), y, N);
    alpha = alpha + da;  phi = pn;
    if abs(da) < 5e-6; break; end
end

%% March to R0
al  = zeros(1,nxp);   al(1)  = alpha;
PHI = zeros(5*N,nxp); PHI(:,1) = phi;
for j = 2:nxp
    dx = xv(j) - xv(j-1);  pp = PHI(:,j-1);
    if j >= 3; alpha = 2*al(j-1) - al(j-2); else; alpha = al(j-1); end
    for it = 1:Opt.pse_maxit
        [E,G] = egcol(bf,j,1i*alpha,omega,beta,D1,D2,prm);
        [A,r] = bcstep(E+G/dx, (G/dx)*pp, N);  pn = A\r;
        da = closure(pn, (pn-pp)/dx, bf.rho(:,j), y, N);
        alpha = alpha + da;
        if abs(da) < Opt.pse_tol; break; end
    end
    al(j) = alpha;  PHI(:,j) = pn;
end
out = struct('Rv',Rv, 'alpha',al, 'phi',PHI);
alpha0 = al(end);

%% Assemble inflow shape at R0 (normalise by max|u|)
q  = PHI(:,end);
rh = q(1:N);  u = q(N+1:2*N);  v = q(2*N+1:3*N);  w = q(3*N+1:4*N);  T = q(4*N+1:5*N);
p  = (rh.*bf.T(:,end) + bf.rho(:,end).*T) / prm.kM2;
un = max(abs(u));  if un == 0; un = 1; end
prof = struct('u',(u/un).', 'v',(v/un).', 'w',(w/un).', ...
              'rho',(rh/un).', 'T',(T/un).', 'p',(p/un).');

end

function bf = interp_bf(SG, R_fine, Rv, D1, D2)
%INTERP_BF  Interpolate the base flow onto the PSE stations, with derivatives.
%   Wall-normal derivatives from D1/D2, streamwise derivatives by FD in x.

ip  = @(F) interp1(R_fine(:), F.', Rv(:), 'pchip', 'extrap').';
xf  = R_fine.^2 / SG.Re;  xv = Rv.^2 / SG.Re;
ipx = @(F) interp1(xf(:), ddx(F,xf).', xv(:), 'pchip', 'extrap').';

[bf.U, bf.V, bf.W, bf.rho, bf.T]  = deal(ip(SG.U), ip(SG.V), ip(SG.W), ip(SG.rho), ip(SG.T));
[bf.dUx, bf.dVx, bf.dTx, bf.dRx]  = deal(ipx(SG.U), ipx(SG.V), ipx(SG.T), ipx(SG.rho));
bf.dUy = D1*bf.U;  bf.dUyy = D2*bf.U;  bf.dVy = D1*bf.V;
bf.dTy = D1*bf.T;  bf.dTyy = D2*bf.T;  bf.dRy = D1*bf.rho;  bf.dWy = D1*bf.W;

% Transport properties (Sutherland) and calorically-perfect closure
[mu,muT,muTT] = sutherland(bf.T, SG.Sstar);
bf.mu = mu;  bf.dmuT = muT;  bf.dmuTT = muTT;
bf.dPx = zeros(size(bf.U));                        % zero pressure gradient
bf.kap = mu;  bf.kapT = muT;  bf.kapTT = muTT;     % k = mu (constant Pr)
bf.cp  = ones(size(bf.U));  bf.cpT = zeros(size(bf.U));

end

function dF = ddx(F, x)
%DDX  Streamwise gradient of each wall-normal row of F w.r.t. x.
dF = zeros(size(F));
for k = 1:size(F,1); dF(k,:) = gradient(F(k,:), x); end
end

function [E,G] = egcol(bf, i, a, om, be, D1, D2, prm)
%EGCOL  Extract station i of the base flow and assemble the operator blocks.
f = {'U','Ux','Uy','Uyy','V','Vx','Vy','T','Tx','Ty','Tyy','R','Rx','Ry','Px','mu','muT','muTT','kap','kapT','kapTT','cp','cpT'};
s = {'U','dUx','dUy','dUyy','V','dVx','dVy','T','dTx','dTy','dTyy','rho','dRx','dRy','dPx','mu','dmuT','dmuTT','kap','kapT','kapTT','cp','cpT'};
for k = 1:numel(f); c.(f{k}) = bf.(s{k})(:,i); end
[E,G] = operator(c, a, om, be, D1, D2, prm);
end

function [E,G] = operator(c, a, omega, beta, D1, D2, prm)
%OPERATOR  E/G blocks of the density-form PSE operator (Chang & Malik 1993).
%   State [rho;u;v;w;T], pressure eliminated through the ideal-gas EOS.
Re = prm.Re;  kM2 = prm.kM2;  Ec = prm.Ec;  sig = prm.sig;
N = prm.N;  I = prm.I;  Z = prm.Z;  iu = prm.iu;
dOmega = diag(prm.Omega(:));
dU = diag(c.U);  dUy = diag(c.Uy);  dUyy = diag(c.Uyy);  dUx = diag(c.Ux);
dV = diag(c.V);  dVy = diag(c.Vy);  dVx = diag(c.Vx);
dT = diag(c.T);  dTy = diag(c.Ty);  dTyy = diag(c.Tyy);  dTx = diag(c.Tx);
dR = diag(c.R);  dRy = diag(c.Ry);  dRx = diag(c.Rx);
dmu = diag(c.mu);  dmuT = diag(c.muT);  dmuTT = diag(c.muTT);  dmuY = diag(D1*c.mu);
dkap = diag(c.kap);  dkapT = diag(c.kapT);  dkapTT = diag(c.kapTT);
dcp = diag(c.cp);  dcpT = diag(c.cpT);  dPx = diag(c.Px);

E11 = -iu*omega*I+dU*a+dUx+dV*D1+dVy;           E12 = dRx+dR*a; E13 = dRy+dR*D1; E14 = iu*beta*dR; E15 = Z;
E21 = dU*dUx+dV*dUy+(1/kM2)*(dT*a+dTx);
E22 = dR*(-iu*omega*I+dU*a+dV*D1+dUx)-(dmu/Re)*((4/3)*a^2*I-beta^2*I+D2)-(1/Re)*dmuY*D1;
E23 = dR*dUy-(a/Re)*(dmuY+(dmu/3)*D1);           E24 = -iu*beta*dmu*a/(3*Re);
E25 = -(dmuT/Re)*(dUyy+dUy*D1)+(1/kM2)*(dRx+dR*a)-(dmuTT/Re)*dTy*dUy;
E31 = dU*dVx+dV*dVy+(1/kM2)*(dTy+dT*D1);         E32 = -(a/(3*Re))*(dmu*D1-2*dmuY);
E33 = dR*(-iu*omega*I+dU*a+dV*D1+dVy)-(dmu/Re)*(a^2*I-beta^2*I+(4/3)*D2)-(4/(3*Re))*dmuY*D1;
E34 = -(iu*beta/(3*Re))*(dmu*D1-2*dmuY);         E35 = (1/kM2)*(dR*D1+dRy)-(dmuT*a/Re)*dUy;
E41 = iu*beta*dT/kM2;                            E42 = -iu*beta*a*dmu/(3*Re);
E43 = -(iu*beta/Re)*((dmu/3)*D1+dmuY);
E44 = dR*(-iu*omega*I+dU*a+dV*D1)-(dmu/Re)*(a^2*I-(4/3)*beta^2*I+D2)-(1/Re)*dmuY*D1;
E45 = iu*beta*dR/kM2;
E51 = (dcp-(Ec/kM2)*I)*(dU*dTx+dV*dTy)-(Ec*dT/kM2)*(-iu*omega*I+dU*a+dV*D1);
E52 = dcp*dR*dTx-2*Ec*(dmu/Re)*dUy*D1-Ec*dPx;    E53 = dcp*dR*dTy-2*Ec*(dmuT/Re)*dUy*a; E54 = Z;
E55 = dcp*dR*(-iu*omega*I+dU*a+dV*D1)-Ec*(dmuT/Re)*dUy.^2+dR*dcpT*(dU*dTx+dV*dTy) ...
     -(Ec/kM2)*(-iu*omega*dR+dR*dU*a+dR*dV*D1+dV*dRy+dU*dRx) ...
     -(1/(Re*sig))*(dkap*(a^2*I-beta^2*I+D2)+dkapT*(dTyy+2*dTy*D1)+dkapTT*dTy.^2);
E = [E11 E12 E13 E14 E15; E21 E22 E23 E24 E25; E31 E32 E33 E34 E35; E41 E42 E43 E44 E45; E51 E52 E53 E54 E55];

G11 = dU; G12 = dR; G13 = Z; G14 = Z; G15 = Z;
G21 = dOmega*(dT/kM2); G22 = dR*dU-(8/3)*(dmu/Re)*a; G23 = -(1/Re)*(dmuY+(dmu/3)*D1);
G24 = -(iu*beta/(3*Re))*dmu; G25 = dOmega*(dR/kM2);
G31 = Z; G32 = -(1/(3*Re))*(dmu*D1-2*dmuY); G33 = dR*dU-2*(dmu/Re)*a; G34 = Z; G35 = Z;
G41 = Z; G42 = -(iu*beta/(3*Re))*dmu; G43 = Z; G44 = dR*dU-2*(dmu/Re)*a; G45 = Z;
G51 = -Ec*dT*dU/kM2; G52 = Z; G53 = Z; G54 = Z; G55 = dR*dU*(dcp-(Ec/kM2)*I)-2*(dkap/(Re*sig))*a;
G = [G11 G12 G13 G14 G15; G21 G22 G23 G24 G25; G31 G32 G33 G34 G35; G41 G42 G43 G44 G45; G51 G52 G53 G54 G55];

end

function [A,r] = bcstep(A, r, N)
%BCSTEP  Marched-step boundary conditions (BCs): u,v,w,T = 0 at wall (row N) and freestream (row 1); rho free.
for blk = 1:4
    for row = [blk*N+1, blk*N+N]
        A(row,:) = 0;  A(row,row) = 1;  r(row) = 0;
    end
end
end

function da = closure(phi, dphi, rho, y, N)
%CLOSURE  Chang & Malik (2.21) normalisation: int rho (q^H q_x) dy = 0 -> dalpha.
w = rho(:);  num = 0;  den = 0;
for blk = 1:3
    q  = phi(blk*N+1:blk*N+N);
    qx = dphi(blk*N+1:blk*N+N);
    num = num + trapz(y, w.*conj(q).*qx);
    den = den + trapz(y, w.*abs(q).^2);
end
if abs(den) < 1e-20; da = 0; else; da = -1i*num/den; end
end

function [L0,L1,L2] = bcevp(L0, L1, L2, N, bc, D1)
%BCEVP  Eigenproblem BCs: Dirichlet u,v,w,T (H_DR) or Neumann freestream (H_NM); rho free.
nm = strcmpi(bc,'H_NM');
for blk = 1:4
    rf = blk*N+1;  rw = blk*N+N;
    if nm
        L0(rf,:) = 0;  L0(rf,blk*N+1:blk*N+N) = D1(1,:);
    else
        L0(rf,:) = 0;  L0(rf,rf) = 1;
    end
    L0(rw,:) = 0;  L0(rw,rw) = 1;
    L1([rf rw],:) = 0;  L2([rf rw],:) = 0;
end
end

function v = getopt(O, f, dv)
%GETOPT  Return O.f if present and non-empty, else default dv.
if isfield(O,f) && ~isempty(O.(f)); v = O.(f); else; v = dv; end
end
