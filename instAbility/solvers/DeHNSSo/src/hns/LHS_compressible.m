function [MLAMFD,MLA1,MLA2,MLA3,MLA4,MLA5] = LHS_compressible(BF,Bound,Grid,Init,Stab,StabRes,Opt)
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
% Assemble the global sparse left-hand-side (LHS) matrices for
%   the compressible harmonic Navier-Stokes (HNS) system, in curvilinear
%   coordinates with variables (u,v,w,rho,T) and block size bsize=5*ny.
%   Rows: 1=x-momentum 2=y-momentum 3=z-momentum 4=continuity 5=energy.
%   Pressure is eliminated using the equation of state (EOS)
%   p'=(T_bar*rho'+rho_bar*T')/(gamma*Ma^2). Outputs MLAMFD, MLA1..MLA5
%   (see LHS.m for the MLA1..4 convention; MLA5 is the extra beta*d/dxi
%   cross-viscous block needed by the compressible operator).

%% Inputs
% Name    size         unit   explanation
% BF      (1,1) struct [-]    base flow on its own grid
% Bound   (1,1) struct [-]    boundary conditions and outflow-buffer coefficients
% Grid    (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Init    (1,1) struct [-]    inflow initialization data
% Stab    (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md

%% Outputs
% Name   size              unit   explanation
% MLAMFD nvar*ny*nx sparse [-]    operator for the mean-flow-distortion mode (omega = 0, beta = 0)
% MLA1   nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2   nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3   nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4   nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5   nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term

nx    = Grid.nx;
ny    = Grid.ny;
bsize = 5*ny;

if ~isfield(Opt,'parfor'); Opt.parfor = 'off'; end

cell_i1MFD = cell(nx,1);  cell_j1MFD = cell(nx,1);  cell_s1MFD = cell(nx,1);
cell_i1    = cell(nx,1);  cell_j1    = cell(nx,1);  cell_s1    = cell(nx,1);
cell_i2    = cell(nx,1);  cell_j2    = cell(nx,1);  cell_s2    = cell(nx,1);
cell_i3    = cell(nx,1);  cell_j3    = cell(nx,1);  cell_s3    = cell(nx,1);
cell_i4    = cell(nx,1);  cell_j4    = cell(nx,1);  cell_s4    = cell(nx,1);

% Row maps come from bc_rows so that RHS writes the inhomogeneous values into
% exactly the rows constrained here. Freestream constrains u, w, rho, T (v is
% free, it carries the displacement); the wall constrains u, v, w, T.
B           = bc_rows(ny, true);
bc_top      = B.top;
bc_top_mfd  = B.top_mfd;
bc_wall     = B.wall;

% init validates the length, so a short cell here means LHS_compressible was
% called directly. Fail rather than silently substituting all-Dirichlet.
if ~isfield(Opt,'bc_top') || isempty(Opt.bc_top)
    bc_top_type = repmat({'H_DR'}, 1, numel(bc_top));
elseif ~iscell(Opt.bc_top) || numel(Opt.bc_top) ~= numel(bc_top)
    error('LHS_compressible:badBCtop', ...
          'Opt.bc_top must be a cell array of %d entries (u, w, rho, T); got %d.', ...
          numel(bc_top), numel(Opt.bc_top));
else
    bc_top_type = Opt.bc_top;
end
wall_thermal = 'isothermal';
if isfield(Opt,'bc_wall_thermal'); wall_thermal = Opt.bc_wall_thermal; end

iu  = Init.iu;
D1  = Grid.D1;
D2  = Grid.D2;
I   = Stab.I;
Z   = Stab.Z;
dxi = Grid.dxi;
Re  = BF.Re;
ga  = 1.4;
Pr  = 0.72;
if isfield(BF,'gamma'); ga = BF.gamma; end
if isfield(BF,'Pr');    Pr = BF.Pr; end
if ~isfield(BF,'Ma')
    error('LHS_compressible requires BF.Ma (base flow Mach number).');
end
Ma  = BF.Ma;
gM2 = ga*Ma^2;

if ~isfield(Grid,'xixy') || ~isfield(Grid,'etaxy')
    error('LHS_compressible:NoMixedMetrics', ...
        'Grid.xixy/etaxy missing; rebuild the grid with build_stab_grid (mixed metrics).');
end

D1_is_sparse = issparse(Grid.D1) && (nnz(Grid.D1) < 0.2 * ny^2);
use_full = (ny >= 60) && ~D1_is_sparse;

fd_order_x = 4;
if isfield(Grid,'FD_eta_order') && Grid.FD_eta_order == 6; fd_order_x = 6; end
if use_full
    D1f = full(D1);  D2f = full(D2);  If = full(I);  Zf = full(Z);
else
    D1f = D1;  D2f = D2;  If = I;  Zf = Z;
end

% Per-station Fornberg weights for xi derivatives (1st: W1, 2nd: W2), matching the
% base flow (compute_bf_derivatives) and NLT_HNS; correct on non-uniform xi grids
% (replaces the hardcoded uniform-dxi stencils; on a uniform grid they are identical).
% Stencils mirror the assembly branches below: i=2 forward-5, i=nx-1 back-4, i=nx
% back-3, 6th-order interior centered-7, else centered-5. Left-aligned to stenr order.
x_xi = Grid.x(:);
W1 = zeros(nx, 7);  W2 = zeros(nx, 7);
for i = 2:nx
    if i == 2
        cols = 1:5;
    elseif i == nx-1
        cols = nx-3:nx;
    elseif i == nx
        cols = nx-2:nx;
    elseif (fd_order_x == 6) && (i >= 4) && (i <= nx-3)
        cols = i-3:i+3;
    else
        cols = i-2:i+2;
    end
    c = fornberg_weights(x_xi(i), x_xi(cols), 2);
    np = numel(cols);
    W1(i, 1:np) = c(:, 2).';
    W2(i, 1:np) = c(:, 3).';
end

%% Station i=1: inflow identity block
[ii,jj,ss] = find(speye(bsize));
cell_i1MFD{1} = reshape(ii,[],1);  cell_j1MFD{1} = reshape(jj,[],1);  cell_s1MFD{1} = reshape(ss,[],1);
cell_i1{1}    = reshape(ii,[],1);  cell_j1{1}    = reshape(jj,[],1);  cell_s1{1}    = reshape(ss,[],1);
cell_i2{1}    = zeros(0,1);  cell_j2{1}    = zeros(0,1);  cell_s2{1}    = zeros(0,1);
cell_i3{1}    = zeros(0,1);  cell_j3{1}    = zeros(0,1);  cell_s3{1}    = zeros(0,1);
cell_i4{1}    = zeros(0,1);  cell_j4{1}    = zeros(0,1);  cell_s4{1}    = zeros(0,1);
cell_i5{1}    = zeros(0,1);  cell_j5{1}    = zeros(0,1);  cell_s5{1}    = zeros(0,1);

%% Stations i=2:nx
if strcmpi(Opt.parfor, 'on')
    fprintf('   LHS_compressible: assembling stations 2:%d in parallel ...\n', nx);
    parfor i = 2:nx
        [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4, ci5,cj5,cs5] = ...
            lhs_station_full(i, nx, bsize, ny, BF, Grid, Bound, use_full, ...
                D1f, D2f, If, Zf, Re, ga, Pr, gM2, iu, dxi, ...
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x, W1(i,:), W2(i,:));
        cell_i1MFD{i}=ci1MFD; cell_j1MFD{i}=cj1MFD; cell_s1MFD{i}=cs1MFD;
        cell_i1{i}=ci1; cell_j1{i}=cj1; cell_s1{i}=cs1;
        cell_i2{i}=ci2; cell_j2{i}=cj2; cell_s2{i}=cs2;
        cell_i3{i}=ci3; cell_j3{i}=cj3; cell_s3{i}=cs3;
        cell_i4{i}=ci4; cell_j4{i}=cj4; cell_s4{i}=cs4;
        cell_i5{i}=ci5; cell_j5{i}=cj5; cell_s5{i}=cs5;
    end
    fprintf('   LHS_compressible: parallel assembly done\n');
else
    t_lhs = tic;
    lastPerc = 0;
    for i = 2:nx
        [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4, ci5,cj5,cs5] = ...
            lhs_station_full(i, nx, bsize, ny, BF, Grid, Bound, use_full, ...
                D1f, D2f, If, Zf, Re, ga, Pr, gM2, iu, dxi, ...
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x, W1(i,:), W2(i,:));
        cell_i1MFD{i}=ci1MFD; cell_j1MFD{i}=cj1MFD; cell_s1MFD{i}=cs1MFD;
        cell_i1{i}=ci1; cell_j1{i}=cj1; cell_s1{i}=cs1;
        cell_i2{i}=ci2; cell_j2{i}=cj2; cell_s2{i}=cs2;
        cell_i3{i}=ci3; cell_j3{i}=cj3; cell_s3{i}=cs3;
        cell_i4{i}=ci4; cell_j4{i}=cj4; cell_s4{i}=cs4;
        cell_i5{i}=ci5; cell_j5{i}=cj5; cell_s5{i}=cs5;

        perc = floor(i * 100 / nx);
        if perc >= lastPerc + 25 && perc < 100
            elapsed = toc(t_lhs);
            eta_t = elapsed / i * (nx - i);
            fprintf('   LHS_compressible: %d%% (%.0fs elapsed, ~%.0fs remaining)\n', perc, elapsed, eta_t);
            lastPerc = perc;
        end
    end
    fprintf('   LHS_compressible: 100%% (%.1fs total)\n', toc(t_lhs));
end

%% Assemble global sparse matrices from triplets
N = bsize * nx;

MLAMFD = sparse(cell2mat(cell_i1MFD), cell2mat(cell_j1MFD), cell2mat(cell_s1MFD), N, N);
clear cell_i1MFD cell_j1MFD cell_s1MFD

MLA1 = sparse(cell2mat(cell_i1), cell2mat(cell_j1), cell2mat(cell_s1), N, N);
clear cell_i1 cell_j1 cell_s1

if any(StabRes.omegavec)
    MLA2 = sparse(cell2mat(cell_i2), cell2mat(cell_j2), cell2mat(cell_s2), N, N);
else
    MLA2 = sparse(N, N);
end
clear cell_i2 cell_j2 cell_s2

if any(StabRes.betavec)
    MLA3 = sparse(cell2mat(cell_i3), cell2mat(cell_j3), cell2mat(cell_s3), N, N);
    MLA4 = sparse(cell2mat(cell_i4), cell2mat(cell_j4), cell2mat(cell_s4), N, N);
    MLA5 = sparse(vertcat(cell_i5{:}), vertcat(cell_j5{:}), vertcat(cell_s5{:}), N, N);
else
    MLA3 = sparse(N, N);
    MLA4 = sparse(N, N);
    MLA5 = sparse(N, N);
end
clear cell_i3 cell_j3 cell_s3 cell_i4 cell_j4 cell_s4 cell_i5 cell_j5 cell_s5

end

%% Local function: assemble triplets for one station i (i >= 2)
function [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4, ci5,cj5,cs5, ...
          Araw,Braw,Craw,Draw,Eraw,Fraw,Kraw] = ...
    lhs_station_full(i, nx, bsize, ny, BF, Grid, Bound, use_full, ...
                D1f, D2f, If, Zf, Re, ga, Pr, gM2, iu, dxi, ...
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x, w1, w2)
if nargin < 24; fd_order_x = 4; end

% Base-flow fields + x/y gradients as diagonals
if use_full
    U   = diag(BF.U(:,i));    dxU = diag(BF.dxU(:,i));  dyU = diag(BF.dyU(:,i));
    V   = diag(BF.V(:,i));    dxV = diag(BF.dxV(:,i));  dyV = diag(BF.dyV(:,i));
    W   = diag(BF.W(:,i));    dxW = diag(BF.dxW(:,i));  dyW = diag(BF.dyW(:,i));
    Rb  = diag(BF.rho(:,i));  Tb  = diag(BF.T(:,i));
    dxR = diag(BF.dxrho(:,i)); dyR = diag(BF.dyrho(:,i));
    dxT = diag(BF.dxT(:,i));   dyT = diag(BF.dyT(:,i));
else
    U   = spdiags(BF.U(:,i),   0,ny,ny);  dxU = spdiags(BF.dxU(:,i),0,ny,ny);  dyU = spdiags(BF.dyU(:,i),0,ny,ny);
    V   = spdiags(BF.V(:,i),   0,ny,ny);  dxV = spdiags(BF.dxV(:,i),0,ny,ny);  dyV = spdiags(BF.dyV(:,i),0,ny,ny);
    W   = spdiags(BF.W(:,i),   0,ny,ny);  dxW = spdiags(BF.dxW(:,i),0,ny,ny);  dyW = spdiags(BF.dyW(:,i),0,ny,ny);
    Rb  = spdiags(BF.rho(:,i), 0,ny,ny);  Tb  = spdiags(BF.T(:,i),  0,ny,ny);
    dxR = spdiags(BF.dxrho(:,i),0,ny,ny); dyR = spdiags(BF.dyrho(:,i),0,ny,ny);
    dxT = spdiags(BF.dxT(:,i),  0,ny,ny); dyT = spdiags(BF.dyT(:,i),  0,ny,ny);
end

% Viscosity mu(T) + T-derivatives (conductivity Kap=mu)
if isfield(BF, 'mu')
    Mu_v   = BF.mu(:,i);
    MuT_v  = BF.dmuT(:,i);
    if isfield(BF, 'd2muT2')
        MuTT_v = BF.d2muT2(:,i);
    else
        MuTT_v = zeros(ny,1);
    end
else
    Mu_v   = ones(ny,1);
    MuT_v  = zeros(ny,1);
    MuTT_v = zeros(ny,1);
end
if use_full
    Mu   = diag(Mu_v);
    MuT  = diag(MuT_v);
    MuTT = diag(MuTT_v);
else
    Mu   = spdiags(Mu_v,   0, ny, ny);
    MuT  = spdiags(MuT_v,  0, ny, ny);
    MuTT = spdiags(MuTT_v, 0, ny, ny);
end
Kap   = Mu;
KapT  = MuT;
KapTT = MuTT;
% Energy equation in c_v form: rho*DT/Dt = -(ga-1)*rho*T*div(u)
%   + ga/(Re*Pr)*div(Kap*grad(T)) + ga*(ga-1)*Ma^2/Re*Phi.
% Coefficients: conduction ga*Kap/(Pr*Re), dissipation 2*gM2*(ga-1)/Re
% (the 2 comes from linearizing the quadratic Phi), pressure-work (ga-1)*Tb*divU.

% Inverse metrics (1st + 2nd order)
xix  = Grid.xix(:,i);   etax = Grid.etax(:,i);
xiy  = Grid.xiy(:,i);   etay = Grid.etay(:,i);
xixx = Grid.xixx(:,i);  xiyy = Grid.xiyy(:,i);
etaxx= Grid.etaxx(:,i); etayy= Grid.etayy(:,i);
xixy = Grid.xixy(:,i);  etaxy= Grid.etaxy(:,i);

% Outflow-buffer weights (bufs damps the solution; bufp fades the streamwise viscous terms)
bufs_i  = Bound.bufs(:,i);
bufsp_i = Bound.bufsp(:,i);
bufc_i  = Bound.bufc(i);
bufp_i  = Bound.bufp(i);

% Base-flow convection U.grad, divergence, Cxi/Ceta
UdU  = U*dxU + V*dyU;
UdV  = U*dxV + V*dyV;
UdW  = U*dxW + V*dyW;
UdT  = U*dxT + V*dyT;
divU = dxU + dyV;
EnT  = UdT + (ga-1)*Tb*divU;

Cxi  = U.*xix  + V.*xiy;
Ceta = U.*etax + V.*etay;

if use_full
    Mxix = diag(xix);  Mxiy = diag(xiy);  Metax = diag(etax);  Metay = diag(etay);
else
    Mxix = spdiags(xix,0,ny,ny);   Mxiy = spdiags(xiy,0,ny,ny);
    Metax= spdiags(etax,0,ny,ny);  Metay= spdiags(etay,0,ny,ny);
end

% Curvilinear viscous operators (second-derivative chain rule); Vcross_* = mixed dilatation
gxx2 = xix.^2;   gyy2 = xiy.^2;
hxx2 = etax.^2;  hyy2 = etay.^2;
cxe  = 2*xix.*etax;  cye = 2*xiy.*etay;
gxx1 = xixx;     gyy1 = xiyy;
hxx1 = etaxx;    hyy1 = etayy;
ViscA   = @(a,b,coef) -coef*((a*hxx2 + b*hyy2).*D2f + (a*hxx1 + b*hyy1).*D1f);
ViscED1 = @(a,b,coef) -coef*((a*cxe  + b*cye ).*D1f);
ViscEI  = @(a,b,coef) -coef*((a*gxx1 + b*gyy1).*If);
ViscF   = @(a,b,coef) -coef*((a*gxx2 + b*gyy2).*If);

pxy_F  = xix.*xiy;
pxy_E1 = xix.*etay + xiy.*etax;
pxy_E0 = xixy;
pxy_A2 = etax.*etay;
pxy_A0 = etaxy;
Vcross_E = -Mu/3/Re*(pxy_E1.*D1f + pxy_E0.*If);
Vcross_A = -Mu/3/Re*(pxy_A2.*D2f + pxy_A0.*D1f);
Vcross_F = -Mu/3/Re*(pxy_F.*If);

% Wall-normal (eta) advection + viscous, per row
Act_u = Rb*Ceta*D1f + ViscA(4/3, 1,   Mu/Re);
Act_v = Rb*Ceta*D1f + ViscA(1,   4/3, Mu/Re);
Act_w = Rb*Ceta*D1f + ViscA(1,   1,   Mu/Re);
Act_T = Rb*Ceta*D1f + ViscA(1,   1,   ga*Kap/(Pr*Re));

% Base-flow stress tau + dissipation Phi (energy source)
tau_xx = 4/3*dxU - 2/3*dyV;
tau_yy = 4/3*dyV - 2/3*dxU;
tau_xy = dyU + dxV;
tau_xz = dxW;
tau_yz = dyW;
tau_zz = -2/3*divU;

Phi_b_norm = tau_xx.*dxU + tau_xy.*(dyU+dxV) + tau_yy.*dyV ...
           + tau_xz.*dxW + tau_yz.*dyW;

Cont_rho_local = Ceta*D1f + divU;

% Pressure-gradient blocks (EOS p'=(Tb*rho'+Rb*T')/gM2)
PresG_Ax_rho = dxT/gM2;   PresG_Ax_T = dxR/gM2;
PresG_Ay_rho = Tb*Metay*D1f/gM2 + dyT/gM2;
PresG_Ay_T   = Rb*Metay*D1f/gM2 + dyR/gM2;

% Ai: MLA1 diagonal (advection, reaction, pressure gradient, eta-viscous)
Ai = [ ...
    (Act_u + Rb*dxU).*bufs_i,       (Rb*dyU).*bufs_i,                 Zf,                   (UdU + PresG_Ax_rho + Tb*Metax*D1f/gM2).*bufs_i,   (PresG_Ax_T + Rb*Metax*D1f/gM2).*bufs_i;
    (Rb*dxV).*bufs_i,               (Act_v + Rb*dyV).*bufs_i,         Zf,                   (UdV + PresG_Ay_rho).*bufs_i,   (PresG_Ay_T).*bufs_i;
    (Rb*dxW).*bufs_i,               (Rb*dyW).*bufs_i,                 (Act_w).*bufs_i,      (UdW).*bufs_i,                   Zf;
    (dxR + Rb*Metax*D1f).*bufsp_i,  (Rb*Metay*D1f + dyR).*bufsp_i,    Zf,                   (Cont_rho_local).*bufsp_i,       Zf;
    (Rb*dxT + (ga-1)*Tb*Rb*Metax*D1f).*bufs_i,  (Rb*dyT + (ga-1)*Tb*Rb*Metay*D1f).*bufs_i, Zf,   (EnT).*bufs_i,   (Act_T + (ga-1)*Rb*divU).*bufs_i ];

ru = 1:ny;  rv = ny+1:2*ny;  rw = 2*ny+1:3*ny;  rT = 4*ny+1:5*ny;

% Variable-mu mu'(T) corrections into Ai (T-column)
Ai(ru,ru) = Ai(ru,ru) - (MuT*dyT/Re*D1f).*bufs_i;
Ai(rv,rv) = Ai(rv,rv) - (4/3*MuT*dyT/Re*D1f).*bufs_i;
Ai(rw,rw) = Ai(rw,rw) - (MuT*dyT/Re*D1f).*bufs_i;
Ai(rT,rT) = Ai(rT,rT) - (2*ga*KapT/(Pr*Re)*dyT*D1f).*bufs_i;

Ai(ru,rT) = Ai(ru,rT) + (-MuT*tau_xy/Re*D1f).*bufs_i;
Ai(rv,rT) = Ai(rv,rT) + (-MuT*tau_yy/Re*D1f).*bufs_i;

d2yU_v = BF.ddyU(:,i);
d2yV_v = BF.ddyV(:,i);
d2yW_v = BF.ddyW(:,i);
dxyV_v = BF.dxyV(:,i);
dxyU_v = BF.dxyU(:,i);
d2yT_v = BF.ddyT(:,i);
d2xU_v = BF.dxxU(:,i);
d2xV_v = BF.dxxV(:,i);
d2xW_v = BF.dxxW(:,i);
d2xT_v = BF.dxxT(:,i);
gradT2_v = BF.dyT(:,i).^2 + BF.dxT(:,i).^2;
KapT_v  = MuT_v;
KapTT_v = MuTT_v;

Lu_v = (4/3)*d2xU_v + d2yU_v + (1/3)*dxyV_v;
Lv_v = d2xV_v + (4/3)*d2yV_v + (1/3)*dxyU_v;
Lw_v = d2xW_v + d2yW_v;
LT_v = -(ga*KapT_v/(Pr*Re)) .* (d2xT_v + d2yT_v) - (ga*KapTT_v/(Pr*Re)) .* gradT2_v;

if use_full
    Ai(ru,rT) = Ai(ru,rT) + (-diag(MuT_v .* Lu_v)/Re).*bufs_i;
    Ai(rv,rT) = Ai(rv,rT) + (-diag(MuT_v .* Lv_v)/Re).*bufs_i;
    Ai(rw,rT) = Ai(rw,rT) + (-diag(MuT_v .* Lw_v)/Re).*bufs_i;
    Ai(rT,rT) = Ai(rT,rT) + diag(LT_v).*bufs_i;
else
    Ai(ru,rT) = Ai(ru,rT) + spdiags(-MuT_v .* Lu_v/Re, 0, ny, ny).*bufs_i;
    Ai(rv,rT) = Ai(rv,rT) + spdiags(-MuT_v .* Lv_v/Re, 0, ny, ny).*bufs_i;
    Ai(rw,rT) = Ai(rw,rT) + spdiags(-MuT_v .* Lw_v/Re, 0, ny, ny).*bufs_i;
    Ai(rT,rT) = Ai(rT,rT) + spdiags(LT_v, 0, ny, ny).*bufs_i;
end

tau_xx_v = (4/3)*BF.dxU(:,i) - (2/3)*BF.dyV(:,i);
tau_yy_v = (4/3)*BF.dyV(:,i) - (2/3)*BF.dxU(:,i);
tau_xy_v = BF.dyU(:,i) + BF.dxV(:,i);
dxT_vec  = BF.dxT(:,i);   dyT_vec = BF.dyT(:,i);
dxW_vec  = BF.dxW(:,i);   dyW_vec = BF.dyW(:,i);
LMTT_u_v = -(MuTT_v/Re) .* (tau_xx_v .* dxT_vec + tau_xy_v .* dyT_vec);
LMTT_v_v = -(MuTT_v/Re) .* (tau_yy_v .* dyT_vec + tau_xy_v .* dxT_vec);
LMTT_w_v = -(MuTT_v/Re) .* (dyW_vec  .* dyT_vec + dxW_vec  .* dxT_vec);

if use_full
    Ai(ru,rT) = Ai(ru,rT) + diag(LMTT_u_v).*bufs_i;
    Ai(rv,rT) = Ai(rv,rT) + diag(LMTT_v_v).*bufs_i;
    Ai(rw,rT) = Ai(rw,rT) + diag(LMTT_w_v).*bufs_i;
else
    Ai(ru,rT) = Ai(ru,rT) + spdiags(LMTT_u_v, 0, ny, ny).*bufs_i;
    Ai(rv,rT) = Ai(rv,rT) + spdiags(LMTT_v_v, 0, ny, ny).*bufs_i;
    Ai(rw,rT) = Ai(rw,rT) + spdiags(LMTT_w_v, 0, ny, ny).*bufs_i;
end

Ai(ru,rv) = Ai(ru,rv) + (2/3*MuT*dxT/Re*D1f).*bufs_i;
Ai(rv,ru) = Ai(rv,ru) - (MuT*dxT/Re*D1f).*bufs_i;
Ai(rw,rT) = Ai(rw,rT) - (MuT*dyW/Re*D1f).*bufs_i;

pref_diss = 2*gM2*(ga-1)/Re;
Ai(rT,ru) = Ai(rT,ru) - (pref_diss * Mu * tau_xy * D1f).*bufs_i;
Ai(rT,rv) = Ai(rT,rv) - (pref_diss * Mu * tau_yy * D1f).*bufs_i;
Ai(rT,rw) = Ai(rT,rw) - (pref_diss * Mu * tau_yz * D1f).*bufs_i;
Ai(rT,rT) = Ai(rT,rT) - (gM2*(ga-1)/Re * MuT * Phi_b_norm).*bufs_i;

% Bi: MLA2 (x omega): d/dt -> -i*omega
Bi = [ ...
    (-iu*Rb*If).*bufs_i,  Zf,                  Zf,                  Zf,                  Zf;
    Zf,                   (-iu*Rb*If).*bufs_i, Zf,                  Zf,                  Zf;
    Zf,                   Zf,                  (-iu*Rb*If).*bufs_i, Zf,                  Zf;
    Zf,                   Zf,                  Zf,                  (-iu*If).*bufsp_i,   Zf;
    Zf,                   Zf,                  Zf,                  Zf,                  (-iu*Rb*If).*bufs_i ];

% Ci: MLA3 (x beta): i*beta convection/pressure/cross-viscous
Ci = [ ...
    (iu*Rb*W*If).*bufs_i, Zf,  (2/3*iu*MuT*dxT/Re*If).*bufs_i,  Zf,  (-iu*MuT*dxW/Re*If).*bufs_i;
    Zf,  (iu*Rb*W*If).*bufs_i, (-iu*Mu/3/Re*Metay*D1f + 2/3*iu*MuT*dyT/Re*If).*bufs_i, Zf, (-iu*MuT*dyW/Re*If).*bufs_i;
    (-iu*MuT*dxT/Re*If).*bufs_i, (-iu*Mu/3/Re*Metay*D1f - iu*MuT*dyT/Re*If).*bufs_i, (iu*Rb*W*If).*bufs_i, (iu*Tb/gM2*If).*bufs_i, (iu*Rb/gM2*If + 2/3*iu*MuT*divU/Re*If).*bufs_i;
    Zf,  Zf,  (iu*Rb).*bufsp_i,  (iu*W*If).*bufsp_i, Zf;
    Zf,  Zf,  (iu*(ga-1)*Tb*Rb).*bufs_i, Zf, (iu*Rb*W*If).*bufs_i ];

Ci(rT,ru) = Ci(rT,ru) - (iu * pref_diss * Mu * tau_xz).*bufs_i;
Ci(rT,rv) = Ci(rT,rv) - (iu * pref_diss * Mu * tau_yz).*bufs_i;
Ci(rT,rw) = Ci(rT,rw) - (iu * pref_diss * Mu * tau_zz).*bufs_i;

% Di: MLA4 (x beta^2): -mu*beta^2
Di = [ ...
    (Mu/Re*If).*bufs_i,    Zf,                     Zf,                    Zf, Zf;
    Zf,                    (Mu/Re*If).*bufs_i,      Zf,                    Zf, Zf;
    Zf,                    Zf,                      (Mu*4/3/Re*If).*bufs_i, Zf, Zf;
    Zf,                    Zf,                      Zf,                    Zf, Zf;
    Zf,                    Zf,                      Zf,                    Zf, (ga*Kap/(Pr*Re)*If).*bufs_i ];

% Ei: MLA1 d/dxi: advection + cross-metric viscous
Ei = [ ...
    (Rb*Cxi - 4/3*MuT*dxT/Re*If).*bufs_i, (Vcross_E - MuT*dyT/Re*If).*bufs_i, Zf, (Tb*Mxix/gM2).*bufs_i, (Rb*Mxix/gM2 - MuT*tau_xx/Re).*bufs_i;
    (Vcross_E + 2/3*MuT*dyT/Re*If).*bufs_i, (Rb*Cxi - MuT*dxT/Re*If).*bufs_i,  Zf, (Tb*Mxiy/gM2).*bufs_i,  (Rb*Mxiy/gM2 - MuT*tau_xy/Re).*bufs_i;
    Zf,  Zf,  (Rb*Cxi - MuT*dxT/Re*If).*bufs_i,  Zf,  (-MuT*dxW/Re*If).*bufs_i;
    (Rb*Mxix).*bufsp_i, (Rb*Mxiy).*bufsp_i, Zf, (Cxi).*bufsp_i, Zf;
    ((ga-1)*Tb*Rb*Mxix).*bufs_i, ((ga-1)*Tb*Rb*Mxiy).*bufs_i, Zf, Zf, (Rb*Cxi - 2*ga*KapT/(Pr*Re)*dxT).*bufs_i ];

Ei(rT,ru) = Ei(rT,ru) - (pref_diss * Mu * tau_xx).*bufs_i;
Ei(rT,rv) = Ei(rT,rv) - (pref_diss * Mu * tau_xy).*bufs_i;
Ei(rT,rw) = Ei(rT,rw) - (pref_diss * Mu * tau_xz).*bufs_i;

Ei(ru,ru) = Ei(ru,ru) + (ViscED1(4/3,1,Mu/Re)         + ViscEI(4/3,1,Mu/Re)).*bufs_i;
Ei(rv,rv) = Ei(rv,rv) + (ViscED1(1,4/3,Mu/Re)         + ViscEI(1,4/3,Mu/Re)).*bufs_i;
Ei(rw,rw) = Ei(rw,rw) + (ViscED1(1,1,Mu/Re)           + ViscEI(1,1,Mu/Re)).*bufs_i;
Ei(rT,rT) = Ei(rT,rT) + (ViscED1(1,1,ga*Kap/(Pr*Re))  + ViscEI(1,1,ga*Kap/(Pr*Re))).*bufs_i;

Ai(ru,rv) = Ai(ru,rv) + Vcross_A.*bufs_i;
Ai(rv,ru) = Ai(rv,ru) + Vcross_A.*bufs_i;


% Fi: MLA1 d2/dxi2: streamwise viscous
Fi = [ ...
    (ViscF(4/3,1,Mu/Re)).*bufs_i.*bufp_i, Zf,                               Zf,                               Zf, Zf;
    Zf,                              (ViscF(1,4/3,Mu/Re)).*bufs_i.*bufp_i,   Zf,                               Zf, Zf;
    Zf,                              Zf,                              (ViscF(1,1,Mu/Re)).*bufs_i.*bufp_i,      Zf, Zf;
    Zf,                              Zf,                              Zf,                              Zf, Zf;
    Zf,                              Zf,                              Zf,                              Zf, (ViscF(1,1,ga*Kap/(Pr*Re))).*bufs_i.*bufp_i ];

Fi(ru,rv) = Fi(ru,rv) + Vcross_F.*bufs_i.*bufp_i;
Fi(rv,ru) = Fi(rv,ru) + Vcross_F.*bufs_i.*bufp_i;

% Curvilinear variable-viscosity cross-metric corrections (added into Ei and Ai)
MXc = Mxix - If;   MYc = Metay - If;
vvE = @(X,Y) (X*MXc + Y*Mxiy).*bufs_i;
vvA = @(X,Y) (X*Metax*D1f + Y*MYc*D1f).*bufs_i;
Xc=-4/3*MuT*dxT/Re;             Yc=-MuT*dyT/Re;            Ei(ru,ru)=Ei(ru,ru)+vvE(Xc,Yc); Ai(ru,ru)=Ai(ru,ru)+vvA(Xc,Yc);
Xc=-MuT*dxT/Re;                 Yc=-4/3*MuT*dyT/Re;        Ei(rv,rv)=Ei(rv,rv)+vvE(Xc,Yc); Ai(rv,rv)=Ai(rv,rv)+vvA(Xc,Yc);
Xc=-MuT*dxT/Re;                 Yc=-MuT*dyT/Re;            Ei(rw,rw)=Ei(rw,rw)+vvE(Xc,Yc); Ai(rw,rw)=Ai(rw,rw)+vvA(Xc,Yc);
Xc=-2*ga*KapT/(Pr*Re)*dxT;      Yc=-2*ga*KapT/(Pr*Re)*dyT; Ei(rT,rT)=Ei(rT,rT)+vvE(Xc,Yc); Ai(rT,rT)=Ai(rT,rT)+vvA(Xc,Yc);
Xc=-MuT*dyT/Re;                 Yc= 2/3*MuT*dxT/Re;        Ei(ru,rv)=Ei(ru,rv)+vvE(Xc,Yc); Ai(ru,rv)=Ai(ru,rv)+vvA(Xc,Yc);
Xc= 2/3*MuT*dyT/Re;             Yc=-MuT*dxT/Re;            Ei(rv,ru)=Ei(rv,ru)+vvE(Xc,Yc); Ai(rv,ru)=Ai(rv,ru)+vvA(Xc,Yc);
Xc=-MuT*tau_xx/Re;              Yc=-MuT*tau_xy/Re;         Ei(ru,rT)=Ei(ru,rT)+vvE(Xc,Yc); Ai(ru,rT)=Ai(ru,rT)+vvA(Xc,Yc);
Xc=-MuT*tau_xy/Re;             Yc=-MuT*tau_yy/Re;          Ei(rv,rT)=Ei(rv,rT)+vvE(Xc,Yc); Ai(rv,rT)=Ai(rv,rT)+vvA(Xc,Yc);
Xc=-MuT*dxW/Re;                 Yc=-MuT*dyW/Re;            Ei(rw,rT)=Ei(rw,rT)+vvE(Xc,Yc); Ai(rw,rT)=Ai(rw,rT)+vvA(Xc,Yc);
Xc=-pref_diss*Mu*tau_xx;        Yc=-pref_diss*Mu*tau_xy;   Ei(rT,ru)=Ei(rT,ru)+vvE(Xc,Yc); Ai(rT,ru)=Ai(rT,ru)+vvA(Xc,Yc);
Xc=-pref_diss*Mu*tau_xy;        Yc=-pref_diss*Mu*tau_yy;   Ei(rT,rv)=Ei(rT,rv)+vvE(Xc,Yc); Ai(rT,rv)=Ai(rT,rv)+vvA(Xc,Yc);
Xc=-pref_diss*Mu*tau_xz;        Yc=-pref_diss*Mu*tau_yz;   Ei(rT,rw)=Ei(rT,rw)+vvE(Xc,Yc); Ai(rT,rw)=Ai(rT,rw)+vvA(Xc,Yc);

% Gi: freestream/buffer identity; Pres = continuity outflow stencil
if use_full
    Pres = diag(1-bufsp_i) + diag(-1*(1-Bound.bufsp(2:end,i)), -1);
else
    Pres = If.*(1-bufsp_i) + spdiags(-1*(1-Bound.bufsp(2:end,i)), -1, ny, ny);
end
Gi_u = If.*(1-bufs_i) + If.*Bound.bufs(:,end)*(1-bufc_i);
Gi_v = If.*(1-bufs_i);
Gi_w = If.*(1-bufs_i) + If.*Bound.bufs(:,end)*(1-bufc_i);
Gi_T = If.*(1-bufs_i);

Gi = [Gi_u  Zf   Zf   Zf   Zf
       Zf  Gi_v  Zf   Zf   Zf
       Zf   Zf  Gi_w  Zf   Zf
       Zf   Zf   Zf  Pres  Zf
       Zf   Zf   Zf   Zf  Gi_T];

% Ki: MLA5 (x beta*d/dxi): cross-viscous
Ki = [Zf, Zf,              (-iu*Mu/3/Re*Mxix).*bufs_i, Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf;
      (-iu*Mu/3/Re*Mxix).*bufs_i, Zf, Zf,              Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf];
Ci(ru,rw) = Ci(ru,rw) + (-iu*Mu/3/Re*Metax*D1f).*bufs_i;
Ci(rw,ru) = Ci(rw,ru) + (-iu*Mu/3/Re*Metax*D1f).*bufs_i;
Ki(rv,rw) = Ki(rv,rw) + (-iu*Mu/3/Re*Mxiy).*bufs_i;
Ki(rw,rv) = Ki(rw,rv) + (-iu*Mu/3/Re*Mxiy).*bufs_i;

sten = bsize*(i-1)+1 : bsize*i;

is_6th_interior = (fd_order_x == 6) && (i >= 4) && (i <= nx-3);

% Streamwise stencil stations + current-station index. Weights come from the
% precomputed Fornberg W1/W2 (passed as w1/w2), correct on non-uniform xi grids;
% on a uniform grid they reduce to the classical 4th/6th-order coefficients.
if i == 2
    stenr   = [sten-bsize, sten, sten+bsize, sten+2*bsize, sten+3*bsize];
    cur_blk = 2;
elseif i == nx-1
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize];
    cur_blk = 3;
elseif i == nx
    stenr   = [sten-2*bsize, sten-bsize, sten];
    cur_blk = 3;
elseif is_6th_interior
    stenr   = [sten-3*bsize, sten-2*bsize, sten-bsize, sten, sten+bsize, sten+2*bsize, sten+3*bsize];
    cur_blk = 4;
else
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize, sten+2*bsize];
    cur_blk = 3;
end

% blk1 (MLA1): Ei*d/dxi + Fi*d2/dxi2 across the stencil, Ai on the current block;
% blk5 (MLA5, spanwise Ki): first-derivative only (Ki*d/dxi).
nblk = numel(stenr) / bsize;         % number of stencil stations (blocks)
blk1MFD = cell(1, nblk);  blk5c = cell(1, nblk);
for kk = 1:nblk
    blk1MFD{kk} = Ei*w1(kk) + Fi*w2(kk);
    if kk == cur_blk; blk1MFD{kk} = blk1MFD{kk} + Ai + Gi/bufc_i; end
    blk5c{kk} = Ki*w1(kk);
end
blk1MFD = cell2mat(blk1MFD) * bufc_i;
blk5    = cell2mat(blk5c)    * bufc_i;

ncols     = numel(stenr);
blk1      = blk1MFD;
cur_off   = (cur_blk-1)*bsize;

blk2 = zeros(bsize, ncols);
blk3 = zeros(bsize, ncols);
blk4 = zeros(bsize, ncols);
blk2(:, cur_off+1:cur_off+bsize) = Bi*bufc_i;
blk3(:, cur_off+1:cur_off+bsize) = Ci*bufc_i;
blk4(:, cur_off+1:cur_off+bsize) = Di*bufc_i;

% 'free' components keep the governing equation (no imposed boundary condition, BC); see LHS.m.
for k = 1:numel(bc_top)
    if strcmpi(bc_top_type{k}, 'free'); continue; end
    r      = bc_top(k);
    is_mfd = ismember(r, bc_top_mfd);
    blk1(r,:) = 0;  blk2(r,:) = 0;  blk3(r,:) = 0;  blk4(r,:) = 0;  blk5(r,:) = 0;
    if is_mfd; blk1MFD(r,:) = 0; end
    if strcmpi(bc_top_type{k}, 'H_NM')
        cols = cur_off + r + (0:ny-1);
        blk1(r, cols) = D1f(1,:);
        if is_mfd; blk1MFD(r, cols) = D1f(1,:); end
    else
        blk1(r, cur_off+r) = 1;
        if is_mfd; blk1MFD(r, cur_off+r) = 1; end
    end
end


blk1MFD(bc_wall,:) = 0;
blk1(bc_wall,:)    = 0;
blk2(bc_wall,:)    = 0;
blk3(bc_wall,:)    = 0;
blk4(bc_wall,:)    = 0;
blk5(bc_wall,:)    = 0;
for r = bc_wall
    blk1MFD(r, cur_off+r) = 1;
    blk1(r,    cur_off+r) = 1;
end

if strcmpi(wall_thermal, 'adiabatic')
    r  = 5*ny;
    cs = cur_off + 4*ny + 1;
    ce = cur_off + 5*ny;
    blk1(r, :) = 0;
    blk1MFD(r, :) = 0;
    blk1(r, cs:ce) = D1f(ny,:);
    blk1MFD(r, cs:ce) = D1f(ny,:);
end

row_off = bsize*(i-1);

[ii,jj,ss] = find(blk1MFD);  ci1MFD=reshape(ii+row_off,[],1); cj1MFD=reshape(stenr(jj),[],1); cs1MFD=reshape(ss,[],1);
[ii,jj,ss] = find(blk1);     ci1   =reshape(ii+row_off,[],1); cj1   =reshape(stenr(jj),[],1); cs1   =reshape(ss,[],1);
[ii,jj,ss] = find(blk2);     ci2   =reshape(ii+row_off,[],1); cj2   =reshape(stenr(jj),[],1); cs2   =reshape(ss,[],1);
[ii,jj,ss] = find(blk3);     ci3   =reshape(ii+row_off,[],1); cj3   =reshape(stenr(jj),[],1); cs3   =reshape(ss,[],1);
[ii,jj,ss] = find(blk4);     ci4   =reshape(ii+row_off,[],1); cj4   =reshape(stenr(jj),[],1); cs4   =reshape(ss,[],1);
[ii,jj,ss] = find(blk5);     ci5   =reshape(ii+row_off,[],1); cj5   =reshape(stenr(jj),[],1); cs5   =reshape(ss,[],1);

% Raw per-station operator blocks (pre-stencil, pre-BC), for the WKB-marched
% (PSE-form) re-solve of DeHNSSo's own operator. P/Q built by the marcher.
if nargout > 18
    Araw = Ai; Braw = Bi; Craw = Ci; Draw = Di; Eraw = Ei; Fraw = Fi; Kraw = Ki;
end

end
