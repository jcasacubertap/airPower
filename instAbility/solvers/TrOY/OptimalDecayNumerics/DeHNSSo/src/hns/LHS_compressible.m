function [MLAMFD,MLA1,MLA2,MLA3,MLA4,MLA5] = LHS_compressible(BF,Bound,Grid,Init,Stab,StabRes,Opt)
%LHS_COMPRESSIBLE  Full compressible HNS LHS in (u,v,w,rho,T) variables.
%
%   State per station: q = (u, v, w, rho, T), bsize = 5*ny.
%   Row order: 1=x-mom, 2=y-mom, 3=z-mom, 4=continuity, 5=energy.
%
%   Pressure perturbation eliminated via EOS:
%       p' = (T_bar * rho' + rho_bar * T') / (gamma * Ma^2)
%
%   Formulation follows Scholten et al. (2024), AIAA J., Appendix A.
%
%   Boundary conditions:
%     Wall      u = v = w = T = 0 (no-slip isothermal), rho free
%     Far-field rho = u = w = T = 0, v free (permits acoustic radiation /
%                                            outgoing mass flux)
%
%   Assumptions:
%     - Variable viscosity (Sutherland): BF.mu, BF.dmuT, BF.d2muT2.
%       Falls back to mu=1, derivatives 0 if BF.mu not provided.
%     - Flat-plate identity metrics (warning if not).
%
%   Outputs (same shape as LHS.m, plus MLA5):
%     MLAMFD  mean-flow-distortion LHS (omega=0, beta=0)
%     MLA1    base LHS
%     MLA2    omega coefficient
%     MLA3    beta coefficient
%     MLA4    beta^2 coefficient
%     MLA5    beta * d_x coefficient (iβ×d_x cross viscous terms)

nx    = Grid.nx;
ny    = Grid.ny;
bsize = 5*ny;

if ~isfield(Opt,'parfor'); Opt.parfor = 'off'; end

% Pre-allocate cell arrays for sparse triplets (avoid O(n^2) concatenation)
cell_i1MFD = cell(nx,1);  cell_j1MFD = cell(nx,1);  cell_s1MFD = cell(nx,1);
cell_i1    = cell(nx,1);  cell_j1    = cell(nx,1);  cell_s1    = cell(nx,1);
cell_i2    = cell(nx,1);  cell_j2    = cell(nx,1);  cell_s2    = cell(nx,1);
cell_i3    = cell(nx,1);  cell_j3    = cell(nx,1);  cell_s3    = cell(nx,1);
cell_i4    = cell(nx,1);  cell_j4    = cell(nx,1);  cell_s4    = cell(nx,1);

% BC indices within one bsize block (state slots u,v,w,rho,T):
%   u: 1..ny, v: ny+1..2ny, w: 2ny+1..3ny, rho: 3ny+1..4ny, T: 4ny+1..5ny
%
% BC physics:
%   Wall     : u = v = w = T = 0  (no-slip isothermal), rho free
%   Far-field: rho = u = w = T = 0  (subsonic decay / hypersonic below shock);
%              v free (permits acoustic radiation / outgoing mass flux)
%
% Cell-array user ordering follows the state-vector slot order (u,v,w,rho,T),
% skipping the free slot:
%   Opt.bc_top  = { u-type, w-type, rho-type, T-type }  (v is free, omitted)
%   bc_wall_type= { u-type, v-type, w-type, T-type }    (rho is free, omitted)
bc_top      = [1, 2*ny+1, 3*ny+1, 4*ny+1];   % (u, w, rho, T) at freestream (v free)
bc_top_mfd  = bc_top;                          % same for MFD (rho pinned at top)
bc_wall     = [ny, 2*ny, 3*ny, 5*ny];          % (u, v, w, T) at wall (rho free)

if isfield(Opt,'bc_top') && iscell(Opt.bc_top) && numel(Opt.bc_top) >= 4
    bc_top_type = Opt.bc_top(1:4);
else
    bc_top_type = {'H_DR','H_DR','H_DR','H_DR'};
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
% 4/3 compressible normal-stress factor (Stokes hypothesis, always active).
% Cross terms in E_i/C_i complete the tensor and recover (1/Re)*nabla^2 at Ma→0.

% Streamwise BF second derivatives (viscous Laplacians) are provided by
% gridgen: Grid.dxxU/V/W/T (see compute_bf_derivatives).

% Sanity check: this function assumes flat-plate metrics
is_flat = all(abs(Grid.xix(:)-1) < 1e-10) && all(abs(Grid.etay(:)-1) < 1e-10) && ...
          all(abs(Grid.xiy(:)) < 1e-10)   && all(abs(Grid.etax(:)) < 1e-10);
if ~is_flat
    warning('LHS_compressible:NotFlat', ['LHS_compressible assumes flat-plate identity metrics ' ...
        '(xix=etay=1, xiy=etax=0). Non-flat metrics detected; results may be wrong.']);
end

% Use full (dense) arithmetic for large Chebyshev operators
D1_is_sparse = issparse(Grid.D1) && (nnz(Grid.D1) < 0.2 * ny^2);
use_full = (ny >= 60) && ~D1_is_sparse;

% Streamwise FD order — derived from Grid.FD_eta_order (single knob: same as wall-normal)
fd_order_x = 4;
if isfield(Grid,'FD_eta_order') && Grid.FD_eta_order == 6; fd_order_x = 6; end
if use_full
    D1f = full(D1);  D2f = full(D2);  If = full(I);  Zf = full(Z);
else
    D1f = D1;  D2f = D2;  If = I;  Zf = Z;
end

%% Station i=1: inflow — identity block
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
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x);
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
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x);
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
%
%  q = (u, v, w, rho, T).   bsize = 5*ny.
%  Full compressible flat-plate operator; pressure perturbation eliminated
%  via linearised ideal-gas EOS  p' = (T_bar*rho' + rho_bar*T')/(gamma*Ma^2).
function [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4, ci5,cj5,cs5] = ...
    lhs_station_full(i, nx, bsize, ny, BF, Grid, Bound, use_full, ...
                D1f, D2f, If, Zf, Re, ga, Pr, gM2, iu, dxi, ...
                bc_top, bc_top_mfd, bc_wall, bc_top_type, wall_thermal, fd_order_x)
if nargin < 24; fd_order_x = 4; end

% Base flow diagonal matrices
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

% Density form: rho is a direct state variable (slot 4); pressure is
% eliminated via EOS. Pressure-gradient terms are written directly in the
% momentum rows using Tb, Rb, dxT, dxR.

% Variable viscosity (Sutherland law): mu_bar(y) and dmu/dT(y)
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
% Thermal conductivity: kappa = mu (ideal gas, constant Pr); see equations.tex
Kap   = Mu;
KapT  = MuT;
KapTT = MuTT;

% Metric coefficients (kept for compatibility with LHS.m structure;
% lhs_full assumes flat-plate values: xix=etay=1, xiy=etax=0, all 2nd zero)
xix  = Grid.xix(:,i);   etax = Grid.etax(:,i);
xiy  = Grid.xiy(:,i);   etay = Grid.etay(:,i);
xixx = Grid.xixx(:,i);  xiyy = Grid.xiyy(:,i);
etaxx= Grid.etaxx(:,i); etayy= Grid.etayy(:,i);

% Buffer
bufs_i  = Bound.bufs(:,i);
bufsp_i = Bound.bufsp(:,i);
bufc_i  = Bound.bufc(i);
bufp_i  = Bound.bufp(i);

% ---- Reusable scalar combinations ----
UdU  = U*dxU + V*dyU;            % ubar . grad U
UdV  = U*dxV + V*dyV;            % ubar . grad V
UdW  = U*dxW + V*dyW;            % ubar . grad W
UdT  = U*dxT + V*dyT;            % ubar . grad T
divU = dxU + dyV;                % div ubar (W-independent of z for 2D base)
EnT  = UdT + (ga-1)*Tb*divU;     % energy ρ'-coupling: u.grad T + (γ-1) T div u

% ---- Per-equation wall-normal Act-like operators ----
% Compressible: rho_bar on inertia, mu_bar=1 on viscous, with the (4/3)
% factor on the normal-stress component of each momentum equation.
%   Act_u: x-mom diagonal (1·D2 + 4/3 in Fi-x)
%   Act_v: y-mom diagonal (4/3·D2 + 1   in Fi-x)
%   Act_w: z-mom diagonal (1·D2 + 1   in Fi-x; 4/3 on β² in Di)
%   Act_T: energy diagonal (γκ/(RePr)·D2)
% (Flat plate: η_x=0, η_y=1.)
visc_fac_v  = 4/3/Re;   % 4/3/Re: normal-stress factor (Stokes hypothesis)
Act_u = Rb*V*D1f -      Mu/Re*D2f;
Act_v = Rb*V*D1f - Mu*visc_fac_v*D2f;
Act_w = Rb*V*D1f -      Mu/Re*D2f;
Act_T = Rb*V*D1f - ga*Kap/(Pr*Re)*D2f;

% ---- Coefficient matrices: full compressible 5x5 in (u,v,w,p,T) ----

% Base-flow stress tensor S_ij = tau_ij / mu_bar (general 3D base flow):
% Used for variable-viscosity T̂ couplings and viscous dissipation Phi'.
tau_xx = 4/3*dxU - 2/3*dyV;   % S_xx diagonal matrix
tau_yy = 4/3*dyV - 2/3*dxU;   % S_yy
tau_xy = dyU + dxV;             % S_xy
tau_xz = dxW;                   % S_xz = d_x W_bar  (zero for flat plate W=0)
tau_yz = dyW;                   % S_yz = d_y W_bar
tau_zz = -2/3*divU;             % S_zz = -2/3 div U_bar (nonzero for compressible)

% Phi_b / mu_bar: base-flow dissipation normalised by mu_bar.
% Phi_b = tau_ij * d(U_i)/d(x_j), sum over i in {u,v,w}, j in {x,y} (d_z(Ubar)=0).
Phi_b_norm = tau_xx.*dxU + tau_xy.*(dyU+dxV) + tau_yy.*dyV ...
           + tau_xz.*dxW + tau_yz.*dyW;     % dxW.*dxW + dyW.*dyW at W!=0

% Continuity row 4: direct rho' evolution
%   rho'_t + ubar.grad rho' + rho' div ubar + rho_bar div u' + u'.grad rho_bar = 0
% Wall-normal operator on rho': V*D1 (advection) + divU (local).
% Streamwise advection Ubar*d_x rho' goes to E_i; spanwise iβ*Wbar*rho' to C_i.
Cont_rho_local = V*D1f + divU;   % A_i row-4, rho-col (local + D1 on rho')

% Pressure-gradient chain-rule helpers:
%   p' = (Tbar*rho' + rhobar*T') / (gamma Ma^2)
%   d_j p' = (Tbar d_j rho' + rho_bar d_j T' + d_j Tbar * rho' + d_j rho_bar * T') / gM2
% We only need the LOCAL (no-derivative-of-pert) part for A_i chain rule,
% plus the D1 (wall-normal derivative-of-pert) for y-mom.  Streamwise (x-mom)
% and spanwise (z-mom) derivatives of p' go into E_i and C_i respectively.
PresG_Ax_rho = dxT/gM2;   PresG_Ax_T = dxR/gM2;  % d_x p' chain-rule -> A_i x-mom (row 1)
PresG_Ay_rho = Tb*D1f/gM2 + dyT/gM2;             % d_y p' -> A_i y-mom (row 2) rho-col
PresG_Ay_T   = Rb*D1f/gM2 + dyR/gM2;             % d_y p' -> A_i y-mom (row 2) T-col

% A_i (no-derivative + ∂_y via D1 + ∂²_y via D2)
Ai = [ ...
    (Act_u + Rb*dxU).*bufs_i,       (Rb*dyU).*bufs_i,                 Zf,                   (UdU + PresG_Ax_rho).*bufs_i,   (PresG_Ax_T).*bufs_i;
    (Rb*dxV).*bufs_i,               (Act_v + Rb*dyV).*bufs_i,         Zf,                   (UdV + PresG_Ay_rho).*bufs_i,   (PresG_Ay_T).*bufs_i;
    (Rb*dxW).*bufs_i,               (Rb*dyW).*bufs_i,                 (Act_w).*bufs_i,      (UdW).*bufs_i,                   Zf;
    (dxR).*bufsp_i,                 (Rb*D1f + dyR).*bufsp_i,          Zf,                   (Cont_rho_local).*bufsp_i,       Zf;
    (Rb*dxT).*bufs_i,               (Rb*dyT + (ga-1)*Tb*Rb*D1f).*bufs_i, Zf,                (EnT).*bufs_i,                   (Act_T + (ga-1)*Rb*divU).*bufs_i ];

% ---- Variable-viscosity additions to Ai ----
% Index ranges within bsize block (state order: u, v, w, rho, T)
% Slot 4 is rho (not p); variable-viscosity terms below only touch
% u, v, w, T blocks so indices ru/rv/rw/rT are unchanged.
ru = 1:ny;  rv = ny+1:2*ny;  rw = 2*ny+1:3*ny;  rT = 4*ny+1:5*ny;

% 1. Product-rule D1 terms: ∂_y(mu_bar * ∂_y q̂) = mu_bar*D2*q̂ + mu_T*(∂_y T̄)*D1*q̂
%    The D2 part is already in Act_*; the new D1 part goes on the LHS with a
%    minus sign (same sign as the existing -mu_bar/Re*D2 viscous term):
%    (u,u): -mu_T*(∂_y T̄)/Re * D1      (coeff 1, wall-normal stress on û)
%    (v,v): -4/3*mu_T*(∂_y T̄)/Re * D1  (4/3 normal stress on v̂)
%    (w,w): -mu_T*(∂_y T̄)/Re * D1
%    (T,T): -2*γκ_T*(∂_y T̄)/Re * D1  (factor 2: two contributions)
Ai(ru,ru) = Ai(ru,ru) - (MuT*dyT/Re*D1f).*bufs_i;
Ai(rv,rv) = Ai(rv,rv) - (4/3*MuT*dyT/Re*D1f).*bufs_i;
Ai(rw,rw) = Ai(rw,rw) - (MuT*dyT/Re*D1f).*bufs_i;
Ai(rT,rT) = Ai(rT,rT) - (2*ga*KapT/(Pr*Re)*dyT*D1f).*bufs_i;

% 2. D1 T̂ couplings: ∂_y(mu' * τ̄_ij) contributes mu_T*τ̄_ij * D1*T̂
%    (u,T): -mu_T*τ̄_xy/Re * D1    (from ∂_y(mu'*τ̄_xy) in x-mom)
%    (v,T): -mu_T*τ̄_yy/Re * D1    (from ∂_y(mu'*τ̄_yy) in y-mom)
Ai(ru,rT) = Ai(ru,rT) + (-MuT*tau_xy/Re*D1f).*bufs_i;
Ai(rv,rT) = Ai(rv,rT) + (-MuT*tau_yy/Re*D1f).*bufs_i;

% 3. Local T̂ couplings: mu_T*T̂ * (base-flow viscous Laplacian).
%    Full form with ∂_xx of base flow precomputed in BF.dxxU/V/W/T.
%    (u,T): -mu_T * [(4/3)*∂_xx Ū + ∂_yy Ū + (1/3)*∂_xy V̄] / Re
%    (v,T): -mu_T * [∂_xx V̄ + (4/3)*∂_yy V̄ + (1/3)*∂_xy Ū] / Re
%    (w,T): -mu_T * [∂_xx W̄ + ∂_yy W̄] / Re
%    (T,T): -γκ_T/(Re Pr) * [∂_xx T̄ + ∂_yy T̄] - γκ_TT/(Re Pr) * |∇T̄|²
d2yU_v = D2f * BF.U(:,i);
d2yV_v = D2f * BF.V(:,i);
d2yW_v = D2f * BF.W(:,i);
dxyV_v = D1f * BF.dxV(:,i);
dxyU_v = D1f * BF.dxU(:,i);
d2yT_v = D2f * BF.T(:,i);
d2xU_v = BF.dxxU(:,i);
d2xV_v = BF.dxxV(:,i);
d2xW_v = BF.dxxW(:,i);
d2xT_v = BF.dxxT(:,i);
gradT2_v = BF.dyT(:,i).^2 + BF.dxT(:,i).^2;
KapT_v  = MuT_v;   % kappa_T = mu_T (ideal gas, const Pr)
KapTT_v = MuTT_v;  % kappa_TT = mu_TT

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

% 3b. μ_TT × T̂ scalar couplings in momentum rows.
%     2D base reduction (∂_z(.)̄=0 → T_ζ=u_ζ=v_ζ=w_ζ=0):
%       (u,T): -(μ_TT/Re) × [τ_xx · ∂_x T̄ + τ_xy · ∂_y T̄]
%       (v,T): -(μ_TT/Re) × [τ_yy · ∂_y T̄ + τ_xy · ∂_x T̄]
%       (w,T): -(μ_TT/Re) × [∂_y W̄ · ∂_y T̄ + ∂_x W̄ · ∂_x T̄]   (W̄-dependent)
%     τ_ij here is the base-flow strain-rate-like quantity (no μ̄ baked in),
%     same definition as the matrix `tau_*` used in the dissipation block below.
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

% 4. Non-parallel variable-viscosity additions to Ai: ∂_x μ̄ × ∂_y q̂ terms.
%    Source: (∂_x μ̄)(∂_y v̂) in x-mom from ∂_x τ'_xx,
%            (∂_x μ̄)(∂_y û) in y-mom from ∂_x τ'_xy.
%    (u,v): +(2/3)*mu_T*(∂_x T̄)/Re * D1
%    (v,u): -mu_T*(∂_x T̄)/Re * D1
%    W̄-dependent:
%    (w,T): -mu_T*(∂_y W̄)/Re * D1  (from ∂_y(mu_T T̂ ∂_y W̄) in z-mom)
Ai(ru,rv) = Ai(ru,rv) + (2/3*MuT*dxT/Re*D1f).*bufs_i;
Ai(rv,ru) = Ai(rv,ru) - (MuT*dxT/Re*D1f).*bufs_i;
Ai(rw,rT) = Ai(rw,rT) - (MuT*dyW/Re*D1f).*bufs_i;

% 5. Viscous dissipation Phi' in energy equation.
%    Phi' = 2 * tau_ij * d(u_hat_i)/d(x_j) + mu_T*T_hat * Phi_b/mu_bar
%    Prefactor from energy eq: -gamma*(gamma-1)*Ma^2/Re (LHS sign, moving RHS to left).
%    Factor 2 comes from linearisation of Phi = tau_ij * dU_i/dx_j:
%      bar_tau_ij * d_xj(u_hat_i)   [Part 1: base stress x perturb grad]
%      mu_bar * S_hat_ij * d_xj(Ubar_i) = bar_tau_ij * d_xj(u_hat_i) [Part 2: same]
%    The mu_T*T_hat part contributes a scalar (T,T) coupling via Phi_b/mu_bar.
pref_diss = 2*gM2*(ga-1)/Re;   % = 2*gamma*(gamma-1)*Ma^2/Re
% d_y u_hat terms -> Ai row (energy):
Ai(rT,ru) = Ai(rT,ru) - (pref_diss * Mu * tau_xy * D1f).*bufs_i;
Ai(rT,rv) = Ai(rT,rv) - (pref_diss * Mu * tau_yy * D1f).*bufs_i;
Ai(rT,rw) = Ai(rT,rw) - (pref_diss * Mu * tau_yz * D1f).*bufs_i;  % zero if W=0
% T_hat scalar coupling from mu_T * Phi_b / mu_bar:
Ai(rT,rT) = Ai(rT,rT) - (gM2*(ga-1)/Re * MuT * Phi_b_norm).*bufs_i;

% B_i: -iω coefficient (diag(1, rho, rho, rho, rho)).
% Continuity slot-4 is direct: rho'_t = -iω rho'.  All other rows carry
% rho_bar on the time derivative (ρ u_t etc.).
Bi = [ ...
    (-iu*Rb*If).*bufs_i,  Zf,                  Zf,                  Zf,                  Zf;
    Zf,                   (-iu*Rb*If).*bufs_i, Zf,                  Zf,                  Zf;
    Zf,                   Zf,                  (-iu*Rb*If).*bufs_i, Zf,                  Zf;
    Zf,                   Zf,                  Zf,                  (-iu*If).*bufsp_i,   Zf;
    Zf,                   Zf,                  Zf,                  Zf,                  (-iu*Rb*If).*bufs_i ];

% C_i: iβ coefficient (G).  Includes spanwise advection (W-dependent),
% pressure gradient ∂_z p' = iβ p' in z-mom, ρ̄ ∂_z w' in continuity, and
% iβ viscous cross-terms:
%   y-mom (row 2, w-col): -(1/3/Re)*iβ*D1f          (1/3 μ̄ ∂_y∂_z ŵ)
%   z-mom (row 3, v-col): -(1/3/Re)*iβ*D1f          (1/3 μ̄ ∂_y∂_z v̂)
% Variable-viscosity iβ scalar terms (∂_y μ̄):
%   y-mom (row 2, w-col): +(2/3)*iβ*mu_T*(∂_y T̄)/Re  (G34)
%   z-mom (row 3, v-col): -iβ*mu_T*(∂_y T̄)/Re         (G43)
% Non-parallel variable-viscosity iβ scalar terms (∂_x μ̄):
%   x-mom (row 1, w-col): +(2/3)*iβ*mu_T*(∂_x T̄)/Re  (from -(2/3)∂_x μ̄ iβŵ in x-mom)
%   z-mom (row 3, u-col): -iβ*mu_T*(∂_x T̄)/Re         (from ∂_x μ̄ iβû in z-mom)
% iβ T̂ scalar terms (variable-viscosity × base-flow):
%   x-mom (row 1, T-col): -iβ*mu_T*(∂_x W̄)/Re         (from ∂_z τ'_xz, W̄-dependent)
%   y-mom (row 2, T-col): -iβ*mu_T*(∂_y W̄)/Re         (from ∂_z τ'_yz, W̄-dependent)
%   z-mom (row 3, T-col): +(2/3)*iβ*mu_T*(∂_x Ū+∂_y V̄)/Re  (from ∂_z τ'_zz, compressible)
% Pressure gradient ∂p'/∂z = iβ*(T̄ρ'+ρ̄T')/gM2
% replaces the direct (iu*If) in the old p-column of z-mom.
% Continuity row: iβ*ρ̄ on w-col (ρ̄·∂_z w'), iβ*W̄ on rho-col (W̄·∂_z ρ'),
% T-col = 0 (no EOS substitution).
Ci = [ ...
    (iu*Rb*W*If).*bufs_i, Zf,  (2/3*iu*MuT*dxT/Re*If).*bufs_i,  Zf,  (-iu*MuT*dxW/Re*If).*bufs_i;
    Zf,  (iu*Rb*W*If).*bufs_i, (-iu*Mu/3/Re*D1f + 2/3*iu*MuT*dyT/Re*If).*bufs_i, Zf, (-iu*MuT*dyW/Re*If).*bufs_i;
    (-iu*MuT*dxT/Re*If).*bufs_i, (-iu*Mu/3/Re*D1f - iu*MuT*dyT/Re*If).*bufs_i, (iu*Rb*W*If).*bufs_i, (iu*Tb/gM2*If).*bufs_i, (iu*Rb/gM2*If + 2/3*iu*MuT*divU/Re*If).*bufs_i;
    Zf,  Zf,  (iu*Rb).*bufsp_i,  (iu*W*If).*bufsp_i, Zf;
    Zf,  Zf,  (iu*(ga-1)*Tb*Rb).*bufs_i, Zf, (iu*Rb*W*If).*bufs_i ];

% Viscous dissipation iβ terms in energy row (Ci additions).
% iβ*u_hat terms: -pref_diss * tau_xz * iβ (zero if W=0)
Ci(rT,ru) = Ci(rT,ru) - (iu * pref_diss * Mu * tau_xz).*bufs_i;
% iβ*v_hat terms: -pref_diss * tau_yz * iβ (zero if W=0)
Ci(rT,rv) = Ci(rT,rv) - (iu * pref_diss * Mu * tau_yz).*bufs_i;
% iβ*w_hat terms: -pref_diss * tau_zz * iβ
Ci(rT,rw) = Ci(rT,rw) - (iu * pref_diss * Mu * tau_zz).*bufs_i;

% D_i: +β² coefficient (M, viscous spanwise diffusion).  4/3 on z-mom.
Di = [ ...
    (Mu/Re*If).*bufs_i,    Zf,                     Zf,                    Zf, Zf;
    Zf,                    (Mu/Re*If).*bufs_i,      Zf,                    Zf, Zf;
    Zf,                    Zf,                      (Mu*4/3/Re*If).*bufs_i, Zf, Zf;
    Zf,                    Zf,                      Zf,                    Zf, Zf;
    Zf,                    Zf,                      Zf,                    Zf, (ga*Kap/(Pr*Re)*If).*bufs_i ];

% E_i: ∂_ξ coefficient.  Streamwise advection (with ρ̄), pressure
% gradient ∂_x p' in x-mom, ρ̄ ∂_x u' in continuity, energy pressure work.
% Viscous cross-derivative terms (∂_y∂_x):
%   x-mom (row 1, v-col): -(1/3)μ̄/Re * D1f
%   y-mom (row 2, u-col): -(1/3)μ̄/Re * D1f
% Variable-viscosity ∂_y T̄ scalar terms:
%   x-mom (row 1, v-col): -mu_T*(∂_y T̄)/Re   (D23)
%   y-mom (row 2, u-col): +(2/3)*mu_T*(∂_y T̄)/Re  (D32)
% Non-parallel variable-viscosity ∂_x T̄ diagonal terms (∂_x μ̄ × ∂_x q̂):
%   x-mom (row 1, u-col): -(4/3)*mu_T*(∂_x T̄)/Re
%   y-mom (row 2, v-col): -mu_T*(∂_x T̄)/Re
%   z-mom (row 3, w-col): -mu_T*(∂_x T̄)/Re
% mu' * τ̄_ij * ∂_x T̂ couplings (T̂-column):
%   x-mom (row 1, T-col): -mu_T*τ̄_xx/Re
%   y-mom (row 2, T-col): -mu_T*τ̄_xy/Re
%   z-mom (row 3, T-col): -mu_T*(∂_x W̄)/Re  (W̄-dependent, from τ̄_xz=μ̄∂_x W̄)
% Pressure gradient ∂p'/∂x streamwise-derivative
% contribution becomes Tb/gM2 on rho-col (x-mom) and Rb/gM2 on T-col (x-mom).
% Continuity (row 4) gets Ubar*∂_x rho' on rho-col (direct advection), no T-col.
Ei = [ ...
    (Rb*U - 4/3*MuT*dxT/Re*If).*bufs_i, (-Mu/3/Re*D1f - MuT*dyT/Re*If).*bufs_i, Zf, (Tb/gM2*If).*bufs_i, (Rb/gM2*If - MuT*tau_xx/Re).*bufs_i;
    (-Mu/3/Re*D1f + 2/3*MuT*dyT/Re*If).*bufs_i, (Rb*U - MuT*dxT/Re*If).*bufs_i,  Zf, Zf,                   (-MuT*tau_xy/Re).*bufs_i;
    Zf,  Zf,  (Rb*U - MuT*dxT/Re*If).*bufs_i,  Zf,  (-MuT*dxW/Re*If).*bufs_i;
    (Rb).*bufsp_i, Zf, Zf, (U).*bufsp_i, Zf;
    ((ga-1)*Tb*Rb).*bufs_i, Zf, Zf, Zf, (Rb*U - 2*ga*KapT/(Pr*Re)*dxT).*bufs_i ];

% Viscous dissipation d_x u_hat terms in energy row (Ei additions):
%   -pref_diss * tau_xx * d_x u_hat  (always nonzero)
%   -pref_diss * tau_xy * d_x v_hat  (always nonzero)
%   -pref_diss * tau_xz * d_x w_hat  (zero if W=0)
Ei(rT,ru) = Ei(rT,ru) - (pref_diss * Mu * tau_xx).*bufs_i;
Ei(rT,rv) = Ei(rT,rv) - (pref_diss * Mu * tau_xy).*bufs_i;
Ei(rT,rw) = Ei(rT,rw) - (pref_diss * Mu * tau_xz).*bufs_i;

% F_i: ∂²_ξ coefficient (streamwise viscous diffusion).
%   x-mom diagonal has the 4/3 factor (normal stress on u).
%   y-mom, z-mom, energy: factor 1.
visc_fac_u_fi = 4/3/Re;   % 4/3/Re: streamwise normal stress on u
Fi = [ ...
    (-Mu*visc_fac_u_fi*If).*bufs_i.*bufp_i, Zf,                               Zf,                               Zf, Zf;
    Zf,                                      (-Mu/Re*If).*bufs_i.*bufp_i,      Zf,                               Zf, Zf;
    Zf,                              Zf,                              (-Mu/Re*If).*bufs_i.*bufp_i,      Zf, Zf;
    Zf,                              Zf,                              Zf,                              Zf, Zf;
    Zf,                              Zf,                              Zf,                              Zf, (-ga*Kap/(Pr*Re)*If).*bufs_i.*bufp_i ];

% G_i: buffer / identity terms (same construction as LHS.m)
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

% K_i: iβ × ∂_x coefficient (cross-viscous block).
% These are the (1/3)/Re cross terms coupling x-mom↔w and z-mom↔u via iβ*∂_x.
% Assembled into MLA5 (β-linear, ∂_x stencil), analogous to Ei for MLA1.
Ki = [Zf, Zf,              (-Mu/3/Re*If).*bufs_i, Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf;
      (-Mu/3/Re*If).*bufs_i, Zf, Zf,              Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf;
      Zf, Zf,              Zf,                     Zf, Zf];

% ---- Assemble FD stencil (same structure as LHS.m) ----
sten = bsize*(i-1)+1 : bsize*i;

% Boundary edges (i=2, nx-1, nx) keep 4th-order biased regardless of fd_order_x.
% For 6th-order interior, also keep i=3 and i=nx-2 at 4th-order central
% (they don't have enough room for the 7-point stencil).
is_6th_interior = (fd_order_x == 6) && (i >= 4) && (i <= nx-3);

if i == 2
    stenr   = [sten-bsize, sten, sten+bsize, sten+2*bsize, sten+3*bsize];
    cur_blk = 2;
    blk1MFD = [-3*Ei/(12*dxi)+11*Fi/(12*dxi^2), ...
               (Ai-10*Ei/(12*dxi)-20*Fi/(12*dxi^2))+Gi/bufc_i, ...
                18*Ei/(12*dxi)+ 6*Fi/(12*dxi^2), ...
               - 6*Ei/(12*dxi)+ 4*Fi/(12*dxi^2), ...
                   Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc_i;
    blk5 = [-3*Ki/(12*dxi), -10*Ki/(12*dxi), 18*Ki/(12*dxi), ...
             -6*Ki/(12*dxi),    Ki/(12*dxi)]*bufc_i;

elseif i == nx-1
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize];
    cur_blk = 3;
    blk1MFD = [   Ei/(6*dxi), ...
               -6*Ei/(6*dxi)+  Fi/(dxi^2), ...
               (Ai+3*Ei/(6*dxi)-2*Fi/(dxi^2))+Gi/bufc_i, ...
                2*Ei/(6*dxi)+  Fi/(dxi^2)]*bufc_i;
    blk5 = [Ki/(6*dxi), -6*Ki/(6*dxi), 3*Ki/(6*dxi), 2*Ki/(6*dxi)]*bufc_i;

elseif i == nx
    stenr   = [sten-2*bsize, sten-bsize, sten];
    cur_blk = 3;
    blk1MFD = [   Ei/(2*dxi)+  Fi/(dxi^2), ...
               -4*Ei/(2*dxi)-2*Fi/(dxi^2), ...
               (Ai+3*Ei/(2*dxi)+  Fi/(dxi^2))+Gi/bufc_i]*bufc_i;
    blk5 = [Ki/(2*dxi), -4*Ki/(2*dxi), 3*Ki/(2*dxi)]*bufc_i;

elseif is_6th_interior
    % 6th-order central, 7-point stencil [j-3, j-2, j-1, j, j+1, j+2, j+3]
    % 1st-deriv coeffs: (-1, 9, -45, 0, 45, -9, 1)/(60*dxi)
    % 2nd-deriv coeffs: ( 2,-27, 270,-490, 270,-27, 2)/(180*dxi^2)
    stenr   = [sten-3*bsize, sten-2*bsize, sten-bsize, sten, sten+bsize, sten+2*bsize, sten+3*bsize];
    cur_blk = 4;
    blk1MFD = [  -Ei/(60*dxi) +   2*Fi/(180*dxi^2), ...   % j-3
                9*Ei/(60*dxi) -  27*Fi/(180*dxi^2), ...   % j-2
              -45*Ei/(60*dxi) + 270*Fi/(180*dxi^2), ...   % j-1
              (Ai            - 490*Fi/(180*dxi^2)) + Gi/bufc_i, ...
               45*Ei/(60*dxi) + 270*Fi/(180*dxi^2), ...   % j+1
               -9*Ei/(60*dxi) -  27*Fi/(180*dxi^2), ...   % j+2
                  Ei/(60*dxi) +   2*Fi/(180*dxi^2)] * bufc_i; % j+3
    blk5 = [  -Ki/(60*dxi),  9*Ki/(60*dxi), -45*Ki/(60*dxi), ...
              zeros(bsize,bsize), ...
              45*Ki/(60*dxi), -9*Ki/(60*dxi),    Ki/(60*dxi)] * bufc_i;

else   % i in {3, nx-2} when fd_order_x=6, OR all interior when fd_order_x=4
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize, sten+2*bsize];
    cur_blk = 3;
    blk1MFD = [ Ei/(12*dxi)- 1*Fi/(12*dxi^2), ...
               -8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
               (Ai             -30*Fi/(12*dxi^2))+Gi/bufc_i, ...
                8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
               -  Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc_i;
    blk5 = [ Ki/(12*dxi), -8*Ki/(12*dxi), zeros(bsize,bsize), ...
              8*Ki/(12*dxi),   -Ki/(12*dxi)]*bufc_i;
end

ncols     = numel(stenr);
blk1      = blk1MFD;
cur_off   = (cur_blk-1)*bsize;

% ML2/ML3/ML4 only write to the diagonal (current station) block
blk2 = zeros(bsize, ncols);
blk3 = zeros(bsize, ncols);
blk4 = zeros(bsize, ncols);
blk2(:, cur_off+1:cur_off+bsize) = Bi*bufc_i;
blk3(:, cur_off+1:cur_off+bsize) = Ci*bufc_i;
blk4(:, cur_off+1:cur_off+bsize) = Di*bufc_i;

% ---- Apply freestream BCs per component ----
% bc_top_type cell-array entries map to (u, w, rho, T) in that order
% -- following the state-vector slot ordering (u,v,w,rho,T) with v omitted.
blk1MFD(bc_top_mfd,:) = 0;
blk1(bc_top,:)   = 0;
blk2(bc_top,:)   = 0;
blk3(bc_top,:)   = 0;
blk4(bc_top,:)   = 0;
blk5(bc_top,:)   = 0;

for k = 1:numel(bc_top)
    r  = bc_top(k);
    is_mfd = ismember(r, bc_top_mfd);
    if strcmpi(bc_top_type{k}, 'H_NM')
        cs = cur_off + r;
        ce = cur_off + r + ny - 1;
        blk1(r, cs:ce) = D1f(1,:);
        if is_mfd; blk1MFD(r, cs:ce) = D1f(1,:); end
    else  % H_DR or IH_DR -> identity row
        blk1(r, cur_off+r) = 1;
        if is_mfd; blk1MFD(r, cur_off+r) = 1; end
    end
end

% ---- Apply wall BCs (u=v=w=0, T isothermal/adiabatic; rho free) ----
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

% Adiabatic wall: replace T Dirichlet with Neumann dT'/dy = 0
if strcmpi(wall_thermal, 'adiabatic')
    r  = 5*ny;
    cs = cur_off + 4*ny + 1;
    ce = cur_off + 5*ny;
    blk1(r, :) = 0;
    blk1MFD(r, :) = 0;
    blk1(r, cs:ce) = D1f(ny,:);
    blk1MFD(r, cs:ce) = D1f(ny,:);
end

% ---- Extract sparse triplets ----
row_off = bsize*(i-1);

[ii,jj,ss] = find(blk1MFD);  ci1MFD=reshape(ii+row_off,[],1); cj1MFD=reshape(stenr(jj),[],1); cs1MFD=reshape(ss,[],1);
[ii,jj,ss] = find(blk1);     ci1   =reshape(ii+row_off,[],1); cj1   =reshape(stenr(jj),[],1); cs1   =reshape(ss,[],1);
[ii,jj,ss] = find(blk2);     ci2   =reshape(ii+row_off,[],1); cj2   =reshape(stenr(jj),[],1); cs2   =reshape(ss,[],1);
[ii,jj,ss] = find(blk3);     ci3   =reshape(ii+row_off,[],1); cj3   =reshape(stenr(jj),[],1); cs3   =reshape(ss,[],1);
[ii,jj,ss] = find(blk4);     ci4   =reshape(ii+row_off,[],1); cj4   =reshape(stenr(jj),[],1); cs4   =reshape(ss,[],1);
[ii,jj,ss] = find(blk5);     ci5   =reshape(ii+row_off,[],1); cj5   =reshape(stenr(jj),[],1); cs5   =reshape(ss,[],1);

end
