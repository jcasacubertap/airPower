function [Mode,StabRes] = inflow(Mode,Stab,StabGrid,StabRes,Opt)
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
% Seed the inflow disturbance for all active modes.
%
%   [Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes, Opt)
%
%   Builds the disturbance at the inlet station, one mode at a time. The
%   method is chosen by Stab.IC: 'ILST'/'CLST' (incompressible/compressible
%   LST, linear stability theory) compute a local eigenmode (via local
%   functions ilst_inflow/clst_inflow), 'PSE' (parabolized stability
%   equations) marches a non-parallel seed instead (via
%   solver_IPSE/solver_CPSE), 'LOAD' warm-starts from a saved file, and
%   'ZERO' leaves the inlet undisturbed so the field is generated entirely
%   by wall/freestream forcing. Each mode is scaled by Stab.A0, given a
%   phase reference (wall pressure or max|u|), and written into
%   StabRes.u/v/w/p(/rho/T)/phi at station 1. Mode.RunJ is then updated to
%   list the active modes (those seeded here plus those forced at the wall
%   or freestream through Stab.bcw/bct).

%% Inputs
% Name     size         unit   explanation
% Mode     (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Stab     (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% StabRes  (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% Opt      (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md

%% Outputs
% Name    size         unit   explanation
% Mode    (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

nonzeromodeswall = [];
nonzeromodesic   = [];

% Modes activated by wall/freestream actuation (Stab.bcw / Stab.bct)
nf_bc = min(size(Stab.bcw,3), Mode.nf);
for j = 1:nf_bc
    if any(Stab.bcw(:,:,j),'all') || any(Stab.bct(:,:,j),'all')
        nonzeromodeswall(end+1) = j;
    end
end

is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + double(is_compr);

% Amplitude convention: Stab.A0 (and Stab.bcw, halved in RHS) is the PHYSICAL peak
% amplitude 2|û| of the real disturbance q' = q̂ e^{-iωt} + c.c. The inflow seeds the
% complex amplitude û = A0/2; post-process then reports StabRes.A = 2|û| = A0, so the
% amplitude you set in equals the amplitude you get out (see docs/DeHNSSo_guide.md §5).
A0s = 0.5 * Stab.A0;

StabRes.phi    = zeros(Mode.nf,nvar*StabGrid.ny,StabGrid.nx);
StabRes.alpha  = zeros(Mode.nf,StabGrid.nx);
StabRes.A      = zeros(Mode.nf,StabGrid.nx);
StabRes.u      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.v      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.w      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.p      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
if is_compr
    StabRes.rho = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
    StabRes.T   = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
end

if is_compr && ~(strcmpi(Stab.IC,'CLST') || strcmpi(Stab.IC,'LOAD') || strcmpi(Stab.IC,'PSE') || strcmpi(Stab.IC,'ZERO'))
    error('inflow:badIC', ...
        'Opt.compressible=''on'' requires Stab.IC=''CLST'', ''PSE'', ''LOAD'' or ''ZERO''; got ''%s''.', Stab.IC);
end
if ~is_compr && strcmpi(Stab.IC,'CLST')
    error('inflow:badIC', 'Stab.IC=''CLST'' requires Opt.compressible=''on''.');
end

% Freestream BC flag for the LST inflow solvers: single 'H_DR'/'H_NM', Neumann
% if any component is. clst_inflow uses it; the ILST seed decays so it ignores it.
if ~exist('Opt','var') || ~isfield(Opt,'bc_top')
    bc_top = 'H_DR';
elseif iscell(Opt.bc_top)
    if any(strcmpi(Opt.bc_top, 'H_NM'))
        bc_top = 'H_NM';
    else
        bc_top = 'H_DR';
    end
else
    bc_top = Opt.bc_top;
end
%% Inflow disturbance, one mode at a time
% Every method ends in the same steps for each seeded mode j: a unit-max|u|
% profile is scaled by A0s(j), assembled into phi and rotated to the phase
% reference (seed_mode); the reality partner of each physical mode is its
% complex conjugate (conj_mode). The branches therefore differ only in where
% the profile comes from. Helpers are local functions at the end of the file.
if is_compr; flds = {'u','v','w','rho','T'}; else; flds = {'u','v','w','p'}; end
ny      = StabGrid.ny;
seeded  = find(Stab.A0(:).' > 0);                 % modes given an inflow amplitude
phys    = @(j) j >= round(Mode.nf/2);             % physical mode (not a conjugate)
partner = @(j) 2*round(Mode.nf/2) - j;            % index of the reality partner

switch upper(char(Stab.IC))

case 'LOAD'
    assert(isfield(Stab,'ICfile') && ~isempty(Stab.ICfile), 'Stab.ICfile must be set when Stab.IC = ''LOAD''.');
    assert(isfile(Stab.ICfile), 'ICfile not found: %s', Stab.ICfile);
    L = load(Stab.ICfile);
    if isfield(L,'phiIC')
        StabRes = load_phiIC(L, StabRes, Stab, Mode, nvar, ny, flds);
    elseif isfield(L,'StabRes')
        StabRes = load_checkpoint(L, StabRes, Stab, StabGrid, Mode, nvar, flds);
    else
        error('ICfile must contain a variable named ''StabRes'' or ''phiIC''.');
    end
    if is_compr     % p' from the equation of state: p' = (rho' T + rho T')/(gamma Ma^2)
        gM2 = StabGrid.gamma*StabGrid.Ma^2;
        for j = 1:Mode.nf
            StabRes.p(j,:,1) = (StabRes.rho(j,:,1).*StabGrid.T(:,1).' + StabGrid.rho(:,1).'.*StabRes.T(j,:,1))/gM2;
        end
    end
    for j = flip(seeded)
        if phys(j)
            StabRes = normalise_umax(StabRes, j, flds, StabGrid.y);   % a loaded shape carries an arbitrary scale
            StabRes = seed_mode(StabRes, j, A0s(j), flds, Stab.phaseRef);
        else
            StabRes = conj_mode(StabRes, j, partner(j), flds);
        end
    end
    nonzeromodesic = seeded;
    fprintf('Initialising with LOAD solution from: %s\n', Stab.ICfile);

case 'ILST'
    for j = flip(seeded)
        if phys(j)
            StabRes = ilst_inflow(j, 1, StabGrid, StabRes, bc_top, Opt);    % unit-max|u| eigenfunction
            StabRes = seed_mode(StabRes, j, A0s(j), flds, Stab.phaseRef);
        else
            StabRes = conj_mode(StabRes, j, partner(j), flds);
        end
    end
    nonzeromodesic = seeded;
    if ~any(Stab.bcw(2:end,:), 'all')
        fprintf('Initialising with ILST solution at inflow. \n')
    else
        fprintf('Initialising with ILST solution at inflow and non-homogeneous b.c. at wall. \n')
    end

case {'CLST','PSE'}
    is_pse = strcmpi(char(Stab.IC),'PSE');       % non-parallel linear-PSE seed instead of the local eigenmode
    ov = real(StabRes.omegavec(:)); bv = real(StabRes.betavec(:)); tol = 1e-9;
    % Opt.sym_inflow: solve only the canonical modes (omega>=0, beta>=0) and build
    % the others as their symmetry images on the 2D base flow, which lands the
    % (-1,1) partner of an oblique pair that the index mirror leaves at zero.
    sym_inflow = isfield(Opt,'sym_inflow') && strcmpi(Opt.sym_inflow,'on');
    if sym_inflow
        is_canon = (ov >= -tol) & (bv >= -tol) & (abs(ov) > tol | abs(bv) > tol);
        runorder = [seeded(is_canon(seeded)), seeded(~is_canon(seeded))];   % canonical first
    else
        runorder = flip(seeded);
    end
    for j = runorder
        if (sym_inflow && is_canon(j)) || (~sym_inflow && phys(j))
            if is_pse
                if is_compr; [prof, al] = solver_CPSE(ov(j), bv(j), StabGrid, Opt);
                else;        [prof, al] = solver_IPSE(ov(j), bv(j), StabGrid, Opt); end
                for f = [flds {'p'}]; if isfield(prof,f{1}); StabRes.(f{1})(j,:,1) = prof.(f{1}); end; end
                StabRes.alpha(j,1) = al;
            else
                StabRes = clst_inflow(j, 1, StabGrid, StabRes, bc_top, Opt);
            end
            StabRes = seed_mode(StabRes, j, A0s(j), flds, Stab.phaseRef);
        elseif sym_inflow
            jc = find(abs(abs(ov)-abs(ov(j))) < tol & abs(abs(bv)-abs(bv(j))) < tol & is_canon, 1);
            do_conj = ov(j) < -tol;                          % reality conjugation (omega sign flip)
            flip_w  = (bv(j) < -tol) ~= do_conj;             % spanwise mirror flips w', unless conjugation did
            StabRes = image_mode(StabRes, j, jc, Stab.A0(j)/Stab.A0(jc), do_conj, flip_w, flds);
        else
            StabRes = conj_mode(StabRes, j, partner(j), flds);
        end
    end
    nonzeromodesic = seeded;
    if is_pse; fprintf('Initialising with linear-PSE (non-parallel) solution at inflow.\n');
    else;      fprintf('Initialising with CLST solution at inflow.\n'); end

case 'ZERO'
    % No inflow disturbance: the field comes from wall/freestream forcing only.
    if isempty(nonzeromodeswall)
        warning('inflow:zeroNoForcing', ...
            'Stab.IC=''ZERO'' with no wall/freestream forcing (Stab.bcw/bct all zero): the disturbance field will be identically zero.');
    end
    fprintf('Initialising with ZERO inflow (disturbance from wall/freestream forcing only).\n')

otherwise
    error('Unknown inflow condition Stab.IC = ''%s''. Valid options: ''ILST'', ''CLST'', ''PSE'', ''LOAD'', ''ZERO''.', Stab.IC);
end

% PSE inflow override (Opt.pse_inflow): impose a developed multi-mode NPSE state
% at the inlet (all modes, not just the linear pair) so that generated modes do
% not re-establish inside the domain. File fields: m,n, eta, u,v,w,rho,T,p (physical u_hat).
if isfield(Opt,'pse_inflow') && ~isempty(Opt.pse_inflow)
    P = load(Opt.pse_inflow);
    ovv = real(StabRes.omegavec(:)); bvv = real(StabRes.betavec(:)); tolp = 1e-9;
    o0 = min(abs(ovv(abs(ovv)>tolp))); b0 = min(abs(bvv(abs(bvv)>tolp)));
    if isempty(o0)||o0==0; o0=1; end; if isempty(b0)||b0==0; b0=1; end
    yq = StabGrid.etaun(:).'; loaded = [];
    for j = 1:Mode.nf
        mj = round(ovv(j)/o0); nj = round(bvv(j)/b0);
        k = find(P.m==mj & P.n==nj, 1); cc = false;
        if isempty(k); k = find(P.m==-mj & P.n==-nj, 1); cc = true; end   % reality conjugate
        if isempty(k); continue; end
        gv = @(F) interp1(P.eta(:).', F(k,:), yq, 'pchip', 0);
        u=gv(P.u); v=gv(P.v); w=gv(P.w); r=gv(P.rho); T=gv(P.T); p=gv(P.p);
        if cc; u=conj(u); v=conj(v); w=conj(w); r=conj(r); T=conj(T); p=conj(p); end
        StabRes.u(j,:,1)=u; StabRes.v(j,:,1)=v; StabRes.w(j,:,1)=w;
        StabRes.rho(j,:,1)=r; StabRes.T(j,:,1)=T; StabRes.p(j,:,1)=p;
        StabRes.phi(j,:,1)=[u,v,w,r,T]; loaded=[loaded j];
    end
    nonzeromodesic = unique([nonzeromodesic, loaded]);
    fprintf('PSE inflow: loaded %d modes from %s\n', numel(loaded), Opt.pse_inflow);
end

% nonzero modes are those initiated at inflow + those forced at wall
nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
Mode.RunJ = nonzeromodes;            % Used in loops to assess only active modes
Mode.Aactive(Mode.RunJ) = 1;      % Flag these modes as active

% Define phiIC
StabRes.phiIC = squeeze(StabRes.phi(:,:,1));

end

%% ========================================================================
%  Per-mode seeding steps shared by every inflow method
%% ========================================================================

function StabRes = seed_mode(StabRes, j, A0, flds, phaseRef)
%SEED_MODE  Scale the unit-max|u| profile of mode j by A0, assemble phi in the
% order of flds, and rotate phi to the phase reference ('pwall' or 'umax').
% p' is scaled with the rest also in the compressible case, where it is not
% part of phi, so that the stored fields stay mutually consistent.
sf = flds; if ~any(strcmp(sf,'p')); sf{end+1} = 'p'; end
for f = sf; StabRes.(f{1})(j,:,1) = A0*StabRes.(f{1})(j,:,1); end
phi = [];
for f = flds; phi = [phi, StabRes.(f{1})(j,:,1)]; end %#ok<AGROW>
switch lower(phaseRef)
    case 'pwall'; theta0 = angle(StabRes.p(j,end,1));
    case 'umax';  [~,iy] = max(abs(StabRes.u(j,:,1))); theta0 = angle(StabRes.u(j,iy,1));
    otherwise;    error('Unknown phaseRef option');
end
StabRes.phi(j,:,1) = phi*exp(-1i*theta0);
end

function StabRes = conj_mode(StabRes, j, j2, flds)
%CONJ_MODE  Mode j as the reality partner of mode j2: conjugate shapes and phi,
% alpha_r flips sign, alpha_i (the growth rate) is unchanged.
sf = flds; if ~any(strcmp(sf,'p')); sf{end+1} = 'p'; end
for f = sf; StabRes.(f{1})(j,:,1) = conj(StabRes.(f{1})(j2,:,1)); end
StabRes.alpha(j,1)  = -conj(StabRes.alpha(j2,1));
StabRes.phi(j,:,1)  = conj(StabRes.phi(j2,:,1));
end

function StabRes = normalise_umax(StabRes, j, flds, y)
%NORMALISE_UMAX  Divide the shape of mode j by max|u|, found on a fine spline.
sf = flds; if ~any(strcmp(sf,'p')); sf{end+1} = 'p'; end
yi = linspace(y(1), y(end), 4000);
amp = max(abs(interp1(y, abs(StabRes.u(j,:,1)), yi, 'spline')));
if amp > 0
    for f = sf; StabRes.(f{1})(j,:,1) = StabRes.(f{1})(j,:,1)/amp; end
end
end

function StabRes = image_mode(StabRes, j, jc, sc, do_conj, flip_w, flds)
%IMAGE_MODE  Mode j as a symmetry image of the canonical mode jc on a 2D base
% flow (W = 0): optional reality conjugation (alpha -> -conj(alpha)), optional
% spanwise mirror (w' -> -w'), and rescaling by sc to this mode's amplitude.
sf = flds; if ~any(strcmp(sf,'p')); sf{end+1} = 'p'; end
for f = sf
    q = StabRes.(f{1})(jc,:,1);
    if do_conj; q = conj(q); end
    if flip_w && strcmp(f{1},'w'); q = -q; end
    StabRes.(f{1})(j,:,1) = sc*q;
end
if do_conj; StabRes.alpha(j,1) = -conj(StabRes.alpha(jc,1)); else; StabRes.alpha(j,1) = StabRes.alpha(jc,1); end
phi = [];
for f = flds; phi = [phi, StabRes.(f{1})(j,:,1)]; end %#ok<AGROW>
StabRes.phi(j,:,1) = phi;
end

function StabRes = load_phiIC(L, StabRes, Stab, Mode, nvar, ny, flds)
%LOAD_PHIIC  Inlet profiles from a phiIC array. Accepted layouts: [ny x nvar]
% or [ny x nvar x k] (one column per variable, u,v,w,p or u,v,w,rho,T), or the
% internal [k x nvar*ny] row stacking. Without phiIC_modes the k rows span the
% whole mode list; with it, phiIC_modes is a k-by-2 list of (m,n) naming each row.
P = L.phiIC;
if size(P,1) == ny && size(P,2) == nvar                    % per-variable layout -> row stacking
    P = reshape(P, nvar*ny, []).';
end
assert(size(P,2) == nvar*ny, 'Loaded phiIC has %d columns but %d are expected (nvar*ny = %d*%d).', size(P,2), nvar*ny, nvar, ny);
if isfield(L,'phiIC_modes')
    pm = L.phiIC_modes;
    assert(size(pm,2) == 2 && size(pm,1) == size(P,1), 'phiIC_modes must be [%d x 2], one (m,n) pair per phiIC row; got [%d x %d].', size(P,1), size(pm,1), size(pm,2));
    rows = zeros(size(pm,1),1);
    for r = 1:size(pm,1)
        jr = find(Stab.Mvec(:).' == pm(r,1) & Stab.Nvec(:).' == pm(r,2));
        assert(isscalar(jr), 'phiIC_modes row %d names mode (%d,%d), which is not in this run (M=%d, N=%d).', r, pm(r,1), pm(r,2), Stab.M, Stab.N);
        rows(r) = jr;
    end
    assert(numel(unique(rows)) == numel(rows), 'phiIC_modes names the same mode twice.');
    full = zeros(Mode.nf, nvar*ny); full(rows,:) = P; P = full;
    fprintf('  LOAD: %d user profile(s) mapped to mode(s)%s\n', numel(rows), sprintf(' (%d,%d)', pm.'));
else
    assert(size(P,1) == Mode.nf, 'Loaded phiIC has %d rows but %d are expected (the full mode list). Supply phiIC_modes, a k-by-2 list of (m,n) pairs, to give only some modes.', size(P,1), Mode.nf);
end
for k = 1:numel(flds)
    StabRes.(flds{k})(:,:,1) = P(:, (k-1)*ny+1:k*ny);
end
end

function StabRes = load_checkpoint(L, StabRes, Stab, StabGrid, Mode, nvar, flds)
%LOAD_CHECKPOINT  Warm start from a checkpoint written by Opt.save_checkpoint
% (a results file from Opt.save_output has been through postprocess, which
% strips q, and cannot restore solver state). Fields are rescaled linearly when
% Stab.A0_fund differs from the value the checkpoint was run at.
SR = L.StabRes;
if ~isfield(SR,'q')
    error('inflow:ICfileNotCheckpoint', ['ICfile "%s" has no StabRes.q: it looks like an Opt.save_output results file. ' ...
          'Stab.IC=''LOAD'' needs a checkpoint written with Opt.save_checkpoint.'], Stab.ICfile);
end
assert(size(SR.u,1) == Mode.nf && size(SR.u,2) == StabGrid.ny, ...
    'Loaded StabRes.u size [%d x %d x ...] does not match [%d x %d x ...]: the checkpoint must come from a run with the same M/N and grid.', ...
    size(SR.u,1), size(SR.u,2), Mode.nf, StabGrid.ny);
assert(all(isfield(SR, flds)), 'Loaded StabRes lacks one of the fields %s.', strjoin(flds, ', '));
scale = 1;
if isfield(L,'Stab') && isfield(L.Stab,'A0_fund') && isfield(Stab,'A0_fund') && L.Stab.A0_fund > 0
    scale = Stab.A0_fund/L.Stab.A0_fund;
end
if size(SR.u,3) == StabGrid.nx            % full-domain warm start
    for f = flds; StabRes.(f{1}) = scale*SR.(f{1}); end
    if size(SR.q,2) == nvar*StabGrid.ny*StabGrid.nx; StabRes.q = scale*SR.q; end
    fprintf('  LOAD: full-domain warm-start (%d modes, %d x %d), A0 scale = %.3f\n', size(SR.u,1), size(SR.u,2), size(SR.u,3), scale);
else                                      % inflow station only
    for f = flds; StabRes.(f{1})(:,:,1) = SR.(f{1})(:,:,1); end
    fprintf('  LOAD: inflow-only (station 1), no downstream warm-start\n');
end
end

%% ========================================================================
%  Per-mode inflow seeding: solve the local eigenvalue problem (solver_ILST /
%  solver_CLST), select the physical mode (EVfilter), normalise by max|u|, and
%  store the shape in StabRes. The eigen-solvers themselves live in init/lst/.
%% ========================================================================

function [StabRes] = ilst_inflow(j,i,StabGrid,StabRes,bc_top,Opt)
%ILST_INFLOW  Incompressible LST seed for mode j (Orr-Sommerfeld/Squire).

Re = StabGrid.Re;
U = StabGrid.U(:,i);
W = StabGrid.W(:,i);
dyU = StabGrid.dyU(:,i);
dyW = StabGrid.dyW(:,i);

omega = StabRes.omegavec(j);
beta = StabRes.betavec(j);
ny = StabGrid.ny;

if ~exist('bc_top','var'); bc_top = 'H_DR'; end % seed mode decays: freestream BC immaterial
if ~exist('Opt','var');   Opt.plot = 'on'; end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

assert(isfield(StabGrid,'H') && isfield(StabGrid,'y_i'), 'ILST:gridMeta', ...
    'StabGrid needs H and y_i (regenerate the grid via gridgen).');

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with ILST: ')
end
% solver_ILST rebuilds the StabGrid.y grid (same H, y_i, ny) and returns the
% eigenfunctions on it, so no re-interpolation is needed on Malik grids.
[EV,u,v,w,p,yOS] = solver_ILST(Re,U,W,dyU,dyW,StabGrid.y,omega,beta,ny,StabGrid.H,StabGrid.y_i);

k = find(U<=0.99*max(U),1,'first');
d99 = StabGrid.y(k);
[alpha, index] = EVfilter(EV, v, yOS, omega, beta, d99);

if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

StabRes.alpha(j,i) = alpha;

% Select the physical mode and put it on StabGrid.y (identity for Malik grids
% where yOS == StabGrid.y; a genuine interpolation for FD grids).
u = interp1(yOS, u(:,index), StabGrid.y, 'spline');
v = interp1(yOS, v(:,index), StabGrid.y, 'spline');
w = interp1(yOS, w(:,index), StabGrid.y, 'spline');
p = interp1(yOS, p(:,index), StabGrid.y, 'spline');

% Amplitude = max|u| after fine spline interpolation
y_mdInt     = linspace(StabGrid.y(1),StabGrid.y(end),4000);
u1int       = interp1(StabGrid.y,abs(u),y_mdInt,'spline');
[ampltd, ~] = max(abs(u1int));

StabRes.u(j,:,i)  = u/ampltd;
StabRes.v(j,:,i)  = v/ampltd;
StabRes.w(j,:,i)  = w/ampltd;
StabRes.p(j,:,i)  = p/ampltd;
StabRes.phi(j,:,i) = [ StabRes.u(j,:,i).'; StabRes.v(j,:,i).'; ...
                       StabRes.w(j,:,i).'; StabRes.p(j,:,i).'];

%% Visual check for mode selection
if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1200 350]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(W, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(dyU, StabGrid.y, 'b--', 'LineWidth', 1.0);
plot(dyW, StabGrid.y, 'r--', 'LineWidth', 1.0);
xlabel('Velocity / derivative', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('U', 'W', 'dU/d\eta', 'dW/d\eta', 'Location', 'best', 'FontSize', 9)
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(real(EV), imag(EV), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
alpha_r_center = omega;  % entropy branch at c=1, so alpha_r = omega
xlim([alpha_r_center - 0.1, alpha_r_center + 0.1])
ylim([-0.1, 0.02])
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.u(j,:,1))), StabGrid.y, 'k-', 'LineWidth', 1.5);
xlabel('$|\hat{u}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Eigenfunction', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

sgtitle(sprintf('Inflow ILST: $\\alpha = %.4f %+.4fi$', real(alpha), imag(alpha)), ...
    'Interpreter', 'latex', 'FontSize', 14);
drawnow;

end

function [StabRes] = clst_inflow(j, i, StabGrid, StabRes, bc_top, Opt)
%CLST_INFLOW  Compressible LST seed for mode j, variables (rho,u,v,w,T).

ny    = StabGrid.ny;
Ma    = StabGrid.Ma;
ga    = StabGrid.gamma;
gM2   = ga * Ma^2;
N5    = 5 * ny;

omega = StabRes.omegavec(j);
beta  = StabRes.betavec(j);

if ~exist('bc_top','var'); bc_top = 'H_DR'; end
if ~exist('Opt','var');    Opt = struct('plot','on'); end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with CLST: ')
end

U  = StabGrid.U(:,i);          % for d99 / plotting
Rb = StabGrid.rho(:,i);        % for p reconstruction (EOS)
Tb = StabGrid.T(:,i);

% Assemble and solve the compressible LST eigenvalue problem
drop_a2 = isfield(Opt,'drop_alpha2') && Opt.drop_alpha2;
[Eval, Evec] = solver_CLST(StabGrid, i, omega, beta, bc_top, drop_a2);

% Select the physical mode
if isfield(Opt, 'alpha_shift') && ~isempty(Opt.alpha_shift)
    % Continuation: nearest eigenvalue to the previous-station alpha, within a
    % physical band, so that spurious and continuous-spectrum modes cannot be
    % picked and a streamwise sweep stays smooth.
    band = real(Eval) > 0.05 & real(Eval) < 0.5 & abs(imag(Eval)) < 0.05;
    d = abs(Eval - Opt.alpha_shift);  d(~band) = inf;
    [~, index] = min(d);
    alpha = Eval(index);
else
    % Full spectrum + EVfilter (station 1: mode identification)
    v_evec = Evec(3:5:N5, :);
    kd99 = find(U <= 0.99*max(U), 1, 'first');
    d99_val = StabGrid.y(kd99);
    % Optional Mack-mode selector (Ma>1): Opt.mode_select = 'auto'|'slow'|'fast'
    mode_select = 'auto';
    if isfield(Opt,'mode_select') && ~isempty(Opt.mode_select)
        mode_select = Opt.mode_select;
    end
    [alpha, index] = EVfilter(Eval, v_evec, StabGrid.y, omega, beta, d99_val, Ma, mode_select);
end

StabRes.alpha(j,i) = alpha;
if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

% Extract and normalise eigenfunction. CLST solves in (rho,u,v,w,T); convert to
% (u,v,w,p,T) via EOS: p' = (rho'*T_bar + rho_bar*T') / (gamma*Ma^2).
q       = Evec(1:N5, index);
rho_hat = q(1:5:N5);
u_hat   = q(2:5:N5);
v_hat   = q(3:5:N5);
w_hat   = q(4:5:N5);
T_hat   = q(5:5:N5);
p_hat   = (rho_hat .* Tb + Rb .* T_hat) / gM2;

y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
u_int  = interp1(StabGrid.y, abs(u_hat), y_int, 'spline');
ampltd = max(abs(u_int));

StabRes.u(j,:,i)   = u_hat.'   / ampltd;
StabRes.v(j,:,i)   = v_hat.'   / ampltd;
StabRes.w(j,:,i)   = w_hat.'   / ampltd;
StabRes.p(j,:,i)   = p_hat.'   / ampltd;
StabRes.rho(j,:,i) = rho_hat.' / ampltd;

% Compressible LHS is density form: slot 4 of phi is rho'.
if isfield(StabRes,'T')
    StabRes.T(j,:,i) = T_hat.' / ampltd;
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),    StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),    StabRes.rho(j,:,i), ...
                           StabRes.T(j,:,i)];
else
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),  StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),  StabRes.rho(j,:,i)];
end

%% Plotting
if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1400 350]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
% delta99: y where U first drops below 0.99*U_max (works ascending OR descending y)
kd99_plot    = find(U >= 0.99*max(U), 1, 'last');
if isempty(kd99_plot); kd99_plot = numel(U); end
delta99_plot = StabGrid.y(kd99_plot);
eta_max_plot = min(2 * delta99_plot, max(StabGrid.y));
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(Tb, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(Rb, StabGrid.y, 'g-', 'LineWidth', 1.5);
ylim([0, eta_max_plot]);
xlabel('Base flow', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$U$', '$T$', '$\rho$', 'Location', 'best', 'Interpreter', 'latex')
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
phys = ~isinf(Eval) & ~isnan(Eval) & abs(Eval) < 10;
plot(real(Eval(phys)), imag(Eval(phys)), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlim([real(alpha)-0.2, real(alpha)+0.2]);
ylim([-0.1, 0.05]);
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.u(j,:,i))), StabGrid.y, 'k-', 'LineWidth', 1.5); hold on;
plot(abs(squeeze(StabRes.v(j,:,i))), StabGrid.y, 'b--', 'LineWidth', 1.2);
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$|\hat{u}|$', '$|\hat{v}|$', 'Interpreter', 'latex', 'Location', 'best')
title('Velocity eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.p(j,:,i))), StabGrid.y, 'g-', 'LineWidth', 1.5); hold on;
plot(abs(T_hat)/ampltd,               StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(abs(rho_hat)/ampltd,             StabGrid.y, 'm-', 'LineWidth', 1.5);
legend({'$|\hat{p}|$', '$|\hat{T}|$', '$|\hat{\rho}|$'}, ...
       'Interpreter', 'latex', 'Location', 'best')
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Thermodynamic eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

end

