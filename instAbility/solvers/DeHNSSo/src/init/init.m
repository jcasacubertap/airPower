function [Grid,Init,Mode,Opt,Stab,StabRes] = init(Grid,Opt,Stab)
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
% Initialise DeHNSSo: set defaults, build operators, allocate arrays.
%
%   [Grid, Init, Mode, Opt, Stab, StabRes] = init(Grid, Opt, Stab)
%
%   Sets defaults for Grid, Opt and Stab, checks the freestream boundary
%   conditions, builds the streamwise FD operators (D_xi, D_xi2), sets up the
%   mode tables and harmonic balancing, and allocates the solution arrays.

%% Inputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Opt  (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes

%% Outputs
% Name    size         unit   explanation
% Grid    (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Init    (1,1) struct [-]    inflow initialization data
% Mode    (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab    (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Input validation: required fields
assert(isfield(Stab,'M') && ~isempty(Stab.M),       'Stab.M (spectral truncation omega) is required.');
assert(isfield(Stab,'N') && ~isempty(Stab.N),       'Stab.N (spectral truncation beta) is required.');
assert(isfield(Stab,'omega_0') && ~isempty(Stab.omega_0), 'Stab.omega_0 (fundamental frequency) is required.');
assert(isfield(Stab,'beta_0') && ~isempty(Stab.beta_0),   'Stab.beta_0 (fundamental spanwise wavenumber) is required.');
assert(isfield(Stab,'IC') && ~isempty(Stab.IC),     'Stab.IC (inflow method: ILST/CLST/PSE/LOAD/ZERO) is required.');
assert(isfield(Grid,'nx') && ~isempty(Grid.nx),      'Grid.nx (streamwise points) is required.');
assert(isfield(Grid,'ny') && ~isempty(Grid.ny),      'Grid.ny (wall-normal points) is required.');

% Steps and blocks on the wall: the geometry travels on the grid in Grid.step, and
% step_geometry derives StepType, StepX, StepH and hasMask from it. A flat wall
% is the default and leaves every mask at one (see setup_boundary).
Grid = step_geometry(Grid);
assert(isfield(Grid,'Re') && ~isempty(Grid.Re),      'Grid.Re (Reynolds number) is required.');

if ~isfield(Opt,'compressible'); Opt.compressible = 'off'; end
is_compr  = strcmpi(Opt.compressible, 'on');
nvar      = 4 + double(is_compr);      % (u,v,w,p) | (u,v,w,rho,T)
n_bc_top  = 3 + double(is_compr);      % freestream BC channels: u,v,w | u,w,rho,T (bc_rows)
if is_compr && ~isfield(Opt,'bc_wall_thermal'); Opt.bc_wall_thermal = 'isothermal'; end

% Stab: homogeneous BCs by default
if ~isfield(Stab,'bcx'); Stab.bcx = linspace(Grid.S,Grid.S+Grid.L,10); end
if ~isfield(Stab,'bcw'); Stab.bcw = zeros(n_bc_top,10,(2*Stab.N+1)*(2*Stab.M+1)); end
if ~isfield(Stab,'bct'); Stab.bct = zeros(n_bc_top,10,(2*Stab.N+1)*(2*Stab.M+1)); end
assert(size(Stab.bcw,1) == n_bc_top, 'init:bcwSize', ...
    'Stab.bcw first dim is %d but n_bc_top is %d (Opt.compressible=%s).', ...
    size(Stab.bcw,1), n_bc_top, Opt.compressible);
assert(size(Stab.bct,1) == n_bc_top, 'init:bctSize', ...
    'Stab.bct first dim is %d but n_bc_top is %d (Opt.compressible=%s).', ...
    size(Stab.bct,1), n_bc_top, Opt.compressible);
% Build A0 from A0_fund if user provided the simplified scalar form.
% Fundamental mode is picked automatically from (M, N, omega_0, beta_0):
%   TS         (M>=1, omega_0~=0, N=0, beta_0=0)      -> (m=1, n=0)
%   CFI        (M=0 or omega_0=0, N>=1, beta_0~=0)    -> (m=0, n=1)
%   Oblique    (M>=1, N>=1, both wavenumbers ~=0)     -> (m=1, n=1)
%   Single-oblique-mode (N=0, beta_0~=0; beta encoded via sign(m) in init)
%                                                     -> (m=1, n=0)
if isfield(Stab,'A0_fund') && ~isfield(Stab,'A0')
    nf = (2*Stab.N+1)*(2*Stab.M+1);
    Stab.A0 = zeros(1, nf);

    m_fund = double( Stab.M >= 1 && (Stab.omega_0 ~= 0 || Stab.N == 0) );
    n_fund = double( Stab.N >= 1 && Stab.beta_0  ~= 0 );
    if m_fund == 0 && n_fund == 0
        error(['A0_fund: no fundamental mode: set omega_0~=0 (with M>=1) ' ...
               'or beta_0~=0 (with N>=1), or set Stab.A0 explicitly.']);
    end

    % A0_fund sets the physical amplitude for both linear and nonlinear runs
    % (linear is scale-invariant, so no rescaling needed).
    A_val = Stab.A0_fund;
    Stab.A0(ModeToModeNumber( m_fund,  n_fund, Stab.M, Stab.N)) = A_val;
    Stab.A0(ModeToModeNumber(-m_fund, -n_fund, Stab.M, Stab.N)) = A_val;
end
if ~isfield(Stab,'A0');  Stab.A0 = []; end

% Opt: Set default solver specifications
if ~isfield(Opt,'xb'),          Opt.xb = 85; end
if ~isfield(Opt,'xbe'),         Opt.xbe = 85; end
if ~isfield(Opt,'AFg'),         Opt.AFg = 1.1; end
if ~isfield(Opt,'Sweep'),       Opt.Sweep = 0; end
if ~isfield(Opt,'nltbufxb');    Opt.nltbufxb = Opt.xb; end
if ~isfield(Opt,'solid') || isempty(Opt.solid); Opt.solid = 'auto'; end   % 'auto' | 'mask' | 'ghost': how a surface between grid lines is carried
if ~isfield(Opt,'solid_stencil') || isempty(Opt.solid_stencil); Opt.solid_stencil = 'wall'; end   % 'wall' | 'fluid': streamwise stencils next to a face reach it, as DeHNSSo-v1, or stop short
if Grid.hasMask && isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on')
    error('init:solid', 'Masked solids on the wall (Grid.StepType ~= FLAT) are implemented for the incompressible operator only.');
end
if ~isfield(Opt,'TH');          Opt.TH = 1e-11; end
if ~isfield(Opt,'Conv');        Opt.Conv = 1e-4 ; end
if ~isfield(Opt,'kappa');       Opt.kappa = 6; end
if ~isfield(Opt,'AMAX');        Opt.AMAX = 0.1; end
if ~isfield(Opt,'ConvF');       Opt.ConvF = 100; end
if ~isfield(Opt,'plot');        Opt.plot = 'on'; end
if ~isfield(Opt,'buffer');     Opt.buffer = 'on'; end     % 'on', 'para', or 'off'
if ~isfield(Opt,'allmodes');   Opt.allmodes = 'off'; end  % 'on' activates the full (m,n) mode set (2D oblique cascade); 'off' = 1D incremental introduction
% Freestream BC: per-component type: 3 entries (u,v,w) incompressible,
% 4 entries (u,w,rho,T) compressible. NB v is not a freestream channel in the
% compressible case: it is left free to carry the displacement (see bc_rows).
%   'H_DR'  = homogeneous Dirichlet  (q̂ = 0)
%   'IH_DR' = inhomogeneous Dirichlet (q̂ = bct value)
%   'H_NM'  = homogeneous Neumann    (∂q̂/∂y = 0)
%   'free'  = unconstrained: keep the governing equation (no imposed BC), as the MFD does for v
% Backward compatibility: scalar 'dirichlet'/'neumann' expands to all channels.
if ~isfield(Opt,'bc_top')
    Opt.bc_top = repmat({'H_DR'}, 1, n_bc_top);
elseif ischar(Opt.bc_top) || isstring(Opt.bc_top)
    switch lower(char(Opt.bc_top))
        case 'dirichlet'; Opt.bc_top = repmat({'H_DR'}, 1, n_bc_top);
        case 'neumann';   Opt.bc_top = repmat({'H_NM'}, 1, n_bc_top);
        otherwise;         error('Unknown bc_top: %s', Opt.bc_top);
    end
end
valid_bc = {'H_DR','IH_DR','H_NM','free'};
assert(all(ismember(lower(Opt.bc_top), lower(valid_bc))), ...
    'init:badBCtop', 'bc_top entries must be one of: %s.', strjoin(valid_bc, ', '));
if ~is_compr && numel(Opt.bc_top) >= 3
    u_pinned = strcmpi(Opt.bc_top{1}, 'H_DR') || strcmpi(Opt.bc_top{1}, 'IH_DR');
    w_pinned = strcmpi(Opt.bc_top{3}, 'H_DR') || strcmpi(Opt.bc_top{3}, 'IH_DR');
    beta0    = ~isfield(Stab,'beta_0') || Stab.beta_0 == 0;
    if strcmpi(Opt.bc_top{2}, 'H_NM') && u_pinned && (beta0 || w_pinned)
        error('init:bcTopSingular', ...
            ['Opt.bc_top = {''%s'',''H_NM'',''%s''} with beta_0 = %g makes the '  ...
             'operator singular: with u'' pinned at the freestream and '          ...
             'dw''/dz = 0 there, continuity already forces dv''/dy = 0, so the '  ...
             'Neumann row on v'' is redundant. Use ''free'' on v'' to leave it '  ...
             'unconstrained, or relax u'' (and w'', if beta_0 ~= 0) to ''H_NM''.'], ...
             Opt.bc_top{1}, Opt.bc_top{3}, Stab.beta_0*(~beta0));
    end
end

B_chan = bc_rows(1, is_compr);
assert(numel(Opt.bc_top) == n_bc_top, 'init:bcTopSize', ...
    ['Opt.bc_top has %d entries but %d are expected (Opt.compressible=%s): ' ...
     'the freestream channels are %s.'], numel(Opt.bc_top), n_bc_top, ...
    Opt.compressible, strjoin(B_chan.top_names, ', '));
if ~isfield(Opt,'lu_mode');    Opt.lu_mode = 'auto'; end  % lapack_band|block_thomas|full|auto|seq|none (see solve_modes)
if ~isfield(Opt,'parfor')
    has_pct = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if has_pct; Opt.parfor = 'on'; else; Opt.parfor = 'off'; end
end

% Define imaginary unit
Init.iu=sqrt(-1);

%% Initialization of the modes

% Find all mode interactions
[Mode.Nmat, Mode.Mmat, Mode.Modevec, Stab.Mvec, Stab.Nvec] = Mint(Stab.M,Stab.N);

% Set up Harmonic Balancing matrix
[Stab.HB] = Hbalancing(Stab.M, Stab.N, Mode.Mmat, Mode.Nmat, Mode.Modevec ); 

% Create vectors of omegas and betas per mode
StabRes.omegavec = Stab.omega_0*Stab.Mvec;
if Stab.N == 0 && Stab.beta_0 ~= 0
    % Single oblique mode: beta follows sign of omega (conjugate symmetry)
    % Mode (m,0): beta = sign(m) * beta_0,  MFD (0,0): beta = 0
    StabRes.betavec = Stab.beta_0 * sign(Stab.Mvec);
else
    StabRes.betavec = Stab.beta_0*Stab.Nvec;
end

% Define mode counter L and count total number of modes nf
Mode.L = 1:(2*Stab.N+1)*(2*Stab.M+1); 
Mode.nf = max(Mode.L); 
Mode.Aactive = zeros(size(Mode.L));       % Active-mode register

%% Initialisation of the matrices

% Initiate identity and zero matrix of order ny
Stab.I      = eye(Grid.ny);
Stab.Z      = zeros(Grid.ny);

% LHS sparse blocks (MLA1..MLA5, MLAMFD) are built from triplet (i,j,s)
% accumulators inside LHS / LHS_compressible, so no preallocation is needed.

StabRes.R = zeros(nvar*Grid.nx*Grid.ny, 1);  % Initiate RHS vector
Stab.f    = zeros(Mode.nf, nvar*Grid.ny*Grid.nx);

% Streamwise FD: dxi is the scalar first-spacing; D_xi is the sparse
% 4th-order FD matrix valid on non-uniform ξ grids.
Stab.dxi   = Grid.x(1,2) - Grid.x(1,1);
Grid.D_xi  = fd1_matrix_nonuniform(Grid.x(:));
% Streamwise 2nd-derivative operator (analog of the wall-normal D2): matrix form
% of gridgen's deriv2_uneven (4th-order at every node, full 5-point boundary stencils),
% the same stencil the base flow is differentiated with in compute_bf_derivatives.
Grid.D_xi2 = sparse(deriv2_uneven(Grid.x(:), eye(numel(Grid.x)), 4));

%% Initialization of nonlinear convergence loop

StabRes.Aold = zeros(Mode.nf,Grid.nx);            % Amplitudes of previous iteration
StabRes.q = zeros(Mode.nf, nvar*Grid.ny*Grid.nx);     % Solution vector

Opt.dal = ones(Mode.nf,1);   % Convergence measure
Opt.dalsave(:,1) = Opt.dal; % Convergence measure register
Opt.mmry = 0;           % memory variable so that LHS is created only once
Opt.k1 = 0;             % Iteration counter
Opt.k2 = 0;             % Converged iteration counter
Opt.AF = 1;             % Amplitude factor = 1 initially. damped later to values <1
Opt.AVAL=0; 

end

