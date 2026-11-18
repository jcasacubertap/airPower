function [Grid,Init,Mode,Opt,Stab,StabRes] = init(Grid,Opt,Stab)
%INIT  Initialise DeHNSSo: set defaults, build operators, allocate arrays.
%
%   [Grid, Init, Mode, Opt, Stab, StabRes] = init(Grid, Opt, Stab)
%
%   Sets default values for all unspecified fields in Grid, Opt, and Stab.
%   Builds Chebyshev differentiation operators (D1, D2) with Malik mapping,
%   sets up mode interaction matrices (HB), and pre-allocates solution arrays.

%% Input validation — required fields
assert(isfield(Stab,'M') && ~isempty(Stab.M),       'Stab.M (spectral truncation omega) is required.');
assert(isfield(Stab,'N') && ~isempty(Stab.N),       'Stab.N (spectral truncation beta) is required.');
assert(isfield(Stab,'omega_0') && ~isempty(Stab.omega_0), 'Stab.omega_0 (fundamental frequency) is required.');
assert(isfield(Stab,'beta_0') && ~isempty(Stab.beta_0),   'Stab.beta_0 (fundamental spanwise wavenumber) is required.');
assert(isfield(Stab,'IC') && ~isempty(Stab.IC),     'Stab.IC (inflow method: ILST or LOAD) is required.');
assert(isfield(Grid,'nx') && ~isempty(Grid.nx),      'Grid.nx (streamwise points) is required.');
assert(isfield(Grid,'ny') && ~isempty(Grid.ny),      'Grid.ny (wall-normal points) is required.');
assert(isfield(Grid,'Re') && ~isempty(Grid.Re),      'Grid.Re (Reynolds number) is required.');

% Stab: Homogeneous BC's by default
if ~isfield(Stab,'bcx'); Stab.bcx = linspace(Grid.S,Grid.S+Grid.L,10); end
if ~isfield(Stab,'bcw'); Stab.bcw = zeros(3,10,(2*Stab.N+1)*(2*Stab.M+1)); end
if ~isfield(Stab,'bct'); Stab.bct = zeros(3,10,(2*Stab.N+1)*(2*Stab.M+1)); end
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
        error(['A0_fund: no fundamental mode — set omega_0~=0 (with M>=1) ' ...
               'or beta_0~=0 (with N>=1), or set Stab.A0 explicitly.']);
    end

    if isfield(Opt,'linear') && strcmpi(Opt.linear,'on')
        A_val = 1e-10;   % small amplitude keeps NLT inactive (1 iteration)
        A_val =  Stab.A0_fund;   % Antonis

    else
        A_val = Stab.A0_fund;
    end
    Stab.A0(ModeToModeNumber( m_fund,  n_fund, Stab.M, Stab.N)) = A_val;
    Stab.A0(ModeToModeNumber(-m_fund, -n_fund, Stab.M, Stab.N)) = A_val;
end
if ~isfield(Stab,'A0');  Stab.A0 = []; end

% % Grid: Straight geometry without streamwise refinement by default
if ~isfield(Grid,'ytype'),      Grid.ytype ="malik"; end
if ~isfield(Grid,'xtype'),      Grid.xtype ="equidistant"; end
if ~isfield(Grid,'mug'),        Grid.mug = 0; end
if ~isfield(Grid,'sig'),        Grid.sig = 1; end
if ~isfield(Grid,'ag'),         Grid.ag = 0; end
if ~isfield(Grid,'ft'),         Grid.ft = 1; end

% Opt: Set default solver specifications
if ~isfield(Opt,'xb'),          Opt.xb = 85; end
if ~isfield(Opt,'xbe'),         Opt.xbe = 85; end
if ~isfield(Opt,'AFg'),         Opt.AFg = 1.1; end
if ~isfield(Opt,'Sweep'),       Opt.Sweep = 0; end
if ~isfield(Opt,'nltbufxb');    Opt.nltbufxb = Opt.xb; end
if ~isfield(Opt,'TH');          Opt.TH = 1e-11; end
if ~isfield(Opt,'Conv');        Opt.Conv = 1e-4 ; end
if ~isfield(Opt,'kappa');       Opt.kappa = 6; end
if ~isfield(Opt,'AMAX');        Opt.AMAX = 0.1; end
if ~isfield(Opt,'ConvF');       Opt.ConvF = 100; end
if ~isfield(Opt,'plot');        Opt.plot = 'on'; end
if ~isfield(Opt,'buffer');     Opt.buffer = 'on'; end     % 'on', 'para', or 'off'
% Freestream BC: per-component type for u, v, w.
%   'H_DR'  = homogeneous Dirichlet  (q̂ = 0)
%   'IH_DR' = inhomogeneous Dirichlet (q̂ = bct value)
%   'H_NM'  = homogeneous Neumann    (∂q̂/∂y = 0)
% Backward compatibility: scalar 'dirichlet'/'neumann' expands to all 3.
if ~isfield(Opt,'bc_top')
    Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
elseif ischar(Opt.bc_top) || isstring(Opt.bc_top)
    switch lower(char(Opt.bc_top))
        case 'dirichlet'; Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
        case 'neumann';   Opt.bc_top = {'H_NM', 'H_NM', 'H_NM'};
        otherwise;         error('Unknown bc_top: %s', Opt.bc_top);
    end
end
if ~isfield(Opt,'lu_mode');    Opt.lu_mode = 'full'; end  % 'full', 'auto', or 'none'
if ~isfield(Opt,'lu_max_cache'); Opt.lu_max_cache = 3; end  % max cached LU for 'auto' mode
if ~isfield(Opt,'parfor')
    has_pct = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if has_pct; Opt.parfor = 'on'; else; Opt.parfor = 'off'; end
end
if ~isfield(Opt,'gpu')
    try; has_gpu = canUseGPU(); catch; has_gpu = false; end
    if has_gpu; Opt.gpu = 'on'; else; Opt.gpu = 'off'; end
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
Mode.Aactive = zeros(size(Mode.L));       % Active mode registrer

%% Initialisation of the matrices

% Initiate identity and zero matrix of order ny
Stab.I      = eye(Grid.ny);
Stab.Z      = zeros(Grid.ny);

% Allocate LHS matrix blocks 1-4 and separately for MFD.
% Sparsity estimates assume 4th-order FD; revisit if FD order changes.
Stab.ML1MFD = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,4*(3*Grid.ny^2+4*Grid.ny)+7*Grid.ny^2+4*Grid.ny);
Stab.ML1    = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,4*(3*Grid.ny^2+4*Grid.ny)+7*Grid.ny^2+4*Grid.ny);
Stab.ML2    = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,3*Grid.ny);
Stab.ML3    = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,5*Grid.ny);
Stab.ML4    = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,3*Grid.ny);

StabRes.R=zeros(4*Grid.nx*Grid.ny,1);  % Initiate RHS vector
Stab.f=zeros(Mode.nf,4*Grid.ny*Grid.nx);

% Streamwise FD: dxi is the scalar first-spacing; D_xi is the sparse
% 4th-order FD matrix valid on non-uniform ξ grids.
Stab.dxi   = Grid.x(1,2) - Grid.x(1,1);
Grid.D_xi  = fd1_matrix_nonuniform(Grid.x(:));

%% Initialization of nonlinear convergence loop

StabRes.Aold = zeros(Mode.nf,Grid.nx);            % Amplitudes of previous iteration
StabRes.q = zeros(Mode.nf,4*Grid.ny*Grid.nx);     % Solution vector

Opt.dal = ones(Mode.nf,1);   % Convergence measure
Opt.dalsave(:,1) = Opt.dal; % Convergence measure register
Opt.mmry = 0;           % memory variable so that LHS is created only once
Opt.k1 = 0;             % Iteration counter
Opt.k2 = 0;             % Converged iteration counter
Opt.AF = 1;             % Amplitude factor = 1 initially. damped later to values <1
Opt.AVAL=0; 

end

