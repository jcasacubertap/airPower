function in = inputsAmplitudeFinisher(airPowerRoot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    AmplitudeFinisher — Inputs                         %
%                                                                       %
%   Find the disturbance amplitude (Stab.A0_fund) that makes the        %
%   DeHNSSo Harmonic Navier-Stokes solution's spanwise-velocity         %
%   perturbation match the experimental PIV, for an OpenFOAM            %
%   flat-plate base flow.                                               %
%                                                                       %
%   Set the MAIN INPUTS below for a new case; the SECONDARY INPUTS      %
%   have sensible defaults and rarely need changing.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in = struct();
in.airPowerRoot   = airPowerRoot;
in.airPowerInputs = fullfile(airPowerRoot, 'inputs.jl');
in.dehnssoRoot    = fullfile(airPowerRoot, 'instAbility', 'DeHNSSo');
ap = local_airpower_scalars(in.airPowerInputs);   % DFP-block scalars from inputs.jl

% =====================================================================
%  MAIN INPUTS  ── set these for a new base flow
% =====================================================================

% (1) Base flow: OpenFOAM flat-plate mid-plane CSV (columns x,y,z,u,v,w,p,omz),
%     which DeHNSSo gridgen resamples onto the stability grid.
%     This CSV is NOT hand-placed here: it is CONVERTED from the reference
%     OpenFOAM output of the DirectFlatPlate module, i.e. the folder
%       airPower/PreProcessing/modules/directFlatPlateModule/postProcessing
%     The path below is just where that conversion writes/caches the result.
in.gridgen.baseFlowCsv = fullfile(in.dehnssoRoot, 'baseflow', 'output', ...
                                  'benchmark', 'bf_sweptwing_flat.csv');

% (2) Spanwise fundamental wavelength of the crossflow instability.
in.Stab.lambda_z = 7.5e-3;             % [m]   (case '7d5mm')

% (3) PIV chordwise window to match, [x/c_min x/c_max] in percent.
in.match.xcWindow = [17 30];           % [% chord]

% (4) Number of spanwise Fourier modes (N).
in.Stab.N = 5;

% (5) Amplitude search multipliers. The linear seed only estimates A0; the
%     true (nonlinear) amplitude differs because saturation shifts the w' peak.
%     So the tool runs one nonlinear solve at each A0 = seed .* these factors,
%     giving a few (A0, peak-ratio) points that BRACKET peak-ratio = 1; the
%     matched A0 is then interpolated from them. Keep values straddling 1.0
%     (need >=2 points, ideally one ratio<1 and one >1) so the root is bracketed.
in.match.sweep = [0.8 1.0 1.25];

% (6) PIV case. By default the generation/case are read from airPower/inputs.jl
%     (VAL block). Override here with integers to force a specific case.
in.VAL.Gen  = [];                      % [] -> read from inputs.jl
in.VAL.Case = [];

% (7) Machine profile — sets the solver memory/speed strategy in one place:
%     'small' : <=16 GB RAM. lu_mode='none' (re-factorises each iter; low RAM,
%               ~100 s/iter). Newton takeover off.
%     'large' : big RAM. lu_mode='full' (cache all mode factorisations; ~10x
%               faster iters but tens of GB), Newton takeover ON (needed for
%               deep-saturation convergence downstream). Use this for high N
%               and multi-station runs.
in.machine = 'small';

% =====================================================================
%  SECONDARY INPUTS  ── defaults are fine; change only if needed
% =====================================================================

%% -- Reference scales (non-dimensionalization), auto from inputs.jl DFP --
in.Ref.U0 = ap.Uinf;                                   % reference velocity [m/s]
in.Ref.nu = ap.nu;                                     % kinematic viscosity [m^2/s]
in.Ref.L0 = sqrt(in.Ref.nu * ap.xInlet / in.Ref.U0);   % Blasius length at inlet [m]

%% -- gridgen (base flow CSV -> StabGrid) --
in.gridgen.outMat = fullfile(in.dehnssoRoot, 'input', 'DeHNSSo_input.mat');
in.gridgen.mode   = 'auto';            % 'auto' | 'always' | 'never'
in.gridgen.params = struct( ...
    'n_eta_new', 60, 'n_xi_new', 1200, ...
    'Uref', ap.Uinf, ...
    'lref', in.Ref.L0, ...
    'Re',   in.Ref.U0 * in.Ref.L0 / in.Ref.nu, ...
    'y_i', [], 'H', [], 'xi_range', [], ...
    'xi_trim_inflow', [], 'xi_trim_outflow', [], ...
    'rescale', true, 'plot', false, ...
    'FD_xi_order_1', 4, 'FD_xi_order_2', 2, ...
    'FD_eta_method', 'cheb', 'FD_eta_order', 4);

%% -- Stability / disturbance (fixed for stationary crossflow) --
in.Stab.M        = 0;                  % temporal (omega) modes
in.Stab.omega_0  = 0;                  % stationary
in.Stab.IC       = 'ILST';             % inflow: local stability eigenmode
in.Stab.phaseRef = 'umax';             % CFI phase reference
in.Stab.A0_fund  = 1.0e-4/2;           % only a starting guess; the search sets it

%% -- Solver options (DeHNSSo main) --
in.Opt.bc_top          = {'H_DR','H_DR','H_DR'};
in.Opt.buffer          = 'on';
in.Opt.xb              = 85;           % outflow buffer start [% of domain]
in.Opt.plot            = 'off';
in.Opt.plot_metric     = 'umax';
in.Opt.parfor          = 'off';
in.Opt.gpu             = 'off';
in.Opt.solver          = 'picard';
in.Opt.anderson        = 5;
% lu_mode + Newton takeover set from the machine profile (MAIN input #7):
switch lower(in.machine)
    case 'large'
        in.Opt.lu_mode         = 'full';   % cache factorisations (fast; big RAM)
        in.Opt.newton_takeover = 'on';     % deep-saturation convergence
    otherwise   % 'small'
        in.Opt.lu_mode         = 'none';   % re-factorise each iter (low RAM)
        in.Opt.newton_takeover = 'off';
end

%% -- Matching details --
% extractWprimeHNS builds the mode amplitude as |shape|*(A/2) (the complex Fourier
% coefficient |a|), then multiplies by rmsFactor to get the spanwise (z) RMS.
% DeHNSSo's A = 2|u_max| = max|u'| (the PEAK), so the physical field is
% u'(z) = 2|a|cos(bz), peak = A*shape, and z-RMS = peak/sqrt(2) = A*shape/sqrt(2).
% From the base |a| = A/2, that RMS is (A/2)*sqrt(2). Hence rmsFactor = sqrt(2)
% (NOT 1/sqrt(2): using 1/sqrt(2) double-counts the halving in A/2 and gives half
% the true RMS, which over-scales the matched A0 by 2x).
in.match.rmsFactor = sqrt(2);          % |a| (=A/2) -> spanwise RMS = sqrt(2)*|a| = A/sqrt(2)
% Streamwise frame: the OpenFOAM/StabGrid x is measured from the INLET, while
% the PIV x/c maps to arc-length S from the LE. They differ by xInlet, exactly
% as in airPower's DFP validation (dfpvalidation.jl: x_DFP = S - xInlet, i.e.
% S = xInlet + x_OF). So HNS npse.x (S-frame) = StabGrid.x*lref + xInlet.
in.match.xOffset   = ap.xInlet;        % [m] = xInlet (align HNS x_OF -> PIV arc-length S)
in.match.bufferFrac= 0.80;             % matched stations kept within this fraction of the run
in.match.stationIdx= [];               % filled from xcWindow in AmplitudeFinisher.m
in.match.xStations = [];

%% -- PIV data mapping --
in.validation.rootDir     = fullfile(airPowerRoot, 'PreProcessing', 'io', 'Validation');
in.validation.matFile     = '';        % '' -> first .mat in the Case dir
in.validation.structVar   = 'output';
in.validation.xField      = 'xc';
in.validation.wFields     = {'w_pert_m_prof_rms_01','w_pert_m_prof_rms_02','w_pert_m_prof_rms_03'};
in.validation.yFields     = {'y_prof_rms_01','y_prof_rms_02','y_prof_rms_03'};
in.validation.totField    = 'w_m_tot';                 % raw plane -> w'-RMS target
in.validation.lengthScale = 1e-3;                      % PIV y [mm]  -> m
in.validation.velScale    = 1.0;                       % PIV w' [m/s]
in.validation.refFlowData = 'M3J_Case_Q_25_AOA_3.017_Re_2171263.mat';

end

% ---------------------------------------------------------------------
function ap = local_airpower_scalars(inputsJlPath)
% Parse the DFP-block scalars needed for the reference scales.
txt = fileread(inputsJlPath);
ap.Uinf   = local_scan_num(txt, 'Uinf');
ap.xInlet = local_scan_num(txt, 'xInlet');
ap.nu     = local_scan_num(txt, 'freeStreamViscosity');
end

function v = local_scan_num(txt, key)
tok = regexp(txt, ['(?m)^\s*' key '\s*=\s*([-\d.eE+]+)'], 'tokens', 'once');
if isempty(tok)
    error('inputsAmplitudeFinisher:scan', 'Could not parse `%s` from inputs.jl', key);
end
v = str2double(tok{1});
end
