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

% (1) Problem geometry — selects the frame offset AND the gridgen caller name:
%       'flat'   -> flat swept plate    (gridgen/benchmark/sweptwing_flat.m)
%       'curved' -> M3J curved airfoil  (gridgen/benchmark/m3j.m)
%     The flat plate omits wall curvature; the curved airfoil includes it (and
%     buildStabGrid auto-computes the arc-length frame offset from the geometry).
%     Run the matching gridgen caller for your problem, then paste the resulting
%     DeHNSSo_input.mat into AmplitudeFinisher/io/ (keep the DeHNSSo name, for
%     consistency). The tool trusts that DeHNSSo_input.mat corresponds to in.problem.
in.problem = 'curved';                 % 'flat' | 'curved'
in.stabGridFile = fullfile(airPowerRoot, 'instAbility', 'AmplitudeFinisher', ...
                           'io', 'DeHNSSo_input.mat');
switch lower(in.problem)
    case 'flat';   in.caller = 'sweptwing_flat';
    case 'curved'; in.caller = 'm3j';
    otherwise
        error('inputsAmplitudeFinisher:problem', ...
              'in.problem must be ''flat'' or ''curved'' (got ''%s'').', in.problem);
end

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

% (6) PIV target. Place the experimental PIV .mat in AmplitudeFinisher/io/ and
%     give its name here. It holds the struct in.validation.structVar ('output')
%     with the w'-RMS profiles and the streamwise stations (xc).
in.pivName = '7d5mm_h0000.mat';        % file in AmplitudeFinisher/io/
in.pivFile = fullfile(airPowerRoot, 'instAbility', 'AmplitudeFinisher', 'io', in.pivName);

% (7) Machine profile — sets the solver memory/speed strategy in one place:
%     'small' : <=16 GB RAM. lu_mode='none' (re-factorises each iter; low RAM,
%               ~100 s/iter). Newton takeover off.
%     'large' : big RAM. lu_mode='full' (cache all mode factorisations; ~10x
%               faster iters but tens of GB), Newton takeover ON (needed for
%               deep-saturation convergence downstream). Use this for high N
%               and multi-station runs.
in.machine = 'small';

% (8) Linear-check anchor. EVERY run produces a linear-check figure
%     (AmplitudeFinisher_linearCheck.png) as a by-product: one linear HNS whose
%     w'-RMS is anchored to PIV at this x/c. The figure title reports the implied
%     A0_fund, and that SAME value seeds the nonlinear search (centre of
%     in.match.sweep). Anchor upstream (near-linear) for the cleanest estimate.
in.linearCheck.xcAnchor = 13;      % [% chord] station to anchor the linear amplitude to PIV

% =====================================================================
%  SECONDARY INPUTS  ── defaults are fine; change only if needed
% =====================================================================

% (Base flow / gridgen inputs removed — the StabGrid is loaded directly from
%  AmplitudeFinisher/io/DeHNSSo_input.mat; see MAIN input (1). Reference scales
%  Uref/lref/Re come baked into that StabGrid.)

%% -- Stability / disturbance (fixed for stationary crossflow) --
in.Stab.M        = 0;                  % temporal (omega) modes
in.Stab.omega_0  = 0;                  % stationary
in.Stab.IC       = 'ILST';             % inflow: local stability eigenmode
in.Stab.phaseRef = 'umax';             % CFI phase reference
in.Stab.A0_fund  = 5.0e-5;           % only a starting guess; the search sets it

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
% Streamwise frame offset: hns.x = StabGrid.x*lref + xOffset aligns the HNS to
% the PIV arc-length S from the LE. Set automatically by buildStabGrid per
% problem -- flat: xOffset = xInlet (S = xInlet + x_OF, as in dfpvalidation.jl);
% curved: xOffset = arc-length(LE -> grid origin), computed from the airfoil
% geometry (constant to machine precision, since the grid wall-arc matches BL.S).
in.match.xOffset   = [];               % [m] filled by buildStabGrid (see in.geom below)
% Geometry for the frame offset. Flat uses xInlet; curved inverts the airfoil
% surface x_phys(x/c) (as in PreProcessing airfoil_surface) to map each grid
% column to x/c, then to S via the reference BL table.
in.geom.xInlet     = ap.xInlet;        % [m]   flat: inlet distance from LE
in.geom.airfoilFile= 'M3J.dat';        % curved: file in PreProcessing/io/airfoilGeometryData/
in.geom.chord      = 0.900;            % [m]   curved: chord (matches inputs.jl TTCP)
in.geom.alphaDeg   = -3.017;           % [deg] curved: AoA   (matches inputs.jl TTCP)
in.match.bufferFrac= 0.80;             % matched stations kept within this fraction of the run
in.match.stationIdx= [];               % filled from xcWindow in AmplitudeFinisher.m
in.match.xStations = [];

%% -- PIV data mapping --  (PIV file itself is MAIN input (6): in.pivFile)
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
