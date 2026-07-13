%% DeHNSSo — Compressible Blasius boundary layer (TS waves)
%  Base flow: BF_compressible.mat (see baseflow/caller/caller_CBL.m)
%  Required: gridgen must produce DeHNSSo_input_compressible.mat that carries
%  rho, T and thermo metadata (Ma, Pr, gamma, S*, mu, dmuT, d2muT2).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))

% Base flow: generate first by running gridgen + baseflow/caller/caller_CBL.m
input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input_compressible.mat'));   % provides StabGrid

% =========================================================================
%% Compressible mode
% =========================================================================

Opt.compressible = 'on';        % activate (u,v,w,rho,T) path

% =========================================================================
%% Spectral content
% =========================================================================
%  total modes: (2M+1)(2N+1), including conjugates

Stab.M = 1;                     % number of omega (temporal) modes
Stab.N = 0;                     % number of beta (spanwise) modes
Stab.F_0     = 40e-6;           % reduced frequency
Stab.omega_0 = Stab.F_0 * StabGrid.Re;
Stab.beta_0  = 0;               % fundamental spanwise wavenumber (0 = 2D)

% =========================================================================
%% Inflow conditions
% =========================================================================

Stab.IC         = 'CLST';       % required for compressible (CLST or LOAD)
Stab.phaseRef   = 'pwall';      % 'pwall' (TS) or 'umax' (CFI)
Opt.mode_select = 'slow';       % 'slow' / 'fast' / 'auto' — Mack F/S pick at Ma>1
% Stab.ICfile = fullfile(input.folder, 'StabRes_previous.mat');  % required for LOAD

Opt.linear   = 'on';            % 'on' linear (NLT off, 1 iter) | 'off' nonlinear
Stab.A0_fund = 1e-6;            % fundamental disturbance amplitude

% =========================================================================
%% Boundary conditions
% =========================================================================
%
% Wall (y=0): always inhomogeneous Dirichlet.
%   Stab.bcw(k,ix,j) — wall BC value.  k=1 u, k=2 v, k=3 w, k=4 T.
%   Default: zero (no-slip + isothermal).
%
% Freestream (y=H): per-component BC type — 4 entries (u,v,w,T).
%   Opt.bc_top{k} — type for component k=1 u, k=2 v, k=3 w, k=4 T.
%     'H_DR'  = homogeneous Dirichlet   (q = 0, default)
%     'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%     'H_NM'  = homogeneous Neumann     (dq/dy = 0; absorbs acoustics at Ma>1)
%   Density: no BC at freestream (determined by the equations).

Opt.bc_top          = {'H_DR', 'H_DR', 'H_DR', 'H_DR'};
Opt.bc_wall_thermal = 'isothermal';   % T'=0 at wall (also valid for adiabatic BF)

% =========================================================================
%% Buffer — outflow treatment
% =========================================================================
%   'on'   = amplitude damping (scales solution -> 0 near outflow)
%   'para' = parabolisation (removes d2/dx2 in buffer -> no upstream reflections)
%   'off'  = no buffer
% At Ma > 0, parabolisation is required: 'on' / 'off' reflect outgoing acoustics.

Opt.buffer = 'para';
Opt.xb     = 90;                % buffer start [% of domain]

% =========================================================================
%% Solver options
% =========================================================================

Opt.plot        = 'on';         % 'on' | 'off'
Opt.plot_metric = 'umax';       % 'umax' | 'urms' | 'energy'
% Opt.Conv = 1e-4;              % final convergence criterion
% Opt.TH   = 1e-11;             % NLT activation threshold

% =========================================================================
%% Performance
% =========================================================================

% LU mode (Opt.lu_mode):
%   'full'        — Cache every LU factorisation. OOM at nx=800 (compressible).
%   'auto'        — LRU cache bounded by Opt.lu_max_cache. OOM at nx=800.
%   'lapack_band' — LAPACK zgbsv/dgbsv MEX. Fastest single solve (~3x vs seq).
%                   No factor cache yet — refactors per call.
%                   Build once per machine: run src/hns/mexbuild_dgbsv.m.
%   'seq'         — symrcm + sparse LU, one mode group at a time. ~210 s on M=0.5.
%                   Robust fallback if MEX unavailable.
%   'none'        — Backslash each time. Reference / debug.
% See docs/user_guide.md §8 for the picker table.
Opt.lu_mode = 'lapack_band';
Opt.parfor  = 'off';            % parallel assembly/solve: 'on' | 'off'
Opt.gpu     = 'off';            % GPU sparse solve: 'on' | 'off'

% =========================================================================
%% Nonlinear solver choice
% =========================================================================
% Compressible only supports Picard. Newton/JFNK is not implemented for the
% compressible path; main.m errors out if Opt.solver='newton' & compressible='on'.

Opt.solver = 'picard';

% =========================================================================
%% Run
% =========================================================================
% Opt.rerun = 'on' reuses LHS from a previous run (skips matrix assembly).

Opt.rerun = 'off';

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
t_total = toc(tic_total);
w = whos;
fprintf('\n========================================\n');
fprintf('  Total execution time: %.2f s\n', t_total);
fprintf('  Workspace memory:     %.1f MB\n', sum([w.bytes])/1e6);
fprintf('========================================\n');

% =========================================================================
%% Post-processing
% =========================================================================
close all

% Zanus/Bertolotti 2017 PSE CPG reference overlay (Ma=0.5, F=40e-6).
% Set ref_csv = [] to skip (different (Ma, F) case).
ref_csv = fullfile(rootdir, 'literature', 'Zanus_PSE_CPG_M05_F40.csv');

if strcmpi(Opt.plot, 'on')
    plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    plot_stability(StabGrid, StabRes, Opt);

    % HNS vs Zanus/Bertolotti PSE comparison
    if ~isempty(ref_csv) && isfile(ref_csv)
        % HNS alpha is in global non-dim (fixed lref at inflow); Zanus/Bertolotti is in
        % local non-dim (per-station lref). Rescale HNS -> local by (R/R_0).
        R_sweep   = sqrt(StabGrid.xun * StabGrid.xun(1));
        alpha_hns = squeeze(StabRes.alpha(end, :)) .* (R_sweep ./ R_sweep(1));
        ZN        = readmatrix(ref_csv);

        figure('Position', [100 100 1000 400]);
        subplot(1, 2, 1);
        plot(R_sweep, -imag(alpha_hns), 'r--', 'LineWidth', 1.5, 'DisplayName', 'HNS'); hold on
        plot(ZN(:,1),  ZN(:,2),         'ms',  'MarkerSize', 5,  'DisplayName', 'Zanus/Bertolotti 2017 PSE');
        xlabel('R = \surd Re_x'); ylabel('-\alpha_i');
        grid on; legend('Location', 'best');
        title(sprintf('Spatial growth rate, Ma=%.2f, F=%.0fe-6', ...
                      StabGrid.Ma, Stab.F_0*1e6));

        subplot(1, 2, 2);
        plot(R_sweep, real(alpha_hns), 'r--', 'LineWidth', 1.5, 'DisplayName', 'HNS');
        xlabel('R = \surd Re_x'); ylabel('\alpha_r');
        grid on; legend('Location', 'best');
        title('Streamwise wavenumber');
    elseif ~isempty(ref_csv)
        fprintf('  [Zanus/Bertolotti overlay skipped — file not found: %s]\n', ref_csv);
    end
end
