%% DeHNSSo — Stationary CFI over a hump in a swept-wing BL
%  Base flow: bf_sweptwing_hump.mat (see gridgen/benchmark/sweptwing_hump.m)

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))

% Base flow: generate first by running gridgen/benchmark/sweptwing_hump.m
input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input.mat'));   % provides StabGrid

% =========================================================================
%% Spectral content
% =========================================================================
%  total modes: (2M+1)(2N+1), including conjugates

Stab.M = 0;                 % number of omega (temporal) modes
Stab.N = 5;                 % number of beta (spanwise) modes
Stab.omega_0 = 0;           % fundamental angular frequency (stationary)
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;   % fundamental spanwise wavenumber (lambda = 7.5 mm)

% =========================================================================
%% Inflow conditions
% =========================================================================

Stab.IC       = 'ILST';     % 'ILST' (local stability) or 'LOAD' (from file)
Stab.phaseRef = 'umax';     % 'pwall' (TS) or 'umax' (CFI)
% Stab.ICfile = fullfile(input.folder, 'StabRes_previous.mat');  % required for LOAD

Opt.linear   = 'on';        % 'on' linear (NLT off, 1 iter) | 'off' nonlinear
Stab.A0_fund = 3.5e-3/2;    % fundamental disturbance amplitude

% =========================================================================
%% Boundary conditions
% =========================================================================
%
% Wall (y=0): always inhomogeneous Dirichlet.
%   Stab.bcw(k,ix,j) — wall BC value.  k=1 u, k=2 v, k=3 w.
%   Default: zero (no-slip).
%
% Freestream (y=H): per-component BC type.
%   Opt.bc_top{k} — type for component k=1 u, k=2 v, k=3 w.
%     'H_DR'  = homogeneous Dirichlet   (q = 0, default)
%     'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%     'H_NM'  = homogeneous Neumann     (dq/dy = 0)
%   Pressure: no BC at freestream (determined by the equations).

Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};

% =========================================================================
%% Buffer — outflow treatment
% =========================================================================
%   'on'   = amplitude damping (scales solution -> 0 near outflow)
%   'para' = parabolisation (removes d2/dx2 in buffer -> no upstream reflections)
%   'off'  = no buffer

Opt.buffer = 'para';
Opt.xb     = 85;            % buffer start [% of domain]

% =========================================================================
%% Solver options
% =========================================================================

Opt.plot        = 'off';     % 'on' | 'off'
Opt.plot_metric = 'umax';   % 'umax' | 'urms' | 'energy'
% Opt.Conv = 1e-4;          % final convergence criterion
% Opt.TH   = 1e-11;         % NLT activation threshold

% =========================================================================
%% Performance
% =========================================================================

% LU mode (Opt.lu_mode):
%   'full'        — Cache every LU factorisation. Fastest for repeated solves,
%                   unbounded memory.
%   'auto'        — LRU cache bounded by Opt.lu_max_cache.
%   'lapack_band' — LAPACK zgbsv/dgbsv MEX. Fastest single solve.
%                   No factor cache yet — refactors per call.
%                   Build once per machine: run src/hns/mexbuild_dgbsv.m.
%   'none'        — Backslash each time. Reference / debug.
% See docs/user_guide.md §8 for the picker table.
Opt.lu_mode = 'full';
Opt.parfor  = 'off';        % parallel assembly/solve: 'on' | 'off'
Opt.gpu     = 'off';        % GPU sparse solve: 'on' | 'off'

% =========================================================================
%% Nonlinear solver choice
% =========================================================================
% 'picard' — Picard with AF ramping + mode activation (default)
% 'newton' — JFNK with per-mode preconditioner

Opt.solver = 'picard';      % 'picard' (+ Newton fallback) | 'newton'
% Opt.newton_tol = 1e-10;

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

% Reference for overlay (set ref_file = [] to skip)
ref_file  = fullfile(rootdir, 'literature', 'StabRes_SweptWing_Hump_AHLNS.mat');
ref_label = 'AHLNS';

if strcmpi(Opt.plot, 'on')
    if ~isempty(ref_file) && isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label), [], Opt.plot_metric);
    else
        if ~isempty(ref_file)
            fprintf('  [ref overlay skipped — file not found: %s]\n', ref_file);
        end
        plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    end
    plot_stability(StabGrid, StabRes, Opt);
end

% =========================================================================
%% Save results in the DeHNSSo folder
% =========================================================================
save(fullfile(rootdir, 'sweptwing_hump_output.mat'), 'StabRes', 'StabGrid');
fprintf('  saved -> %s\n', fullfile(rootdir, 'sweptwing_hump_output.mat'));