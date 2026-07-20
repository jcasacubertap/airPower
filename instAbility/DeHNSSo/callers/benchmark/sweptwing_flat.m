%% DeHNSSo — Stationary CFI in a swept-wing flat-plate BL
%  Base flow: bf_sweptwing_flat.mat (see gridgen/benchmark/sweptwing_flat.m)
%  Nonlinear N=5 = Fig. 9; needs Opt.anderson=5 (set below).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))

input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input.mat'));   % provides StabGrid

% =========================================================================
%% Spectral content
% =========================================================================
Stab.M       = 0;                               % omega (temporal) modes
Stab.N       = 5;                               % beta (spanwise) modes
Stab.omega_0 = 0;                               % stationary
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;   % lambda = 7.5 mm

% =========================================================================
%% Inflow conditions
% =========================================================================
Stab.IC       = 'ILST';      % 'ILST' | 'LOAD' (needs Stab.ICfile)
Stab.phaseRef = 'umax';      % 'pwall' (TS) | 'umax' (CFI)
Opt.linear    = 'on';        % 'on' linear (1 iter) | 'off' nonlinear
Stab.A0_fund  = 0.7524e-4/2;    % fundamental amplitude

% =========================================================================
%% Boundary conditions
% =========================================================================
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % per component: H_DR | IH_DR | H_NM

% =========================================================================
%% Buffer — outflow treatment
% =========================================================================
Opt.buffer = 'on';           % 'on' damping | 'para' parabolisation | 'off'
Opt.xb     = 85;             % buffer start [% of domain]

% =========================================================================
%% Solver options
% =========================================================================
Opt.plot        = 'on';      % 'on' | 'off'
Opt.plot_metric = 'umax';    % 'umax' | 'urms' | 'energy'

% =========================================================================
%% Performance
% =========================================================================
Opt.lu_mode = 'none';        % 'full' | 'auto' | 'lapack_band' | 'none' (docs/user_guide.md §8)
Opt.parfor  = 'off';
Opt.gpu     = 'off';

% =========================================================================
%% Nonlinear solver
% =========================================================================
Opt.solver   = 'picard';     % 'picard' | 'newton'
Opt.anderson = 5;            % required for nonlinear (no-op for linear)
Opt.newton_takeover = 'off';

% =========================================================================
%% Run
% =========================================================================
Opt.rerun = 'off';           % 'on' reuses LHS from a previous run

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
ref_file  = fullfile(rootdir, 'literature', 'StabRes_SweptWing_flat_NPSE.mat');
ref_label = 'NPSE';

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
save(fullfile(rootdir, 'sweptwing_flat_output.mat'), 'StabRes', 'StabGrid');
fprintf('  saved -> %s\n', fullfile(rootdir, 'sweptwing_flat_output.mat'));