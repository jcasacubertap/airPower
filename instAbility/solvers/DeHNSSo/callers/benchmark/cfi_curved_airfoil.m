%% DeHNSSo — Stationary CFI on the M3J curved-wall BL
%  Base flow: DeHNSSo_input_cfi_curved_airfoil.mat (run gridgen/benchmark/cfi_curved_airfoil.m first).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_cfi_curved_airfoil.mat'));   % StabGrid

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 0;                 % temporal modes
Stab.N = 1;                 % spanwise modes
Stab.omega_0 = 0;           % stationary
Stab.beta_0  = 0.40;        % spanwise wavenumber (most-amplified, lambda ~ 7 mm)

%% Inflow
Stab.IC       = 'ILST';     % 'ILST' or 'LOAD'
Stab.phaseRef = 'umax';     % 'pwall' (TS) | 'umax' (CFI)
Opt.linear    = 'on';       % 'on' linear (1 iter) | 'off' nonlinear
Stab.A0_fund  = 3.5e-3;     % fundamental amplitude (physical 2|u|)

%% Boundary conditions (freestream type per component u,v,w; wall Dirichlet)
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % 'H_DR' | 'IH_DR' | 'H_NM' | 'free' (v only)

%% Buffer — 'on' damping | 'para' parabolisation | 'off'
Opt.buffer = 'on';
Opt.xb     = 85;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';   % 'umax' | 'urms' | 'energy'
Opt.lu_mode     = 'full';   % 'full' | 'auto' | 'lapack_band' | 'none' (DeHNSSo_guide §8)
Opt.parfor      = 'off';
Opt.solver      = 'picard'; % 'picard' | 'newton'

%% Run
Opt.rerun       = 'off';    % 'on' reuses LHS from a previous run
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_cfi_curved_airfoil.mat');   % [] to skip

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing (no published M3J reference yet)
close all
if strcmpi(Opt.plot, 'on')
    plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    plot_stability(StabGrid, StabRes, Opt);
end
