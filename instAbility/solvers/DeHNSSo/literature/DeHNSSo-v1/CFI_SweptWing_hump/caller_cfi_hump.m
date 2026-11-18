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
% DeHNSSo-v1 test case, reproduced with DeHNSSo-v2:
%   crossflow instability over a hump in a swept-wing boundary layer
%
% Original v1 caller: Callers/CFI over a hump in a swept-wing BL/Example_Caller_CFI_hump_in_SweptWing_BL_paper.m
% Reference data in this folder: StabRes_SweptWing_Hump_AHLNS.mat (AHLNS)
%
% Run gridgen_cfi_hump.m (in this folder) first to build the stability grid.
% This is the same calculation as callers/benchmark/cfi_hump.m; it is repeated here so
% the v1 cases are self-contained, each beside the data it is compared against.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_v1_cfi_hump.mat'));   % StabGrid

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 0;                 % temporal modes
Stab.N = 1;                 % spanwise modes
Stab.omega_0 = 0;           % stationary
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;   % lambda = 7.5 mm

%% Inflow
Stab.IC       = 'ILST';     % 'ILST' or 'LOAD'
Stab.phaseRef = 'umax';     % 'pwall' (TS) | 'umax' (CFI)
Opt.linear    = 'on';       % 'on' linear (1 iter) | 'off' nonlinear
Stab.A0_fund  = 3.5e-3;     % fundamental amplitude (physical 2|u|)

%% Boundary conditions (freestream type per component u,v,w; wall Dirichlet)
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % 'H_DR' | 'IH_DR' | 'H_NM' | 'free' (v only)

%% Buffer — 'on' damping | 'para' parabolisation | 'off'
Opt.buffer = 'para';
Opt.xb     = 85;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';   % 'umax' | 'urms' | 'energy'
Opt.lu_mode     = 'auto';   % 'full' | 'auto' | 'lapack_band' | 'none' (user_guide §8)
                            % 'full' caches every factor and exceeds 32 GB on this grid
Opt.parfor      = 'off';
Opt.solver      = 'picard'; % 'picard' | 'newton'

%% Run
Opt.rerun       = 'off';    % 'on' reuses LHS from a previous run
Opt.save_output = fullfile(fileparts(mfilename('fullpath')), 'output', 'StabRes_cfi_hump.mat');   % [] to skip

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing
close all
ref_file  = fullfile(fileparts(mfilename('fullpath')), 'StabRes_SweptWing_Hump_AHLNS.mat');
ref_label = 'AHLNS';
if strcmpi(Opt.plot, 'on')
    if isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label), [], Opt.plot_metric);
    else
        plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    end
    plot_stability(StabGrid, StabRes, Opt);
end
