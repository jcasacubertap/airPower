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
% DeHNSSo — Blasius boundary layer (TS waves)
%  Base flow: DeHNSSo_input_ts_blasius.mat (run gridgen/benchmark/ts_blasius.m first).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_ts_blasius.mat'));   % StabGrid

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 5;                 % temporal modes
Stab.N = 0;                 % spanwise modes (0 = 2D)
Stab.omega_0 = 0.034576;    % fundamental frequency
Stab.beta_0  = 0;           % fundamental spanwise wavenumber

%% Inflow
Stab.IC       = 'ILST';     % 'ILST' (local stability) or 'LOAD'
Stab.phaseRef = 'pwall';    % 'pwall' (TS) or 'umax' (CFI)
Opt.linear    = 'off';      % 'on' linear (1 iter, NLT off) | 'off' nonlinear
Stab.A0_fund  = 0.0025*sqrt(2);   % fundamental amplitude (physical 2|u|)

%% Boundary conditions (freestream type per component u,v,w; wall is Dirichlet)
% 'H_DR' Dirichlet | 'IH_DR' inhomogeneous | 'H_NM' Neumann | 'free' (v only)
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};

%% Buffer — 'on' amplitude damping | 'para' parabolisation | 'off'
Opt.buffer = 'on';
Opt.xb     = 70;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';   % 'umax' | 'urms' | 'energy'
Opt.lu_mode     = 'full';   % 'full' | 'auto' | 'lapack_band' | 'none' (DeHNSSo_guide §8)
Opt.parfor      = 'off';
Opt.solver      = 'picard'; % 'picard' | 'newton'

%% Run
Opt.rerun       = 'off';    % 'on' reuses LHS from a previous run
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_ts_blasius.mat');   % [] to skip

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing
close all
ref_file  = fullfile(rootdir, 'literature', 'StabRes_Blasius_NPSE.mat');
ref_label = 'NPSE';
if strcmpi(Opt.plot, 'on')
    if isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label), [], Opt.plot_metric);
    else
        plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    end
    plot_stability(StabGrid, StabRes, Opt);
end
