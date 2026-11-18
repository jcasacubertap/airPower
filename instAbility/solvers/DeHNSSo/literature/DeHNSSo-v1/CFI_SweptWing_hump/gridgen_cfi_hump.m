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
% Stability grid for the DeHNSSo-v1 swept-wing hump test case.
%
% Wall-normal recipe of the v1 caller
%   Callers/CFI over a hump in a swept-wing BL/Example_Caller_CFI_hump_in_SweptWing_BL_paper.m
% Grid.nx = 1000, Grid.ny = 50 Chebyshev, Grid.y_i = 4.45, Grid.H = max(BF.Y(:,1)).
% Streamwise clustering about the hump as in gridgen/benchmark/cfi_hump.m,
% whose spread was set to match the heritage grid. The benchmark itself uses
% 1200 x 80 finite differences; this file keeps v1's numbers and writes a
% separate input so the benchmark is left untouched.
%
% Writes input/DeHNSSo_input_v1_cfi_hump.mat.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(rootdir, 'gridgen'), genpath(fullfile(rootdir, 'gridgen', 'src')))

%% Input
input.folder     = fullfile(rootdir, 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_cfi_hump.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — v1 values
bf  = load(fullfile(input.folder, input.filename));
fn  = fieldnames(bf);
BF  = bf.(fn{1});
params.n_eta_new = 50;
params.n_xi_new  = 1000;
params.H         = max(BF.Y(:,1));   % v1: Grid.H = max(BF.Y(:,1))
params.y_i       = 4.45;             % v1: Grid.y_i = 4.45
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;
params.plot      = false;

%% Streamwise clustering about the hump (as gridgen/benchmark/cfi_hump.m)
x_m = 0.1837410879 / 2.1394e-4;      % hump centre ~859
params.mug    = x_m;
params.spread = 55.137096;
params.refine = 2;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, 'input');
output.filename = 'DeHNSSo_input_v1_cfi_hump.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
