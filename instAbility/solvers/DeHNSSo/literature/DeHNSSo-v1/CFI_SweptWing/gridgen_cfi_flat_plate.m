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
% Stability grid for the DeHNSSo-v1 swept-wing crossflow test case.
%
% Same recipe as the v1 caller
%   Callers/CFI in a swept-wing BL/Example_Caller_CFI_in_SweptWing_BL.m
% which used Grid.nx = 1200, Grid.ny = 40 Chebyshev, Grid.H = 0.99*max(BF.Y) and
% Grid.y_i = Grid.H/10. Those are also gridgen's defaults for H and y_i, so this
% is the gridgen/benchmark/cfi_flat_plate.m grid written to a separate file, so
% the v1 case is self-contained and the benchmark input is left untouched.
%
% Writes input/DeHNSSo_input_v1_cfi_flat_plate.mat.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(rootdir, 'gridgen'), genpath(fullfile(rootdir, 'gridgen', 'src')))

%% Input
input.folder     = fullfile(rootdir, 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_cfi_flat_plate.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — 1200 equidistant stations, 40 Chebyshev points (v1 values)
params.n_eta_new = 40;
params.n_xi_new  = 1200;
params.y_i       = [];        % default H/10 = v1's Grid.y_i
params.H         = [];        % default 0.99*max(BF.Y) = v1's Grid.H
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, 'input');
output.filename = 'DeHNSSo_input_v1_cfi_flat_plate.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
