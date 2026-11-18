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
% Gridgen benchmark — Compressible Blasius (FSC) base flow
%  Cartesian pass-through from baseflow/caller/caller_CBL.m. build_stab_grid
%  auto-detects rho/T, runs Sutherland, and forwards the thermo metadata.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_ts_compressible.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — pass-through
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators (FD_xi_order_2 = 4: BF d2/dx2 viscous terms in LHS_compressible)
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_ts_compressible.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
