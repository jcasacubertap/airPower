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

%% Description
% Gridgen benchmark — swept-wing CFI over a forward-facing step

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_cfi_step.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
% xi_range is a pair of fractions of the streamwise span, not coordinates.
BF       = load(fullfile(input.folder, input.filename)).BF;
x_window = [129.5, 215.1];
params.n_eta_new = 80;
params.n_xi_new  = 1500;
params.y_i       = 1.0;
params.H         = [];
params.xi_range  = (x_window - BF.X(1)) / (BF.X(end) - BF.X(1));
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_cfi_step.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
