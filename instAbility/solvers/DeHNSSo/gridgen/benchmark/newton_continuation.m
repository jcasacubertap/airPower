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
% Gridgen benchmark — Blasius BL for the Newton-continuation demo
%  The cheapest grid on which a Newton pseudo-arclength run is worth doing:
%  the ts_blasius base flow resampled to 400 x 40, equidistant in x over the
%  full domain (R = 400 to 1067), Chebyshev in eta. Laptop-sized: the demo in
%  callers/benchmark/newton_continuation.m needs ~8.6 GB and a few minutes.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_ts_blasius.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — equidistant x, Chebyshev eta; H and y_i take gridgen's defaults
params.n_eta_new = 40;
params.n_xi_new  = 400;
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_newton_continuation.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
