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
% Stability grid for the DeHNSSo-v1 TS Blasius test case.
%
% This is the gridgen benchmark grid with one change: the median collocation
% height is halved to y_i = H/20, which is what the v1 caller
%   Callers/TS-Waves in a Blasius BL/Example_Caller_Blasius_BL_TS.m
% prescribed via Grid.y_i = Grid.H/20. The shared gridgen/benchmark/ts_blasius.m
% leaves params.y_i empty and so falls back to the gridgen default H/10, giving
% half the near-wall resolution v1 used.
%
% Everything else matches v1: 800 equidistant streamwise stations and 40
% Chebyshev collocation points, taken straight from the base flow.
%
% Writes input/DeHNSSo_input_v1_blasius.mat, a separate file from the benchmark
% grid, so the benchmark and verification cases that rely on
% DeHNSSo_input_ts_blasius.mat are left untouched.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(rootdir, 'gridgen'), genpath(fullfile(rootdir, 'gridgen', 'src')))

%% Input
input.folder     = fullfile(rootdir, 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_ts_blasius.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — pass-through (empty counts keep the source 800 x 40 grid)
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;
params.plot      = false;

%% Wall-normal clustering — v1 used y_i = H/20
% resample_cartesian takes H as 0.99 x the base-flow height when params.H is
% empty; reproduce that here so y_i follows the same rule v1 applied.
bf  = load(fullfile(input.folder, input.filename));
fn  = fieldnames(bf);
BF  = bf.(fn{1});
H   = abs(max(BF.Y(:)) - min(BF.Y(:))) * 0.99;
params.H   = H;
params.y_i = H / 20;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, 'input');
output.filename = 'DeHNSSo_input_v1_blasius.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
