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
% Gridgen benchmark — Blasius BL, mesh clustered at the synthetic jet
%  Same base flow as blasius.m, refined around the jet orifice. The forcing is
%  defined ONCE here (params.forcing); the DeHNSSo caller only sets amplitude/phase.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_ts_blasius.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — pass-through (clustering redistributes the source stations)
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;
params.plot      = false;

%% Wall forcing — gridgen refines the mesh for it and stores StabGrid.forcing.
% Actuator (raised cosine) at 40% of the domain, half-width 18 (bf_ts_blasius x ~ [400, 2848]).
% For a measured distribution instead: params.forcing = struct('x', x_data, 'v', v_data).
% Several actuators: give xc a vector (w scalar or matching). Mixed types: a
% cell array. The caller selects each with spec.which.
%   params.forcing = struct('xc', [1379.1, 2000, 3000], 'w', 18);
%   params.forcing = {struct('xc',1379.1,'w',18), struct('x',x_data,'v',v_data)};
params.forcing = struct('xc', 1379.1, 'w', 18, 'shape', 'cosine');
params.refine  = 10/3;           % ag = 0.3

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_wall_jet.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
