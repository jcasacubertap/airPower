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
% Gridgen benchmark — Blasius boundary layer with a Helmholtz-resonator neck (Sec. 6.4)
%  Michelis et al. (2023) conditions: U0 = 20.9 m/s, resonance 295 Hz, a 1 mm
%  neck opening at x0 = 0.55 m. The 1 mm mouth has to carry enough stations for
%  the impedance row to be resolved, so the streamwise mesh is clustered on it
%  with a Gaussian point density: the target ratio is 4 at the neck, falling off
%  over 60 displacement thicknesses, which puts about 20 stations in the mouth.
%
%  Base flow: baseflow/caller/caller_CBL_resonator.m (writes bf_compressible_m006_fine.mat
%  at high resolution, so this step is a clean downsampling).
%  Writes input/DeHNSSo_input_resonator.mat  (Table C.7: 5000 x 40, FD4)

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_compressible_m006_fine.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
params.n_eta_new = 40;
params.n_xi_new  = 4000;
params.rescale   = false;
params.plot      = false;

%% Neck clustering — declare the mouth, gridgen sizes the mesh to it
% params.forcing is in the units of the base flow's X field: displacement
% thicknesses here, so the 1 mm neck becomes 1e-3/lref. Declaring it also stores
% the mouth on the grid as StabGrid.forcing, which the caller reads back, so the
% two cannot drift apart.
lref = load(fullfile(input.folder, input.filename)).BF.lref;
% x_r = 0.55 m is the neck CENTRE, not its leading edge: the COMSOL reference has
% wall |v| over 0.549545..0.550445 m. 'Xc = 0.55/lref + s_nd/2' put the mouth half
% a width downstream of it.
s_nd = 1e-3/lref;                       % 1 mm neck, non-dimensional
Xc   = 0.55/lref;                       % neck centre

params.forcing = struct('xc', Xc, 'w', s_nd/2, 'shape', 'cosine');
params.spread  = 30;                    % refinement reach, ~16x the mouth
params.refine  = 5;                     % dxi at the neck vs equidistant

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_resonator.mat';

%% Run
StabGrid = main_gridgen(input, params, output);

x = StabGrid.xun(1,:);
fprintf('\n%d x %d, %d stations in the 1 mm neck at R = %.0f\n', ...
        StabGrid.nx, StabGrid.ny, ...
        nnz(x >= Xc-s_nd/2 & x <= Xc+s_nd/2), sqrt(StabGrid.Re*Xc));
