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
% Gridgen benchmark — SweptWing CFI over a smooth hump
%  Curvilinear structured input (80 x 1200, Y follows a Gaussian hump),
%  resampled onto Malik-Chebyshev eta x wall-arc-length xi.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_cfi_hump.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling
params.n_eta_new       = 80;
params.n_xi_new        = 5000; %2500;
params.y_i             = 8;
params.H               = [];
params.xi_range        = [];
params.xi_trim_inflow  = 0.10;
params.xi_trim_outflow = [];
params.plot            = false;

% CSV is dimensional. PLACEHOLDER reference scales (revise from the CFD setup):
% Uref ~ max|u|, lref = Blasius delta at inflow, Re = Uref*lref/nu_air.
params.rescale = true;
params.Uref    = 12.41792766;
params.lref    = 2.3403348e-04;
params.Re      = 199.518703662;

%% Streamwise refinement — Gaussian cluster at the hump centre
%params.mug    = -1356.38;   % bump crest, x = -317.4381 mm / lref (global coordinates)
%params.spread =    11.73;   % sigma/3 = 2.7459 mm / lref; sigma = bump 'representative' width
%params.refine =        1;   % 1 = uniform; raise only if sigma < 2*lambda_x/N

%% Grid operators (fd eta: sparse pentadiagonal, avoids cheb OOM on laptop)
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_cfi_hump.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
