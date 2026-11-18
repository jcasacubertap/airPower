%% Gridgen benchmark — flat-plate SweptWing CFI base flow
%  Curvilinear structured input (80 × 1200 body-fitted, Y constant per column).
%  Resampled onto Malik-Chebyshev η × equidistant ξ arc-length.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_flat.csv'; %'bf_sweptwing_flat.mat';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters
params.n_eta_new = 60;           % validated config uses 60 × 1200
params.n_xi_new  = 1200;
params.Uref      = 12.41792766;
params.lref      = 2.3403348e-04;   % Blasius length at R=400 (matches bf_blasius.mat)
params.Re        = 199.518703662;
params.y_i       = [];           % auto (H/10)
params.H         = [];           % auto
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = true;        % BF already non-dim
params.plot      = false;

%% Grid-generation options
params.FD_xi_order_1    = 4;
params.FD_xi_order_2    = 2;
params.FD_eta_method   = 'cheb';
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
