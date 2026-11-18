%% Gridgen benchmark — flat-plate SweptWing CFI base flow
%  Curvilinear structured input (80 × 1200 body-fitted, Y constant per column).
%  Resampled onto Malik-Chebyshev η × equidistant ξ arc-length.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_flat.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling parameters
params.n_eta_new = 60;           % validated config uses 60 × 1200
params.n_xi_new  = 1200;
params.y_i       = [];           % auto (H/10)
params.H         = [];           % auto
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;        % BF already non-dim
params.plot      = false;

%% Grid-generation options
params.FD1_order    = 4;
params.FD2_order    = 2;
params.eta_method   = 'cheb';
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
