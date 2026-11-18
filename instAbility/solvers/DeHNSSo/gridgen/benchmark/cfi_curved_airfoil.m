%% Gridgen benchmark — M3J curved surface (unstructured CSV from CFD)
%  Wall detected from the boundary polygon + velocity threshold, then a
%  body-fitted Malik-Chebyshev mesh is launched by column interpolation.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_cfi_curved_airfoil.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling (unstructured path)
params.n_eta_new = 80;
params.n_xi_new  = 5000;
params.y_i       = 8;            % auto from delta_99
params.H         = [];           % auto from projection
params.xi_range  = [];
params.xi_trim_inflow  = 0.10;   % drop noisy inflow
params.xi_trim_outflow = [];
params.elliptic  = false;        % true for strongly concave walls
params.plot      = true;
params.plot_opts.xi_slices  = linspace(0, 1, 8);
params.plot_opts.eta_slices = [0.0, 0.1, 0.3];

% CSV is dimensional. PLACEHOLDER reference scales (revise from the CFD setup):
% Uref ~ max|u|, lref = Blasius delta at inflow, Re = Uref*lref/nu_air.
params.rescale = true;
params.Uref    = 12.41792766;
params.lref    = 2.3403348e-04;
params.Re      = 199.518703662;

%% Grid operators (fd eta: cheb + curvilinear metrics too large for UMFPACK on 32 GB)
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_cfi_curved_airfoil.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
