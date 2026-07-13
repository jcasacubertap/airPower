%% Gridgen benchmark — M3J curved surface (unstructured CSV from CFD)
%  Wall is detected via boundary polygon + velocity threshold, then a body-
%  fitted Malik-Chebyshev mesh is launched from column-by-column interpolation.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_m3j.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters (unstructured path)
params.n_eta_new = 60;
params.n_xi_new  = 1200;
params.y_i       = [];           % auto from delta_99
params.H         = [];           % auto from projection
params.xi_range  = [];
params.xi_trim_inflow  = 0.10;   % drop first 10% of stations (noisy inflow)
params.xi_trim_outflow = [];


% M3J CSV exports are dimensional — rescale explicitly.
% *** PLACEHOLDER SCALES — revise once the CFD setup's reference values are
%     known. Derived from m3j.csv ranges + standard air at 15°C:
%       max|u| ≈ 29.5 m/s, max|w| ≈ 20.6 m/s    → |Ue| ~ 36 m/s; use u-max for now
%       x-domain ∈ [-0.44, -0.0046] m           → |x₀| ≈ 0.44 m at inflow
%       ν_air ≈ 1.47e-5 m²/s
%       δ̄₀ = sqrt(|x₀|·ν/Uref) ≈ 4.6e-4 m
params.rescale = true;
params.Uref    = 30;             % m/s       (≈ max|u|, tentative)
params.lref    = 4.6e-4;         % m         (Blasius at inflow, tentative)
params.Re      = 940;            % = Uref·lref/ν_air ≈ 30·4.6e-4/1.47e-5

params.elliptic = false;         % set true for strongly concave walls
params.plot     = true;
params.plot_opts.xi_slices  = linspace(0, 1, 8);
params.plot_opts.eta_slices = [0.0, 0.1, 0.3];

%% Grid-generation options
params.FD1_order    = 4;
params.FD2_order    = 2;
% Use sparse FD for the wall-normal direction — same reason as hump caller:
% dense Chebyshev + curvilinear metric cross-coupling produces a LHS too large
% for UMFPACK to factor on 32 GB. FD is ~indistinguishable for smooth BL.
params.eta_method   = 'fd';
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
