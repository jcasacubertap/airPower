%% Gridgen benchmark — SweptWing CFI over a smooth hump
%  Curvilinear structured input (80 × 1200, Y follows a Gaussian hump).
%  Resampled onto Malik-Chebyshev η × equidistant wall-arc-length ξ.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_hump.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters
params.n_eta_new = 60;              % heritage hump caller uses 50 × 1000
params.n_xi_new  = 1200;
params.Uref      = 12.41792766;
params.lref      = 2.3403348e-04;   % Blasius length at R=400 (matches bf_blasius.mat)
params.Re        = 199.518703662;
params.y_i       = [];              % auto (H/10)
params.H         = [];              % auto
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = true;            % BF already non-dim
params.plot      = false;

%% Streamwise refinement — Gaussian cluster around the hump centre.
% Same parameters as heritage Example_Caller_CFI_hump_in_SweptWing_BL
% (Grid.mug = x_m, Grid.sig = 0.2, Grid.ag = 0.5). Accepts vectors for
% multiple humps: e.g. params.mug = [x_m1, x_m2], sig/ag scalar or vector.
x_m = 0.08833 / 2.1394e-4;  % hump centre ≈ 859 (matches heritage)
params.mug = [x_m];
params.sig = [0.2];
params.ag  = [0.5];

%% Grid-generation options
params.FD_xi_order_1    = 4;
params.FD_xi_order_2    = 2;
params.FD_eta_method   = 'fd';   % sparse FD pentadiagonal — ~80× fewer nnz than cheb;
                              % avoids OOM in UMFPACK factorisation on laptop
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
