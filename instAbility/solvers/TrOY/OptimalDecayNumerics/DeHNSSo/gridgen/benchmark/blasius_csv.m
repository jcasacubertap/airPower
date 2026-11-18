%% Gridgen benchmark — Blasius BL base flow (cartesian CSV from OpenFOAM)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_blasius.csv';
input.format     = 'csv';
input.structured = false;

% base flow is stored gzipped in the repo; unpack on first use
if ~isfile(fullfile(input.folder, input.filename))
    gunzip(fullfile(input.folder, [input.filename '.gz']));
end

%% Resampling parameters
params.n_eta_new = 60;
params.n_xi_new  = 1200;
params.rescale   = true;
params.Uref      = 20.9;
params.lref      = 2.891866e-04;   % Blasius length at R=400 (matches bf_blasius.mat)
params.Re        = 400.03309;
params.plot      = true;

% R window [400, 1067], as fractions of the x-span (matches bf_blasius.mat)
params.xi_range = [0.08877, 0.63344];

%% Grid-generation options
params.FD_xi_order_1  = 4;
params.FD_xi_order_2  = 2;
params.FD_eta_method  = 'cheb';
params.FD_eta_order   = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
