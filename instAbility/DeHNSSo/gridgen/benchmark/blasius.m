%% Gridgen benchmark — Blasius BL base flow
%  Cartesian structured input, already at DeHNSSo grid resolution.
%  Passed through unchanged (no resampling).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_blasius.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling parameters
% Pass-through: leave n_eta_new empty so the Cartesian branch keeps the
% source grid as-is (bf_blasius already matches DeHNSSo's target).
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;    % already non-dim
params.plot      = false;

%% Optional: streamwise trim & interpolation method (cartesian path)
% Drop fractions of the x-span at inflow / outflow before resampling.
% params.xi_trim_inflow  = 0.02;     % drop first 2% of x
% params.xi_trim_outflow = 0.0;      % keep full outflow

% interp2 method used by resample_cartesian.
% Default 'makima' (handles non-uniform CFD spacing, no overshoot).
% 'spline' is smoother in the interior but can overshoot at boundaries.
% 'cubic' silently falls back to 'spline' on non-uniform grids.
% params.interp_method = 'makima';  % 'makima' | 'spline' | 'linear' | 'cubic'

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
