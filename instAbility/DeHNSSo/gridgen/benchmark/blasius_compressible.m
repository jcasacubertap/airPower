%% Gridgen benchmark — Compressible Blasius (FSC) base flow
%  Cartesian structured input from baseflow/caller/caller_CBL.m
%  (already at DeHNSSo grid resolution — passed through unchanged).
%  build_stab_grid auto-detects ρ/T and runs Sutherland to populate
%  μ, ∂μ/∂T, ∂²μ/∂T².  Thermo metadata (Ma, Pr, γ, S*) is forwarded.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output');
input.filename   = 'BF_compressible.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling parameters (pass-through: BF already on target DeHNSSo grid)
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;       % BF is already non-dim
params.plot      = false;

%% Optional: streamwise trim & interpolation method (cartesian path)
% params.xi_trim_inflow  = 0.0;
% params.xi_trim_outflow = 0.0;
% params.interp_method   = 'makima';  % 'makima' | 'spline' | 'linear' | 'cubic'

%% Grid-generation options
% FD_xi_order_2 = 4 needed by LHS_compressible for the BF ∂²/∂x² viscous terms.
params.FD_xi_order_1   = 4;
params.FD_xi_order_2   = 4;
params.FD_eta_method   = 'fd';
params.FD_eta_order    = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_compressible.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
