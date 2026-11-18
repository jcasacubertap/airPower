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

%% Grid-generation options
params.FD1_order    = 4;
params.FD2_order    = 2;
params.eta_method   = 'cheb';
params.FD_eta_order = 4;
%params.plot =1 ; Antonis for plotting. 

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);



