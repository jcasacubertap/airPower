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
% Generate an incompressible boundary-layer base flow.
%
%   Solves the Falkner-Skan-Cooke similarity problem for the case set in the
%   USER INPUTS block and writes baseflow/output/bf_ts_blasius.mat, ready for
%   gridgen. BL.model selects the incompressible (0) or compressible (1) path.
%
%   Run this before gridgen/benchmark/ts_blasius.m.
clc; clear all; format long
caller_dir = fileparts(mfilename('fullpath'));
rootdir    = fullfile(caller_dir, '..', '..');
addpath(genpath(rootdir));

BL.model = 0;                % 0 = incompressible, 1 = compressible

%% Inputs

% Mesh
BL.S   = 0.2418;             % streamwise domain start [m]
BL.L   = 1.7213;             % streamwise domain end [m]
BL.H   = 0.0585;             % wall-normal domain height [m]
BL.y_i = BL.H/20;            % Malik median node height [m]
BL.nx  = 800;                % streamwise gridpoints
BL.ny  = 40;                 % wall-normal gridpoints

% Freestream
V      = 10;                 % edge velocity at inlet [m/s]
BL.Ue  = ones(1,BL.nx)*V;    % streamwise edge velocity (ZPG: constant)
BL.We  = 0;                  % spanwise edge velocity (unswept: 0)
BL.nu  = 1.511e-5;           % kinematic viscosity [m^2/s]

plt.mainFig = 1;             % plot after solving (0 = off)

%% Boundary-layer fields
[BL] = incompressibleBL_B1(BL);

%% Reference scales and non-dimensionalisation
BF.Uref = V;                          % reference velocity [m/s]
BF.lref = sqrt(BL.S*BL.nu/BF.Uref);   % Blasius length scale [m]
BF.nu   = BL.nu;
BF.Re   = BF.Uref*BF.lref/BF.nu;      % Re_0 = sqrt(Re_x)

BF.X = BL.x/BF.lref;
BF.Y = BL.y/BF.lref;
BF.U = BL.u/BF.Uref;
BF.V = BL.v/BF.Uref;
BF.W = BL.w/BF.Uref;
BF.H   = BL.H   / BF.lref;
BF.y_i = BL.y_i / BF.lref;
% Derivatives are produced by gridgen (compute_bf_derivatives).

%% Postprocessing
[BL] = post_BL(BL);

%% Save
% Written to baseflow/output/, which is scratch and gitignored. The gridgen
% benchmarks read the tracked copies in baseflow/output/benchmark/; copy a
% regenerated file there to make it the benchmark input.
outdir = fullfile(caller_dir, '..', 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end
save(fullfile(outdir, 'bf_ts_blasius.mat'), 'BF');
