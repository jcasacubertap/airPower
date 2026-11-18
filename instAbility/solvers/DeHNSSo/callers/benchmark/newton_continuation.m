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
% DeHNSSo — Blasius K-type breakdown via Newton pseudo-arclength continuation
%  Continuation in the inflow amplitude from A0_fund = 2.40e-3 to 3.00e-3, on a
%  400 x 40 grid (R = 400-1067) at M=1. The Jacobian is 2N = 256k, so this runs
%  in minutes on a laptop. A direct solve at the target does not converge from a
%  cold start, which is why this is a continuation and not a single Newton solve.
%
%  Meant as a first Newton run and as a check that MUMPS is working. It is NOT a
%  benchmark: the grid is coarse and the domain short, so the amplitudes are not
%  converged. The full case (1050 x 60, M=2, Jacobian 2N = 1.512M) needs ~80 GB
%  and a workstation; see docs/DeHNSSo_guide.md for the sizing rule.
%
%  REFERENCE (this grid, these settings). If your run matches these, Newton and
%  MUMPS are working:
%      converged in 6 steps, 1 factorisation
%      max A(1,0) = 2.641824e-02
%      N-factor   = 2.175482
%  Took 4.0 min and ~8.6 GB on an M-series laptop. The step count and timing may
%  shift with a different MUMPS build; the two numbers should not.
%
%  MUMPS is required and is a separate dependency. If you have not built it, run
%  src/external/mumps_5.3.5/install/install_mumps_{macos,linux}.sh, and start MATLAB with
%  OMP_NUM_THREADS=1 (see docs/DeHNSSo_guide.md, "Running it: MUMPS, memory and where").

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))

% Your own MUMPS build goes on the path AFTER the genpath calls above, whose
% src/external/mumps_5.3.5/dmumps.m would otherwise shadow it.
mumpsdir = fullfile(getenv('HOME'), 'mumps_build', 'MUMPS_5.3.5', 'MATLAB');
if isfolder(mumpsdir); addpath(mumpsdir, '-begin'); end
fprintf('dmumps resolves to: %s\n', which('dmumps'));

input.folder = fullfile(rootdir, 'input');
gridname     = 'DeHNSSo_input_newton_continuation.mat';   % 400 x 40 Cheb, R = 400-1067: gridgen/benchmark/newton_continuation.m
G = load(fullfile(input.folder, gridname)); StabGrid = G.StabGrid;

%% Spectral content
Stab.M = 1;                 % fundamental + MFD only: the cheapest useful Newton run
Stab.N = 0;
Stab.omega_0 = 0.034576;    % same fundamental as newton_continuation.m
Stab.beta_0  = 0;

%% Inflow
Stab.IC       = 'ILST';
Stab.phaseRef = 'pwall';
Opt.linear    = 'off';      % nonlinear
Stab.A0_fund  = 2.40e-3;    % continuation START amplitude

%% Boundary conditions and buffer
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
Opt.buffer = 'para';
Opt.xb     = 85;

%% Solver
Opt.plot   = 'off';
Opt.rerun  = 'off';
Opt.solver = 'newton';

%% Pseudo-arclength continuation (src/newton/solve_newton.m)
Opt.newton_target = 3.00e-3;   % target A0_fund; a short continuation on purpose
Opt.newton_start  = 2.40e-3;
Opt.newton_seed2  = 2.50e-3;

%% Run
t = tic;
[StabRes, info] = solve_newton(Stab, StabGrid, Opt);
fprintf('\n%s in %d steps, %d factorisation(s), %.1f min\n', ...
        info.status, info.steps, info.nfac_total, toc(t)/60);

%% Check — compare these three against the reference in the header
jF = find(StabRes.omegavec(:).' == Stab.omega_0, 1);   % fundamental
fprintf('\ncheck:  steps = %d,  max A(1,0) = %.6e,  N-factor = %.6f\n', ...
        info.steps, max(StabRes.A(jF,:)), ...
        log(max(StabRes.A(jF,:)) / StabRes.A(jF,1)));
