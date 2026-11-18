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
% DeHNSSo — Synthetic jet actuator on a Blasius TS wave (benchmark)
%  Reference caller for the wall-actuator workflow. A linear TS wave (M=1) in a
%  Blasius BL is actuated by a wall-mounted synthetic jet: a wall-normal velocity
%  forcing v'|_wall over the orifice, at the TS frequency with prescribed
%  amplitude and phase. Models an active TS-control actuator (cf. the TS-control
%  study).
%
%  New actuator workflow (see docs/DeHNSSo_guide.md "Wall actuation"):
%    1. gridgen/benchmark/wall_jet.m clusters the mesh at the jet centre
%       (params.mug) and stores it as StabGrid.cluster.xc.
%    2. add_wall_forcing reads that centre and builds a raised-cosine forcing of a
%       PHYSICAL half-width w (not a grid-point count), so the refinement
%       resolves the same orifice better instead of shrinking it. Calls are
%       additive, so several actuators just stack.
%
%  Base flow / grid: run baseflow/.../bf_ts_blasius.mat, then
%  gridgen/benchmark/wall_jet.m (writes DeHNSSo_input_wall_jet.mat).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()

input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input_wall_jet.mat'));   % clustered grid + cluster.xc

%% Spectral content — single TS mode
Stab.M = 1;                 % omega (temporal) modes
Stab.N = 0;                 % beta  (spanwise) modes
Stab.omega_0 = 0.034576;    % fundamental angular frequency (F = 86e-6)
Stab.beta_0  = 0;           % 2D

%% Inflow — linear TS from local stability
Stab.IC       = 'ILST';
Stab.phaseRef = 'pwall';
Opt.linear    = 'on';       % linear: no NLT, single solve
Stab.A0_fund  = 0.0025*sqrt(2);

%% Freestream BC — free (unconstrained) v; u,w Dirichlet
Opt.bc_top = {'H_DR', 'free', 'H_DR'};

%% Buffer — 'on' (amplitude damping) | 'para' (parabolisation)
Opt.buffer = 'on';
Opt.xb     = 80;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';
Opt.lu_mode     = 'full';
Opt.parfor      = 'off';
Opt.solver      = 'picard';

%% Baseline — incoming TS wave only (assembles and returns LHS for reuse)
Opt.rerun       = 'off';
Opt.save_output = [];
[StabRes0, StabGrid, LHS] = main(Stab, StabGrid, Opt);

%% Wall forcing: apply the profile gridgen stored (add_wall_forcing) ----------
% The shape, location and width are defined in gridgen/benchmark/wall_jet.m
% (StabGrid.forcing). Here we only set the strength. The conjugate mode is added
% automatically for nonlinear runs (here linear, so just the fundamental).
A_jet   = 2e-3;          % jet velocity amplitude (physical 2|u|, as Stab.A0)
phi_jet = 0;             % jet phase [rad] (control knob)

[Stab, jet] = add_wall_forcing(Stab, StabGrid, Opt, struct('A', A_jet, 'phi', phi_jet));
% To add a second forcing, call again — it superposes:
%   Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('xc',x2,'w',12,'A',A2));

%% Actuated run — reuse LHS (only the RHS changes)
Opt.rerun = 'on';
[StabResA, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);

%% Jet response at the orifice — the wall-normal velocity actually delivered
vr = squeeze(StabResA.v(2,:,:)) .* StabResA.A(2,:);

% average over the orifice stations (jet.idx from add_wall_forcing)
vav  = mean(vr(end,jet.idx));
vmax = max(vr(end,jet.idx));

fprintf('\n orifice: mean |v''| = %.3e, peak |v''| = %.3e, arg = %+.1f deg  (set A = %.1e)\n', ...
        abs(vav), abs(vmax), rad2deg(angle(vav)), A_jet);

%% Effect of the jet — TS amplitude with and without forcing
x     = StabGrid.xun(1, :);   nx = numel(x);
Re0   = StabGrid.Re;          Rx = sqrt(Re0 * x);
xc    = StabGrid.forcing(1).xc;
jA    = find(abs(StabRes0.omegavec(:).' - Stab.omega_0) < 1e-9, 1);
imeas = max(1, round(Opt.xb*nx/100) - 1);

Au0 = StabRes0.A(jA, :);   AuA = StabResA.A(jA, :);
dNu = log(AuA(imeas) / Au0(imeas));
fprintf(' jet Rx = %.0f,  measure Rx = %.0f,  dN_u'' = %+.3f\n', ...
        sqrt(Re0*xc), Rx(imeas), dNu);

xr = 1:imeas;
figure('Position', [100 100 760 460]);
plot(Rx(xr), log(Au0(xr)./Au0(1)), 'k-',  'LineWidth', 1.5, 'DisplayName', 'baseline (no jet)'); hold on
plot(Rx(xr), log(AuA(xr)./AuA(1)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'jet on');
xline(sqrt(Re0*xc), 'b:', 'jet', 'HandleVisibility', 'off');
xlabel('$R_x = \sqrt{Re_0\,x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N_{u''}$', 'Interpreter', 'latex', 'FontSize', 12)
title(sprintf('TS N-factor, synthetic jet (dN_{u''} = %+.2f)', dNu))
legend('Location', 'best'); grid on; box on
