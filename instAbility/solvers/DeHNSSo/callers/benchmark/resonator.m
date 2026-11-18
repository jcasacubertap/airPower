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
% DeHNSSo — Tollmien-Schlichting wave over a Helmholtz resonator (paper Sec. 6.4)
%  A 1 mm resonator neck is embedded in the wall of a Ma = 0.061 Blasius boundary
%  layer (Michelis et al. 2023: U0 = 20.9 m/s, f_r = 295 Hz, x0 = 0.55 m). Over
%  the mouth the wall condition is the reactive relation v' - Y p' = 0, with Y the
%  cavity admittance at the wave frequency; everywhere else the wall is no-slip
%  and isothermal. The case is run twice, rigid wall and resonator, and the
%  N-factor increment between them is compared with the reference.
%
%  The compressible formulation is required even at this Mach number: the wave
%  and the cavity couple acoustically, and the incompressible limit over-
%  constrains the near-wall pressure.
%
%  Grid: gridgen/benchmark/resonator.m (writes DeHNSSo_input_resonator.mat)

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_resonator.mat'));   % StabGrid

Opt.compressible = 'on';        % (u,v,w,rho,T) path

U0   = StabGrid.Uref;   lref = StabGrid.lref;
f_r  = 295;                     % cavity resonance [Hz], fixed by the geometry
f_TS = f_r;                     % on resonance, where the admittance is largest

%% Spectral content — one TS wave
Stab.M = 1;
Stab.N = 0;                     % 0 = 2D
Stab.omega_0 = 2*pi*f_TS*lref/U0;
Stab.beta_0  = 0;

%% Inflow
Stab.IC         = 'CLST';       % compressible local stability solution
Stab.phaseRef   = 'pwall';
Opt.mode_select = 'slow';
Opt.linear      = 'on';         % the amplitude is arbitrary in a linear run
Stab.A0_fund    = 1e-4;

%% Boundary conditions — no-slip isothermal wall; the neck row is set by Opt.hr
Opt.bc_top          = {'H_DR', 'H_DR', 'H_DR', 'H_DR'};
Opt.bc_wall_thermal = 'isothermal';

%% Buffer — amplitude damping, cleaner than parabolisation at Ma = 0.061:
% the same dN, but a flat field to the outflow, where 'para' rings at the boundary
Opt.buffer = 'on';
Opt.xb     = 90;                % buffer start [% of domain]

%% Solver
Opt.plot    = 'off';
Opt.lu_mode = 'lapack_band';    % banded LU (compressible); build via mexbuild_dgbsv.m
Opt.parfor  = 'off';
Opt.gpu     = 'off';
Opt.solver  = 'picard';

%% Resonator — admittance at the wave frequency (Michelis 2023, on resonance)
Y    = 12.7688*exp(1i*pi);      % non-dimensional |Y| exp(i arg Y) at f_r
x    = StabGrid.xun(1,:);
% x_r = 0.55 m is the neck CENTRE
s_nd = 1e-3/lref;   Xc = 0.55/lref;
mask = (x >= Xc - s_nd/2) & (x <= Xc + s_nd/2);    % mouth stations
fprintf('f = %.0f Hz | Y = |%.3f| exp(%+.0f deg) | %d stations in the neck\n', ...
        f_TS, abs(Y), rad2deg(angle(Y)), nnz(mask));

%% Run — rigid wall first, then the resonator, each a fresh operator
Opt.rerun       = 'off';
Opt.save_output = [];

tic_total = tic;
StabRes0 = main(Stab, StabGrid, Opt);                                  % baseline
Opt.hr   = struct('Y', Y, 'mask', mask, 'uslip', false, 'taper', 0, 'x', x);
StabResR = main(Stab, StabGrid, Opt);                                  % resonator
fprintf('\nTotal execution time: %.1f s\n', toc(tic_total));

%% Post-processing — N-factor increment across the neck
% N is measured from the baseline branch-I neutral point and the resonator curve
% is referenced to the SAME station, so the two are directly comparable.
j  = find(abs(StabRes0.omegavec - Stab.omega_0) < 1e-9, 1);
A0 = abs(StabRes0.A(j,:));   AR = abs(StabResR.A(j,:));
[~, i1] = min(A0);                                  % branch I of the baseline
N0 = log(A0/A0(i1));   NR = log(AR/A0(i1));
fprintf('dN over the flat-plate value: %.3f\n', NR(end) - N0(end));

figure; plot(x, N0, x, NR, 'LineWidth', 1.4); grid on
xline(Xc, ':'); legend('rigid wall', 'resonator', 'Location', 'northwest');
xlabel('$x/\delta^{*}$', 'Interpreter', 'latex'); ylabel('$N$', 'Interpreter', 'latex');
