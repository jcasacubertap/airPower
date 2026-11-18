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
%% Description
% DeHNSSo — Stationary CFI over a forward-facing step in a swept-wing BL
%  Base flow: DeHNSSo_input_cfi_step.mat (run gridgen/benchmark/cfi_step.m first).
%
%  The step is an embedded boundary: gridgen finds the solid in the base flow,
%  snaps a station to its face and a node to its top, and the solver imposes
%  no-slip on the surface nodes (Opt.solid = 'auto'). Linear, one stationary
%  crossflow mode introduced by a blowing/suction strip upstream of the step
%  (DeHNSSo-v1 Eq. 28). Compared with the DNS on the shape of u'_max(x): both
%  are linear, so the HNS curve is scaled to the DNS at the crossflow peak.
%
%  Memory: 1500 x 80 FD stations. Keep Opt.lu_ordering = 'symrcm'; COLAMD fill
%  grows with nx and needs far more than a 32 GB laptop for one factorisation.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_cfi_step.mat'));   % StabGrid

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 0;                 % temporal modes (stationary)
Stab.N = 1;                 % spanwise modes: MFD and the +-beta pair
Stab.omega_0 = 0;
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;   % lambda = 7.5 mm

%% Inflow — none: the mode is forced at the wall below
Stab.IC       = 'ILST';
Stab.phaseRef = 'umax';     % 'pwall' (TS) | 'umax' (CFI)
Stab.A0       = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Opt.linear    = 'on';

%% Wall actuator — blowing/suction strip on mode (0,1), x = 135.5 .. 143.6
nf    = (2*Stab.N+1)*(2*Stab.M+1);
j_act = ModeToModeNumber(0, 1, Stab.M, Stab.N);
xbcw  = [135.5, 143.6];
Abs   = 1.414*0.5e-5;       % peak wall-normal velocity (DeHNSSo-v1 value)
Stab.bcx = linspace(StabGrid.xun(1), StabGrid.xun(end), 5000);
Stab.bcw = zeros(3, numel(Stab.bcx), nf);
Stab.bct = zeros(3, numel(Stab.bcx), nf);
in_win = Stab.bcx >= xbcw(1) & Stab.bcx <= xbcw(2);
env    = 4*(Stab.bcx - xbcw(1)).*(xbcw(2) - Stab.bcx) / (xbcw(2)-xbcw(1))^2;
Stab.bcw(2, in_win, j_act) = Abs * env(in_win).^3;

%% Boundary conditions (freestream type per component u,v,w; wall Dirichlet)
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
Opt.solid  = 'auto';        % 'auto' masks on a snapped grid, ghost row otherwise | 'mask' | 'ghost'

%% Buffer — 'on' damping | 'para' parabolisation | 'off'
Opt.buffer   = 'para';
Opt.xb       = 85;          % buffer start [% of domain]
Opt.nltbufxb = 80;

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';   % 'umax' | 'urms' | 'energy'
Opt.lu_mode     = 'auto';   % 'full' | 'auto' | 'lapack_band' | 'none' (DeHNSSo_guide §8)
Opt.lu_ordering = 'symrcm'; % banded ordering: fill stays flat in nx (see header)
Opt.parfor      = 'off';
Opt.solver      = 'picard'; % 'picard' | 'newton'

%% Run
Opt.rerun       = 'off';    % 'on' reuses LHS from a previous run
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_cfi_step.mat');   % [] to skip
tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing — overlay on the DNS in paper coordinates
close all
ref_file = fullfile(rootdir, 'literature', 'StabRes_SweptWing_step_DNS.mat');
x_paper  = 3.604*StabGrid.xun + 220.6;          % DNS coordinates
A_hns    = StabRes.A(end, :);                   % mode (0,1)
if isfile(ref_file)
    D = load(ref_file);  X_dns = D.StabResDNS.X(:);  A_dns = D.StabResDNS.A(:);
    [~, ip] = max(A_dns);
    scale = A_dns(ip) / interp1(x_paper, A_hns, X_dns(ip), 'linear', 'extrap');
else
    X_dns = [];  A_dns = [];  scale = 1;
end
if strcmpi(Opt.plot, 'on')
    figure(1); clf
    plot(x_paper, scale*A_hns, 'k-', 'LineWidth', 1.5); hold on
    if ~isempty(X_dns)
        plot(X_dns, A_dns, 'm:', 'LineWidth', 1.5);
        legend({'DeHNSSo-v2', 'DNS'}, 'Location', 'best');
    end
    xline(3.604*StabGrid.step.x + 220.6, 'k--');
    grid on; xlim([700 920]); xlabel('x (paper)'); ylabel('u''_{max}');
end
