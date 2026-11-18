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
% DeHNSSo — Compressible Blasius boundary layer (TS waves)
%  Base flow: DeHNSSo_input_ts_compressible.mat (gridgen + baseflow/caller/caller_CBL.m;
%  must carry rho, T and thermo metadata Ma, Pr, gamma, S*, mu, dmuT, d2muT2).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(rootdir, 'input', 'DeHNSSo_input_ts_compressible.mat'));   % StabGrid

Opt.compressible = 'on';        % (u,v,w,rho,T) path

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 1;
Stab.N = 0;                     % 0 = 2D
Stab.F_0     = 40e-6;           % reduced frequency
Stab.omega_0 = Stab.F_0 * StabGrid.Re;
Stab.beta_0  = 0;

%% Inflow
Stab.IC         = 'CLST';       % compressible: CLST or LOAD
Stab.phaseRef   = 'pwall';
Opt.mode_select = 'slow';       % Mack F/S pick at Ma>1: 'slow'/'fast'/'auto'
Opt.linear      = 'on';         % 'on' linear | 'off' nonlinear
Stab.A0_fund    = 2e-6;         % fundamental amplitude (physical 2|u|)

%% Boundary conditions (freestream type per component u,v,w,T; wall Dirichlet)
% 'H_DR' Dirichlet | 'IH_DR' inhomogeneous | 'H_NM' Neumann (absorbs acoustics)
Opt.bc_top          = {'H_DR', 'H_DR', 'H_DR', 'H_DR'};
Opt.bc_wall_thermal = 'isothermal';   % T'=0 at wall

%% Buffer — at Ma>0 parabolisation is required ('on'/'off' reflect acoustics)
Opt.buffer = 'para';
Opt.xb     = 90;                % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';       % 'umax' | 'urms' | 'energy'
Opt.lu_mode     = 'lapack_band';  % 'full'/'auto' OOM at nx=800 (DeHNSSo_guide §8)
Opt.parfor      = 'off';
Opt.solver      = 'picard';     % 'picard' | 'newton' (pseudo-arclength continuation, needs Opt.newton_target)

%% Run
Opt.rerun       = 'off';
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_ts_compressible.mat');   % [] to skip

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing — HNS vs Zanus/Bertolotti 2017 PSE (Ma=0.5, F=40e-6)
close all
ref_mat = fullfile(rootdir, 'literature', 'Zanus_PSE_CPG_M05_F40.mat');   % [] to skip

if strcmpi(Opt.plot, 'on')
    plot_amplitudes(StabGrid, StabRes, [], [], Opt.plot_metric);
    plot_stability(StabGrid, StabRes, Opt);

    if ~isempty(ref_mat) && isfile(ref_mat)
        % HNS alpha is global non-dim (fixed lref at inflow); Zanus is local
        % (per-station lref). Rescale HNS -> local by R/R_0.
        R_sweep   = sqrt(StabGrid.xun * StabGrid.xun(1));
        alpha_hns = squeeze(StabRes.alpha(end, :)) .* (R_sweep ./ R_sweep(1));
        ZN        = load(ref_mat).ZN;   % [R, -alpha_i]

        figure('Position', [100 100 1000 400]);
        subplot(1, 2, 1);
        plot(R_sweep, -imag(alpha_hns), 'r--', 'LineWidth', 1.5, 'DisplayName', 'HNS'); hold on
        plot(ZN(:,1),  ZN(:,2),         'ms',  'MarkerSize', 5,  'DisplayName', 'Zanus/Bertolotti 2017 PSE');
        xlabel('R = \surd Re_x'); ylabel('-\alpha_i');
        grid on; legend('Location', 'best');
        title(sprintf('Spatial growth rate, Ma=%.2f, F=%.0fe-6', StabGrid.Ma, Stab.F_0*1e6));

        subplot(1, 2, 2);
        plot(R_sweep, real(alpha_hns), 'r--', 'LineWidth', 1.5, 'DisplayName', 'HNS');
        xlabel('R = \surd Re_x'); ylabel('\alpha_r');
        grid on; legend('Location', 'best');
        title('Streamwise wavenumber');
    end
end
