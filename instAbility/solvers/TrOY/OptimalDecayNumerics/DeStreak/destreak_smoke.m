%% DESTREAK_SMOKE  End-to-end smoke test of the DeHNSSo drive + streak extraction.
%
%  Purpose: de-risk the three unknowns before building the operator assembly —
%    (1) does main() run headless on a freshly built flat grid,
%    (2) does the ILST inflow converge for the swept-wing CFI streak,
%    (3) can we pull out the fundamental beta0 streak u_hat(x,y), v_hat(x,y).
%
%  Reduced resolution on purpose (fast). Not a physics result.

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here);
srcroot = fullfile(here, '..', 'DeHNSSo', 'src');
addpath(genpath(srcroot));

%% 1. Flat base-flow grid (small)
gridfile = fullfile(here, 'flat_grid_smoke.mat');
if isfile(gridfile)
    fprintf('[smoke] loading cached grid %s\n', gridfile);
    S = load(gridfile); StabGrid = S.StabGrid;
else
    StabGrid = destreak_make_flatgrid(40, 300, gridfile);
end
fprintf('[smoke] grid: ny=%d nx=%d  x=[%.2f %.2f]\n', ...
        StabGrid.ny, StabGrid.nx, StabGrid.x(1), StabGrid.x(end));

%% 2. Stability setup (swept-wing CFI fundamental, linear)
Stab.M = 0;  Stab.N = 1;  Stab.omega_0 = 0;
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;
Stab.IC       = 'ILST';
Stab.phaseRef = 'umax';
Stab.A0_fund  = 3.5e-3/2;

Opt.linear   = 'on';
Opt.bc_top   = {'H_DR','H_DR','H_DR'};
Opt.buffer   = 'para';
Opt.xb       = 85;
Opt.lu_mode  = 'none';     % backslash — no MEX needed for the small smoke grid
Opt.parfor   = 'off';
Opt.gpu      = 'off';
Opt.plot     = 'off';

%% 3. Run
t0 = tic;
[StabRes, StabGrid] = main(Stab, StabGrid, Opt);
fprintf('[smoke] main() done in %.1f s\n', toc(t0));

%% 4. Extract the fundamental beta0 streak
tol   = 1e-9 * max(1, abs(Stab.beta_0));
jfund = find(abs(StabRes.betavec - Stab.beta_0) < tol & abs(StabRes.omegavec) < 1e-12, 1);
assert(~isempty(jfund), 'fundamental beta0 mode not found in StabRes');
fprintf('[smoke] fundamental mode index = %d  (beta=%.4f, omega=%.4f)\n', ...
        jfund, StabRes.betavec(jfund), StabRes.omegavec(jfund));

uhat = squeeze(StabRes.u(jfund,:,:));   % [ny x nx] complex streak
vhat = squeeze(StabRes.v(jfund,:,:));   % [ny x nx] complex roller (wall-normal)

%% 5. Sanity diagnostics
%  Wall / freestream orientation from base flow (no-slip => U~0 at wall).
[~, iwall] = min(abs(StabGrid.U(:, round(end/2))));
[~, itop ] = max(abs(StabGrid.U(:, round(end/2))));
fprintf('[smoke] wall row=%d (U=%.3e), top row=%d (U=%.3e)\n', ...
        iwall, StabGrid.U(iwall,round(end/2)), itop, StabGrid.U(itop,round(end/2)));

Umax_x  = max(abs(uhat), [], 1);        % |u| peak vs x
fprintf('[smoke] |u_hat| peak: inflow=%.3e  mid=%.3e  outflow(pre-buffer)=%.3e\n', ...
        Umax_x(1), Umax_x(round(end/2)), Umax_x(round(0.8*end)));
fprintf('[smoke] max|v_hat|/max|u_hat| = %.3e (roller/streak ratio)\n', ...
        max(abs(vhat(:)))/max(abs(uhat(:))));

save(fullfile(here,'smoke_result.mat'), 'StabGrid','StabRes','Stab','Opt','jfund','-v7.3');
fprintf('[smoke] OK — saved smoke_result.mat\n');
