function [uhat, vhat, what, StabRes] = destreak_run_case(G, beta0, opt)
%DESTREAK_RUN_CASE  Linear HNS solve on grid G; return the fundamental streak.
%
%   [uhat, vhat, what, StabRes] = destreak_run_case(G, beta0, opt)
%
%   Runs main() in LINEAR mode with the swept-wing CFI settings and an ILST
%   inflow at beta0, and extracts the fundamental (m=0,n=1) mode fields as
%   [ny x nx] complex arrays. ILST is used on every grid: since the relief is
%   zero at the inflow station, the inflow streak is identical across grids.
%
%   opt (optional): .lu_mode (default 'none'), .xb (default 85), .quiet (true).

if nargin < 3; opt = struct(); end
if ~isfield(opt,'lu_mode'); opt.lu_mode = 'none'; end
if ~isfield(opt,'xb');      opt.xb = 85;          end
if ~isfield(opt,'quiet');   opt.quiet = true;     end

Stab.M = 0;  Stab.N = 1;  Stab.omega_0 = 0;  Stab.beta_0 = beta0;
Stab.IC = 'ILST';  Stab.phaseRef = 'umax';  Stab.A0_fund = 3.5e-3/2;

Opt.linear = 'on';
Opt.bc_top = {'H_DR','H_DR','H_DR'};
Opt.buffer = 'para';  Opt.xb = opt.xb;
Opt.lu_mode = opt.lu_mode;  Opt.parfor = 'off';  Opt.gpu = 'off';  Opt.plot = 'off';

if opt.quiet
    [~, StabRes] = evalc('main(Stab, G, Opt)');
else
    StabRes = main(Stab, G, Opt);
end

tol   = 1e-9 * max(1, abs(beta0));
jf = find(abs(StabRes.betavec - beta0) < tol & abs(StabRes.omegavec) < 1e-12, 1);
assert(~isempty(jf), 'destreak_run_case: fundamental beta0 mode not found');

uhat = squeeze(StabRes.u(jf,:,:));
vhat = squeeze(StabRes.v(jf,:,:));
what = squeeze(StabRes.w(jf,:,:));
end
