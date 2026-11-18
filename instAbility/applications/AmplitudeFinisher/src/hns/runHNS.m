function out = runHNS(StabGrid, A0_fund, N, in, mode, LHS)
% RUNHNS  One DeHNSSo Harmonic-NS solve for a given amplitude / mode count.
%
%   out = runHNS(StabGrid, A0_fund, N, in, mode)
%   out = runHNS(StabGrid, A0_fund, N, in, mode, LHS)   % reuse factorised LHS
%
% Assembles Stab/Opt from `in`, sets the amplitude knob Stab.A0_fund = A0_fund
% and spanwise mode count Stab.N = N, and calls DeHNSSo main().
%
% mode: 'linear'    -> Opt.linear = 'on'  (single solve, fundamental only)
%       'nonlinear' -> Opt.linear = 'off' (Picard/Anderson convergence loop)
%
% Pass LHS (from a previous out.LHS at the SAME N and base flow) to reuse the
% assembled + factorised operators across A0 sweeps (Opt.rerun = 'on'); the LHS
% is A0-independent so this skips assembly/factorisation.
%
% Returns a bundle:
%   out.res      StabRes  (post-processed: u/v/w/p [nmode x ny x nx], A, betavec)
%   out.grid     StabGrid returned by main (adds computed fields)
%   out.LHS      factorised operators (feed back for reuse at same N)
%   out.lref     reference length [m]  (re-dimensionalization)
%   out.A0_fund  amplitude used
%   out.N        mode count used

Stab = in.Stab;
Opt  = in.Opt;

Stab.N       = N;
Stab.A0_fund = A0_fund;
Stab.beta_0  = 2*pi * StabGrid.lref / Stab.lambda_z;

switch lower(mode)
    case 'linear';    Opt.linear = 'on';
    case 'nonlinear'; Opt.linear = 'off';
    otherwise; error('runHNS:mode','mode must be ''linear'' or ''nonlinear''');
end

if nargin >= 6 && ~isempty(LHS) && isstruct(LHS) && ~isempty(fieldnames(LHS))
    Opt.rerun = 'on';
    [StabRes, StabGridOut, LHSout] = main(Stab, StabGrid, Opt, LHS);
else
    Opt.rerun = 'off';
    [StabRes, StabGridOut, LHSout] = main(Stab, StabGrid, Opt);
end

out.res     = StabRes;
out.grid    = StabGridOut;
out.LHS     = LHSout;
out.lref    = StabGrid.lref;
out.Uref    = StabGrid.Uref;          % [m/s] for re-dimensionalizing w'
out.A0_fund = A0_fund;
out.N       = N;
end
