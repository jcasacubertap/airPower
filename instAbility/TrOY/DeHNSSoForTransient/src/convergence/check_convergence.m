function [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%CHECK_CONVERGENCE  Evaluate convergence, ramp amplitude, extrapolate forcing.

%% Convergence measure
ib = min(Bound.ibuf, StabGrid.nx);  % 'off' buffer sets ibuf = nx+1 as sentinel
for j = Mode.RunJ
    numerator = abs(StabRes.A(j,1:ib) - StabRes.Aold(j,1:ib));
    numerator(numerator < 1e-14) = 0;   % absolute floor (paper Eq. 21)
    Opt.dal(j) = sum(numerator) / max(StabRes.A(j,1:ib)) / StabGrid.nx;
end

% Linear simulation: force convergence after first iteration
if length(Mode.RunJ) == 1
    Opt.dal = zeros(size(Opt.dal));
end

% Save amplitudes
StabRes.Aold = StabRes.A;

%% Amplitude ramping
if Opt.k1 ~= 0
    if max(abs(Opt.dal)) <= (1-Opt.Sweep)*Opt.Conv*Opt.ConvF + Opt.Sweep*Opt.Conv
        Opt.AF = Opt.AF * Opt.AFg;
        if ge(Opt.AF, 1); Opt.AF = 1; end

        % Keep a rolling window of the last 3 AF values and forcing fields
        % (only 3 entries needed for the 2nd-order extrapolation; avoids unbounded growth)
        if ~isfield(Opt, 'AFlist')
            Opt.AFlist     = [Opt.AF, Opt.AF, Opt.AF];
            Stab.flist     = cat(3, Stab.f, Stab.f, Stab.f);
        end
        Opt.AFlist = [Opt.AFlist(2:3), Opt.AF];
        Stab.flist = cat(3, Stab.flist(:,:,2:3), Stab.f);

        % Linearly (in A^2) extrapolate the forcing term
        if Opt.AFlist(end) ~= Opt.AFlist(end-1) && Opt.AFlist(end-1) ~= Opt.AFlist(end-2)
            Stab.f = Stab.flist(:,:,end) + ...
                (Stab.flist(:,:,end) - Stab.flist(:,:,end-1)) / ...
                (Opt.AFlist(end-1)^2 - Opt.AFlist(end-2)^2) * ...
                (Opt.AFlist(end)^2 - Opt.AFlist(end-1)^2);
        end
    end
end

%% Update iteration counter and convergence log
Opt.k1 = Opt.k1 + 1;
Opt.dalsave(:, Opt.k1) = abs(Opt.dal);
Opt.dalsave(Opt.dalsave <= 0) = NaN;
Opt.dalsave(Opt.dalsave == 1) = NaN;

end
