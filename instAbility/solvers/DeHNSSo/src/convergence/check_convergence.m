function [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes)
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
% Evaluate convergence, ramp amplitude, extrapolate forcing.
%
%   [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%
%   Measures the relative change in amplitude between successive iterations
%   (paper Eq. 21), and once it drops below threshold ramps the amplitude
%   factor Opt.AF towards one. The nonlinear forcing is extrapolated in A^2
%   from the last three amplitude levels to accelerate convergence.
%

%% Inputs
% Name     size unit   explanation
% StabGrid      [-]    grid structure (nx, D_xi, ...)
% Bound         [-]    boundary/buffer structure (ibuf sets convergence window)
% Mode          [-]    mode structure (RunJ: active modes)
% Opt           [-]    options and iteration state (AF, dal, k1, ...)
% Stab          [-]    stability data (forcing field f, extrapolation window flist)
% StabRes       [-]    stability results (amplitudes A, Aold)

%% Outputs
% Name    size unit   explanation
% Opt          [-]    options and iteration state (AF, dal, k1, ...)
% Stab         [-]    stability data (forcing field f, extrapolation window flist)
% StabRes      [-]    stability results (amplitudes A, Aold)

%% Convergence measure
ib = min(Bound.ibuf, StabGrid.nx);  % when buffering is off, ibuf = nx+1 (sentinel for "no limit"); clamp to nx
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
