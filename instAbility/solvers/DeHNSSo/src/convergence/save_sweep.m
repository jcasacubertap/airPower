function [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes)
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
% Save intermediate results during convergence iterations.
%
%   [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes)
%
%   When Opt.Sweep is enabled, saves the current amplitude and wavenumber
%   data to disk at each convergence iteration for post-hoc analysis.

%% Inputs
% Name     size         unit   explanation
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Mode     (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt      (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab     (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes  (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Outputs
% Name    size         unit   explanation
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Advance the sweep counter and store the current amplitude factor
iu = sqrt(-1);
Opt.k2 = Opt.k2 + 1;
Opt.AF(Opt.k2) = Opt.AF;

%% Streamwise derivative of u on the (possibly non-uniform) xi grid
D_xi = StabGrid.D_xi;
for j = Mode.RunJ
    dux(j,:,:) = squeeze(StabRes.u(j,:,:)) * D_xi.';
end

%% Store wavenumber, eigenfunctions and amplitude at each station
yinteg = linspace(0, max(StabGrid.y), 4000);

for j = Mode.RunJ
    for i = 1:StabGrid.nx
        [~, yimax] = max(abs(interp1(StabGrid.y, StabRes.u(j,:,i), yinteg, 'spline')));
        dumax = interp1(StabGrid.y, dux(j,:,i), yinteg(yimax), 'spline');
        umax  = interp1(StabGrid.y, StabRes.u(j,:,i), yinteg(yimax), 'spline');

        StabRes.alphasweep(j, i, Opt.k2) = (dumax / umax) / iu;
        StabRes.usweep(j,:,i, Opt.k2) = StabRes.u(j,:,i) ./ abs(umax);
        StabRes.vsweep(j,:,i, Opt.k2) = StabRes.v(j,:,i) ./ abs(umax);
        StabRes.wsweep(j,:,i, Opt.k2) = StabRes.w(j,:,i) ./ abs(umax);
        StabRes.psweep(j,:,i, Opt.k2) = StabRes.p(j,:,i) ./ abs(umax);

        if j == round(Mode.nf/2)
            StabRes.Asweep(j, i, Opt.k2) = abs(umax);
        else
            StabRes.Asweep(j, i, Opt.k2) = 2 * abs(umax);
        end
    end
end

end
