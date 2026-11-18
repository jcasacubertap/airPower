function plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt)
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
% Plot amplitude evolution during the convergence loop.
%
%   plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt)
%
%   Live monitor (figure 2) showing the amplitude of every mode with nonzero
%   amplitude against R = sqrt(Re_x), together with the buffer location, so the
%   convergence of the current iteration can be followed at a glance.
%

%% Inputs
% Name     size unit   explanation
% StabGrid      [-]    grid structure (xi1D, nx)
% StabRes       [-]    stability results (u, A)
% Mode          [-]    mode structure (nf)
% Stab          [-]    stability data (Mvec, Nvec for legend labels)
% Opt           [-]    options (k1, AF, xb)

%% Outputs
% None.

%% Target figure 2 (reuse if it already exists)
if isempty(findobj('Type','figure','Number',2))
    figure(2); set(gcf, 'Visible', 'on');
else
    set(0, 'CurrentFigure', 2); clf;
end

% R = sqrt(Re_x) = sqrt(x_norm) * R0, with R0 = sqrt(x_norm(1))
R0    = sqrt(StabGrid.xi1D(1));
R_vec = sqrt(StabGrid.xi1D) * R0;

% Plot all modes that have nonzero amplitude (not just RunJ)
half = round(Mode.nf/2);
legs = {};
for j = half:Mode.nf
    A_j = max(abs(StabRes.u(j,:,:)), [], 2);
    if max(A_j(:)) < 1e-10; continue; end
    semilogy(R_vec, A_j(:), '-', 'LineWidth', 1.5); hold on;
    legs{end+1} = sprintf('(%d,%d)', Stab.Mvec(j), Stab.Nvec(j));
end
xline(R_vec(round(StabGrid.nx*Opt.xb/100)), 'k--'); legs{end+1} = 'Buffer';
set(gca,'YScale','log'); xlim(R_vec([1 end]));
A_vals = StabRes.A(half:end, :); A_pos = A_vals(A_vals > 0);
if ~isempty(A_pos)
    ylim([max(1e-12, 10^floor(log10(min(A_pos)))), min(1, 10^ceil(log10(max(A_pos))))]);
end
xlabel('R = \surd Re_x'); ylabel('A');
title(sprintf('Iteration %d, AF = %.3f (%d modes)', Opt.k1, Opt.AF, numel(legs)-1));
legend(legs, 'Location','best'); grid on; box on; drawnow;

end
