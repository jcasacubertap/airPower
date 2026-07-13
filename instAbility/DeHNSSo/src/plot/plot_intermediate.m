function plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt)
%PLOT_INTERMEDIATE  Plot amplitude evolution during convergence loop.

if isempty(findobj('Type','figure','Number',2))
    figure(2); set(gcf, 'Visible', 'on');
else
    set(0, 'CurrentFigure', 2); clf;
end

% R = sqrt(Re_x) = sqrt(x_norm * R0), with R0 = sqrt(x_norm(1))
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
