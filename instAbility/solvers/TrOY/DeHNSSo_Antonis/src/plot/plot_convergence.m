function plot_convergence(StabGrid, Mode, Stab, Opt)
%PLOT_CONVERGENCE  Convergence history with optional Newton subplot.

has_newton = isfield(Opt, 'newton_F_history') && ~isempty(Opt.newton_F_history);

if isempty(findobj('Type','figure','Number',3))
    figure(3); set(gcf, 'Visible', 'on');
else
    set(0, 'CurrentFigure', 3); clf;
end

if has_newton
    n_sub = 2;
else
    n_sub = 1;
end

%% Subplot 1: Picard fixed-point residual (dal)
subplot(n_sub, 1, 1);
set(gca, 'YScale', 'log');
hold on;

half = round(Mode.nf/2);
all_phys = half:Mode.nf;
colors = lines(numel(all_phys));
phys = Mode.RunJ(Mode.RunJ >= half);
modeslegend = {};
iters = 1:Opt.k1;

for j = phys
    cidx = j - half + 1;
    semilogy(iters, Opt.dalsave(j, iters), '-o', 'Color', colors(cidx,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4);
    modeslegend{end+1} = sprintf('(m,n) = (%d,%d)', Stab.Mvec(j), Stab.Nvec(j));
end

if Opt.ConvF ~= 1
    semilogy(iters, Opt.Conv*Opt.ConvF .* ones(1,Opt.k1), 'k--', 'LineWidth', 0.8);
    modeslegend{end+1} = 'Intermediate threshold';
end
semilogy(iters, Opt.Conv .* ones(1,Opt.k1), 'k-', 'LineWidth', 1.0);
modeslegend{end+1} = 'Convergence threshold';

% AF change markers
if isfield(Opt, 'AF_history') && numel(Opt.AF_history) >= Opt.k1
    AF_h = Opt.AF_history(1:Opt.k1);
    AF_changes = find(diff(AF_h) ~= 0);
    for ii = 1:numel(AF_changes)
        xline(AF_changes(ii) + 0.5, '--', sprintf('AF=%.2f', AF_h(AF_changes(ii)+1)), ...
            'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'center', 'FontSize', 7);
    end
end

xlim([0 Opt.k1+1])
xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Relative error $\delta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Picard fixed-point residual history', 'Interpreter', 'latex', 'FontSize', 13)
legend(modeslegend, 'Location', 'best', 'FontSize', 9)
grid on; box on;

%% Subplot 2: Newton convergence (||F|| per mode)
if has_newton
    subplot(n_sub, 1, 2);
    set(gca, 'YScale', 'log');
    hold on;

    F_hist  = Opt.newton_F_history;
    AF_hist = Opt.newton_AF_history;
    [nphys, n_pts] = size(F_hist);
    iters_n = 1:n_pts;

    modeslegend_n = {};
    for kk = 1:min(nphys, numel(phys))
        j = phys(kk);
        cidx = j - half + 1;
        semilogy(iters_n, F_hist(kk,:), '-o', 'Color', colors(cidx,:), ...
            'LineWidth', 1.2, 'MarkerSize', 4);
        modeslegend_n{end+1} = sprintf('(m,n) = (%d,%d)', Stab.Mvec(j), Stab.Nvec(j));
    end

    tol = 1e-10;
    if isfield(Opt, 'newton_tol'); tol = Opt.newton_tol; end
    semilogy(iters_n, tol .* ones(1, n_pts), 'k-', 'LineWidth', 1.0);
    modeslegend_n{end+1} = 'Newton tolerance';

    % AF change markers
    AF_changes = find(diff(AF_hist) ~= 0);
    for ii = 1:numel(AF_changes)
        xline(AF_changes(ii) + 0.5, '--', sprintf('AF=%.3f', AF_hist(AF_changes(ii)+1)), ...
            'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'center', 'FontSize', 7);
    end

    xlim([0.5 n_pts+0.5])
    xlabel('Newton sub-iteration', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$\|F_k\|$', 'Interpreter', 'latex', 'FontSize', 12)
    title('Newton convergence history', 'Interpreter', 'latex', 'FontSize', 13)
    legend(modeslegend_n, 'Location', 'best', 'FontSize', 9)
    grid on; box on;
end

set(gcf, 'Position', [100 100 700 250 + 250*n_sub]);
drawnow;

end
