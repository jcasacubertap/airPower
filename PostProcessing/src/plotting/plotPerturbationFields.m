function figs = plotPerturbationFields(sPert, sBF, inp, savedir)
% plotPerturbationFields  Perturbation shape functions + amplitude (DeHNSSo).
%
%   figs = plotPerturbationFields(sPert, sBF, inp)           % show only
%   figs = plotPerturbationFields(sPert, sBF, inp, savedir)  % also save PNGs
%
% For the fundamental mode (0,1), produces:
%   * three shape-function figures (uhat, vhat, what): y-x plane contours of the
%     amplitude-scaled magnitude |q_hat|(x,y) (NOT peak-normalized -- the stored
%     amplitude A is folded in by importData);
%   * one amplitude figure with two stacked panels: max_y|uhat|(x) (top) and the
%     perturbation-energy amplitude max_y E(x), E = 0.5(|uhat|^2+|vhat|^2+|what|^2)
%     (bottom).
% Quantities are non-dimensional (x by delta_0, velocities by u_infty). The x-axis
% is cut at the outflow buffer and the contour y-axis is a near-wall window (both
% via plotWindow / inp.ro).
%
% Save behaviour: savedir non-empty -> PNGs written there (terminal/dispatcher);
% empty/omitted -> figures left open for interactive inspection (direct MATLAB).

    if nargin < 4; savedir = ''; end
    figs = gobjects(0);

    X = sBF.x;  Y = sBF.y;                  % Ny x Nx, non-dimensional (x/delta_0, y/delta_0)
    x = sBF.x(1, :);                        % 1 x Nx
    [xl, yl] = plotWindow(X, Y, inp);       % buffer-cut in x, near-wall window in y

    % fundamental mode (0,1): second entry of betavec = [0, beta_0, ...]
    Nmode = size(sPert.u, 1);
    m  = min(2, Nmode);
    ms = sprintf('(0,%d)', m - 1);          % mode label for the amplitude-function subscript

    % ===== (1-3) shape-function contours, one figure per component ==========
    comps = { 'u', '\tilde{u}', 'uhat'; ...
              'v', '\tilde{v}', 'vhat'; ...
              'w', '\tilde{w}', 'what' };
    for c = 1:size(comps, 1)
        fld = comps{c,1};  sym = comps{c,2};  tag = comps{c,3};
        if ~isfield(sPert, fld) || isempty(sPert.(fld)); continue; end
        Qm = abs(squeeze(sPert.(fld)(m, :, :)));   % Ny x Nx

        f = figure('Name', sprintf('Shape %s', tag), 'Color', 'w', ...
                   'Position', [60 60 900 320]);
        % Colour limits from ONLY the plotted window (bufferFrac x-cut x near-wall
        % y-cut), so the visible structure is resolved rather than compressed by a
        % peak sitting in the outflow buffer; |q| is a magnitude, so single-signed.
        inwin = (X >= xl(1) & X <= xl(2) & Y >= yl(1) & Y <= yl(2) & isfinite(Qm));
        Qw = Qm(inwin); if isempty(Qw); Qw = Qm(isfinite(Qm)); end
        lo = min(Qw); hi = max(Qw);
        if ~(hi > lo); hi = lo + eps(max(abs(lo),1)); end
        Qc = min(max(Qm, lo), hi);                 % clamp; out-of-window saturates
        contourf(X, Y, Qc, 30, 'LineStyle', 'none');
        clim([lo hi]);
        colormap(gca, parula);
        cb = colorbar; cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        local_sciColorbar(cb, sprintf('$|%s|_{%s} / u_\\infty$', sym, ms));
        set(gca, 'TickLabelInterpreter', 'latex');
        xlabel('$x / \delta_0$', 'Interpreter', 'latex');
        ylabel('$y / \delta_0$', 'Interpreter', 'latex');
        title(sprintf('\\textrm{Perturbation shape}\\ \\ $|%s|_{%s}(x,y)$', sym, ms), ...
              'Interpreter', 'latex');
        xlim(xl);  ylim(yl);
        figs(end+1) = f; %#ok<AGROW>
        if ~isempty(savedir)
            if ~exist(savedir, 'dir'); mkdir(savedir); end
            out = fullfile(savedir, sprintf('perturbShape_%s.png', tag));
            exportgraphics(f, out, 'Resolution', 150);
            fprintf('plotPerturbationFields: saved %s\n', out);
        end
    end

    % ===== (4) amplitude: max_y|uhat| (top) and max_y Energy (bottom) =======
    fA = figure('Name', 'Perturbation amplitude', 'Color', 'w', 'Position', [60 60 720 620]);

    Au = max(abs(squeeze(sPert.u(m,:,:))), [], 1, 'omitnan');
    ax1 = subplot(2,1,1);
    plot(ax1, x, Au, '-', 'LineWidth', 1.5);
    set(ax1,'YScale','log');   % after plot: with hold off, plot() resets YScale
    grid(ax1,'on'); box(ax1,'on'); set(ax1,'TickLabelInterpreter','latex');
    ylabel(ax1, sprintf('$\\max_y\\,|\\tilde{u}|_{%s} / u_\\infty$', ms), 'Interpreter', 'latex');
    title(ax1, sprintf('\\textrm{Chordwise amplitude of}\\ \\ $\\tilde{u}_{%s}$', ms), 'Interpreter', 'latex');
    xlim(ax1, xl);

    E = 0.5 * ( abs(squeeze(sPert.u(m,:,:))).^2 ...
              + abs(squeeze(sPert.v(m,:,:))).^2 ...
              + abs(squeeze(sPert.w(m,:,:))).^2 );     % Ny x Nx
    AE = max(E, [], 1, 'omitnan');
    ax2 = subplot(2,1,2);
    plot(ax2, x, AE, '-', 'LineWidth', 1.5);
    set(ax2,'YScale','log');   % after plot: with hold off, plot() resets YScale
    grid(ax2,'on'); box(ax2,'on'); set(ax2,'TickLabelInterpreter','latex');
    xlabel(ax2, '$x / \delta_0$', 'Interpreter', 'latex');
    ylabel(ax2, sprintf(['$\\max_y\\,\\frac{1}{2}\\left(|\\tilde{u}|_{%s}^2+|\\tilde{v}|_{%s}^2' ...
           '+|\\tilde{w}|_{%s}^2\\right) / u_\\infty^2$'], ms, ms, ms), 'Interpreter', 'latex');
    title(ax2, 'Chordwise amplitude of the perturbation energy', 'Interpreter', 'latex');
    xlim(ax2, xl);

    figs(end+1) = fA;
    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, 'perturbAmplitude.png');
        exportgraphics(fA, out, 'Resolution', 150);
        fprintf('plotPerturbationFields: saved %s\n', out);
    end
end

% --- factor the order of magnitude out of the tick labels into the label ---
function local_sciColorbar(cb, baseLabel)
    ticks = cb.Ticks;
    m = max(abs(ticks));
    if ~(m > 0); cb.Label.String = baseLabel; return; end
    e = floor(log10(m));
    if e ~= 0
        cb.Ruler.Exponent = 0;
        cb.TickLabels = arrayfun(@(t) sprintf('$%.2g$', t/10^e), cb.Ticks, ...
                                 'UniformOutput', false);
        cb.Label.String = sprintf('%s $\\;(\\times 10^{%d})$', baseLabel, e);
    else
        cb.Label.String = baseLabel;
    end
end
