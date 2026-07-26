function figs = plotPerturbationFields(sPert, sBF, inp, savedir)
% plotPerturbationFields  Perturbation shape functions + amplitude (DeHNSSo).
%
%   figs = plotPerturbationFields(sPert, sBF, inp[, savedir])
%
% Modes to plot are selected by inp.modeIdx (1-based index into the mode array
% betavec = [0, beta0, 2*beta0, ...]; index 2 = fundamental (0,1)). [] -> all.
% Produces:
%   * shape-function contours |q_hat|(x,y) for uhat/vhat/what, one figure per
%     (component, selected mode);
%   * one amplitude figure overlaying the selected modes: max_y|uhat|(x) (top)
%     and the perturbation-energy amplitude max_y 0.5(|uhat|^2+|vhat|^2+|what|^2)
%     (bottom).
% Fields are amplitude-scaled (A folded in by importData). Non-dimensional
% (x by delta_0, velocities by u_infty); plot window via plotWindow / inp.ro.
%
% Save behaviour: savedir non-empty -> PNGs written there; empty/omitted ->
% figures left open for interactive inspection.

    if nargin < 4; savedir = ''; end
    figs = gobjects(0);

    [X, Y] = plotCoords(sBF, inp);  x = X(1, :);   % wall-fitted rectangle (unwraps curved TTCP wall)
    [xl, yl] = plotWindow(X, Y, inp);          % buffer-cut in x, near-wall window in y

    % modes to plot (inp.modeIdx; [] -> all)
    Nmode = size(sPert.u, 1);
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx)
        modes = inp.modeIdx(:).';
    else
        modes = 1:Nmode;
    end
    modes = modes(modes >= 1 & modes <= Nmode);
    multi = numel(modes) > 1;

    % ===== shape-function contours: one figure per (component, mode) =========
    comps = { 'u', '\tilde{u}', 'uhat'; ...
              'v', '\tilde{v}', 'vhat'; ...
              'w', '\tilde{w}', 'what' };
    for c = 1:size(comps, 1)
        fld = comps{c,1};  sym = comps{c,2};  tag = comps{c,3};
        if ~isfield(sPert, fld) || isempty(sPert.(fld)); continue; end
        for mm = modes
            ms = sprintf('(0,%d)', mm - 1);
            Qm = abs(squeeze(sPert.(fld)(mm, :, :)));   % Ny x Nx

            f = figure('Name', sprintf('Shape %s %s', tag, ms), 'Color', 'w', ...
                       'Position', [60 60 900 320]);
            % colour limits from ONLY the plotted window (single-signed magnitude)
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
            xlabel('$S^{*} / \delta_0$', 'Interpreter', 'latex');
            ylabel('$y / \delta_0$', 'Interpreter', 'latex');
            title(sprintf('\\textrm{Perturbation shape}\\ \\ $|%s|_{%s}(x,y)$', sym, ms), ...
                  'Interpreter', 'latex');
            xlim(xl);  ylim(yl);
            figs(end+1) = f; %#ok<AGROW>
            if ~isempty(savedir)
                if ~exist(savedir, 'dir'); mkdir(savedir); end
                if multi
                    out = fullfile(savedir, sprintf('perturbShape_%s_m%d.png', tag, mm-1));
                else
                    out = fullfile(savedir, sprintf('perturbShape_%s.png', tag));
                end
                exportgraphics(f, out, 'Resolution', 150);
                fprintf('plotPerturbationFields: saved %s\n', out);
            end
        end
    end

    % ===== amplitude: max_y|uhat| (top) + max_y Energy (bottom), modes overlaid ==
    fA = figure('Name', 'Perturbation amplitude', 'Color', 'w', 'Position', [60 60 720 620]);
    colors = lines(max(numel(modes), 1));
    leg = arrayfun(@(mm) sprintf('$(0,%d)$', mm-1), modes, 'UniformOutput', false);

    ax1 = subplot(2,1,1); hold(ax1, 'on');
    for i = 1:numel(modes)
        Au = max(abs(squeeze(sPert.u(modes(i),:,:))), [], 1, 'omitnan');
        plot(ax1, x, Au, '-', 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    set(ax1,'YScale','log'); grid(ax1,'on'); box(ax1,'on'); set(ax1,'TickLabelInterpreter','latex');
    ylabel(ax1, '$\max_y\,|\tilde{u}| / u_\infty$', 'Interpreter', 'latex');
    title(ax1, '\textrm{Chordwise amplitude of}\ \ $\tilde{u}$', 'Interpreter', 'latex');
    xlim(ax1, xl);
    if multi; legend(ax1, leg, 'Interpreter', 'latex', 'Location', 'best'); end

    ax2 = subplot(2,1,2); hold(ax2, 'on');
    for i = 1:numel(modes)
        mm = modes(i);
        E = 0.5 * ( abs(squeeze(sPert.u(mm,:,:))).^2 ...
                  + abs(squeeze(sPert.v(mm,:,:))).^2 ...
                  + abs(squeeze(sPert.w(mm,:,:))).^2 );
        AE = max(E, [], 1, 'omitnan');
        plot(ax2, x, AE, '-', 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    set(ax2,'YScale','log'); grid(ax2,'on'); box(ax2,'on'); set(ax2,'TickLabelInterpreter','latex');
    xlabel(ax2, '$S^{*} / \delta_0$', 'Interpreter', 'latex');
    ylabel(ax2, '$\max_y\,\frac{1}{2}(|\tilde{u}|^2+|\tilde{v}|^2+|\tilde{w}|^2) / u_\infty^2$', ...
           'Interpreter', 'latex');
    title(ax2, 'Chordwise amplitude of the perturbation energy', 'Interpreter', 'latex');
    xlim(ax2, xl);
    if multi; legend(ax2, leg, 'Interpreter', 'latex', 'Location', 'best'); end

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
