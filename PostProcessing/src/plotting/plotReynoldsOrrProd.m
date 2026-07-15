function plotReynoldsOrrProd(sRO, sBF, inp, savedir)
% PLOTREYNOLDSORRPROD  Contour of Reynolds-Orr production P(x,y) per mode.
%
%   plotReynoldsOrrProd(sRO, sBF, inp, savedir)
%
% One panel per processed mode: filled contour of the production field
% sRO.P(:,:,m) on the (x,y) stability grid. Diverging colormap centred at 0
% (red = production of perturbation energy, blue = dissipation). A base-flow
% U contour is overlaid for reference. Saves reynoldsOrrProd.png in savedir.

    if nargin < 4 || isempty(savedir)
        savedir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'plotting');
    end
    if ~exist(savedir, 'dir'); mkdir(savedir); end

    X = sRO.x;  Y = sRO.y;  P = sRO.P;
    M = size(P, 3);
    nc = min(M, 3);  nr = ceil(M / nc);

    % near-wall focus: rows run free-stream(1,:) -> wall(end,:), so small y is
    % the wall; show the lower part of the domain where production concentrates.
    y0 = min(Y(:));  y1 = max(Y(:));
    yTop = y0 + 0.30 * (y1 - y0);

    % x-axis ends where the DeHNSSo outflow buffer begins (default 85% of the
    % domain, matching Opt.xb=85), so the damped tail is not shown. Override
    % via inp.ro.bufferFrac (1 -> show the full domain).
    bufFrac = 0.85;
    if isfield(inp, 'ro') && isfield(inp.ro, 'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end
    Nx = size(X, 2);
    ib = min(Nx, max(2, round(bufFrac * Nx)));
    x0 = min(X(:));  xEnd = X(1, ib);

    fig = figure('Position', [60 60 500*nc 380*nr], 'Color', 'w');
    for m = 1:M
        subplot(nr, nc, m);
        Pm = P(:, :, m);

        % robust symmetric colour limit (manual 99th pct -> no toolbox needed)
        a = sort(abs(Pm(isfinite(Pm))));
        if isempty(a) || a(end) == 0
            cl = 1;
        else
            cl = a(max(1, round(0.99 * numel(a))));
            if ~(cl > 0); cl = a(end); end
        end

        contourf(X, Y, Pm, 30, 'LineStyle', 'none'); hold on;
        clim([-cl cl]);
        colormap(gca, local_diverging());
        cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.Interpreter    = 'latex';
        % clean LaTeX scientific notation: factor the order of magnitude out of
        % the tick labels and into the colorbar label ($P\;(\times 10^{e})$).
        e = 0;
        if cl > 0; e = floor(log10(cl)); end
        if e ~= 0
            cb.Ruler.Exponent = 0;                       % suppress MATLAB's own x10^n
            cb.TickLabels = arrayfun(@(t) sprintf('$%.2g$', t/10^e), cb.Ticks, ...
                                     'UniformOutput', false);
            cb.Label.String = sprintf('$P \\; (\\times 10^{%d})$', e);
        else
            cb.Label.String = '$P$';
        end

        % base-flow U reference (thin grey lines)
        if isfield(sBF, 'u')
            contour(X, Y, sBF.u, 6, 'LineColor', [0.4 0.4 0.4], 'LineWidth', 0.4);
        end

        set(gca, 'TickLabelInterpreter', 'latex');
        xlabel('$x$', 'Interpreter', 'latex');
        ylabel('$y$', 'Interpreter', 'latex');
        xlim([x0 xEnd]);
        ylim([y0 yTop]);

    end
    sgtitle('\textrm{Total production}', 'Interpreter', 'latex');

    out = fullfile(savedir, 'reynoldsOrrProd.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('plotReynoldsOrrProd: saved %s\n', out);
end

% --- blue-white-red diverging colormap (no toolbox dependency) ---
function cmap = local_diverging()
    n = 128;
    t = linspace(0, 1, n)';
    b2w = [t,        t,        ones(n,1)];   % blue  -> white
    w2r = [ones(n,1), flipud(t), flipud(t)]; % white -> red
    cmap = [b2w; w2r];
end
