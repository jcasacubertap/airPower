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
    % which computed modes to plot: inp.modeIdx (original sPert indices); [] -> all
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx) && isfield(sRO,'modeNo')
        sel = find(ismember(sRO.modeNo, inp.modeIdx(:).'));
    else
        sel = 1:size(P, 3);
    end
    if isempty(sel); sel = 1:size(P, 3); end
    M = numel(sel);
    nc = min(M, 3);  nr = ceil(M / nc);

    % near-wall focus: rows run free-stream(1,:) -> wall(end,:), so small y is
    % the wall; show the lower part of the domain where production concentrates.
    % y-axis top: inp.ro.yMax (absolute) if set, else 30% of the domain height.
    y0 = min(Y(:));  y1 = max(Y(:));
    if isfield(inp,'ro') && isfield(inp.ro,'yMax') && ~isempty(inp.ro.yMax)
        yTop = y0 + inp.ro.yMax;
    else
        yTop = y0 + 0.30 * (y1 - y0);
    end

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
    for p = 1:M
        m = sel(p);                  % index into sRO.P (computed mode)
        subplot(nr, nc, p);
        Pm = P(:, :, m);

        % Robust symmetric colour limit (manual 99th pct -> no toolbox needed),
        % from ONLY the plotted window (x up to the buffer, near-wall y) so the
        % colours fit the chosen region instead of being set by out-of-window peaks.
        inwin = (X >= x0 & X <= xEnd & Y >= y0 & Y <= yTop & isfinite(Pm));
        a = sort(abs(Pm(inwin)));
        if isempty(a); a = sort(abs(Pm(isfinite(Pm)))); end
        if isempty(a) || a(end) == 0
            cl = 1;
        else
            cl = a(max(1, round(0.99 * numel(a))));
            if ~(cl > 0); cl = a(end); end
        end

        % Odd number of bands with symmetric edges -> a central band straddles 0
        % and maps to WHITE, so very small (near-zero) values read as white rather
        % than a faint blue/lavender tint. Fewer partitions than before (15 vs 30).
        Pm = min(max(Pm, -cl), cl);    % clamp so the bands resolve the window
        nb = 15;
        contourf(X, Y, Pm, linspace(-cl, cl, nb+1), 'LineStyle', 'none'); hold on;
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
            cb.Label.String = sprintf('$P / u_\\infty^3 \\; (\\times 10^{%d})$', e);
        else
            cb.Label.String = '$P / u_\infty^3$';
        end

        % base-flow U reference (thin grey lines)
        if isfield(sBF, 'u')
            contour(X, Y, sBF.u, 6, 'LineColor', [0.4 0.4 0.4], 'LineWidth', 0.4);
        end

        set(gca, 'TickLabelInterpreter', 'latex');
        xlabel('$x / \delta_0$', 'Interpreter', 'latex');
        ylabel('$y / \delta_0$', 'Interpreter', 'latex');
        xlim([x0 xEnd]);
        ylim([y0 yTop]);
        if isfield(sRO, 'modeNo')
            title(sprintf('$(0,%d)$', sRO.modeNo(m) - 1), 'Interpreter', 'latex');
        end
    end
    sgtitle('\textrm{Total production}', 'Interpreter', 'latex');

    out = fullfile(savedir, 'reynoldsOrrProd.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('plotReynoldsOrrProd: saved %s\n', out);
end

% --- blue-grey-red diverging colormap with a central light-GREY plateau ---
% Zero maps to light grey (not white) so it is distinct from the unfilled hump
% interior (below the wall), which shows the white axes background. The plateau
% (central `w` fraction) also defeats contourf's lower-edge band colouring, so
% the zero-straddling band renders grey rather than a faint tint.
function cmap = local_diverging()
    n = 256;  w = 0.16;  g = 0.93;           % g = grey level for zero
    h = round((1 - w) * n / 2);              % coloured rows on each side
    t = linspace(0, 1, h)';
    cmap = g * ones(n, 3);                    % middle rows = light grey
    cmap(1:h, :)     = [t*g, t*g, 1 - t*(1-g)];             % blue [0 0 1] -> grey
    cmap(n-h+1:n, :) = [g + t*(1-g), g*(1-t), g*(1-t)];     % grey -> red [1 0 0]
end
