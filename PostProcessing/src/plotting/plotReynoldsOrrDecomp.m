function plotReynoldsOrrDecomp(sRO, sBF, inp, savedir)
% PLOTREYNOLDSORRDECOMP  Tangential/normal production decomposition I1..I4.
%
%   plotReynoldsOrrDecomp(sRO, sBF, inp, savedir)
%
% 2x2 layout per mode (I1 I2 / I3 I4; I1 normal-normal, I2 tang-normal,
% I3 normal-tang, I4 tang-tang, with I1+I2+I3+I4 = P). All four panels share a
% single colorbar whose limit is led by I2 (the dominant term), so relative
% magnitudes are visually honest. Same styling as plotReynoldsOrrProd (diverging
% map, LaTeX, buffer-limited x, near-wall y, base-flow U overlay).
% Saves reynoldsOrrDecomp.png (one file per mode when several are processed).

    if nargin < 4 || isempty(savedir)
        savedir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'plotting');
    end
    if ~exist(savedir, 'dir'); mkdir(savedir); end

    X = sRO.x;  Y = sRO.y;
    % which computed modes to plot: inp.modeIdx (original sPert indices); [] -> all
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx) && isfield(sRO,'modeNo')
        sel = find(ismember(sRO.modeNo, inp.modeIdx(:).'));
    else
        sel = 1:size(sRO.I1, 3);
    end
    if isempty(sel); sel = 1:size(sRO.I1, 3); end

    % near-wall y-focus + x-axis ending at the buffer start (as in plotReynoldsOrrProd)
    y0 = min(Y(:));
    if isfield(inp,'ro') && isfield(inp.ro,'yMax') && ~isempty(inp.ro.yMax)
        yTop = y0 + inp.ro.yMax;
    else
        yTop = y0 + 0.30 * (max(Y(:)) - y0);
    end
    bufFrac = 0.85;
    if isfield(inp, 'ro') && isfield(inp.ro, 'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end
    Nx = size(X, 2);  ib = min(Nx, max(2, round(bufFrac * Nx)));
    x0 = min(X(:));   xEnd = X(1, ib);

    flds = {'I1', 'I2', 'I3', 'I4'};      % tiles: I1 I2 / I3 I4
    labs = {'I_1', 'I_2', 'I_3', 'I_4'};

    inwin = (X >= x0 & X <= xEnd & Y >= y0 & Y <= yTop);   % plotted window

    for m = sel
        % single shared colour limit, led by I2 (dominant term); robust 99th pct
        % from ONLY the plotted window so the colours fit the chosen region.
        I2m = sRO.I2(:, :, m);
        a = sort(abs(I2m(inwin & isfinite(I2m))));
        if isempty(a); a = sort(abs(I2m(isfinite(I2m)))); end
        if isempty(a) || a(end) == 0
            cl = 1;
        else
            cl = a(max(1, round(0.99 * numel(a))));
            if ~(cl > 0); cl = a(end); end
        end
        lvls = linspace(-cl, cl, 16);   % 15 bands (odd) -> central band = white at 0
        e = 0;  if cl > 0; e = floor(log10(cl)); end

        fig = figure('Position', [40 40 720 540], 'Color', 'w');
        tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        axLast = [];
        for c = 1:4
            ax = nexttile(tl);  axLast = ax;
            Fc = min(max(sRO.(flds{c})(:, :, m), -cl), cl);   % clamp -> window colours
            contourf(ax, X, Y, Fc, lvls, 'LineStyle', 'none');
            hold(ax, 'on');
            clim(ax, [-cl cl]);
            colormap(ax, local_div());
            if isfield(sBF, 'u')
                contour(ax, X, Y, sBF.u, 6, 'LineColor', [0.4 0.4 0.4], 'LineWidth', 0.4);
            end
            set(ax, 'TickLabelInterpreter', 'latex');
            xlim(ax, [x0 xEnd]);  ylim(ax, [y0 yTop]);
            title(ax, sprintf('$%s$', labs{c}), 'Interpreter', 'latex');
            if c == 1 || c == 3; ylabel(ax, '$y / \delta_0$', 'Interpreter', 'latex'); end
            if c == 3 || c == 4; xlabel(ax, '$x / \delta_0$', 'Interpreter', 'latex'); end
        end

        % one shared colorbar (LaTeX scientific notation), scale led by I2
        cb = colorbar(axLast);  cb.Layout.Tile = 'east';
        cb.TickLabelInterpreter = 'latex';  cb.Label.Interpreter = 'latex';
        if e ~= 0
            cb.Ruler.Exponent = 0;
            cb.TickLabels = arrayfun(@(t) sprintf('$%.2g$', t/10^e), cb.Ticks, ...
                                     'UniformOutput', false);
            cb.Label.String = sprintf('$I_n / u_\\infty^3 \\; (\\times 10^{%d})$', e);
        else
            cb.Label.String = '$I_n / u_\infty^3$';
        end

        if isfield(sRO,'modeNo'); nlab = sRO.modeNo(m) - 1; else; nlab = m; end
        title(tl, sprintf(['\\textrm{Production decomposition (tangential / normal)}' ...
              '\\ \\ mode $(0,%d)$'], nlab), 'Interpreter', 'latex');

        if numel(sel) == 1
            name = 'reynoldsOrrDecomp.png';
        else
            name = sprintf('reynoldsOrrDecomp_m%d.png', nlab);
        end
        out = fullfile(savedir, name);
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('plotReynoldsOrrDecomp: saved %s\n', out);
    end
end

% --- blue-grey-red diverging colormap with a central light-GREY plateau ---
% Zero = light grey (distinct from the unfilled hump interior = white axes bg);
% plateau also keeps the zero-straddling contourf band grey (not tinted).
function cmap = local_div()
    n = 256;  w = 0.16;  g = 0.93;
    h = round((1 - w) * n / 2);
    t = linspace(0, 1, h)';
    cmap = g * ones(n, 3);
    cmap(1:h, :)     = [t*g, t*g, 1 - t*(1-g)];
    cmap(n-h+1:n, :) = [g + t*(1-g), g*(1-t), g*(1-t)];
end
