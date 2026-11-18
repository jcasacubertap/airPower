function plotReynoldsOrrProd(sRO, sBF, inp, savedir)
% PLOTREYNOLDSORRPROD  Contour of Reynolds-Orr production P(x,y) per mode.
%
%   plotReynoldsOrrProd(sRO, sBF, inp, savedir)
%
% One wide panel per selected mode (inp.modeIdx), stacked: filled contour of the
% production field sRO.P(:,:,m) on the plotCoords grid (x,y or s,n), like every
% other field figure. Diverging map with white at zero (red = production of
% perturbation energy, blue = its removal). Colour range from roLimits: automatic (symmetric,
% shared by the panels) or fixed in settings.jl (PostProcessing.ro.climManual /
% clim). Saves reynoldsOrrProd.png in savedir.

    if nargin < 4 || isempty(savedir)
        savedir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'plotting');
    end
    if ~exist(savedir, 'dir'); mkdir(savedir); end

    [X, Y, lab] = plotCoords(sBF, inp);   % same grid as sRO.s/n; x,y when Cartesian
    P = sRO.P;
    % which computed modes to plot: inp.modeIdx (original sPert indices); [] -> all
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx) && isfield(sRO,'modeNo')
        sel = find(ismember(sRO.modeNo, inp.modeIdx(:).'));
    else
        sel = 1:size(P, 3);
    end
    if isempty(sel); sel = 1:size(P, 3); end
    M = numel(sel);

    % Window from the shared definition (near-wall y, x cut at the buffer start),
    % so this figure frames the same region as every other one.
    [xl, yl] = plotWindow(X, Y, inp);
    inwin = X >= xl(1) & X <= xl(2) & Y >= yl(1) & Y <= yl(2);
    [cl, e] = roLimits(arrayfun(@(m) P(:,:,m), sel, 'UniformOutput', false), inwin, inp);

    fig = figure('Position', [60 60 900 70 + 230*M], 'Color', 'w');
    tl  = tiledlayout(fig, M, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:M
        m  = sel(p);
        ax = nexttile(tl);
        roPanel(ax, X, Y, P(:,:,m), cl, lab, xl, yl, false);
        if isfield(sRO, 'modeNo')
            title(ax, sprintf('$P$, \\ mode $(0,%d)$', sRO.modeNo(m) - 1), 'Interpreter', 'latex');
        end
        if p < M; xlabel(ax, ''); end
    end
    cb = colorbar(ax);  cb.Layout.Tile = 'east';
    cb.TickLabelInterpreter = 'latex';  cb.Label.Interpreter = 'latex';  cb.TickDirection = 'out';
    sciColorbar(cb, '$P / u_\infty^3$', e);

    out = fullfile(savedir, 'reynoldsOrrProd.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('plotReynoldsOrrProd: saved %s\n', out);
end

% factor the power of ten out of the tick labels into the colorbar label
function sciColorbar(cb, base, e)
    if e ~= 0
        cb.Ruler.Exponent = 0;
        cb.TickLabels = arrayfun(@(t) sprintf('$%.3g$', t/10^e), cb.Ticks, 'UniformOutput', false);
        cb.Label.String = sprintf('%s $\\;(\\times 10^{%d})$', base, e);
    else
        cb.Label.String = base;
    end
end
