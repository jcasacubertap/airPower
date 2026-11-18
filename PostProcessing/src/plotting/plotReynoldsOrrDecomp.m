function plotReynoldsOrrDecomp(sRO, sBF, inp, savedir)
% PLOTREYNOLDSORRDECOMP  Tangential/normal production decomposition I1..I4.
%
%   plotReynoldsOrrDecomp(sRO, sBF, inp, savedir)
%
% Per selected mode, four wide panels stacked: I1 normal-normal, I2
% tangential-normal, I3 normal-tangential, I4 tangential-tangential, with
% I1+I2+I3+I4 = P. One colour scale for all four (roLimits over the four
% together, or fixed in settings.jl), so relative magnitudes are visually
% honest. Same style as plotReynoldsOrrProd (roPanel).
% Saves reynoldsOrrDecomp.png (one file per mode when several are selected).

    if nargin < 4 || isempty(savedir)
        savedir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'plotting');
    end
    if ~exist(savedir, 'dir'); mkdir(savedir); end

    [X, Y, lab] = plotCoords(sBF, inp);   % same grid as sRO.s/n; x,y when Cartesian
    % which computed modes to plot: inp.modeIdx (original sPert indices); [] -> all
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx) && isfield(sRO,'modeNo')
        sel = find(ismember(sRO.modeNo, inp.modeIdx(:).'));
    else
        sel = 1:size(sRO.I1, 3);
    end
    if isempty(sel); sel = 1:size(sRO.I1, 3); end

    % Same shared window as plotReynoldsOrrProd and every other figure.
    [xl, yl] = plotWindow(X, Y, inp);
    inwin = X >= xl(1) & X <= xl(2) & Y >= yl(1) & Y <= yl(2);

    flds = {'I1', 'I2', 'I3', 'I4'};
    ttl  = {'$I_1$ \ (normal--normal)', '$I_2$ \ (tangential--normal)', ...
            '$I_3$ \ (normal--tangential)', '$I_4$ \ (tangential--tangential)'};

    for m = sel
        F = cellfun(@(f) sRO.(f)(:, :, m), flds, 'UniformOutput', false);
        [cl, e] = roLimits(F, inwin, inp);

        fig = figure('Position', [40 40 900 760], 'Color', 'w');
        tl  = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        for c = 1:4
            ax = nexttile(tl);
            roPanel(ax, X, Y, F{c}, cl, lab, xl, yl, false);
            title(ax, ttl{c}, 'Interpreter', 'latex');
            if c < 4; xlabel(ax, ''); end
        end
        cb = colorbar(ax);  cb.Layout.Tile = 'east';
        cb.TickLabelInterpreter = 'latex';  cb.Label.Interpreter = 'latex';  cb.TickDirection = 'out';
        if e ~= 0
            cb.Ruler.Exponent = 0;
            cb.TickLabels = arrayfun(@(t) sprintf('$%.3g$', t/10^e), cb.Ticks, 'UniformOutput', false);
            cb.Label.String = sprintf('$I_n / u_\\infty^3 \\;(\\times 10^{%d})$', e);
        else
            cb.Label.String = '$I_n / u_\infty^3$';
        end

        if isfield(sRO,'modeNo'); nlab = sRO.modeNo(m) - 1; else; nlab = m; end
        title(tl, sprintf('\\textrm{Production decomposition, mode}~$(0,%d)$', nlab), ...
              'Interpreter', 'latex');

        if numel(sel) == 1
            name = 'reynoldsOrrDecomp.png';
        else
            name = sprintf('reynoldsOrrDecomp_m%d.png', nlab);
        end
        out = fullfile(savedir, name);
        exportgraphics(fig, out, 'Resolution', 200);
        fprintf('plotReynoldsOrrDecomp: saved %s\n', out);
    end
end
