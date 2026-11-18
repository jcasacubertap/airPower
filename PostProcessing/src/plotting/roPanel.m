function cb = roPanel(ax, X, Y, F, cl, lab, xl, yl, showCb)
% roPanel  One Reynolds-Orr field panel, in the style shared by
% plotReynoldsOrrProd and plotReynoldsOrrDecomp.
%
%   cb = roPanel(ax, X, Y, F, cl, lab, xl, yl, showCb)
%
% F is drawn as a smooth filled contour clamped to the colour range cl = [lo hi],
% with a diverging blue-white-red map whose WHITE is always at zero (also for an
% asymmetric or single-signed cl). Only the loaded data is drawn: nothing derived
% (no boundary-layer edge or other estimate) is overlaid. Returns the colorbar handle ([] when showCb is false).

    lo = cl(1);  hi = cl(2);
    nLev = 64;
    Fc = min(max(F, lo), hi);                      % clamp so the levels resolve cl
    contourf(ax, X, Y, Fc, linspace(lo, hi, nLev), 'LineStyle', 'none');
    hold(ax, 'on');
    clim(ax, [lo hi]);
    colormap(ax, zeroWhiteMap(lo, hi));

    set(ax, 'TickLabelInterpreter', 'latex', 'Layer', 'top', 'TickDir', 'in', ...
            'FontSize', 11, 'Box', 'on');
    xlim(ax, xl);  ylim(ax, yl);
    xlabel(ax, lab.x, 'Interpreter', 'latex');
    ylabel(ax, lab.y, 'Interpreter', 'latex');

    cb = [];
    if showCb
        cb = colorbar(ax);
        cb.TickLabelInterpreter = 'latex';
        cb.Label.Interpreter    = 'latex';
        cb.TickDirection        = 'out';
    end
end

% Blue-white-red with white at zero: the share of each side follows where zero
% falls in [lo hi], so the zero contour is white whatever the range.
function cmap = zeroWhiteMap(lo, hi)
    n  = 256;
    f0 = min(max(-lo / (hi - lo), 0), 1);          % position of zero in [0 1]
    nb = round(f0 * n);  nr = n - nb;
    blue = [0.13 0.30 0.62];  red = [0.70 0.09 0.11];  white = [1 1 1];
    tb = linspace(0, 1, nb)';  tr = linspace(0, 1, nr)';
    cmap = [blue + tb .* (white - blue); white + tr .* (red - white)];
end
