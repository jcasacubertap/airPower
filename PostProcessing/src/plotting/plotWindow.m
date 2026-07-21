function [xl, yl] = plotWindow(X, Y, inp)
% plotWindow  Shared plotting window for base-flow / perturbation plots.
%
%   [xl, yl] = plotWindow(X, Y, inp)
%
% x: cut where the DeHNSSo outflow buffer begins, so the damped tail is not
%    shown -- inp.ro.bufferFrac (default 0.85; 1 -> full domain).
% y: from the wall (min y) up to inp.ro.yMax (absolute, in the plot's y units);
%    if yMax is unset, a near-wall fraction inp.ro.yWallFrac (default 0.30) of the
%    domain height is used. Rows may run either way; the wall is min(Y).

    bufFrac = 0.85;
    if isfield(inp, 'ro') && isfield(inp.ro, 'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end

    Nx = size(X, 2);
    ib = min(Nx, max(2, round(bufFrac * Nx)));
    xl = [min(X(:)), X(1, ib)];

    y0 = min(Y(:));  y1 = max(Y(:));
    if isfield(inp, 'ro') && isfield(inp.ro, 'yMax') && ~isempty(inp.ro.yMax)
        yl = [y0, y0 + inp.ro.yMax];              % absolute y-axis max
    else
        yFrac = 0.30;
        if isfield(inp, 'ro') && isfield(inp.ro, 'yWallFrac') && ~isempty(inp.ro.yWallFrac)
            yFrac = inp.ro.yWallFrac;
        end
        yl = [y0, y0 + yFrac * (y1 - y0)];
    end
end
