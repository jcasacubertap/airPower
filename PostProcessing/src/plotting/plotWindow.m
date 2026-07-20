function [xl, yl] = plotWindow(X, Y, inp)
% plotWindow  Shared plotting window for base-flow / perturbation plots.
%
%   [xl, yl] = plotWindow(X, Y, inp)
%
% x: cut where the DeHNSSo outflow buffer begins, so the damped tail is not
%    shown -- inp.ro.bufferFrac (default 0.85; 1 -> full domain).
% y: near-wall window measured from the wall (min y) -- inp.ro.yWallFrac
%    (default 0.30; 1 -> full height).
% Rows may run either wall->freestream or freestream->wall; the wall is the
% smallest y value in the grid, so the window is anchored at min(Y).

    bufFrac = 0.85;
    if isfield(inp, 'ro') && isfield(inp.ro, 'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end
    yFrac = 0.30;
    if isfield(inp, 'ro') && isfield(inp.ro, 'yWallFrac') && ~isempty(inp.ro.yWallFrac)
        yFrac = inp.ro.yWallFrac;
    end

    Nx = size(X, 2);
    ib = min(Nx, max(2, round(bufFrac * Nx)));
    xl = [min(X(:)), X(1, ib)];

    y0 = min(Y(:));  y1 = max(Y(:));
    yl = [y0, y0 + yFrac * (y1 - y0)];
end
