function [X, Y, lab] = plotCoords(sBF, inp) %#ok<INUSD>
% plotCoords  Plotting coordinates (and their axis labels) for field / profile plots.
%
%   [X, Y]      = plotCoords(sBF, inp)
%   [X, Y, lab] = plotCoords(sBF, inp)
%
% Returns Ny x Nx coordinate matrices to plot fields against, chosen from the
% grid itself (no case or source tag), so it works for any dataset importData
% returns:
%
%   Cartesian grid (flat wall along X, every grid column along Y):
%       X = sBF.geom.X           the true streamwise position
%       Y = sBF.geom.Y - Y_wall  height above the wall (= sBF.n)
%   anything else (curved / leaning grid):
%       X = sBF.s   streamwise arc length along the wall, from the base-flow source
%       Y = sBF.n   distance up each grid column from the wall
%   The wall-fitted pair unwraps a curved wall into a rectangle, so the boundary
%   layer fills the plot rather than a thin curved sliver. On a Cartesian grid it
%   would carry the same lines, but with an arbitrary origin for s.
%
% Y is a wall distance in both cases (the wall row is Y = 0), so windows and
% wall-distance quantities (PIV profiles, valYTop) mean the same either way.
%
% lab (all LaTeX, units from sBF.units: delta_0 for 'nondim', m for 'SI'):
%   .cartesian  true for the first case
%   .xs, .ys    the bare symbols: 'x','y' or 'S^{*}','n'
%   .unit       the unit suffix, '/ \delta_0' or '\; [\mathrm{m}]'
%   .x, .y      complete axis labels, e.g. '$x / \delta_0$'
%
% Row convention (set by importData): (1,:) = free-stream, (end,:) = wall.

    if ~(isfield(sBF,'s') && isfield(sBF,'n') && isequal(size(sBF.s), size(sBF.n)))
        error('plotCoords:noWallFitted', ...
              'sBF has no wall-fitted (s,n) pair — load it through importData.');
    end

    lab.cartesian = isCartesian(sBF);
    if lab.cartesian
        X = sBF.geom.X;
        Y = sBF.geom.Y - sBF.geom.Y(end, :);
        lab.xs = 'x';      lab.ys = 'y';
    else
        X = sBF.s;  Y = sBF.n;
        lab.xs = 'S^{*}';  lab.ys = 'n';
    end

    if isfield(sBF, 'units') && strcmpi(sBF.units, 'SI')
        lab.unit = '\; [\mathrm{m}]';
    else
        lab.unit = '/ \delta_0';
    end
    lab.x = sprintf('$%s %s$', lab.xs, lab.unit);
    lab.y = sprintf('$%s %s$', lab.ys, lab.unit);
end

% Flat wall along X and straight grid columns along Y, to round-off relative to
% the domain length (the OpenFOAM midPlane carries ~1e-10 of noise).
function tf = isCartesian(sBF)
    tf = false;
    if ~(isfield(sBF, 'geom') && all(isfield(sBF.geom, {'X', 'Y'}))); return; end
    X = sBF.geom.X;  Y = sBF.geom.Y;
    if ~isequal(size(X), size(sBF.s)); return; end
    L = max(X(:)) - min(X(:));
    if ~(L > 0); return; end
    tol = 1e-6 * L;
    tf = max(Y(end, :)) - min(Y(end, :)) < tol ...       % wall row at constant Y
      && max(max(X, [], 1) - min(X, [], 1)) < tol;       % each column at constant X
end
