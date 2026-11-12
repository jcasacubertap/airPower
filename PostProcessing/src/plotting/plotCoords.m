function [X, Y] = plotCoords(sBF, inp) %#ok<INUSD>
% plotCoords  Rectangular wall-fitted coordinates for 2D field / profile plots.
%
%   [X, Y] = plotCoords(sBF, inp)
%
% Returns Ny x Nx coordinate matrices to plot fields against, chosen so the
% boundary layer fills a rectangle rather than a thin curved sliver:
%   X = sBF.s   streamwise arc length along the wall, from the base-flow source
%   Y = sBF.n   distance up each grid column from the wall
% Both loaders return this pair (see importData), in their own units. The
% physical Cartesian coordinates stay in sBF.geom for whatever needs them (the
% PIV validation).
%
% Row convention (set by importData): (1,:) = free-stream, (end,:) = wall.

    if ~(isfield(sBF,'s') && isfield(sBF,'n') && isequal(size(sBF.s), size(sBF.n)))
        error('plotCoords:noWallFitted', ...
              'sBF has no wall-fitted (s,n) pair — load it through importData.');
    end
    X = sBF.s;  Y = sBF.n;
end
