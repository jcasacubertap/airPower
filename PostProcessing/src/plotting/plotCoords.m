function [X, Y] = plotCoords(sBF, inp)
% plotCoords  Rectangular wall-fitted coordinates for 2D field / profile plots.
%
%   [X, Y] = plotCoords(sBF, inp)
%
% Returns Ny x Nx coordinate matrices to plot fields against, chosen so the
% boundary layer fills a rectangle rather than a thin curved sliver:
%
%   DFP  — sBF.x/sBF.y are ALREADY wall-fitted (x = arc-length from the inlet,
%          y = wall-normal distance, both in delta_0). Passed through unchanged.
%
%   TTCP — sBF.x/sBF.y are PHYSICAL Cartesian coordinates on a body-fitted curved
%          grid (in delta_0). Plotting them directly crams the BL into a diagonal
%          wedge. Here they are unwrapped into
%              X = streamwise arc-length s along the wall  (delta_0)
%              Y = wall-normal distance n up each grid column from the wall (delta_0)
%          so the curved BL maps onto a clean rectangle. The physical arrays in
%          sBF are left untouched (the PIV validation still needs them).
%
% Row convention (set by importData): (1,:) = free-stream, (end,:) = wall.

    if ~isfield(inp,'caseType') || ~strcmpi(inp.caseType, 'TTCP')
        X = sBF.x;  Y = sBF.y;   % DFP (and any already-wall-fitted grid)
        return;
    end

    Xc = sBF.x;  Yc = sBF.y;      % physical, delta_0 units, wall at row end
    [Ny, Nx] = size(Xc);

    % wall-normal distance up each column, measured from the wall (row end)
    seg = hypot(diff(Xc, 1, 1), diff(Yc, 1, 1));      % (Ny-1) x Nx, between-row lengths
    Y   = [flipud(cumsum(flipud(seg), 1)); zeros(1, Nx)];   % Ny x Nx, 0 at wall (row end)

    % streamwise arc-length along the wall (row end), broadcast down the rows
    sw = [0, cumsum(hypot(diff(Xc(end,:)), diff(Yc(end,:))))];   % 1 x Nx
    X  = repmat(sw, Ny, 1);
end
