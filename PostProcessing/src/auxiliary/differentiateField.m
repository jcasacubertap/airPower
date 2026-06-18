function so = differentiateField(si, x, dim, order)
% differentiateField  Numerical derivative of a field along one grid direction.
%
%   so = differentiateField(si, x, dim)
%   so = differentiateField(si, x, dim, order)
%
% Differentiates `si` with respect to the coordinate `x` along dimension `dim`,
% so it applies directly to the Ny x Nx base-flow fields returned by importData
% (loadBF), which follow the convention row 1 = free-stream, row end = wall:
%   du/dx :  differentiateField(sBF.u, sBF.x, 2)   % across columns (x)
%   du/dy :  differentiateField(sBF.u, sBF.y, 1)   % down rows      (y)
%
% Inputs:
%   si     field to differentiate (vector, or Ny x Nx matrix)
%   x      coordinates of si along `dim`. Either a vector along that dimension,
%          or a full matrix the same size as si (e.g. the meshgrid sBF.x/sBF.y),
%          in which case each line uses its own coordinates.
%   dim    1 = down rows, 2 = across columns. Default: first non-singleton dim.
%   order  finite-difference order. Default 2 (only 2 is implemented).
% Output:
%   so     derivative, same size as si.
%
% The scheme (diffOrder2NonUniform) is 2nd-order on non-uniform grids. NaNs
% (e.g. outside a non-rectangular grid) propagate locally through the stencil.

    if nargin < 4 || isempty(order), order = 2; end
    if nargin < 3 || isempty(dim)
        dim = find(size(si) ~= 1, 1);
        if isempty(dim), dim = 1; end
    end

    if order ~= 2
        error('differentiateField:orderNotSupported', ...
              'Only order = 2 is implemented (got %g).', order);
    end
    if ~ismember(dim, [1, 2])
        error('differentiateField:badDim', 'dim must be 1 or 2 (got %g).', dim);
    end
    if isvector(x)
        if numel(x) ~= size(si, dim)
            error('differentiateField:coordLength', ...
                  'coordinate vector length (%d) must equal size(si,%d) = %d.', ...
                  numel(x), dim, size(si, dim));
        end
    elseif ~isequal(size(x), size(si))
        error('differentiateField:coordSize', ...
              'x must be a vector along dim %d or the same size as si.', dim);
    end

    % Orient so the derivative always runs down columns, then restore.
    byColumn = (dim == 2);
    F = si;  C = x;
    if byColumn
        F = F.';
        if ~isvector(C), C = C.'; end
    end

    nLine = size(F, 1);
    if nLine < 3
        error('differentiateField:tooFewPoints', ...
              'Need at least 3 points along dim %d (got %d).', dim, nLine);
    end

    so = F * 0;                  % preserve size, class and complexity
    coordIsVector = isvector(C);
    cVec = C(:);
    for k = 1:size(F, 2)
        if coordIsVector
            xk = cVec;
        else
            xk = C(:, k);
        end
        so(:, k) = diffOrder2NonUniform(F(:, k), xk);
    end

    if byColumn, so = so.'; end
end
