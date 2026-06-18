function so = diffOrder2NonUniform(si, x)
% diffOrder2NonUniform  Second-order first derivative on a non-uniform 1-D grid.
%
%   so = diffOrder2NonUniform(si, x)
%
% First derivative d(si)/dx of the vector `si` sampled at coordinates `x`, with
% a uniform second-order truncation error on arbitrarily-spaced grids:
%   - 2nd-order central  difference at interior points
%   - 2nd-order forward  difference at the first point
%   - 2nd-order backward difference at the last  point
%
% Inputs:
%   si   values, vector of length N >= 3
%   x    grid coordinates matching si (need not be uniform)
% Output:
%   so   derivative at the coordinates x, same shape as si

    n = numel(si);
    if numel(x) ~= n
        error('diffOrder2NonUniform:sizeMismatch', ...
              'si and x must have the same length (got %d and %d).', n, numel(x));
    end
    if n < 3
        error('diffOrder2NonUniform:tooFewPoints', ...
              'Need at least 3 points for a 2nd-order scheme (got %d).', n);
    end

    so = si * 0;   % preserve shape, class and complexity

    % --- first point: 2nd-order forward ---
    num   = -si(3) * (x(2) - x(1))^2 + si(2) * (x(3) - x(1))^2 ...
            - si(1) * ((x(3) - x(1))^2 - (x(2) - x(1))^2);
    den   = (x(2) - x(1)) * (x(3) - x(1)) * (x(3) - x(2));
    so(1) = num / den;

    % --- last point: 2nd-order backward ---
    num     = -si(end-2) * (x(end-1) - x(end))^2 + si(end-1) * (x(end-2) - x(end))^2 ...
              - si(end) * ((x(end-2) - x(end))^2 - (x(end-1) - x(end))^2);
    den     = (x(end-1) - x(end)) * (x(end-2) - x(end)) * (x(end-2) - x(end-1));
    so(end) = num / den;

    % --- interior points: 2nd-order central ---
    for j = 2:n-1
        hMinus = x(j)   - x(j-1);
        hPlus  = x(j+1) - x(j);
        num    = si(j+1) * hMinus^2 - si(j-1) * hPlus^2 + si(j) * (hPlus^2 - hMinus^2);
        den    = hPlus * hMinus * (hMinus + hPlus);
        so(j)  = num / den;
    end
end
