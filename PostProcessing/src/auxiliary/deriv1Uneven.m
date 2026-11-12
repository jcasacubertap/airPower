function D = deriv1Uneven(x, F, order)
% deriv1Uneven  First derivative on a non-uniform 1-D grid, column by column.
%
%   D = deriv1Uneven(x, F)
%   D = deriv1Uneven(x, F, order)
%
% Differentiates F with respect to x down its columns: x is the n x 1
% coordinate, F an n x m matrix whose columns share that coordinate (a vector
% is treated as one column). D has the size of F.
%
% order = 2, 4 or 6 (default 4). Fornberg-type weights of that order: centred
% (order+1)-point stencils in the interior, the nearest (order+1) points at the
% ends. Order 2 is the same 3-point scheme the old diffOrder2NonUniform used.
% x may increase or decrease, and need not be uniformly spaced. Complex F is
% fine.

    if nargin < 3 || isempty(order); order = 4; end
    if ~ismember(order, [2 4 6])
        error('deriv1Uneven:badOrder', 'order must be 2, 4 or 6 (got %g).', order);
    end
    isRow = isrow(F);
    if isRow; F = F.'; end
    x = x(:);
    n = numel(x);
    if size(F, 1) ~= n
        error('deriv1Uneven:sizeMismatch', ...
              'x has %d points but F has %d rows.', n, size(F, 1));
    end
    p = order + 1;                           % stencil width
    if n < p
        error('deriv1Uneven:tooFewPoints', ...
              'Order %d needs at least %d points (got %d).', order, p, n);
    end

    h = order / 2;
    D = zeros(size(F), 'like', F);
    for i = 1:n
        i0 = min(max(1, i - h), n - p + 1);  % centred, shifted inwards at the ends
        idx = i0:(i0 + p - 1);
        D(i, :) = weights(x(i), x(idx)) * F(idx, :);
    end
    if isRow; D = D.'; end
end

% first-derivative weights at z from the nodes xs (Taylor/Vandermonde solve,
% scaled by the local spacing for conditioning)
function w = weights(z, xs)
    d  = xs(:).' - z;
    sc = max(abs(d));
    k  = (0:numel(d) - 1).';
    A  = (d / sc) .^ k;
    b  = zeros(numel(d), 1);  b(2) = 1;
    w  = (A \ b).' / sc;
end
