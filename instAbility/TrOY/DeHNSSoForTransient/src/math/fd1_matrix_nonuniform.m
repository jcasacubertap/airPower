function D = fd1_matrix_nonuniform(x)
%FD1_MATRIX_NONUNIFORM  Sparse 4th-order 1st-derivative matrix for a non-uniform 1-D grid.
%
%   D = fd1_matrix_nonuniform(x)   with x a column vector of n stations.
%
%   D is an (n × n) sparse matrix such that y = D*f gives a 4th-order-accurate
%   estimate of df/dx at each station, including near the boundaries where
%   the stencil becomes biased. On a uniform grid D reduces to the classical
%   4th-order finite-difference coefficients. Stencils mirror those used by
%   the HNS LHS assembly:
%       i = 1        : 5-point forward, x(1:5)
%       i = 2        : 5-point forward-biased, x(1:5)
%       3 ≤ i ≤ n-2  : 5-point centred
%       i = n-1      : 4-point backward-biased, x(n-3:n)
%       i = n        : 3-point backward, x(n-2:n)

x = x(:);
n = length(x);
if n < 5
    error('fd1_matrix_nonuniform: need at least 5 stations (got %d).', n);
end

nnz_est = 5 * n;
row = zeros(nnz_est, 1);
col = zeros(nnz_est, 1);
val = zeros(nnz_est, 1);
k   = 0;

for i = 1:n
    if i == 1
        cols  = 1:5;
    elseif i == 2
        cols  = 1:5;
    elseif i <= n-2
        cols  = (i-2) : (i+2);
    elseif i == n-1
        cols  = (n-3) : n;
    else
        cols  = (n-2) : n;
    end
    c  = fornberg_weights(x(i), x(cols), 1);
    w1 = c(:, 2).';
    idx = k + (1:numel(cols));
    row(idx) = i;
    col(idx) = cols;
    val(idx) = w1;
    k = k + numel(cols);
end

row = row(1:k); col = col(1:k); val = val(1:k);
D = sparse(row, col, val, n, n);
end
