function c = fornberg_weights(z, x, m)
%FORNBERG_WEIGHTS  Finite-difference weights on non-uniform grids.
%
%   c = fornberg_weights(z, x, m)
%
%   Returns an (n × (m+1)) weight matrix c such that, for any function f
%   sampled at stencil points x(1..n),
%       f^(k)(z) ≈ sum_j c(j, k+1) * f(x(j))
%   for k = 0, 1, ..., m.
%
%   Implements B. Fornberg, "Generation of finite difference formulas on
%   arbitrarily spaced grids", Math. Comp. 51 (1988) 699-706.

x = x(:);
n = length(x) - 1;
c = zeros(n+1, m+1);

c1 = 1;
c4 = x(1) - z;
c(1, 1) = 1;

for i = 1:n
    mn = min(i, m);
    c2 = 1;
    c5 = c4;
    c4 = x(i+1) - z;
    for j = 0:i-1
        c3 = x(i+1) - x(j+1);
        c2 = c2 * c3;
        if j == i-1
            for k = mn:-1:1
                c(i+1, k+1) = c1 * (k * c(i, k) - c5 * c(i, k+1)) / c2;
            end
            c(i+1, 1) = -c1 * c5 * c(i, 1) / c2;
        end
        for k = mn:-1:1
            c(j+1, k+1) = (c4 * c(j+1, k+1) - k * c(j+1, k)) / c3;
        end
        c(j+1, 1) = c4 * c(j+1, 1) / c3;
    end
    c1 = c2;
end
end
