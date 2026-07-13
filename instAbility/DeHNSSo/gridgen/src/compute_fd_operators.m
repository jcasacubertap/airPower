function Grid = compute_fd_operators(Grid)
%COMPUTE_FD_OPERATORS  Sparse FD D1, D2 on a non-uniform wall-normal grid.
%
%   Builds sparse differentiation matrices using Fornberg weights on the
%   non-uniform Malik-Chebyshev grid. Order is set by Grid.FD_eta_order:
%     2 : tridiagonal  (bandwidth 3) — fast, lower accuracy
%     4 : pentadiagonal (bandwidth 5) — good balance (recommended)
%     6 : heptadiagonal (bandwidth 7) — higher accuracy, wider stencil
%   Interior rows use local stencils of the specified order.
%   Boundary rows (i=1 freestream, i=N wall) use all N grid points,
%   giving spectral-quality (global polynomial) accuracy at the boundaries.
%   This is required for the (0,0) MFD mode whose v-freestream BC row
%   retains the differential operator in MLAMFD.
%
%   Storage: O(N) sparse for interior rows; boundary rows are dense.

y   = Grid.y(:);
N   = Grid.ny;
p   = Grid.FD_eta_order;   % derivative order (2, 4, or 6)
hw  = p / 2;               % half-width of stencil (interior)

%% Malik mapping derivatives (needed by ILST_solver even in FD mode)
a = Grid.H * Grid.y_i / (Grid.H - 2.0 * Grid.y_i);
b = 1.0 + 2.0 * a / Grid.H;
Grid.dzdy   = a * (1 + b) ./ (a + y).^2;
Grid.d2zdy2 = -2 * a * (1 + b) ./ (a + y).^3;

D1 = spalloc(N, N, N * (p + 1));
D2 = spalloc(N, N, N * (p + 1));

for i = 1:N
    % Interior rows: local stencil of order p (one-sided near boundaries)
    i0 = max(1,   i - hw);
    i1 = min(N,   i0 + p);
    i0 = max(1,   i1 - p);
    idx = i0:i1;

    nodes = y(idx);

    w1 = fornberg_weights(y(i), nodes, 1);
    w2 = fornberg_weights(y(i), nodes, 2);

    D1(i, idx) = w1;   %#ok<SPRIX>
    D2(i, idx) = w2;   %#ok<SPRIX>
end

% Override boundary rows (i=1 freestream, i=N wall) with Chebyshev values.
% The one-sided FD stencil at boundaries is inaccurate on the Malik grid;
% Chebyshev rows are spectrally accurate and critical for the (0,0) MFD mode
% whose v-freestream BC row retains the differential operator in MLAMFD.
Dn = -chebDiff(N);
D1cheb = diag(Grid.dzdy)   * Dn;
D2cheb = diag(Grid.d2zdy2) * Dn + diag(Grid.dzdy.^2) * (Dn * Dn);
D1(1,:) = D1cheb(1,:);
D1(N,:) = D1cheb(N,:);
D2(1,:) = D2cheb(1,:);
D2(N,:) = D2cheb(N,:);

Grid.D1 = D1;
Grid.D2 = D2;

end


function w = fornberg_weights(z, x, m)
%FORNBERG_WEIGHTS  Finite difference weights for derivative order m at z
%   on nodes x, using Fornberg (1998) algorithm.
%   Returns row vector of weights (length = numel(x)).

n  = numel(x) - 1;
c  = zeros(n+1, m+1);
c1 = 1;
c4 = x(1) - z;
c(1,1) = 1;

for i = 1:n
    mn  = min(i, m);
    c2  = 1;
    c5  = c4;
    c4  = x(i+1) - z;

    for j = 0:i-1
        c3 = x(i+1) - x(j+1);
        c2 = c2 * c3;

        for k = mn:-1:1
            c(i+1, k+1) = c1 * (k * c(i, k) - c5 * c(i, k+1)) / c2;
        end
        c(i+1, 1) = -c1 * c5 * c(i, 1) / c2;

        for k = mn:-1:1
            c(j+1, k+1) = (c4 * c(j+1, k+1) - k * c(j+1, k)) / c3;
        end
        c(j+1, 1) = c4 * c(j+1, 1) / c3;
    end

    c1 = c2;
end

w = c(:, m+1)';

end
