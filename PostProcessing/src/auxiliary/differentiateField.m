function [Fx, Fy] = differentiateField(F, sBF, order)
% differentiateField  Cartesian gradient of a field on the importData grid.
%
%   [Fx, Fy] = differentiateField(F, sBF)
%   [Fx, Fy] = differentiateField(F, sBF, order)
%
% F is Ny x Nx (or Ny x Nx x K, K fields differentiated together) on the grid of
% sBF, as returned by importData from EITHER loader. Fx = dF/dX and Fy = dF/dY in
% the Cartesian frame of sBF.geom.X/Y, per unit length of sBF (m for loadBF,
% delta_0 for loadFields). Real or complex.
%
% Works on any structured grid, curved or leaning, by the chain rule through the
% grid metrics:
%   1. derivatives along the grid lines, with xi = s (wall arc length, along the
%      rows) and eta = n (distance up each column), by finite differences of
%      the given order (deriv1Uneven, default 4): F_xi, F_eta, and the same
%      for X and Y;
%   2. the Jacobian J = X_xi Y_eta - X_eta Y_xi and the inverse metrics
%         xi_x = Y_eta/J,  eta_x = -Y_xi/J,  xi_y = -X_eta/J,  eta_y = X_xi/J;
%   3. Fx = xi_x F_xi + eta_x F_eta,   Fy = xi_y F_xi + eta_y F_eta.
% On a Cartesian grid this reduces to plain differences along x and y. Every
% row, wall and free-stream included, uses the same FD stencils (one-sided at
% the ends), so nothing depends on how the grid was built.

    if nargin < 3 || isempty(order); order = 4; end
    if ~all(isfield(sBF, {'s', 'n'})) || ~isfield(sBF, 'geom')
        error('differentiateField:grid', ...
              'sBF needs .s, .n and .geom.X/.Y (load it through importData).');
    end
    X = sBF.geom.X;  Y = sBF.geom.Y;
    [Ny, Nx] = size(X);
    if size(F, 1) ~= Ny || size(F, 2) ~= Nx
        error('differentiateField:size', ...
              'F is %dx%d but the grid is %dx%d.', size(F, 1), size(F, 2), Ny, Nx);
    end
    K  = size(F, 3);
    xi = sBF.s(end, :).';                          % one wall station per column

    % stack geometry and fields so each line is differentiated in one pass:
    % pages 1,2 = X,Y; pages 3.. = the fields
    Q = cat(3, X, Y, F);
    Q_xi  = zeros(size(Q), 'like', Q);
    Q_eta = zeros(size(Q), 'like', Q);

    % along the rows: xi is the same for every row
    for k = 1:size(Q, 3)
        Q_xi(:, :, k) = deriv1Uneven(xi, Q(:, :, k).', order).';
    end

    % up the columns: on a body-fitted mesh n differs per column, so each column
    % is taken with its own n; on a grid whose columns share n, in one pass
    nSame = all(sBF.n == sBF.n(:, 1), 'all');
    if nSame
        for k = 1:size(Q, 3)
            Q_eta(:, :, k) = deriv1Uneven(sBF.n(:, 1), Q(:, :, k), order);
        end
    else
        for c = 1:Nx
            col = reshape(Q(:, c, :), Ny, []);
            Q_eta(:, c, :) = reshape(deriv1Uneven(sBF.n(:, c), col, order), Ny, 1, []);
        end
    end

    % inverse metrics
    X_xi = real(Q_xi(:, :, 1));   Y_xi = real(Q_xi(:, :, 2));
    X_eta = real(Q_eta(:, :, 1)); Y_eta = real(Q_eta(:, :, 2));
    J = X_xi .* Y_eta - X_eta .* Y_xi;
    if any(abs(J(:)) < eps * max(abs(J(:))))
        warning('differentiateField:singular', 'Near-zero grid Jacobian.');
    end
    xi_x  =  Y_eta ./ J;   eta_x = -Y_xi ./ J;
    xi_y  = -X_eta ./ J;   eta_y =  X_xi ./ J;

    F_xi  = Q_xi(:, :, 3:end);
    F_eta = Q_eta(:, :, 3:end);
    Fx = xi_x .* F_xi + eta_x .* F_eta;
    Fy = xi_y .* F_xi + eta_y .* F_eta;
    if K == 1; Fx = Fx(:, :, 1); Fy = Fy(:, :, 1); end
end
