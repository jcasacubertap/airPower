function plot_raw_bf(BF, class)
%PLOT_RAW_BF  Visualise the raw (pre-resample) base-flow mesh and components.
%
%   plot_raw_bf(BF, class)
%
%   One figure with four stacked panels: mesh, then U, V, W contours on the
%   raw (pre-resample) grid. Works for all three input classes:
%     'cartesian'              — 1-D axes; mesh via meshgrid, contourf on (Xg, Yg, field)
%     'curvilinear_structured' — 2-D X, Y arrays; contourf directly
%     'unstructured'           — scatter points; coloured scatter per component

U = BF.U; V = BF.V; W = BF.W;
comps = {U, V, W};
names = {'U','V','W'};

figure('Name','Raw BF');

switch class
    case 'cartesian'
        % Flat-lattice scatter (paired vectors) → 1-D axes + 2-D fields.
        if isvector(BF.X) && isvector(BF.Y) && numel(BF.X) == numel(BF.Y)
            BF = reshape_scatter(BF);
            U = BF.U; V = BF.V; W = BF.W;
            comps = {U, V, W};
        end
        x = BF.X(:).';
        y = BF.Y(:);
        [Xg, Yg] = meshgrid(x, y);
        sx = max(1, round(numel(x) / 60));
        sy = max(1, round(numel(y) / 20));

        % Orient U/V/W to match meshgrid (ny × nx)
        nx = numel(x); ny = numel(y);
        for c = 1:3
            if size(comps{c},1) == nx && size(comps{c},2) == ny
                comps{c} = comps{c}.';
            end
        end

        % Identify the wall row from velocity: where |U|+|V|+|W| is minimal.
        speed_per_row = mean(abs(comps{1}) + abs(comps{2}) + abs(comps{3}), 2);
        [~, k_wall] = min(speed_per_row);

        subplot(4,1,1);
        plot(Xg(:, 1:sx:end),   Yg(:, 1:sx:end),   'k-', 'LineWidth', 0.3); hold on;
        plot(Xg(1:sy:end, :)', Yg(1:sy:end, :)', 'k-', 'LineWidth', 0.3);
        plot(Xg(k_wall, :), Yg(k_wall, :), 'r-', 'LineWidth', 1.5);
        daspect([1 1 1]); xlabel('x'); ylabel('y'); axis tight;
        title(sprintf('Raw BF mesh (Cartesian): %d y × %d x  (wall = row %d)', ny, nx, k_wall));

        for c = 1:3
            subplot(4,1,1+c);
            if range(comps{c}(:)) > eps
                contourf(Xg, Yg, comps{c}, 40, 'LineColor', 'none');
            else
                imagesc(x, y, comps{c}); axis xy;   % constant field: flat colour map
            end
            colorbar; daspect([1 1 1]);
            xlabel('x'); ylabel('y'); title(names{c});
        end

        % ----- Diagnostic figure: wall-row values & dy* on the RAW source -----
        x_wall = Xg(k_wall, :);

        figure('Name','Raw BF — wall row diagnostics');
        for c = 1:3
            subplot(3, 2, 2*(c-1) + 1);
            plot(x_wall, comps{c}(k_wall, :), 'k-', 'LineWidth', 1.0);
            xlabel('x'); ylabel(names{c});
            title(sprintf('%s at wall (raw, row %d)', names{c}, k_wall));
            grid on;

            [~, dy_c] = gradient(comps{c}, 1, y(:));
            subplot(3, 2, 2*(c-1) + 2);
            plot(x_wall, dy_c(k_wall, :), 'k-', 'LineWidth', 1.0);
            xlabel('x'); ylabel(['\partial ' names{c} ' / \partial y']);
            title(sprintf('\\partial %s/\\partial y at wall (raw, row %d)', names{c}, k_wall));
            grid on;
        end

    case 'curvilinear_structured'
        X = BF.X; Y = BF.Y;
        transposed = range(X(1, :)) < range(X(:, 1));
        if transposed
            X = X.'; Y = Y.';
            for c = 1:3;  comps{c} = comps{c}.';  end
        end
        [ny, nx] = size(X);
        sx = max(1, round(nx / 60));
        sy = max(1, round(ny / 20));

        % Identify the wall row from velocity (|U|+|V|+|W| minimal).
        speed_per_row = mean(abs(comps{1}) + abs(comps{2}) + abs(comps{3}), 2);
        [~, k_wall] = min(speed_per_row);

        subplot(4,1,1);
        plot(X(:, 1:sx:end),   Y(:, 1:sx:end),   'k-', 'LineWidth', 0.3); hold on;
        plot(X(1:sy:end, :)', Y(1:sy:end, :)', 'k-', 'LineWidth', 0.3);
        plot(X(k_wall, :), Y(k_wall, :), 'r-', 'LineWidth', 1.5);
        daspect([1 1 1]); xlabel('x'); ylabel('y'); axis tight;
        title(sprintf('Raw BF mesh (curvilinear structured): %d η × %d ξ  (wall = row %d)', ny, nx, k_wall));

        for c = 1:3
            subplot(4,1,1+c);
            if range(comps{c}(:)) > eps
                contourf(X, Y, comps{c}, 40, 'LineColor', 'none');
            else
                pcolor(X, Y, comps{c}); shading flat;
            end
            colorbar; daspect([1 1 1]);
            xlabel('x'); ylabel('y'); title(names{c});
        end

    case 'unstructured'
        xs = BF.X(:); ys = BF.Y(:);

        subplot(4,1,1);
        plot(xs, ys, 'k.', 'MarkerSize', 2);
        daspect([1 1 1]); xlabel('x'); ylabel('y'); axis tight;
        title(sprintf('Raw BF mesh (unstructured): %d scattered points', numel(xs)));

        for c = 1:3
            subplot(4,1,1+c);
            scatter(xs, ys, 3, comps{c}(:), 'filled');
            colorbar; daspect([1 1 1]);
            xlabel('x'); ylabel('y'); title(names{c});
        end

    otherwise
        error('plot_raw_bf: unknown class ''%s''.', class);
end
end


function BF = reshape_scatter(BF)
%RESHAPE_SCATTER  Flat-lattice scatter → canonical (ny × nx) gridded form.
%   X, Y are paired column vectors lying on a rectangular lattice; this
%   helper turns them into 1-D axes (X = 1×nx, Y = ny×1) and reshapes any
%   per-point numeric field of length numel(X) into a (ny × nx) array.
%   No-op when BF is already 2-D or when X/Y lengths differ.
if ~isvector(BF.X) || ~isvector(BF.Y) || numel(BF.X) ~= numel(BF.Y)
    return;
end

tolX = max(1e-6 * range(BF.X), 1e-12);
tolY = max(1e-6 * range(BF.Y), 1e-12);
[Xu, ~, ix] = uniquetol(BF.X(:), tolX, 'DataScale', 1);
[Yu, ~, iy] = uniquetol(BF.Y(:), tolY, 'DataScale', 1);
nx = numel(Xu);
ny = numel(Yu);

fprintf('Reshaping flat scatter → (ny=%d × nx=%d)\n', ny, nx);

fnames = setdiff(fieldnames(BF), {'X','Y'});
for k = 1:numel(fnames)
    f   = fnames{k};
    val = BF.(f);
    if ~isnumeric(val) || numel(val) ~= numel(ix)
        continue
    end
    M   = nan(ny, nx);
    lin = sub2ind([ny, nx], iy, ix);
    M(lin) = val(:);
    if any(isnan(M(:)))
        warning('plot_raw_bf:gaps', ...
            'Field %s has %d missing lattice points after reshape.', f, nnz(isnan(M)));
    end
    BF.(f) = M;
end

BF.X = Xu(:).';
BF.Y = Yu(:);
end
