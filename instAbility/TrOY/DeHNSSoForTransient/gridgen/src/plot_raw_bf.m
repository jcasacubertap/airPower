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

        subplot(4,1,1);
        plot(Xg(:, 1:sx:end),   Yg(:, 1:sx:end),   'k-', 'LineWidth', 0.3); hold on;
        plot(Xg(1:sy:end, :)', Yg(1:sy:end, :)', 'k-', 'LineWidth', 0.3);
        plot(Xg(end, :), Yg(end, :), 'r-', 'LineWidth', 1.5);
        daspect([1 1 1]); xlabel('x'); ylabel('y'); axis tight;
        title(sprintf('Raw BF mesh (Cartesian): %d y × %d x', ny, nx));

        for c = 1:3
            subplot(4,1,1+c);
            contourf(Xg, Yg, comps{c}, 40, 'LineColor', 'none');
            colorbar; daspect([1 1 1]);
            xlabel('x'); ylabel('y'); title(names{c});
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

        subplot(4,1,1);
        plot(X(:, 1:sx:end),   Y(:, 1:sx:end),   'k-', 'LineWidth', 0.3); hold on;
        plot(X(1:sy:end, :)', Y(1:sy:end, :)', 'k-', 'LineWidth', 0.3);
        plot(X(end, :), Y(end, :), 'r-', 'LineWidth', 1.5);
        daspect([1 1 1]); xlabel('x'); ylabel('y'); axis tight;
        title(sprintf('Raw BF mesh (curvilinear structured): %d η × %d ξ', ny, nx));

        for c = 1:3
            subplot(4,1,1+c);
            contourf(X, Y, comps{c}, 40, 'LineColor', 'none');
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
