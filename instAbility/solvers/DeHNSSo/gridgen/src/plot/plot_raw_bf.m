function plot_raw_bf(BF, class)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Visualise the raw (pre-resample) base-flow mesh and components.
%
%   plot_raw_bf(BF, class)
%
%   One figure, four stacked panels: mesh, then U, V, W contours on the raw grid.
%   Handles all three input classes: 'cartesian' via meshgrid and contourf,
%   'curvilinear_structured' by contouring the 2-D X, Y directly, and
%   'unstructured' as a coloured scatter per component.

%% Inputs
% Name  size         unit   explanation
% BF    (1,1) struct [-]    base flow on its own grid
% class char         [-]    input layout from classify_input, selecting how the raw base flow is drawn

%% Outputs
% None.

U = BF.U; V = BF.V; W = BF.W;
comps = {U, V, W};
names = {'U','V','W'};

figure('Name','Raw BF');

% Each class needs a different call: contourf on a meshgrid, contourf on the
% 2-D arrays directly, or a coloured scatter.
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
            if (max(comps{c}(:)) - min(comps{c}(:))) > eps
                contourf(Xg, Yg, comps{c}, 40, 'LineColor', 'none');
            else
                imagesc(x, y, comps{c}); axis xy;   % constant field: flat colour map
            end
            colorbar; daspect([1 1 1]);
            xlabel('x'); ylabel('y'); title(names{c});
        end

        % Wall-row diagnostics on the RAW source, before any interpolation, so a
        % bad wall row can be told apart from a resampling artefact.
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
        transposed = (max(X(1, :)) - min(X(1, :))) < (max(X(:, 1)) - min(X(:, 1)));
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
            if (max(comps{c}(:)) - min(comps{c}(:))) > eps
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

tolX = max(1e-6 * (max(BF.X) - min(BF.X)), 1e-12);
tolY = max(1e-6 * (max(BF.Y) - min(BF.Y)), 1e-12);
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
