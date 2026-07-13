function plot_resampled_bf(BF, opts)
%PLOT_RESAMPLED_BF  Visualise the interpolated base-flow on the new curvilinear mesh.
%
% Produces up to five figures:
%   1.  Mesh and U, V, W contours in physical (x, y)
%   1b. Same mesh and contours in computational (ξ, η) — body-fitted view
%   2.  Wall-normal U, V, W profiles at several ξ stations
%   3.  (optional) Streamwise U, V, W profiles at several η stations
%   4.  Wall-normal ∂/∂η profiles at the same ξ stations as Figure 2
%
% Optional opts fields:
%   .xi_slices      : xi indices or fractions [0,1] to plot wall-normal profiles
%                     (default: 8 equally spaced stations)
%   .eta_slices     : eta indices or fractions [0,1] to plot streamwise profiles
%   .mesh_skip_xi   : draw every N-th ξ-line (default: max(1, round(n_xi/60)))
%   .mesh_skip_eta  : draw every N-th η-line (default: max(1, round(n_eta/20)))

if nargin < 2;  opts = struct();  end

X = BF.X;
Y = BF.Y;
U = BF.U;
V = BF.V;
W = BF.W;
comps = {U, V, W};
names = {'U','V','W'};

% Cartesian inputs have 1-D X, Y; promote to 2-D so the rest of this
% function (which assumes (n_eta × n_xi) arrays) works uniformly.
if isvector(X) && isvector(Y)
    [X, Y] = meshgrid(X(:).', Y(:));
end

[n_eta, n_xi] = size(U);

%% Figure 1: Mesh grid + U, V, W contours on the physical mesh
skip_xi  = get_opt(opts, 'mesh_skip_xi',  max(1, round(n_xi  / 60)));
skip_eta = get_opt(opts, 'mesh_skip_eta', max(1, round(n_eta / 20)));

figure('Name','Resampled mesh and BF — physical (x, y)');

% Identify the wall row from velocity: where |U|+|V|+|W| is minimal.
speed_per_row = mean(abs(comps{1}) + abs(comps{2}) + abs(comps{3}), 2);
[~, k_wall]   = min(speed_per_row);

subplot(4,1,1);
plot(X(:, 1:skip_xi:end), Y(:, 1:skip_xi:end), 'k-', 'LineWidth', 0.3); hold on;
plot(X(1:skip_eta:end, :)', Y(1:skip_eta:end, :)', 'k-', 'LineWidth', 0.3);
plot(X(k_wall, :), Y(k_wall, :), 'r-', 'LineWidth', 1.5);
daspect([1 1 1]);
xlabel('x'); ylabel('y');
title(sprintf('Resampled mesh — physical (x, y): %d \\eta \\times %d \\xi (every %dth \\xi, %dth \\eta; wall = row %d)', ...
      n_eta, n_xi, skip_xi, skip_eta, k_wall));
axis tight;

for c = 1:3
    subplot(4, 1, 1 + c);
    if range(comps{c}(:)) > eps
        contourf(X, Y, comps{c}, 40, 'LineColor', 'none');
    else
        pcolor(X, Y, comps{c}); shading flat;
    end
    colorbar;
    daspect([1 1 1]);
    xlabel('x'); ylabel('y');
    title(sprintf('%s — physical (x, y)', names{c}));
end

%% Figure 1b: Computational (ξ, η) view — unrolls a curvilinear wall flat
[Xi_g, Eta_g] = meshgrid(BF.xi1D(:).', BF.eta1D(:));
figure('Name','Resampled mesh — computational (ξ, η)');
subplot(4,1,1);
plot(Xi_g(:, 1:skip_xi:end), Eta_g(:, 1:skip_xi:end), 'k-', 'LineWidth', 0.3); hold on;
plot(Xi_g(1:skip_eta:end, :)', Eta_g(1:skip_eta:end, :)', 'k-', 'LineWidth', 0.3);
plot(Xi_g(k_wall, :), Eta_g(k_wall, :), 'r-', 'LineWidth', 1.5);
xlabel('\xi'); ylabel('\eta');
title(sprintf('Resampled mesh — computational (\\xi, \\eta): %d \\eta \\times %d \\xi (wall = row %d)', ...
      n_eta, n_xi, k_wall));
axis tight;

for c = 1:3
    subplot(4, 1, 1 + c);
    if range(comps{c}(:)) > eps
        contourf(Xi_g, Eta_g, comps{c}, 40, 'LineColor', 'none');
    else
        pcolor(Xi_g, Eta_g, comps{c}); shading flat;
    end
    colorbar;
    xlabel('\xi'); ylabel('\eta');
    title(sprintf('%s — computational (\\xi, \\eta)', names{c}));
end

%% Figure 2: Wall-normal profiles at several xi stations
if isfield(opts, 'xi_slices') && ~isempty(opts.xi_slices)
    xi_slices = opts.xi_slices;
    if all(xi_slices <= 1) && all(xi_slices >= 0)
        slice_idx = max(1, round(xi_slices * n_xi));
    else
        slice_idx = xi_slices;
    end
else
    slice_idx = round(linspace(1, n_xi, 8));
end

figure('Name','Wall-normal profiles');
colors = lines(numel(slice_idx));
for c = 1:3
    subplot(1,3,c); hold on;
    for k = 1:numel(slice_idx)
        j = slice_idx(k);
        displayname = '';
        if c == 1
            displayname = sprintf('\\xi = %.3f', BF.xi1D(j));
        end
        plot(comps{c}(:, j), BF.eta1D, '-', 'Color', colors(k,:), ...
             'DisplayName', displayname, 'HandleVisibility', ternary(c==1,'on','off'));
    end
    xlabel(names{c});  ylabel('\eta');
    title(sprintf('%s profiles', names{c}));
    grid on;
    if c == 1;  legend('Location','best');  end
end

%% Figure 3: Streamwise profiles at several eta stations (optional)
if isfield(opts, 'eta_slices') && ~isempty(opts.eta_slices)
    eta_slices = opts.eta_slices;
    if all(eta_slices <= 1) && all(eta_slices >= 0)
        eta_idx = max(1, round(eta_slices * n_eta));
    else
        eta_idx = eta_slices;
    end

    figure('Name','Streamwise profiles');
    colors = lines(numel(eta_idx));
    for c = 1:3
        subplot(3,1,c); hold on;
        for k = 1:numel(eta_idx)
            i = eta_idx(k);
            displayname = '';
            if c == 1
                displayname = sprintf('\\eta = %.3f', BF.eta1D(i));
            end
            plot(BF.xi1D, comps{c}(i, :), '-', 'Color', colors(k,:), ...
                 'DisplayName', displayname, 'HandleVisibility', ternary(c==1,'on','off'));
        end
        xlabel('\xi');  ylabel(names{c});
        title(sprintf('%s profiles', names{c}));
        grid on;
        if c == 1;  legend('Location','best');  end
    end
end

%% Figure: wall-normal ∂/∂η profiles at the same ξ stations as Figure 2
eta_vec = BF.eta1D(:);
[~, dyU] = gradient(BF.U, 1, eta_vec);
[~, dyV] = gradient(BF.V, 1, eta_vec);
[~, dyW] = gradient(BF.W, 1, eta_vec);
dys    = {dyU, dyV, dyW};
dnames = {'\partial U / \partial \eta', '\partial V / \partial \eta', '\partial W / \partial \eta'};

figure('Name','Resampled BF — wall-normal derivative profiles');
colors = lines(numel(slice_idx));
for c = 1:3
    subplot(1,3,c); hold on;
    for k = 1:numel(slice_idx)
        j = slice_idx(k);
        displayname = '';
        if c == 1
            displayname = sprintf('\\xi = %.3f', BF.xi1D(j));
        end
        plot(dys{c}(:, j), BF.eta1D, '-', 'Color', colors(k,:), ...
             'DisplayName', displayname, 'HandleVisibility', ternary(c==1,'on','off'));
    end
    xlabel(dnames{c});  ylabel('\eta');
    title(sprintf('%s profiles', dnames{c}));
    grid on;
    if c == 1;  legend('Location','best');  end
end

end


function s = ternary(cond, a, b)
if cond;  s = a;  else;  s = b;  end
end
