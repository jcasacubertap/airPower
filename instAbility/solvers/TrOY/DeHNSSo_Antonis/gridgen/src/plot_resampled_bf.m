function plot_resampled_bf(BF, opts)
%PLOT_RESAMPLED_BF  Visualise the interpolated base-flow on the new curvilinear mesh.
%
% Produces two or three figures:
%   1. Mesh grid lines (top panel) and U, V, W contours (panels 2-4) on the
%      physical (X,Y) mesh
%   2. Wall-normal profiles of U, V, W at several xi stations (tiled)
%   3. (optional) Streamwise profiles of U, V, W at several eta stations (tiled)
%
% Optional opts fields:
%   .xi_slices      : xi indices or fractions [0,1] to plot wall-normal profiles
%                     (default: 8 equally spaced stations)
%   .eta_slices     : eta indices or fractions [0,1] to plot streamwise profiles
%   .mesh_skip_xi   : draw every N-th ξ-line (default: max(1, round(n_xi/60)))
%   .mesh_skip_eta  : draw every N-th η-line (default: max(1, round(n_eta/20)))

if nargin < 2;  opts = struct();  end

X = BF.X;   % [n_eta x n_xi]
Y = BF.Y;
U = BF.U;
V = BF.V;
W = BF.W;
comps = {U, V, W};
names = {'U','V','W'};

[n_eta, n_xi] = size(U);

%% Figure 1: Mesh grid + U, V, W contours on the physical mesh
skip_xi  = get_opt(opts, 'mesh_skip_xi',  max(1, round(n_xi  / 60)));
skip_eta = get_opt(opts, 'mesh_skip_eta', max(1, round(n_eta / 20)));

figure('Name','Resampled mesh and BF components');

subplot(4,1,1);
plot(X(:, 1:skip_xi:end), Y(:, 1:skip_xi:end), 'k-', 'LineWidth', 0.3); hold on;
plot(X(1:skip_eta:end, :)', Y(1:skip_eta:end, :)', 'k-', 'LineWidth', 0.3);
plot(X(end, :), Y(end, :), 'r-', 'LineWidth', 1.5);
daspect([1 1 1]);
xlabel('x'); ylabel('y');
title(sprintf('Resampled mesh (%d \\eta \\times %d \\xi, every %dth \\xi and %dth \\eta line)', ...
      n_eta, n_xi, skip_xi, skip_eta));
axis tight;

for c = 1:3
    subplot(4, 1, 1 + c);
    contourf(X, Y, comps{c}, 40, 'LineColor', 'none');
    colorbar;
    daspect([1 1 1]);
    xlabel('x'); ylabel('y');
    title(names{c});
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

end


function s = ternary(cond, a, b)
if cond;  s = a;  else;  s = b;  end
end


function v = get_opt(opts, name, default)
% Fetch opts.(name) if set and non-empty, otherwise return default.
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
