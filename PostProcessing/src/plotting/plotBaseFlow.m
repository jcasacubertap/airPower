function fig = plotBaseFlow(sBF, inp, savedir)
% plotBaseFlow  Contour plots of the base-flow fields on the (x,y) grid.
%
%   fig = plotBaseFlow(sBF, inp)            % show figure only
%   fig = plotBaseFlow(sBF, inp, savedir)   % also save baseFlow.png in savedir
%
% One stacked panel per available field among {u, v, w, p} (loadFields provides
% u,v,w; loadBF adds p). Base-flow quantities carry a subscript B. Data are
% plotted as-is; only the labels reflect the units:
%   loadFields (DeHNSSo) -> non-dimensional: x/delta_0, u_B/u_infty, ...
%   loadBF     (OpenFOAM) -> dimensional:    x [m], u_B [m/s], p_B [m^2/s^2]
% (OpenFOAM pressure is kinematic, p/rho, hence m^2/s^2.)
%
% Save behaviour: savedir non-empty -> PNG written there (terminal/dispatcher);
% empty/omitted -> figure left open for interactive inspection (direct MATLAB).

    if nargin < 3; savedir = ''; end

    nd = isfield(inp, 'loadMode') && strcmpi(inp.loadMode, 'loadFields');

    if nd
        xlab = '$x / \delta_0$';  ylab = '$y / \delta_0$';
        catalog = { ...
            'u', '$u_{\mathrm{B}} / u_\infty$'; ...
            'v', '$v_{\mathrm{B}} / u_\infty$'; ...
            'w', '$w_{\mathrm{B}} / u_\infty$'; ...
            'p', '$p_{\mathrm{B}} / (\rho\, u_\infty^2)$'};
    else
        xlab = '$x \; [\mathrm{m}]$';  ylab = '$y \; [\mathrm{m}]$';
        catalog = { ...
            'u', '$u_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'v', '$v_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'w', '$w_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'p', '$p_{\mathrm{B}} \; [\mathrm{m^2/s^2}]$'};
    end

    have = false(size(catalog,1),1);
    for k = 1:size(catalog,1)
        f = catalog{k,1};
        have(k) = isfield(sBF, f) && ~isempty(sBF.(f)) && isequal(size(sBF.(f)), size(sBF.x));
    end
    catalog = catalog(have,:);
    nF = size(catalog,1);
    if nF == 0
        warning('plotBaseFlow:noFields', 'No base-flow fields to plot.'); fig = []; return;
    end

    X = sBF.x;  Y = sBF.y;
    [xl, yl] = plotWindow(X, Y, inp);   % buffer-cut in x, near-wall window in y

    fig = figure('Name', 'Base flow', 'Color', 'w', ...
                 'Position', [60 60 900 220*nF]);
    for k = 1:nF
        f   = catalog{k,1};
        lbl = catalog{k,2};
        F   = sBF.(f);

        subplot(nF, 1, k);
        contourf(X, Y, F, 30, 'LineStyle', 'none'); hold on;

        % Diverging map centred at 0 only for fields that actually change sign;
        % single-signed fields use parula (autoscaled) so their structure shows.
        Fmin = min(F(isfinite(F)));  Fmax = max(F(isfinite(F)));
        if ~isempty(Fmin) && Fmin < 0 && Fmax > 0
            a = max(abs([Fmin Fmax]));
            clim([-a a]);
            colormap(gca, local_diverging());
        else
            colormap(gca, parula);
        end

        cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.Interpreter    = 'latex';
        local_sciColorbar(cb, lbl);

        set(gca, 'TickLabelInterpreter', 'latex');
        xlabel(xlab, 'Interpreter', 'latex');
        ylabel(ylab, 'Interpreter', 'latex');
        xlim(xl);  ylim(yl);
    end

    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, 'baseFlow.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('plotBaseFlow: saved %s\n', out);
    end
end

% --- factor the order of magnitude out of the tick labels into the label ---
function local_sciColorbar(cb, baseLabel)
    ticks = cb.Ticks;
    m = max(abs(ticks));
    if ~(m > 0)
        cb.Label.String = baseLabel; return;
    end
    e = floor(log10(m));
    if e ~= 0
        cb.Ruler.Exponent = 0;
        cb.TickLabels = arrayfun(@(t) sprintf('$%.2g$', t/10^e), cb.Ticks, ...
                                 'UniformOutput', false);
        cb.Label.String = sprintf('%s $\\;(\\times 10^{%d})$', baseLabel, e);
    else
        cb.Label.String = baseLabel;
    end
end

% --- blue-white-red diverging colormap (no toolbox dependency) ---
function cmap = local_diverging()
    n = 128;
    t = linspace(0, 1, n)';
    b2w = [t,        t,        ones(n,1)];
    w2r = [ones(n,1), flipud(t), flipud(t)];
    cmap = [b2w; w2r];
end
