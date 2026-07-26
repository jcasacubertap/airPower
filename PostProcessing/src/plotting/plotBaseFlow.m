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
        xlab = '$S^{*} / \delta_0$';  ylab = '$y / \delta_0$';
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

    [X, Y] = plotCoords(sBF, inp);      % wall-fitted rectangle (unwraps curved TTCP wall)
    [xl, yl] = plotWindow(X, Y, inp);   % buffer-cut in x, near-wall window in y

    fig = figure('Name', 'Base flow', 'Color', 'w', ...
                 'Position', [60 60 900 220*nF]);
    for k = 1:nF
        f   = catalog{k,1};
        lbl = catalog{k,2};
        F   = sBF.(f);

        subplot(nF, 1, k);

        % --- colour limits from ONLY the plotted window (bufferFrac x-cut x
        %     near-wall y-cut). Changing bufferFrac re-derives the scale from the
        %     visible region, so its structure is isolated rather than washed out
        %     by out-of-window extremes (e.g. the outflow buffer). ---
        inwin = (X >= xl(1) & X <= xl(2) & Y >= yl(1) & Y <= yl(2) & isfinite(F));
        Fw = F(inwin);
        if isempty(Fw); Fw = F(isfinite(F)); end
        Fmin = min(Fw);  Fmax = max(Fw);
        pos = max(Fmax, 0);  neg = max(-Fmin, 0);

        % Diverging map (centred at 0) ONLY when genuinely two-signed: the
        % minority-sign magnitude is >= 5% of the majority. A tiny recirculation
        % negative then reads as single-signed -> parula fills the full range so
        % the boundary layer shows clearly.
        if pos > 0 && neg > 0 && min(pos,neg) >= 0.05*max(pos,neg)
            a = max(pos, neg);  lo = -a;  hi = a;  cmap = local_diverging();
        else
            lo = Fmin;  hi = Fmax;  cmap = parula;
        end
        if ~(hi > lo); hi = lo + eps(max(abs(lo),1)); end     % guard constant field

        % Clamp to the window range so the 30 bands + colourbar resolve the
        % visible region; out-of-window values saturate to the end colour.
        Fc = min(max(F, lo), hi);
        contourf(X, Y, Fc, 30, 'LineStyle', 'none'); hold on;
        clim([lo hi]);
        colormap(gca, cmap);

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
