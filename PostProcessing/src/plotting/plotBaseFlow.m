function fig = plotBaseFlow(sBF, inp, savedir, fieldSet)
% plotBaseFlow  Contour plots of the base-flow fields on the (x,y) grid.
%
%   fig = plotBaseFlow(sBF, inp)                         % show figure only
%   fig = plotBaseFlow(sBF, inp, savedir)                % also save baseFlow.png
%   fig = plotBaseFlow(sBF, inp, savedir, 'gradients')   % baseFlowGradients.png
%
% fieldSet = 'fields' (default): one stacked panel per available field among
% {u, v, w, p} (loadFields provides u,v,w; loadBF adds p).
% fieldSet = 'gradients': du/dx, du/dy, dw/dx, dw/dy, the Cartesian gradients both
% loaders return (loadBF computes them, loadFields reads DeHNSSo's), along
% the global X, Y; on a curved grid they are labelled X, Y, not x, y.
% Base-flow quantities carry a subscript B. Data are plotted as-is; only the
% labels reflect the units, read from sBF.units:
%   'nondim' (loadFields, DeHNSSo) -> lengths /delta_0, u_B/u_infty,
%                                     gradients delta_0/u_infty, ...
%   'SI'     (loadBF, OpenFOAM)    -> lengths [m], u_B [m/s], p_B [m^2/s^2],
%                                     gradients [1/s]
% Axes are x, y on a Cartesian grid and S*, n otherwise (see plotCoords).
% (OpenFOAM pressure is kinematic, p/rho, hence m^2/s^2.)
%
% Save behaviour: savedir non-empty -> PNG written there (terminal/dispatcher);
% empty/omitted -> figure left open for interactive inspection (direct MATLAB).

    if nargin < 3; savedir = ''; end
    if nargin < 4 || isempty(fieldSet); fieldSet = 'fields'; end

    nd = isfield(sBF, 'units') && strcmpi(sBF.units, 'nondim');

    if strcmpi(fieldSet, 'gradients')
        name = 'Base-flow gradients';  file = 'baseFlowGradients.png';
        % 4th column: map for a single-signed panel; du/dy and dw/dy look alike
        % (both wall-peaked shear), so they get different maps
        grads = {'ux','u','x','parula'; 'uy','u','y','parula'; ...
                 'wx','w','x','parula'; 'wy','w','y','turbo'};
        % The gradients are along the global X, Y. On a curved grid the axes are
        % the wall-fitted S*, n, so the labels say X, Y (not wall-aligned).
        [~, ~, lab0] = plotCoords(sBF, inp);
        if ~lab0.cartesian; grads(:,3) = upper(grads(:,3)); end
        catalog = cell(size(grads, 1), 3);
        for k = 1:size(grads, 1)
            d = sprintf('\\partial %s_{\\mathrm{B}} / \\partial %s', grads{k,2}, grads{k,3});
            if nd
                catalog(k,:) = {grads{k,1}, sprintf('$(%s)\\, \\delta_0 / u_\\infty$', d), grads{k,4}};
            else
                catalog(k,:) = {grads{k,1}, sprintf('$%s \\; [\\mathrm{1/s}]$', d), grads{k,4}};
            end
        end
    elseif nd
        name = 'Base flow';  file = 'baseFlow.png';
        catalog = { ...
            'u', '$u_{\mathrm{B}} / u_\infty$'; ...
            'v', '$v_{\mathrm{B}} / u_\infty$'; ...
            'w', '$w_{\mathrm{B}} / u_\infty$'; ...
            'p', '$p_{\mathrm{B}} / (\rho\, u_\infty^2)$'};
    else
        name = 'Base flow';  file = 'baseFlow.png';
        catalog = { ...
            'u', '$u_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'v', '$v_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'w', '$w_{\mathrm{B}} \; [\mathrm{m/s}]$'; ...
            'p', '$p_{\mathrm{B}} \; [\mathrm{m^2/s^2}]$'};
    end

    % Grid-size reference: both loaders supply the wall-fitted pair.
    if isfield(sBF, 's'); gsz = size(sBF.s); else; gsz = [0 0]; end

    have = false(size(catalog,1),1);
    for k = 1:size(catalog,1)
        f = catalog{k,1};
        have(k) = isfield(sBF, f) && ~isempty(sBF.(f)) && isequal(size(sBF.(f)), gsz);
    end
    catalog = catalog(have,:);
    nF = size(catalog,1);
    if nF == 0
        warning('plotBaseFlow:noFields', 'No base-flow %s to plot.', fieldSet); fig = []; return;
    end

    [X, Y, lab] = plotCoords(sBF, inp); % x,y on a Cartesian grid, else wall-fitted s,n
    [xl, yl] = plotWindow(X, Y, inp);   % buffer-cut in x, near-wall window in y

    fig = figure('Name', name, 'Color', 'w', ...
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
            if size(catalog, 2) >= 3 && ~isempty(catalog{k,3}); cmap = feval(catalog{k,3}); end
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
        xlabel(lab.x, 'Interpreter', 'latex');
        ylabel(lab.y, 'Interpreter', 'latex');
        xlim(xl);  ylim(yl);
    end

    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, file);
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
