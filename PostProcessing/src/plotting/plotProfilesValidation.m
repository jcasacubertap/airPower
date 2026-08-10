function figs = plotProfilesValidation(sBF, sPert, inp, savedir)
% plotProfilesValidation  Dimensional w-profile comparison vs PIV (DeHNSSo).
%
%   figs = plotProfilesValidation(sBF, sPert, inp, savedir)
%
% Produces TWO figures, both with the same layout and styling (DeHNSSo = black
% line with dot markers, Experimental (PIV) = red filled circles, w [m/s] vs wall
% distance [m], LaTeX serif; a 2-row grid, one column per station -- top row the
% base-flow w, bottom row the fundamental perturbation w RMS (|w~|/sqrt(2), the
% z-RMS of a single mode, the same quantity the PIV reports)):
%
%   profiles_w_validation.png        ZOOM  -- the x/c window inp.valXcZoom
%                                    (default [10 25] %). Fixed 4 mm y-range.
%
%   profiles_w_validation_broad.png  BROAD -- the whole valid domain, inlet to
%                                    outflow-buffer start. The y-range is shared
%                                    and sized on the thickest delta_99 among the
%                                    stations, since the BL grows over the domain.
%
% Both draw their stations the same way (see `stations`), which is what
% inp.valPIV switches:
%   valPIV = true   the PIV stations guide the plot -- the stations ARE PIV
%                   stations, so every panel has an exact experimental overlay.
%   valPIV = false  no experimental data; the stations are grid columns spread
%                   evenly over the window, DeHNSSo alone.
% Both figures are cut at the outflow buffer (inp.ro.bufferFrac), as everywhere
% else in PostProcessing. Panel titles carry x/c and, in brackets, S*/delta_0 --
% the same S*/delta_0 the plotProfiles u-figure uses.
%
% The x/c <-> column mapping is case-specific (see gridMap):
%   DFP  (flat reference plate):  column x_DFP -> S = x_DFP + xInlet -> x/c via the
%        airfoilFlowData BL.S/BL.x/BL.c table. Wall distance is the grid's
%        wall-normal coordinate.
%   TTCP (curved airfoil wall):   the stability-grid wall row IS the airfoil
%        surface, so each column's x/c is recovered by the inverse rigid transform
%        (undo AoA rotation + centre offset) of its physical wall point. Wall
%        distance is the arc-length up the (body-fitted) grid column.

    figs = gobjects(0);
    if nargin < 4; savedir = ''; end

    if ~isfield(sPert,'uref') || ~isfield(sPert,'lref') || isempty(sPert.uref) || isempty(sPert.lref)
        warning('plotProfilesValidation:noRef', ...
                'StabGrid Uref/lref missing; cannot dimensionalize DeHNSSo data.'); return;
    end
    uref = sPert.uref;  lref = sPert.lref;

    % ---- column -> x/c and column -> wall-distance maps (case-specific) ----
    G = gridMap(sBF, inp, lref);
    if isempty(G); return; end

    % ---- PIV Gen/Case data (optional: the broad figure works without it) ----
    o = loadPIV(inp);

    % ---- fundamental mode (or the first in modeIdx) ----
    Nmode = size(sPert.w, 1);
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx)
        m = inp.modeIdx(1);
    else
        m = min(2, Nmode);
    end
    m = max(2, min(m, Nmode));                 % (0,0) has no PIV counterpart

    % ---- figure 1: zoom on the requested x/c window ----
    Sz = stations(G, o, inp, xcZoom(inp));
    if ~isempty(Sz)
        f = plotStations(Sz, sBF, sPert, o, inp, m, uref, 0.004);
        figs(end+1) = f;
        saveFig(f, savedir, 'profiles_w_validation.png', m, Sz);
    end

    % ---- figure 2: broad view over the whole valid domain ----
    Sb = stations(G, o, inp, []);
    if ~isempty(Sb)
        f = plotStations(Sb, sBF, sPert, o, inp, m, uref, autoYTop(Sb, sBF));
        figs(end+1) = f;
        saveFig(f, savedir, 'profiles_w_validation_broad.png', m, Sb);
    end
end

% ======================================================================
%  PIV loader (shared) -- non-fatal: [] means "no experimental data"
% ======================================================================
function o = loadPIV(inp)
    o = [];
    valdir = fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'Validation', ...
                      sprintf('Gen%d', inp.valGen), 'Experimental', sprintf('Case%d', inp.valCase));
    mats = dir(fullfile(valdir, '*.mat'));
    if isempty(mats)
        warning('plotProfilesValidation:noPIV', ...
                'No PIV .mat in %s — DeHNSSo profiles will be plotted alone.', valdir); return;
    end
    o = load(fullfile(valdir, mats(1).name)).output;
end

% ======================================================================
%  Column maps: x/c per column, and wall distance [m] up a given column
%
%  G.xic       1 x Nx, chord fraction of each grid column [-]
%  G.sd0       1 x Nx, streamwise wall arc-length from the inlet, in delta_0
%              (same S*/delta_0 the plotProfiles u-figure titles use)
%  G.wallDist  @(c) -> Ny x 1 wall distance [m] for column c
%  G.Nx        number of columns
% ======================================================================
function G = gridMap(sBF, inp, lref)
    G = [];
    switch upper(inp.caseType)
        case 'DFP'
            % x/c from the airfoilFlowData arc-length table: column x_DFP [m]
            % -> S = x_DFP + xInlet -> x (chordwise) -> x/c. This is the inverse
            % of the forward x/c -> x_DFP map, so nearest-column lookups done in
            % x/c space land on the same columns as before.
            fd = dir(fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'airfoilFlowData', '*.mat'));
            if isempty(fd)
                warning('plotProfilesValidation:noBL', 'No airfoilFlowData .mat for x/c<->S.'); return;
            end
            BL = load(fullfile(fd(1).folder, fd(1).name)).BL;
            x_ref = double(BL.x(:));  S_ref = double(BL.S(:));  c_ref = double(BL.c);

            xg    = sBF.x(1,:) * lref;                 % columns in physical x_DFP [m]
            xchord = interp1(S_ref, x_ref, xg + inp.xInlet, 'linear', 'extrap');
            G.xic = xchord(:).' / c_ref;

            G.sd0 = sBF.x(1,:);                        % flat wall: x from the inlet IS S*, in delta_0

            y_m = sBF.y(:,1) * lref;                   % flat wall: same for every column
            G.wallDist = @(c) y_m;

        case 'TTCP'
            req = {'airfoilChord','airfoilAlphaDeg','airfoilXCenter','airfoilYCenter'};
            for r = req
                if ~isfield(inp, r{1})
                    warning('plotProfilesValidation:noGeom', ...
                            ['TTCP validation needs airfoil geometry (%s) from inputs.jl. ', ...
                             'Re-generate inputs_gen.m (julia run.jl PostProcessing config).'], r{1});
                    return;
                end
            end
            chord = inp.airfoilChord;  xC = inp.airfoilXCenter;  yC = inp.airfoilYCenter;
            ca = cosd(inp.airfoilAlphaDeg);  sa = sind(inp.airfoilAlphaDeg);

            X = sBF.x * lref;  Y = sBF.y * lref;       % physical grid [m], wall at row end
            Xw = X(end,:);  Yw = Y(end,:);             % airfoil surface (wall row)

            % Forward map is x = ca*xb + sa*yb + xC with xb = (x/c - 0.5)*chord.
            % Undo the rigid part (rotation + offset) to recover xb, hence x/c; the
            % airfoil thickness (yb) drops out — it is baked into the wall point.
            G.xic = (ca*(Xw - xC) - sa*(Yw - yC)) / chord + 0.5;

            % S*/delta_0: arc-length along the wall from the inlet column, in
            % delta_0 -- sBF.x/sBF.y are already normalized, so no lref here.
            % Same construction as plotCoords, so the two figures agree.
            G.sd0 = [0, cumsum(hypot(diff(sBF.x(end,:)), diff(sBF.y(end,:))))];

            % wall-normal distance up a column: arc-length from the wall (row end)
            G.wallDist = @(c) wallArc(X(:,c), Y(:,c));

        otherwise
            warning('plotProfilesValidation:notWired', ...
                    'w-vs-PIV station mapping is only wired for DFP/TTCP (got ''%s'').', inp.caseType);
            return;
    end
    G.Nx = numel(G.xic);
end

% arc-length up one grid column, 0 at the wall (row end)
function n = wallArc(Xc, Yc)
    seg = hypot(diff(Xc), diff(Yc));               % (Ny-1) x 1, between-row lengths
    n   = [flipud(cumsum(flipud(seg))); 0];        % Ny x 1
end

% ======================================================================
%  Station selection — shared by both figures, they differ only in the x/c
%  window they are given (the zoom gets inp.valXcZoom, the broad gets []).
%
%  valPIV = true   the PIV stations GUIDE the plot: the stations ARE PIV
%                  stations, the NSTA of them nearest to evenly-spaced x/c
%                  targets across the window. Every panel then carries an
%                  exact, like-for-like experimental overlay — no station is
%                  ever compared against PIV taken somewhere else.
%  valPIV = false  no experimental data: NSTA grid columns spread evenly over
%                  the window, DeHNSSo alone.
%
%  Either way the window is clipped to the valid domain, i.e. up to the
%  outflow-buffer start (inp.ro.bufferFrac).
% ======================================================================
function S = stations(G, o, inp, xcWin)
    S = emptyStations();

    % valid domain: inlet -> buffer start, in % chord
    ib = colAt(G, bufferFrac(inp));
    xcLo = min(G.xic(1:ib))*100;  xcHi = max(G.xic(1:ib))*100;
    if ~isempty(xcWin)
        xcLo = max(xcLo, xcWin(1));  xcHi = min(xcHi, xcWin(2));
    end
    if xcHi < xcLo
        warning('plotProfilesValidation:emptyWindow', ...
                'Requested x/c window is outside the valid DeHNSSo domain.'); return;
    end

    if usePIV(inp) && ~isempty(o)
        % ---- PIV stations guide the plot ----
        xc_all = cellfun(@(v) double(v(1)), o.xc(:));
        [xc_all, ord] = sort(xc_all);
        valid = xc_all >= xcLo & xc_all <= xcHi;
        if ~any(valid)
            warning('plotProfilesValidation:noValidStation', ...
                    'No PIV station in x/c %.1f–%.1f %%.', xcLo, xcHi); return;
        end
        xc_v = xc_all(valid);  ord_v = ord(valid);

        % the PIV stations nearest to NSTA evenly-spaced x/c targets. On the
        % dense part of the PIV list this is the same as taking every k-th
        % station; where the list is sparse it keeps the spread even instead
        % of clustering the panels where PIV happens to be dense.
        tgt = linspace(min(xc_v), max(xc_v), NSTA);
        idx = arrayfun(@(t) findNearest(xc_v, t), tgt);
        idx = unique(idx, 'stable');

        for j = 1:numel(idx)
            k = idx(j);
            c = findNearest(G.xic, xc_v(k)/100);
            S(j) = mkStation(G, c, xc_v(k), '%.0f', ord_v(k));
        end
    else
        % ---- no experimental data: evenly-spaced grid columns ----
        cLo  = findNearest(G.xic, xcLo/100);
        cHi  = findNearest(G.xic, xcHi/100);
        cols = unique(round(linspace(cLo, cHi, NSTA)));
        for j = 1:numel(cols)
            c = cols(j);
            S(j) = mkStation(G, c, G.xic(c)*100, '%.1f', []);
        end
    end
end

% one station record; `label` is the x/c text, `sd0` the S*/delta_0 shown
% alongside it in the panel title
function s = mkStation(G, c, xcPct, fmt, korig)
    s = struct('xcPct', xcPct, 'label', sprintf(fmt, xcPct), ...
               'sd0', G.sd0(c), 'col', c, 'n', G.wallDist(c), 'korig', korig);
end

function i = findNearest(v, t)
    [~, i] = min(abs(v - t));
end

% ======================================================================
%  Plotting (shared by both figures)
% ======================================================================
function fig = plotStations(S, sBF, sPert, o, inp, m, uref, yTop)

    ns = numel(S);
    cSolver = [0 0 0];  cExp = [0.75 0.08 0.08];

    % put the legend on the first station that actually has an experimental
    % overlay (in the broad figure that is usually not the first column)
    jLeg = find(~cellfun(@isempty, {S.korig}), 1);
    if isempty(jLeg); jLeg = 1; end

    fig = figure('Color', 'w', 'Position', [40 60 330*ns 650]);
    for j = 1:ns
        c = S(j).col;  k = S(j).korig;  yprof = S(j).n;   % wall distance [m], per column

        % ---- top: base-flow w ----
        ax = subplot(2, ns, j); hold(ax, 'on');
        h1 = plot(ax, sBF.w(:,c)*uref, yprof, '-', 'Color', cSolver, ...
                  'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 7);
        h2 = [];
        if ~isempty(k)
            wavg = local_rowmean(o.w_m_mean{k});  ye = double(o.y{k});  ye = ye(:,1)*1e-3;
            h2 = plot(ax, wavg, ye, 'o', 'MarkerSize', 3.5, ...
                      'MarkerFaceColor', cExp, 'MarkerEdgeColor', cExp);
        end
        local_style(ax); ylim(ax, [0 yTop]);
        xlabel(ax, '$w_{\mathrm{B}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        if j == 1; ylabel(ax, '$\mathrm{Wall\ distance\ [m]}$', 'Interpreter', 'latex'); end
        % two lines: a single line is wide enough to run into the axes' 1e-3
        % exponent label at the top left
        title(ax, {sprintf('$x/c = %s\\%%$', S(j).label), ...
                   sprintf('$(S^{*}/\\delta_0 = %.0f)$', S(j).sd0)}, ...
              'Interpreter', 'latex', 'FontSize', 11);
        if j == jLeg
            local_legend(ax, h1, h2, 'northwest');
        end

        % ---- bottom: fundamental perturbation w RMS ----
        ax = subplot(2, ns, ns + j); hold(ax, 'on');
        wrms = abs(squeeze(sPert.w(m,:,c))) / sqrt(2) * uref;
        g1 = plot(ax, wrms, yprof, '-', 'Color', cSolver, ...
                  'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 7);
        fn = sprintf('w_pert_m_prof_rms_%02d', m-1);  yn = sprintf('y_prof_rms_%02d', m-1);
        g2 = [];
        if ~isempty(k) && isfield(o, fn) && isfield(o, yn)
            g2 = plot(ax, double(o.(fn){k}), double(o.(yn){k})*1e-3, 'o', 'MarkerSize', 3.5, ...
                      'MarkerFaceColor', cExp, 'MarkerEdgeColor', cExp);
        end
        local_style(ax); ylim(ax, [0 yTop]);
        xlabel(ax, '$w_{\mathrm{rms}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        if j == 1; ylabel(ax, '$\mathrm{Wall\ distance\ [m]}$', 'Interpreter', 'latex'); end
        if j == jLeg
            local_legend(ax, g1, g2, 'northeast');
        end
    end
end

function saveFig(fig, savedir, fname, m, S)
    if isempty(savedir); return; end
    if ~exist(savedir, 'dir'); mkdir(savedir); end
    out = fullfile(savedir, fname);
    exportgraphics(fig, out, 'Resolution', 150);
    nPIV = nnz(~cellfun(@isempty, {S.korig}));
    fprintf(['plotProfilesValidation: saved %s (fundamental mode (0,%d), ', ...
             '%d stations, %d with PIV)\n'], out, m-1, numel(S), nPIV);
end

% ======================================================================
%  Settings
%  Only the zoom window comes from inputs.jl (inp.valXcZoom); the broad
%  figure has no knobs -- it deliberately spans the whole numerical domain.
% ======================================================================
% stations (columns) per figure
function n = NSTA
    n = 6;
end

% whether experimental data drives the station choice (inp.valPIV, from the
% VAL block of inputs.jl). False -> DeHNSSo-only figures on grid columns.
function t = usePIV(inp)
    t = true;
    if isfield(inp, 'valPIV') && ~isempty(inp.valPIV)
        t = logical(inp.valPIV);
    end
end

% streamwise fraction of the domain where the outflow buffer starts
function f = bufferFrac(inp)
    f = 0.85;
    if isfield(inp,'ro') && isfield(inp.ro,'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        f = inp.ro.bufferFrac;
    end
end

% zoom figure: PIV-station x/c window [%] as [lo hi]; [] -> no window.
% Also accepts the legacy scalar inp.valXcMin (lower bound only).
function w = xcZoom(inp)
    w = [];
    if isfield(inp, 'valXcZoom') && ~isempty(inp.valXcZoom)
        w = double(inp.valXcZoom(:)).';
    elseif isfield(inp, 'valXcMin') && ~isempty(inp.valXcMin)
        w = double(inp.valXcMin);
    end
    if numel(w) == 1; w = [w, Inf]; end
end

% ======================================================================
%  Shared helpers
% ======================================================================
function S = emptyStations()
    % same field order as mkStation
    S = struct('xcPct', {}, 'label', {}, 'sd0', {}, 'col', {}, 'n', {}, 'korig', {});
end

% column index at a given fraction of the streamwise extent
function ib = colAt(G, frac)
    ib = min(G.Nx, max(2, round(frac * G.Nx)));
end

% shared y-limit for the broad figure: 1.5 x the thickest boundary layer among
% the plotted stations (u-based delta_99), capped by the grid height.
function yTop = autoYTop(S, sBF)
    d99 = 0;  nMax = 0;
    for j = 1:numel(S)
        n = S(j).n;  u = abs(sBF.u(:, S(j).col));
        nMax = max(nMax, max(n));
        ue = max(u);
        if ~isfinite(ue) || ue <= 0; continue; end
        % rows run free-stream (1) -> wall (end), so the LAST row still at
        % 0.99*ue is the boundary-layer edge
        i99 = find(u >= 0.99*ue, 1, 'last');
        if ~isempty(i99); d99 = max(d99, n(i99)); end
    end
    yTop = min(1.5*d99, nMax);
    if ~isfinite(yTop) || yTop <= 0; yTop = nMax; end
end

% --- z-average a (ny x nz) PIV field to a 1-D profile, ignoring NaN/0 ---
function wavg = local_rowmean(W)
    W = double(W);  W(W == 0) = NaN;
    wavg = mean(W, 2, 'omitnan');
end

% --- legend with the experimental entry only when there is one ---
function local_legend(ax, hNum, hExp, loc)
    if isempty(hExp)
        legend(ax, hNum, {'\textrm{DeHNSSo}'}, ...
               'Interpreter', 'latex', 'Location', loc);
    else
        legend(ax, [hNum hExp], {'\textrm{DeHNSSo}','\textrm{Experimental}'}, ...
               'Interpreter', 'latex', 'Location', loc);
    end
end

% --- common axis styling (serif/LaTeX, grid, box) ---
function local_style(ax)
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 11, ...
            'GridAlpha', 0.3, 'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
end
