function figs = plotProfilesValidation(sBF, sPert, inp, savedir)
% plotProfilesValidation  Dimensional w-profile comparison vs PIV (DeHNSSo).
%
%   figs = plotProfilesValidation(sBF, sPert, inp, savedir)
%
% Produces TWO figures, both with the same layout and styling (DeHNSSo = black
% line, Experimental (PIV) = red line; BOTH carry one white-filled open circle
% per sample point, so each curve shows its own wall-normal resolution -- the
% black circles are the stability grid's y-nodes, the red ones the PIV
% measurement points. w [m/s] vs wall distance in MILLIMETRES, LaTeX serif;
% a 2 x nStation tiledlayout with one legend under
% the whole figure -- top row the MEAN flow w (base flow + the (0,0) mean-flow
% distortion of the stability solution, see meanW; the base flow alone is kept as
% a grey dashed reference), bottom row the fundamental perturbation w RMS
% (|w~|/sqrt(2), the z-RMS of a single mode, the same quantity the PIV reports)):
%
%   profiles_w_validation.png        ZOOM  -- the x/c window inp.valXcZoom
%                                    (default [10 25] %).
%
%   profiles_w_validation_broad.png  BROAD -- the whole valid domain, inlet to
%                                    outflow-buffer start.
%
% Both figures auto-size their y-range (autoYTop: 1.5 x the perturbation extent
% over their stations, capped by the grid; override with inp.valYTop [mm]) and
% SHARE it across both rows, so the eye carries vertically from the mean flow to
% the perturbation it drives.
% They also share an x-limit along each row (inp.valShareX, default true), which
% is what makes the streamwise growth readable station-to-station; set it false
% for per-panel autoscaling if a single station's shape needs resolving.
%
% These are DIAGNOSTIC figures, so each perturbation panel carries its own
% numbers: the peak w_rms of each, their ratio N/P (coloured amber past 10% and
% red past 25% error) and the wall distance of each peak -- enough to tell an
% amplitude deficit from a mode-shape error without leaving the figure.
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
    % The y-range is auto-sized here too (it used to be a hard-coded 4 mm, which
    % left two thirds of every panel empty); see autoYTop.
    Sz = stations(G, o, inp, xcZoom(inp));
    if ~isempty(Sz)
        f = plotStations(Sz, sBF, sPert, o, inp, m, uref, autoYTop(Sz, sPert, m, uref, inp));
        figs(end+1) = f;
        saveFig(f, savedir, 'profiles_w_validation.png', m, Sz);
    end

    % ---- figure 2: broad view over the whole valid domain ----
    Sb = stations(G, o, inp, []);
    if ~isempty(Sb)
        f = plotStations(Sb, sBF, sPert, o, inp, m, uref, autoYTop(Sb, sPert, m, uref, inp));
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
%
%  Layout: a 2 x ns tiledlayout, top row the mean flow, bottom row the
%  fundamental w RMS, ONE legend under the whole figure. Wall distance is in
%  MILLIMETRES and both rows share the same y-limit (yTop, passed in metres),
%  so the eye carries vertically from the mean flow to the perturbation it
%  drives; only the left column keeps its y tick labels.
%
%  Both rows also share their x-limit ACROSS the stations (inp.valShareX,
%  default true). Per-panel autoscaling made every station look alike and hid
%  the streamwise growth, which is the thing being validated.
%
%  Diagnostic annotations (this is a diagnostic figure, not a paper figure):
%    top row     max|dw| between HNS and PIV over the overlapping y-range
%    bottom row  peak w_rms, HNS / PIV, their ratio N/P (coloured by how far
%                it is from 1), and the peak wall distance of each. Together
%                these separate an amplitude deficit from a mode-shape error
%                without leaving the figure.
% ======================================================================
function fig = plotStations(S, sBF, sPert, o, inp, m, uref, yTopM)

    ns   = numel(S);
    yTop = yTopM * 1e3;                              % plot in mm throughout
    cSolver = [0 0 0];  cExp = [0.75 0.08 0.08];  cBase = [0.55 0.55 0.55];

    % ---- pass 1: collect every curve, so the shared limits and the
    %      annotations can be built before anything is drawn ----
    D = collectStations(S, sBF, sPert, o, m, uref, yTop);

    % shared x-limits per row (0 -> max over all stations, padded)
    shareX = shareXLim(inp);
    xlM = rowXLim(D, {'wMean','wBase','wMeanExp'}, D(1).yTop);
    xlR = rowXLim(D, {'wRms','wRmsExp'},           D(1).yTop);

    fig = figure('Color', 'w', 'Position', [40 60 300*ns 700]);
    t = tiledlayout(fig, 2, ns, 'TileSpacing', 'compact', 'Padding', 'compact');
    L = struct('solver', [], 'base', [], 'piv', []);

    for j = 1:ns
        d = D(j);

        % ---- top: mean flow w = base flow + mean-flow distortion (0,0) ----
        % The PIV overlay is a MEAN velocity, so the like-for-like numerical
        % quantity is w_B + w_(0,0), not w_B alone (see meanW).
        ax = nexttile(t, j); hold(ax, 'on');
        h = profCurve(ax, d.wMean, d.y, cSolver, 1.3, MS);
        L.solver = firstOf(L.solver, h);
        if ~isempty(d.wMeanExp)
            L.piv = firstOf(L.piv, profCurve(ax, d.wMeanExp, d.yMeanExp, cExp, 1.1, MS));
        end
        % base flow drawn LAST, so the dashed reference sits on top instead of
        % being buried under the two marker trains it is meant to be read against
        if d.hasMFD
            h = plot(ax, d.wBase, d.y, '--', 'Color', cBase, 'LineWidth', 1.1);
            L.base = firstOf(L.base, h);
        end
        finishPanel(ax, j, yTop, shareX, xlM);
        if d.hasMFD
            xlabel(ax, '$w_{\mathrm{B}} + w_{(0,0)} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        else
            xlabel(ax, '$w_{\mathrm{B}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        end
        % single line now: in mm there is no 1e-3 exponent label to collide with
        title(ax, sprintf('$x/c = %s\\%%$ \\ \\ $(S^{*}/\\delta_0 = %.0f)$', ...
                          S(j).label, S(j).sd0), 'Interpreter', 'latex', 'FontSize', 11);

        % ---- bottom: fundamental perturbation w RMS ----
        ax = nexttile(t, ns + j); hold(ax, 'on');
        profCurve(ax, d.wRms, d.y, cSolver, 1.3, MS);
        if ~isempty(d.wRmsExp)
            profCurve(ax, d.wRmsExp, d.yRmsExp, cExp, 1.1, MS);
        end
        finishPanel(ax, j, yTop, shareX, xlR);
        xlabel(ax, '$w_{\mathrm{rms}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        annot(ax, peakText(d), ratioColor(d.ratio));
    end

    ylabel(t, '$\mathrm{Wall\ distance\ [mm]}$', 'Interpreter', 'latex', 'FontSize', 12);
    sharedLegend(t, L);
end

% ---- one pass over the stations: everything needed to draw and to annotate --
function D = collectStations(S, sBF, sPert, o, m, uref, yTop)
    fn = sprintf('w_pert_m_prof_rms_%02d', m-1);
    yn = sprintf('y_prof_rms_%02d', m-1);
    D  = repmat(emptyCurves(yTop), 1, numel(S));

    for j = 1:numel(S)
        c = S(j).col;  k = S(j).korig;
        d = emptyCurves(yTop);
        d.y = S(j).n(:) * 1e3;                       % wall distance [mm]
        [d.wMean, d.wBase, d.hasMFD] = meanW(sBF, sPert, c, uref);
        wr = abs(squeeze(sPert.w(m,:,c))) / sqrt(2) * uref;
        d.wRms = wr(:);

        if ~isempty(k)
            ye = double(o.y{k});                     % PIV y is already in mm
            d.yMeanExp = ye(:,1);
            d.wMeanExp = local_rowmean(o.w_m_mean{k});
            if isfield(o, fn) && isfield(o, yn)
                we = double(o.(fn){k});  d.wRmsExp = we(:);
                ye = double(o.(yn){k});  d.yRmsExp = ye(:);
            end
        end

        % diagnostics, all restricted to the plotted window
        [d.pkN, d.ypkN] = peakIn(d.wRms,    d.y,       yTop);
        [d.pkP, d.ypkP] = peakIn(d.wRmsExp, d.yRmsExp, yTop);
        if isfinite(d.pkN) && isfinite(d.pkP) && d.pkP > 0; d.ratio = d.pkN / d.pkP; end
        D(j) = d;
    end
end

function d = emptyCurves(yTop)
    d = struct('y', [], 'yTop', yTop, 'hasMFD', false, ...
               'wMean', [], 'wBase', [], 'wRms', [], ...
               'wMeanExp', [], 'yMeanExp', [], 'wRmsExp', [], 'yRmsExp', [], ...
               'pkN', NaN, 'ypkN', NaN, 'pkP', NaN, 'ypkP', NaN, 'ratio', NaN);
end

% ---- shared x-limit for one row: [0, max over all stations] + 5% headroom ----
function xl = rowXLim(D, fields, yTop)
    hi = 0;
    for j = 1:numel(D)
        yy = D(j).y;
        for f = fields
            v = D(j).(f{1});
            if isempty(v); continue; end
            % pair each curve with its own y (PIV curves have their own grid)
            switch f{1}
                case 'wMeanExp'; yv = D(j).yMeanExp;
                case 'wRmsExp';  yv = D(j).yRmsExp;
                otherwise;       yv = yy;
            end
            in = yv(:) <= yTop & isfinite(v(:));
            if any(in); hi = max(hi, max(v(in))); end
        end
    end
    if ~isfinite(hi) || hi <= 0; xl = []; return; end
    xl = [0, 1.05*hi];
end

% ---- a profile curve: line + one open circle per SAMPLE POINT ----
% Both HNS and PIV are drawn this way, so each curve shows its own wall-normal
% resolution: the HNS circles are the stability grid's y-nodes, the PIV circles
% its measurement points. The markers are white-filled with a coloured edge so a
% dense cluster stays readable as individual points instead of a solid blob, and
% so the two sets remain distinguishable where they overlap.
function h = profCurve(ax, x, y, c, lw, ms)
    ok = isfinite(x(:)) & isfinite(y(:));
    h  = plot(ax, x(ok), y(ok), '-o', 'Color', c, 'LineWidth', lw, ...
              'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c);
end

% ---- common per-panel finishing: style, limits, tick-label suppression ----
function finishPanel(ax, j, yTop, shareX, xl)
    local_style(ax);
    ylim(ax, [0 yTop]);
    if shareX && ~isempty(xl); xlim(ax, xl); end
    if j > 1; set(ax, 'YTickLabel', []); end
end

% ---- top-right annotation box (perturbation row only) ----
function annot(ax, lines, col)
    if isempty(lines); return; end
    text(ax, 0.96, 0.97, lines, 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'Interpreter', 'latex', 'FontSize', 9, 'Color', col, ...
         'BackgroundColor', 'w', 'Margin', 1.5);
end

% ---- bottom-row diagnostic text: peaks, ratio, peak heights ----
function lines = peakText(d)
    lines = {};
    if ~isfinite(d.pkN); return; end
    if isfinite(d.pkP)
        lines = {sprintf('$\\hat{w}_{\\max}:\\ %.2f\\,/\\,%.2f$', d.pkN, d.pkP), ...
                 sprintf('$N/P = %.2f$', d.ratio), ...
                 sprintf('$y_{\\max}:\\ %.2f\\,/\\,%.2f$', d.ypkN, d.ypkP)};
    else
        lines = {sprintf('$\\hat{w}_{\\max} = %.2f$', d.pkN), ...
                 sprintf('$y_{\\max} = %.2f$', d.ypkN)};
    end
end

% ratio colouring: black within 10%, amber to 25%, red beyond
function c = ratioColor(r)
    c = [0 0 0];
    if ~isfinite(r); return; end
    e = abs(1 - r);
    if     e > 0.25; c = [0.75 0.08 0.08];
    elseif e > 0.10; c = [0.80 0.45 0.00];
    end
end

% ---- one legend for the whole figure, under the tiles ----
function sharedLegend(t, L)
    hs   = {L.solver, L.base, L.piv};
    labs = {'\textrm{HNS (DeHNSSo)}', '$w_{\mathrm{B}}$ \textrm{(base flow)}', ...
            '\textrm{PIV}'};
    keep = ~cellfun(@isempty, hs);
    if ~any(keep); return; end
    lg = legend([hs{keep}], labs(keep), 'Interpreter', 'latex', ...
                'Orientation', 'horizontal', 'FontSize', 11);
    lg.Layout.Tile = 'south';
end

function h = firstOf(h, cand)
    if isempty(h); h = cand; end
end

function saveFig(fig, savedir, fname, m, S)
    if isempty(savedir); return; end
    if ~exist(savedir, 'dir'); mkdir(savedir); end
    out = fullfile(savedir, fname);
    % drop the interactive axes toolbars: exportgraphics otherwise warns that
    % they may appear in the image (they do not here, but the warning is noise)
    set(findall(fig, 'Type', 'axes'), 'Toolbar', []);
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

% sample-point marker size, the SAME in both rows: the circles carry a meaning
% (one per grid node / measurement point), so a size difference between the rows
% would read as a difference in the data rather than in the styling.
function s = MS
    s = 3;
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

% Shared y-limit (BOTH rows, BOTH figures), as a wall distance [m].
%
% Sized on the PERTURBATION extent: 1.5 x the highest point at which w_rms is
% still 5% of its station peak, maxed over the plotted stations and capped by
% the grid height. inp.valYTop [mm] overrides it outright.
%
% Deliberately NOT a delta_99 criterion. The obvious one (deepest row still at
% 99% of the column's max u) is wrong for this geometry: near the leading edge
% the outer flow is still accelerating with height, so max(u) sits at the top of
% the domain and the criterion returns where the INVISCID profile flattens, not
% the boundary-layer edge -- it gave 4 mm at x/c = 12% against a profile that
% plateaus by 1.3 mm, which is what used to inflate the broad figure to 6 mm.
% The perturbation extent needs no edge detection and is what both rows are
% actually about.
function yTop = autoYTop(S, sPert, m, uref, inp)
    if isfield(inp, 'valYTop') && ~isempty(inp.valYTop)
        yTop = double(inp.valYTop) * 1e-3;  return;       % inputs.jl gives mm
    end
    ext = 0;  nMax = 0;
    for j = 1:numel(S)
        n = S(j).n;
        nMax = max(nMax, max(n));
        wr = abs(squeeze(sPert.w(m,:,S(j).col))) / sqrt(2) * uref;
        wr = wr(:);
        pk = max(wr);
        if ~isfinite(pk) || pk <= 0; continue; end
        % rows run free-stream (1) -> wall (end), so n DESCENDS with the row
        % index: the FIRST row above the threshold is the highest such point.
        i = find(wr >= 0.05*pk, 1, 'first');
        if ~isempty(i); ext = max(ext, n(i)); end
    end
    yTop = min(1.5*ext, nMax);
    if ~isfinite(yTop) || yTop <= 0; yTop = nMax; end
end

% peak of a profile and the wall distance where it occurs, restricted to the
% plotted window (so free-stream noise cannot win the max)
function [pk, ypk] = peakIn(v, y, yTop)
    pk = NaN;  ypk = NaN;
    if isempty(v); return; end
    v = v(:);  y = y(:);
    in = y <= yTop & isfinite(v);
    if ~any(in); return; end
    [pk, i] = max(v(in));
    yy = y(in);  ypk = yy(i);
end

% whether the stations share one x-limit per row (inp.valShareX, default true).
% Shared limits make the streamwise growth readable; false restores per-panel
% autoscaling, which resolves each station's shape at the cost of that.
function t = shareXLim(inp)
    t = true;
    if isfield(inp, 'valShareX') && ~isempty(inp.valShareX)
        t = logical(inp.valShareX);
    end
end

% ======================================================================
%  Mean flow at one column: base flow + mean-flow distortion (mode (0,0)).
%
%  Why: the PIV overlay (o.w_m_mean) is a MEASURED mean velocity, so it
%  already contains whatever mean-flow distortion the disturbance drove. The
%  like-for-like numerical quantity is therefore w_B + w_(0,0), not w_B alone
%  -- the base flow the HNS was linearized about knows nothing about the MFD.
%
%  Amplitude convention -- no factor of 2 here, and that is deliberate:
%  DeHNSSo (src/io/postprocess.m) stores A = 2|u_max| for the oscillatory
%  modes, because u'(z) = 2|a|cos(bz) peaks at 2|a|, but A = |u_max| for the
%  MFD, which is constant in z. importData folds the FULL A into sPert for
%  every mode alike, so row (0,0) already IS the physical distortion field and
%  is added straight to the base flow.
%
%  Row 1 is the (0,0) mode: postprocess.m deletes the conjugate half, leaving
%  (0,0) first. We still verify it via beta = omega = 0 when those vectors are
%  present, and fall back to the base flow alone if they say otherwise.
%
%  Returns wTot (Ny x 1, [m/s]), wBase (same, base flow alone) and hasMFD.
% ======================================================================
function [wTot, wBase, hasMFD] = meanW(sBF, sPert, c, uref)
    wBase  = sBF.w(:,c) * uref;
    wTot   = wBase;
    hasMFD = false;
    if isempty(sPert) || ~isfield(sPert,'w') || size(sPert.w,1) < 1; return; end

    % (0,0) identification: zero spanwise wavenumber AND zero frequency.
    for f = {'beta','omega'}
        if isfield(sPert, f{1}) && ~isempty(sPert.(f{1})) && sPert.(f{1})(1) ~= 0
            warnOnce('plotProfilesValidation:noMFD', ...
                     ['Mode row 1 has %s = %g, not 0 — it is not the (0,0) ', ...
                      'mean-flow distortion. Plotting the base flow alone.'], ...
                     f{1}, sPert.(f{1})(1));
            return;
        end
    end

    wm = squeeze(sPert.w(1,:,c));  wm = wm(:);          % Ny x 1, physical MFD
    % A zero-beta, zero-omega mode of a real operator is real; keep only the
    % real part and flag it if the imaginary part is not round-off.
    aIm = max(abs(imag(wm)));  aRe = max(abs(real(wm)));
    if aIm > 1e-6 * max(aRe, eps)
        warnOnce('plotProfilesValidation:complexMFD', ...
                 ['(0,0) mean-flow distortion has a non-negligible imaginary ', ...
                  'part (|Im|/|Re| = %.2g); using real() — check the mode ordering.'], ...
                 aIm/max(aRe, eps));
    end
    wTot   = wBase + real(wm) * uref;
    hasMFD = true;
end

% --- warn at most once per identifier per MATLAB session (12 panels/figure) ---
function warnOnce(id, fmt, varargin)
    persistent seen
    if isempty(seen); seen = containers.Map('KeyType','char','ValueType','logical'); end
    if isKey(seen, id); return; end
    seen(id) = true;
    warning(id, fmt, varargin{:});
end

% --- z-average a (ny x nz) PIV field to a 1-D profile, ignoring NaN/0 ---
function wavg = local_rowmean(W)
    W = double(W);  W(W == 0) = NaN;
    wavg = mean(W, 2, 'omitnan');
end

% --- common axis styling (serif/LaTeX, grid, box) ---
function local_style(ax)
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 11, ...
            'GridAlpha', 0.3, 'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
end
