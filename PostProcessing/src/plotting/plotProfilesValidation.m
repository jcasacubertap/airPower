function fig = plotProfilesValidation(sBF, sPert, inp, savedir)
% plotProfilesValidation  Dimensional w-profile comparison vs PIV (DeHNSSo).
%
%   fig = plotProfilesValidation(sBF, sPert, inp, savedir)
%
% Styled to match the PreProcessing figure
% (wFieldValidationExperimental<module>): DeHNSSo = black line with dot markers,
% Experimental (PIV) = red filled circles, w [m/s] vs wall distance [m], LaTeX
% serif. A 2-row grid, one column per station: top row base-flow w, bottom row
% the fundamental perturbation w RMS (|w~|/sqrt(2), the z-RMS of a single mode --
% the same quantity the PIV reports).
%
% Stations are a subset (<=6) of the PIV Gen/Case xc list that map INTO the valid
% DeHNSSo domain (up to the outflow-buffer start). The x/c -> column mapping is
% case-specific (see mapStationsDFP / mapStationsTTCP):
%   DFP  (flat reference plate):  x/c -> S (arc-length, airfoilFlowData BL.S/x/c)
%        -> x_DFP = S - xInlet -> nearest DeHNSSo column (x_phys = x/delta0*lref).
%   TTCP (curved airfoil wall):   the stability-grid wall row IS the airfoil
%        surface, so each column's x/c is recovered by the inverse rigid
%        transform (undo AoA rotation + centre offset) of its physical wall
%        point; the station picks the nearest column. Wall-normal distance is the
%        arc-length up the (body-fitted) grid column from the wall.

    fig = [];
    if nargin < 4; savedir = ''; end

    if ~isfield(sPert,'uref') || ~isfield(sPert,'lref') || isempty(sPert.uref) || isempty(sPert.lref)
        warning('plotProfilesValidation:noRef', ...
                'StabGrid Uref/lref missing; cannot dimensionalize DeHNSSo data.'); return;
    end
    uref = sPert.uref;  lref = sPert.lref;

    % ---- PIV Gen/Case data (shared by both cases) ----
    o = loadPIV(inp);
    if isempty(o); return; end

    % ---- station -> column + wall-distance mapping (case-specific) ----
    switch upper(inp.caseType)
        case 'DFP'
            S = mapStationsDFP(sBF, inp, o, lref);
        case 'TTCP'
            S = mapStationsTTCP(sBF, inp, o, lref);
        otherwise
            warning('plotProfilesValidation:notWired', ...
                    'w-vs-PIV station mapping is only wired for DFP/TTCP (got ''%s'').', inp.caseType);
            return;
    end
    if isempty(S)
        warning('plotProfilesValidation:noValidStation', ...
                'No PIV station maps inside the valid DeHNSSo domain.'); return;
    end

    % ---- fundamental mode (or the first in modeIdx) ----
    Nmode = size(sPert.w, 1);
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx)
        m = inp.modeIdx(1);
    else
        m = min(2, Nmode);
    end
    m = max(2, min(m, Nmode));                 % (0,0) has no PIV counterpart

    % ---- plot: one column per station (shared across cases) ----
    ns   = numel(S);
    yTop = 0.004;                              % fixed range, matching the reference
    cSolver = [0 0 0];  cExp = [0.75 0.08 0.08];

    fig = figure('Color', 'w', 'Position', [40 60 330*ns 620]);
    for j = 1:ns
        c = S(j).col;  k = S(j).korig;  yprof = S(j).n;   % wall distance [m], per column

        % ---- top: base-flow w ----
        ax = subplot(2, ns, j); hold(ax, 'on');
        h1 = plot(ax, sBF.w(:,c)*uref, yprof, '-', 'Color', cSolver, ...
                  'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 7);
        wavg = local_rowmean(o.w_m_mean{k});  ye = double(o.y{k});  ye = ye(:,1)*1e-3;
        h2 = plot(ax, wavg, ye, 'o', 'MarkerSize', 3.5, ...
                  'MarkerFaceColor', cExp, 'MarkerEdgeColor', cExp);
        local_style(ax); ylim(ax, [0 yTop]);
        xlabel(ax, '$w_{\mathrm{B}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        if j == 1; ylabel(ax, '$\mathrm{Wall\ distance\ [m]}$', 'Interpreter', 'latex'); end
        title(ax, sprintf('$x/c = %.0f\\%%$', S(j).xcPct), 'Interpreter', 'latex');
        if j == 1
            legend(ax, [h1 h2], {'\textrm{DeHNSSo}','\textrm{Experimental}'}, ...
                   'Interpreter', 'latex', 'Location', 'northwest');
        end

        % ---- bottom: fundamental perturbation w RMS ----
        ax = subplot(2, ns, ns + j); hold(ax, 'on');
        wrms = abs(squeeze(sPert.w(m,:,c))) / sqrt(2) * uref;
        g1 = plot(ax, wrms, yprof, '-', 'Color', cSolver, ...
                  'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 7);
        fn = sprintf('w_pert_m_prof_rms_%02d', m-1);  yn = sprintf('y_prof_rms_%02d', m-1);
        g2 = [];
        if isfield(o, fn) && isfield(o, yn)
            g2 = plot(ax, double(o.(fn){k}), double(o.(yn){k})*1e-3, 'o', 'MarkerSize', 3.5, ...
                      'MarkerFaceColor', cExp, 'MarkerEdgeColor', cExp);
        end
        local_style(ax); ylim(ax, [0 yTop]);
        xlabel(ax, '$w_{\mathrm{rms}} \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        if j == 1; ylabel(ax, '$\mathrm{Wall\ distance\ [m]}$', 'Interpreter', 'latex'); end
        if j == 1 && ~isempty(g2)
            legend(ax, [g1 g2], {'\textrm{DeHNSSo}','\textrm{Experimental}'}, ...
                   'Interpreter', 'latex', 'Location', 'northeast');
        end
    end

    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, 'profiles_w_validation.png');
        exportgraphics(fig, out, 'Resolution', 150);
        fprintf('plotProfilesValidation: saved %s (fundamental mode (0,%d))\n', out, m-1);
    end
end

% ======================================================================
%  PIV loader (shared)
% ======================================================================
function o = loadPIV(inp)
    o = [];
    valdir = fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'Validation', ...
                      sprintf('Gen%d', inp.valGen), 'Experimental', sprintf('Case%d', inp.valCase));
    mats = dir(fullfile(valdir, '*.mat'));
    if isempty(mats)
        warning('plotProfilesValidation:noPIV', 'No PIV .mat in %s', valdir); return;
    end
    o = load(fullfile(valdir, mats(1).name)).output;
end

% ======================================================================
%  DFP station mapping: x/c -> arc-length x_DFP -> nearest column
%  (flat reference plate; wall distance is the grid's wall-normal coord)
% ======================================================================
function S = mapStationsDFP(sBF, inp, o, lref)
    S = struct('xcPct', {}, 'col', {}, 'n', {}, 'korig', {});

    fd = dir(fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'airfoilFlowData', '*.mat'));
    if isempty(fd)
        warning('plotProfilesValidation:noBL', 'No airfoilFlowData .mat for x/c->S.'); return;
    end
    BL = load(fullfile(fd(1).folder, fd(1).name)).BL;
    x_ref = double(BL.x(:));  S_ref = double(BL.S(:));  c_ref = double(BL.c);
    xc2xdfp = @(xic) interp1(x_ref, S_ref, xic*c_ref, 'linear', 'extrap') - inp.xInlet;

    xg = sBF.x(1,:) * lref;           % DeHNSSo columns in physical x_DFP [m]
    Nx = numel(xg);
    bufFrac = bufferFrac(inp);
    xbuf = xg(min(Nx, max(2, round(bufFrac*Nx))));   % buffer start in physical x [m]

    % stations: valid (0 < x_DFP <= buffer) -> subset of <=6
    xc_all = cellfun(@(v) double(v(1)), o.xc(:));
    [xc_all, ord] = sort(xc_all);
    xdfp   = arrayfun(@(p) xc2xdfp(p/100), xc_all);
    valid  = xdfp > 0 & xdfp <= xbuf & xc_all >= xcMin(inp);
    if ~any(valid)
        warning('plotProfilesValidation:noValidStation', ...
                'No PIV station maps inside the valid DeHNSSo domain (0, %.3f m].', xbuf); return;
    end
    [xc, korig] = pickStations(xc_all, ord, valid);

    y_m = sBF.y(:,1) * lref;                   % model wall distance [m] (flat wall)
    for j = 1:numel(xc)
        [~, c] = min(abs(xg - xc2xdfp(xc(j)/100)));
        S(j) = struct('xcPct', xc(j), 'col', c, 'n', y_m, 'korig', korig(j));
    end
end

% ======================================================================
%  TTCP station mapping: airfoil coordinates (inverse rigid transform)
%  The stability-grid wall row IS the airfoil surface, so each column's
%  x/c follows from undoing the AoA rotation + centre offset of its
%  physical wall point. Wall distance is the arc-length up the grid column.
% ======================================================================
function S = mapStationsTTCP(sBF, inp, o, lref)
    S = struct('xcPct', {}, 'col', {}, 'n', {}, 'korig', {});

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

    X = sBF.x * lref;  Y = sBF.y * lref;       % physical grid [m], Ny x Nx (wall at row end)
    Nx = size(X, 2);
    Xw = X(end,:);  Yw = Y(end,:);             % airfoil surface (wall row)

    % x/c per column: forward map is x = ca*xb + sa*yb + xC, xb = (x/c - 0.5)*chord.
    % Undo the rigid part (rotation + offset) to recover xb, hence x/c; the airfoil
    % thickness (yb) drops out because it is already baked into the wall point.
    xic_col = (ca*(Xw - xC) - sa*(Yw - yC)) / chord + 0.5;

    % buffer cut: drop stations beyond the outflow-buffer start (in x/c)
    ib   = min(Nx, max(2, round(bufferFrac(inp)*Nx)));
    xcLo = min(xic_col);  xcHi = xic_col(ib);
    if xcHi < xcLo; tmp = xcLo; xcLo = xcHi; xcHi = tmp; end

    xc_all = cellfun(@(v) double(v(1)), o.xc(:));
    [xc_all, ord] = sort(xc_all);
    frac  = xc_all / 100;
    valid = frac >= xcLo & frac <= xcHi & xc_all >= xcMin(inp);
    if ~any(valid)
        warning('plotProfilesValidation:noValidStation', ...
                'No PIV station maps inside the TTCP domain (x/c in [%.3f, %.3f]).', xcLo, xcHi); return;
    end
    [xc, korig] = pickStations(xc_all, ord, valid);

    for j = 1:numel(xc)
        [~, c] = min(abs(xic_col - xc(j)/100));
        % wall-normal distance up column c: arc-length from the wall (row end)
        seg = hypot(diff(X(:,c)), diff(Y(:,c)));       % (Ny-1) x 1, between-row lengths
        n   = [flipud(cumsum(flipud(seg))); 0];        % Ny x 1, 0 at wall (row end)
        S(j) = struct('xcPct', xc(j), 'col', c, 'n', n, 'korig', korig(j));
    end
end

% ======================================================================
%  shared helpers
% ======================================================================
% streamwise fraction of the domain where the outflow buffer starts
function f = bufferFrac(inp)
    f = 0.85;
    if isfield(inp,'ro') && isfield(inp.ro,'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        f = inp.ro.bufferFrac;
    end
end

% lowest PIV station x/c [%] to include (inp.valXcMin; 0/absent -> all)
function v = xcMin(inp)
    v = 0;
    if isfield(inp, 'valXcMin') && ~isempty(inp.valXcMin)
        v = inp.valXcMin;
    end
end

% pick an ordered subset (<=6) of valid stations, evenly spaced
function [xc, korig] = pickStations(xc_all, ord, valid)
    xc_v = xc_all(valid);  ord_v = ord(valid);
    pick = unique(round(linspace(1, numel(xc_v), min(6, numel(xc_v)))));
    xc   = xc_v(pick);  korig = ord_v(pick);
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
