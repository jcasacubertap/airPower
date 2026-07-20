function fig = plotProfilesValidation(sBF, sPert, inp, savedir)
% plotProfilesValidation  Dimensional w-profile comparison vs PIV (DeHNSSo).
%
%   fig = plotProfilesValidation(sBF, sPert, inp, savedir)
%
% Styled to match the PreProcessing DFP figure
% (wFieldValidationExperimentaldirectFlatPlateModule): DeHNSSo = black line with
% dot markers, Experimental (PIV) = red filled circles, w [m/s] vs wall distance
% [m], LaTeX serif. A 2-row grid, one column per station: top row base-flow w,
% bottom row the fundamental perturbation w RMS (|w~|/sqrt(2), the z-RMS of a
% single mode -- the same quantity the PIV reports).
%
% Stations are a subset (<=6) of the PIV Gen/Case xc list that map INTO the valid
% DeHNSSo domain (0 < x_DFP <= buffer start); stations landing in the outflow
% buffer are dropped. Mapping is DFP-only for now (reference airfoil data):
%   x/c -> S (arc-length, airfoilFlowData BL.S/x/c) -> x_DFP = S - xInlet
%   -> nearest DeHNSSo column (x_phys = x/delta_0 * l_ref).

    fig = [];
    if nargin < 4; savedir = ''; end

    if ~strcmpi(inp.caseType, 'DFP')
        warning('plotProfilesValidation:notWired', ...
                'w-vs-PIV station mapping is only wired for caseType=DFP (got ''%s'').', inp.caseType);
        return;
    end
    if ~isfield(sPert,'uref') || ~isfield(sPert,'lref') || isempty(sPert.uref) || isempty(sPert.lref)
        warning('plotProfilesValidation:noRef', ...
                'StabGrid Uref/lref missing; cannot dimensionalize DeHNSSo data.'); return;
    end
    uref = sPert.uref;  lref = sPert.lref;

    % ---- PIV Gen/Case data ----
    valdir = fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'Validation', ...
                      sprintf('Gen%d', inp.valGen), 'Experimental', sprintf('Case%d', inp.valCase));
    mats = dir(fullfile(valdir, '*.mat'));
    if isempty(mats); warning('plotProfilesValidation:noPIV','No PIV .mat in %s', valdir); return; end
    o = load(fullfile(valdir, mats(1).name)).output;

    % ---- DFP station mapping x/c -> x_DFP (reference airfoil data) ----
    fd = dir(fullfile(inp.airPowerRoot, 'PreProcessing', 'io', 'airfoilFlowData', '*.mat'));
    if isempty(fd); warning('plotProfilesValidation:noBL','No airfoilFlowData .mat for x/c->S.'); return; end
    BL = load(fullfile(fd(1).folder, fd(1).name)).BL;
    x_ref = double(BL.x(:));  S_ref = double(BL.S(:));  c_ref = double(BL.c);
    xc2xdfp = @(xic) interp1(x_ref, S_ref, xic*c_ref, 'linear', 'extrap') - inp.xInlet;

    xg = sBF.x(1,:) * lref;           % DeHNSSo columns in physical x_DFP [m]
    Nx = numel(xg);
    bufFrac = 0.85;
    if isfield(inp,'ro') && isfield(inp.ro,'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end
    xbuf = xg(min(Nx, max(2, round(bufFrac*Nx))));   % buffer start in physical x [m]

    % stations: valid (0 < x_DFP <= buffer) -> subset of <=6
    xc_all = cellfun(@(v) double(v(1)), o.xc(:));
    [xc_all, ord] = sort(xc_all);
    xdfp   = arrayfun(@(p) xc2xdfp(p/100), xc_all);
    valid  = xdfp > 0 & xdfp <= xbuf;
    if ~any(valid)
        warning('plotProfilesValidation:noValidStation', ...
                'No PIV station maps inside the valid DeHNSSo domain (0, %.3f m].', xbuf); return;
    end
    xc_v = xc_all(valid);  ord_v = ord(valid);
    pick = unique(round(linspace(1, numel(xc_v), min(6, numel(xc_v)))));
    xc = xc_v(pick);  korig = ord_v(pick);  ns = numel(xc);

    % fundamental mode (or the first in modeIdx)
    Nmode = size(sPert.w, 1);
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx)
        m = inp.modeIdx(1);
    else
        m = min(2, Nmode);
    end
    m = max(2, min(m, Nmode));                 % (0,0) has no PIV counterpart

    y_m  = sBF.y(:,1) * lref;                   % model wall distance [m]
    yTop = 0.004;                              % fixed range, matching the reference
    cSolver = [0 0 0];  cExp = [0.75 0.08 0.08];

    fig = figure('Color', 'w', 'Position', [40 60 330*ns 620]);
    for j = 1:ns
        xic = xc(j)/100;  k = korig(j);
        [~, c] = min(abs(xg - xc2xdfp(xic)));

        % ---- top: base-flow w ----
        ax = subplot(2, ns, j); hold(ax, 'on');
        h1 = plot(ax, sBF.w(:,c)*uref, y_m, '-', 'Color', cSolver, ...
                  'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 7);
        wavg = local_rowmean(o.w_m_mean{k});  ye = double(o.y{k});  ye = ye(:,1)*1e-3;
        h2 = plot(ax, wavg, ye, 'o', 'MarkerSize', 3.5, ...
                  'MarkerFaceColor', cExp, 'MarkerEdgeColor', cExp);
        local_style(ax); ylim(ax, [0 yTop]);
        xlabel(ax, '$w \ \mathrm{[m/s]}$', 'Interpreter', 'latex');
        if j == 1; ylabel(ax, '$\mathrm{Wall\ distance\ [m]}$', 'Interpreter', 'latex'); end
        title(ax, sprintf('$x/c = %.0f\\%%$', xc(j)), 'Interpreter', 'latex');
        if j == 1
            legend(ax, [h1 h2], {'\textrm{DeHNSSo}','\textrm{Experimental}'}, ...
                   'Interpreter', 'latex', 'Location', 'northwest');
        end

        % ---- bottom: fundamental perturbation w RMS ----
        ax = subplot(2, ns, ns + j); hold(ax, 'on');
        wrms = abs(squeeze(sPert.w(m,:,c))) / sqrt(2) * uref;
        g1 = plot(ax, wrms, y_m, '-', 'Color', cSolver, ...
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
