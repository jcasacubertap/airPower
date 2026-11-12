function figs = plotProfiles(sBF, sPert, inp, savedir)
% plotProfiles  Wall-normal profiles at streamwise stations (DeHNSSo, loadFields).
%
%   figs = plotProfiles(sBF, sPert, inp)            % show only
%   figs = plotProfiles(sBF, sPert, inp, savedir)   % also save PNGs
%
% u-profiles (non-dimensional), only when inp.plotUProfiles is true — OFF by
% default for now: a 2-row grid, one column per station — top row the base flow
% u_B, bottom row the perturbation |u~| for the selected modes (inp.modeIdx;
% [] -> all). Stations are the inlet, 1/2, 3/4 and the buffer start (fractions
% of the streamwise extent). The y-axis is a near-wall window.
%
% w-profiles are ALWAYS drawn (see plotProfilesValidation): a broad view of
% w_B + w_(0,0) and w_rms at stations across the domain, plus their peak against
% the streamwise station. inp.validation only decides how:
%   true   stations in x/c, with the PIV Gen/Case overlay, plus a zoom on the
%          PIV stations inside inp.valXcZoom (needs the airfoil geometry);
%   false  stations in S*/delta_0, DeHNSSo alone -- nothing airfoil-related.
%
% Save behaviour: savedir non-empty -> PNGs written there (terminal/dispatcher);
% empty/omitted -> figures left open for interactive inspection (direct MATLAB).

    if nargin < 4; savedir = ''; end
    figs = gobjects(0);

    if showU(inp)
        figs(end+1) = uFigure(sBF, sPert, inp, savedir);
    end

    % w-profiles; against PIV in x/c only when inp.validation is on
    fw = plotProfilesValidation(sBF, sPert, inp, savedir);
    if ~isempty(fw); figs = [figs, fw]; end
end

% whether to draw the non-dimensional u-profile figure (inp.plotUProfiles).
% Default FALSE: profiles_u.png is off for now, the w-vs-PIV comparison is what
% is being worked on. Set inp.plotUProfiles = true in inputs.jl to bring it back.
function t = showU(inp)
    t = false;
    if isfield(inp, 'plotUProfiles') && ~isempty(inp.plotUProfiles)
        t = logical(inp.plotUProfiles);
    end
end

function f = uFigure(sBF, sPert, inp, savedir)

    [X, Y] = plotCoords(sBF, inp);             % wall-fitted rectangle (unwraps curved TTCP wall)
    Nx = size(X, 2);
    [~, yl, bufFrac] = plotWindow(X, Y, inp);  % near-wall y window + buffer cut

    % stations: inlet, 2/4, 3/4, buffer start (fractions of the streamwise extent)
    % four stations evenly spaced from the inlet to the buffer start, so they
    % stay ordered and inside the plotted (buffer-cut) domain for any bufFrac.
    fracs = linspace(0, bufFrac, 4);
    cols  = min(Nx, max(1, round(fracs * Nx)));  cols(1) = 1;
    xs    = X(1, cols);
    ns    = numel(cols);

    % modes to show (inp.modeIdx; [] -> all)
    Nmode = size(sPert.u, 1);
    if isfield(inp,'modeIdx') && ~isempty(inp.modeIdx)
        modes = inp.modeIdx(:).';
    else
        modes = 1:Nmode;
    end
    modes  = modes(modes >= 1 & modes <= Nmode);
    colors = lines(max(numel(modes), 1));

    f = figure('Name', 'Profiles u', 'Color', 'w', 'Position', [50 50 260*ns 560]);
    for j = 1:ns
        c  = cols(j);
        yc = Y(:, c);

        % top row: base-flow u_B
        subplot(2, ns, j); hold on;
        plot(sBF.u(:, c), yc, '-', 'LineWidth', 1.5, 'Color', [0 0 0]);
        set(gca, 'TickLabelInterpreter', 'latex'); box on; grid on;
        xlabel('$u_{\mathrm{B}} / u_\infty$', 'Interpreter', 'latex');
        if j == 1; ylabel('$n / \delta_0$', 'Interpreter', 'latex'); end
        title(sprintf('$S^{*}/\\delta_0 = %.0f$', xs(j)), 'Interpreter', 'latex');
        ylim(yl);

        % bottom row: perturbation |u~| for the selected modes
        subplot(2, ns, ns + j); hold on;
        for i = 1:numel(modes)
            m = modes(i);
            plot(abs(squeeze(sPert.u(m, :, c))), yc, '-', 'LineWidth', 1.3, 'Color', colors(i,:));
        end
        set(gca, 'TickLabelInterpreter', 'latex'); box on; grid on;
        xlabel('$|\tilde{u}| / u_\infty$', 'Interpreter', 'latex');
        if j == 1; ylabel('$n / \delta_0$', 'Interpreter', 'latex'); end
        ylim(yl);
        if j == ns
            legend(arrayfun(@(m) sprintf('$(0,%d)$', m-1), modes, 'UniformOutput', false), ...
                   'Interpreter', 'latex', 'Location', 'best');
        end
    end
    sgtitle('\textrm{Wall-normal profiles of}\ \ $u$', 'Interpreter', 'latex');

    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, 'profiles_u.png');
        exportgraphics(f, out, 'Resolution', 150);
        fprintf('plotProfiles: saved %s\n', out);
    end
end
