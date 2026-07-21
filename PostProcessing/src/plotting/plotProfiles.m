function figs = plotProfiles(sBF, sPert, inp, savedir)
% plotProfiles  Wall-normal profiles at streamwise stations (DeHNSSo, loadFields).
%
%   figs = plotProfiles(sBF, sPert, inp)            % show only
%   figs = plotProfiles(sBF, sPert, inp, savedir)   % also save PNGs
%
% u-profiles (non-dimensional), ALWAYS: a 2-row grid, one column per station —
% top row the base flow u_B, bottom row the perturbation |u~| for the selected
% modes (inp.modeIdx; [] -> all). Stations are the inlet, 1/2, 3/4 and the buffer
% start (fractions of the streamwise extent). The y-axis is a near-wall window.
%
% When inp.validation is true a second, DIMENSIONAL w-profile figure comparing
% against the PIV Gen/Case data is added (see plotProfilesValidation).
%
% Save behaviour: savedir non-empty -> PNGs written there (terminal/dispatcher);
% empty/omitted -> figures left open for interactive inspection (direct MATLAB).

    if nargin < 4; savedir = ''; end
    figs = gobjects(0);

    X = sBF.x;  Y = sBF.y;
    Nx = size(X, 2);
    [~, yl] = plotWindow(X, Y, inp);           % near-wall y window

    % stations: inlet, 2/4, 3/4, buffer start (fractions of the streamwise extent)
    bufFrac = 0.85;
    if isfield(inp,'ro') && isfield(inp.ro,'bufferFrac') && ~isempty(inp.ro.bufferFrac)
        bufFrac = inp.ro.bufferFrac;
    end
    fracs = [0, 0.5, 0.75, bufFrac];
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
        if j == 1; ylabel('$y / \delta_0$', 'Interpreter', 'latex'); end
        title(sprintf('$x/\\delta_0 = %.0f$', xs(j)), 'Interpreter', 'latex');
        ylim(yl);

        % bottom row: perturbation |u~| for the selected modes
        subplot(2, ns, ns + j); hold on;
        for i = 1:numel(modes)
            m = modes(i);
            plot(abs(squeeze(sPert.u(m, :, c))), yc, '-', 'LineWidth', 1.3, 'Color', colors(i,:));
        end
        set(gca, 'TickLabelInterpreter', 'latex'); box on; grid on;
        xlabel('$|\tilde{u}| / u_\infty$', 'Interpreter', 'latex');
        if j == 1; ylabel('$y / \delta_0$', 'Interpreter', 'latex'); end
        ylim(yl);
        if j == ns
            legend(arrayfun(@(m) sprintf('$(0,%d)$', m-1), modes, 'UniformOutput', false), ...
                   'Interpreter', 'latex', 'Location', 'best');
        end
    end
    sgtitle('\textrm{Wall-normal profiles of}\ \ $u$', 'Interpreter', 'latex');
    figs(end+1) = f;

    if ~isempty(savedir)
        if ~exist(savedir, 'dir'); mkdir(savedir); end
        out = fullfile(savedir, 'profiles_u.png');
        exportgraphics(f, out, 'Resolution', 150);
        fprintf('plotProfiles: saved %s\n', out);
    end

    % dimensional w-profile comparison vs PIV (only when validation is on)
    if isfield(inp, 'validation') && inp.validation
        fw = plotProfilesValidation(sBF, sPert, inp, savedir);
        if ~isempty(fw); figs(end+1) = fw; end
    end
end
