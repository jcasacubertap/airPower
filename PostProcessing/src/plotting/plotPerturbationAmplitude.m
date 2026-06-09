function fig = plotPerturbationAmplitude(sBF, sPert, inp, sVal)
% plotPerturbationAmplitude  Chordwise perturbation amplitude of u'.
%
%   fig = plotPerturbationAmplitude(sBF, sPert, inp)
%   fig = plotPerturbationAmplitude(sBF, sPert, inp, sVal)
%
% Plots the perturbation amplitude A_u(x) = max over y of |u'(x,y)|, one
% curve per spanwise mode, versus the chordwise coordinate. 
%
% The amplitude is the stored modal amplitude sPert.A: u/v/w are peak-
% normalized shape functions, so max_y|u'| = A. (sPert.A is preferred over
% re-deriving max_y|sPert.u| because the shape functions carry NaNs for
% several modes.) If sPert.A is absent, falls back to max_y|sPert.u|.
%
% Inputs:
%   sBF    base flow (uses sBF.x; columns index x)
%   sPert  perturbation: sPert.A (Nmode x Nx) amplitude, sPert.u
%          (Nmode x Ny x Nx) shape; optional sPert.beta for the legend
%   inp    config
%   sVal   (optional) validation reference from importValidation; when given,
%          sVal.x and sVal.AMaxU are overlaid as dashed per-mode curves

    % streamwise axis in loaded grid units (x varies along columns)
    x = sBF.x(1, :);                      % 1 x Nx

    % amplitude per mode: max over y (wall-normal) of |u'| == modal amplitude A.
    % Prefer sPert.A (the stored amplitude); fall back to max_y|sPert.u| if absent.
    if isfield(sPert, 'A')
        Au = sPert.A;                                  % Nmode x Nx
    else
        [Nmode, ~, Nx] = size(sPert.u);
        Au = zeros(Nmode, Nx);
        for k = 1:Nmode
            Au(k, :) = max(abs(squeeze(sPert.u(k, :, :))), [], 1, 'omitnan');
        end
    end
    Nmode = size(Au, 1);

    fig = figure('Name', 'Perturbation amplitude |u''|', 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');
    set(ax, 'YScale', 'log');     
    colors = lines(Nmode);
    leg = strings(1, Nmode);
    for k = 1:Nmode
        plot(ax, x, Au(k, :), '-', 'LineWidth', 1.5, 'Color', colors(k, :));
        % (frequency harmonic, spanwise harmonic): steady modes -> (0, k-1)
        leg(k) = sprintf('$(0,%d)$', k - 1);
    end

    % --- validation overlay: reference curves from sVal -----------------
    % Reference amplitude max_y|u'| per mode, drawn dashed in the same
    % per-mode colours (kept out of the legend). The scaling between the
    % reference frame and the loaded grid units is set here:
    if nargin >= 4 && ~isempty(sVal)
        xScaleVal = 7.71E-04;     % <-- chordwise scaling (reference x -> loaded x/l)
        yScaleVal = 15.10;        % <-- amplitude scaling (reference |u'| -> /u_inf)
        nRef = size(sVal.AMaxU, 1);
        for k = 1:min(nRef, Nmode)
            plot(ax, sVal.x * xScaleVal / sPert.lref, sVal.AMaxU(k, :) * yScaleVal / sPert.uref, ...
                '--', 'LineWidth', 1.5, 'Color', colors(k, :), 'HandleVisibility', 'off');
        end
    end
    % --------------------------------------------------------------------

    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [200 1400]);
    ylim(ax, [1e-12 1]);
    xlabel(ax, '$x / \ell$','interpreter','latex');
    ylabel(ax, '$\mathrm{max_y} |u^{\prime}| / u_{\infty}$','interpreter','latex');
    title(ax, 'Chordwise perturbation amplitude of $u^{\prime}$','interpreter','latex');
    legend(ax, leg, 'Location', 'best', 'interpreter', 'latex');
end
