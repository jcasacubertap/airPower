function [xl, yl, bufFrac] = plotWindow(X, Y, inp)
% plotWindow  The plotting window, defined once for every PostProcessing figure.
%
%   [xl, yl]          = plotWindow(X, Y, inp)
%   [xl, yl, bufFrac] = plotWindow(X, Y, inp)
%
% Every figure -- base flow, perturbation shapes, profiles, Reynolds-Orr -- takes
% its limits from here, so they always frame the same region. This is the only
% place a window is computed; no plotter defines its own.
%
% x: cut where the DeHNSSo outflow buffer begins, so the damped tail is not
%    shown. inp.plot.bufferFrac is a fraction of the streamwise extent
%    (default 0.85; 1 -> full domain), so it needs no units and serves both modes.
%
% y: from the wall up to a ceiling given in the units of the loaded grid. The two
%    load modes are deliberately NOT reconciled to a common unit, because each is
%    plotted in the system it is computed in -- axes and legends included:
%
%       loadFields -> inp.plot.yMaxFields, in y/delta_0
%                     (the stability analysis is run non-dimensionally)
%       loadBF     -> inp.plot.yMaxBF,     in metres
%                     (the OpenFOAM midPlane is dimensional)
%
%    [] falls back to a near-wall fraction inp.plot.yWallFrac (default 0.30) of
%    the domain height, which is unit-free and so behaves the same either way.
%
% A ceiling above the available data is clamped to the domain, with a warning: an
% out-of-range y-limit otherwise squashes the boundary layer into an invisible
% sliver, which reads as an empty plot rather than as a bad setting.
%
% Rows may run either way; the wall is min(Y).

    % --- options, with defaults for anything the config does not carry ---
    opt = struct('bufferFrac', 0.85, 'yMaxFields', [], 'yMaxBF', [], 'yWallFrac', 0.30);
    if isfield(inp, 'plot')
        f = fieldnames(opt);
        for k = 1:numel(f)
            if isfield(inp.plot, f{k}) && ~isempty(inp.plot.(f{k}))
                opt.(f{k}) = inp.plot.(f{k});
            end
        end
    end
    bufFrac = opt.bufferFrac;

    % --- x: up to the buffer start ---
    Nx = size(X, 2);
    ib = min(Nx, max(2, round(bufFrac * Nx)));
    xl = [min(X(:)), X(1, ib)];

    % --- y: from the wall up, in this load mode's own units ---
    if isfield(inp, 'loadMode') && strcmpi(inp.loadMode, 'loadBF')
        key = 'yMaxBF';      unitLbl = 'm';
    else
        key = 'yMaxFields';  unitLbl = 'delta_0';
    end
    yMax = opt.(key);

    y0 = min(Y(:));  y1 = max(Y(:));
    if isempty(yMax)
        yl = [y0, y0 + opt.yWallFrac * (y1 - y0)];
        return;
    end

    if yMax > (y1 - y0)
        warning('plotWindow:yMaxOutOfRange', ...
                ['inp.plot.%s = %g %s is above the data, which reaches %.4g %s ' ...
                 'from the wall; clamping to the domain.'], ...
                key, yMax, unitLbl, y1 - y0, unitLbl);
        yMax = y1 - y0;
    end
    yl = [y0, y0 + yMax];
end
