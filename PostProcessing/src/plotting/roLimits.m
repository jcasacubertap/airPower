function [cl, e] = roLimits(fields, inwin, inp)
% roLimits  Colour range of the Reynolds-Orr figures, and its power of ten.
%
%   [cl, e] = roLimits(fields, inwin, inp)
%
% inp.ro.climManual true  -> cl = inp.ro.clim, as given ([lower upper], in the
%                            figures' units, P/u_inf^3 for loadFields).
% otherwise (default)     -> symmetric about zero, +-(99th percentile of |F|)
%                            over the plotted window inwin, taken over every
%                            field in the cell array `fields` together, so the
%                            panels of one figure share a scale.
% e is the power of ten factored out of the colorbar tick labels.

    manual = isfield(inp, 'ro') && isfield(inp.ro, 'climManual') && logical(inp.ro.climManual);
    if manual
        cl = double(inp.ro.clim(:)).';
        if numel(cl) ~= 2 || ~(cl(2) > cl(1))
            error('roLimits:clim', ...
                  'settings.jl PostProcessing.ro.clim must be [lower, upper] with lower < upper.');
        end
    else
        a = [];
        for k = 1:numel(fields)
            F = fields{k};
            a = [a; abs(F(inwin & isfinite(F)))]; %#ok<AGROW>
        end
        if isempty(a)
            for k = 1:numel(fields); a = [a; abs(fields{k}(isfinite(fields{k})))]; end %#ok<AGROW>
        end
        a = sort(a);
        if isempty(a) || a(end) == 0
            c = 1;
        else
            c = a(max(1, round(0.99 * numel(a))));
            if ~(c > 0); c = a(end); end
        end
        cl = [-c c];
    end
    m = max(abs(cl));
    e = 0;  if m > 0; e = floor(log10(m)); end
end
