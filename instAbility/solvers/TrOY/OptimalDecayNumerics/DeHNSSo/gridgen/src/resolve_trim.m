function n = resolve_trim(params, fname, n_total)
%RESOLVE_TRIM  Convert params.<fname> to an integer station count.
%   []          -> 0
%   integer ≥ 1 -> that many stations
%   0 < x < 1   -> round(x * n_total)
n = 0;
if ~isfield(params, fname) || isempty(params.(fname));  return;  end
v = params.(fname);
if ~isscalar(v) || ~isnumeric(v) || v < 0
    error('%s must be a non-negative scalar.', fname);
end
if v < 1
    n = round(v * n_total);
else
    n = round(v);
end
end
