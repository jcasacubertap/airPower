function [J, diag] = wPrimeObjective(npse, piv, in)
% WPRIMEOBJECTIVE  Peak w'-RMS misfit between HNS and PIV.
%
%   [J, diag] = wPrimeObjective(npse, piv, in)
%
% Scalar diagnostic for the amplitude search. At each scored streamwise
% station, take the PEAK (global max over the profile) of the w'-RMS profile
% for both HNS and PIV, and sum their mismatch over the stations. Within the
% matched window (x/c above the near-wall-artifact region) that main peak is
% the CFI-peak amplitude.
%
% Note: this J is only a diagnostic. In the workflow the REPORTED amplitude
% comes from the peak-ratio -> 1 interpolation in findAmplitude, not from
% minimizing J.
%
% J is the sum over stations of the relative peak deviation,
%   | peak_HNS - peak_PIV | / peak_PIV  =  sum |ratio - 1|.
%
% `diag` carries per-station detail for post-hoc inspection:
%   diag.stations  scored PIV station indices
%   diag.x         [m] their streamwise locations
%   diag.res       per-station residual contribution

stations = local_scored_stations(piv, in);

J = 0;
diag.stations = stations;
diag.x   = piv.x(stations);
diag.res = zeros(1, numel(stations));

for s = 1:numel(stations)
    i  = stations(s);
    [~, ix] = min(abs(npse.x - piv.x(i)));      % nearest HNS station

    yq    = piv.y{i};
    rms_n = interp1(npse.y, npse.rmsFull(:, ix), yq, 'linear', 'extrap');
    rms_p = local_target(piv, i);
    pn = max(rms_n);  pp = max(rms_p);          % main (CFI) peak over the profile

    res = abs(pn - pp) / max(pp, eps);          % relative peak deviation |ratio - 1|
    diag.res(s) = res;  J = J + res;
end
end

% -------------------------------------------------------------------------
function g = local_target(piv, i)
% PIV target profile: raw full-plane RMS if available, else the mode-sum.
if isfield(piv, 'rmsFull') && ~isempty(piv.rmsFull{i})
    g = piv.rmsFull{i};
else
    g = piv.rms{i};
end
end

% -------------------------------------------------------------------------
function stations = local_scored_stations(piv, in)
% Resolve which PIV stations to score from in.match.{xStations,stationIdx}.
ns = numel(piv.x);
if ~isempty(in.match.xStations)
    stations = zeros(1, numel(in.match.xStations));
    for q = 1:numel(in.match.xStations)
        [~, stations(q)] = min(abs(piv.x - in.match.xStations(q)));
    end
elseif ~isempty(in.match.stationIdx)
    stations = in.match.stationIdx(:).';
else
    stations = 1:ns;
end
end
