function lc = linear_check(StabGrid, piv, in)
% LINEAR_CHECK  One linear HNS, amplitude-anchored to the PIV at a single x/c.
%
%   lc = linear_check(StabGrid, piv, in)
%
% Runs a SINGLE linear HNS solve, then scales its w'-RMS so the peak at the
% anchor station (in.linearCheck.xcAnchor) matches the PIV peak there. Because a
% linear solution scales uniformly with the inlet amplitude, that one scale
% applies at every station — so comparing the anchored linear curve to the PIV
% up- and down-stream reveals whether the LINEAR growth rate matches the
% experiment (the linear-vs-nonlinear diagnostic). It is immune to the amplitude
% convention: the anchor is a ratio, so any overall factor cancels.
%
% Returns:
%   lc.hns        extracted linear HNS w' bundle (UNSCALED; x, y, rmsFull, ...)
%   lc.c          anchor scale = PIV_peak / HNS_peak at the anchor station
%   lc.A0implied  implied inlet amplitude = c * A(fund,inlet)/2  (the linear seed
%                 that would give the PIV amplitude at the anchor via linear growth)
%   lc.ianchor    PIV station index used as the anchor
%   lc.N          spanwise mode count

    N = in.Stab.N;

    % one linear solve; A0_fund is irrelevant in linear mode (DeHNSSo uses a tiny
    % seed and the result scales linearly), so pass the reference value.
    out = runHNS(StabGrid, in.Stab.A0_fund, N, in, 'linear');
    hns = extractWprimeHNS(out, in);

    % anchor = PIV station nearest xcAnchor
    [~, ia] = min(abs(piv.xc - in.linearCheck.xcAnchor));

    % peak w'-RMS at the anchor: PIV target and (interpolated) linear HNS
    [~, ix] = min(abs(hns.x - piv.x(ia)));
    hn = interp1(hns.y, hns.rmsFull(:, ix), piv.y{ia}, 'linear', 'extrap');
    if isfield(piv,'rmsFull') && ~isempty(piv.rmsFull{ia})
        g = piv.rmsFull{ia};
    else
        g = piv.rms{ia};
    end
    c = max(g) / max(max(hn), eps);

    % implied inlet A0_fund: linear scales, so c times this run's inlet coefficient
    % A(fund,inlet)/2 (same mapping as findAmplitude's linear seed).
    b0   = 2*pi * out.lref / in.Stab.lambda_z;
    fidx = find(abs(out.res.betavec - b0) <= 1e-6*max(1,b0), 1, 'first');
    A0implied = c * out.res.A(fidx, 1) / 2;

    lc = struct('hns', hns, 'c', c, 'A0implied', A0implied, 'ianchor', ia, 'N', N);
end
