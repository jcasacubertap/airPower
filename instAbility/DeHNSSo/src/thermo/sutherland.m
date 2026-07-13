function [mu, dmuT, d2muT2] = sutherland(T, Sstar)
%SUTHERLAND  Sutherland's viscosity law (non-dimensional).
%
%   [mu, dmuT, d2muT2] = sutherland(T, Sstar)
%
%   Non-dimensional form: mu = T^(3/2) (1 + S*) / (T + S*)
%   where S* = S / T_inf  (S = 110.4 K for air).
%
%   Inputs:
%     T      — non-dimensional temperature T/T_inf  (array)
%     Sstar  — S/T_inf  (scalar, e.g. 110.4/288 for air at T_inf = 288 K)
%
%   Outputs:
%     mu     — dynamic viscosity mu(T)
%     dmuT   — first derivative d(mu)/dT
%     d2muT2 — second derivative d^2(mu)/dT^2

mu = T.^(3/2) .* (1 + Sstar) ./ (T + Sstar);

if nargout >= 2
    % dmu/dT = mu/T * (3/2 - T/(T + S*))
    dmuT = mu ./ T .* (3/2 - T ./ (T + Sstar));
end

if nargout >= 3
    % d2mu/dT2 from differentiating dmu/dT
    TpS = T + Sstar;
    d2muT2 = mu ./ T.^2 .* (3/4 + Sstar * (3*Sstar - T) ./ (2 * TpS.^2));
end

end
