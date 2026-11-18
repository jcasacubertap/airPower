function [applyP, Vopt] = destreak_Popt(G, beta0, stations, uhat)
%DESTREAK_POPT  Sealed cross-plane Poisson operator P (the de-streaking ceiling).
%
%   [applyP, Vopt] = destreak_Popt(G, beta0, stations, uhat)
%
%   Per streamwise station i, solves the local sealed Poisson problem for the
%   cross-plane streamfunction psi,
%        d^2 psi/dy^2 - beta0^2 psi = i*beta0 * (dU/dy) * uhat ,   psi=0 at both walls,
%   and returns the optimal roller  v_opt = i*beta0*psi. This is the operator P
%   of the handoff brief: the local optimal (de-streaking) roller for a given
%   streak, and its transfer norm is the ceiling.
%
%   Discretisation is consistent with the solver's own wall-normal operator:
%   physical d/dy is taken as diag(etay(:,i)) * G.D1 (local cross-plane, the
%   small xiy*d/dxi coupling is dropped — exact for the flat grid, O(slope) on
%   perturbed grids, which is consistent with a local ceiling target).
%
%   Returns applyP: a function handle mapping a full [ny x nx] streak field to
%   the optimal roller [ny x numel(stations)] over the window; and Vopt = applyP(uhat)
%   if uhat is supplied.

ny  = G.ny;
ns  = numel(stations);
I   = eye(ny);
decs = cell(1,ns);
dyU  = cell(1,ns);
for idx = 1:ns
    i  = stations(idx);
    Dy = diag(G.etay(:,i)) * G.D1;         % physical wall-normal derivative
    L  = Dy*Dy - beta0^2 * I;
    L(1,:)  = 0;  L(1,1)   = 1;             % psi = 0 at freestream (row 1)
    L(ny,:) = 0;  L(ny,ny) = 1;            % psi = 0 at wall (row ny)
    decs{idx} = decomposition(L, 'lu');
    dyU{idx}  = G.dyU(:,i);
end

applyP = @(U) apply_impl(U, decs, dyU, beta0, stations, ny);

if nargin >= 4 && ~isempty(uhat)
    Vopt = applyP(uhat);
else
    Vopt = [];
end
end

function Vout = apply_impl(U, decs, dyU, beta0, stations, ny)
ns = numel(stations);
Vout = zeros(ny, ns);
for idx = 1:ns
    i   = stations(idx);
    rhs = 1i*beta0 * (dyU{idx} .* U(:,i));
    rhs(1) = 0;  rhs(ny) = 0;
    psi = decs{idx} \ rhs;
    Vout(:,idx) = 1i*beta0 * psi;
end
end
