function ibl_driver(infile, outfile)
% IBL_DRIVER  Batch driver for the spectral integral boundary-layer solver.
%
%   ibl_driver(infile, outfile)
%
% Reads solver inputs from the MAT-file `infile` (written by the Julia bridge,
% runIBL.jl), runs solver_IBL_spectral, and writes the velocity fields to
% `outfile`. Invoked headless:  matlab -batch "ibl_driver('in.mat','out.mat')"
%
% The streamwise virtual origin S is set one of two ways:
%   • given directly  (in.S)            — the DFP path (S = xInlet).
%   • by δ*-calibration (in.calibrate_dstar > 0) — the TTCP path: choose S so
%     the Falkner-Skan-Cooke inflow displacement thickness equals the target,
%     using the SAME β convention as the solver. The marching span L−S is fixed
%     by in.Lspan, so the Ue array (sampled over the physical range) is
%     unaffected by S; only the inflow scaling/β depend on it.
%
% Expected fields in `infile`:
%   nx, ny    streamwise points / wall-normal Chebyshev nodes
%   Lspan     streamwise span  [m]   (L = S + Lspan)
%   H         wall-normal height [m]
%   y_i       Chebyshev node median [m]
%   nu        kinematic viscosity [m^2/s]
%   We        spanwise freestream velocity [m/s]
%   Ue        edge velocity, 1 x nx over the uniform grid
%   S         virtual origin [m]            (used when calibrate_dstar <= 0)
%   calibrate_dstar  target δ* [m]          (>0 activates calibration)
%
% Output fields in `outfile`:
%   u,v,w   ny x nxout (row 1 = freestream, row ny = wall)
%   Y       ny x 1 Chebyshev wall-normal coordinate [m] (H -> 0)
%   Xgrid   1 x nxout streamwise coordinate measured from the inlet
%           (0 -> Lspan), truncated if the solver stopped at separation.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', 'external', 'IBL')));
addpath(genpath(fullfile(here, '..', 'external', 'SimilarFlow')));

in = load(infile);

nx    = double(in.nx);
ny    = double(in.ny);
Lspan = double(in.Lspan);
H     = double(in.H);
y_i   = double(in.y_i);
nu    = double(in.nu);
We    = double(in.We);
Ue    = double(in.Ue(:)).';          % 1 x nx row vector

if numel(Ue) ~= nx
    error('ibl_driver:UeSize', 'Ue has %d entries but nx = %d.', numel(Ue), nx);
end

calib = 0;
if isfield(in, 'calibrate_dstar'); calib = double(in.calibrate_dstar); end

if calib > 0
    % --- δ*-calibration of the virtual origin (option B) ---
    % Replicate the solver's β: dUe = -FD1d4o(Ue,dx), m = dUe(1)*S/Ue(1).
    dx   = Lspan / (nx - 1);
    dUe  = -FD1d4o(Ue, dx);
    dUe1 = dUe(1);
    Ue1  = Ue(1);
    fsc_ds = @(S) local_fsc_dstar(S, dUe1, Ue1, nu, We);
    % Blasius-based initial guess (δ* ≈ 1.72·sqrt(ν S / Ue)).
    S0 = (calib / 1.72)^2 * Ue1 / nu;
    S  = fzero(@(S) fsc_ds(S) - calib, S0);
    fprintf('Calibrated virtual origin S = %.6g m (target δ* = %.4g, FSC δ* = %.4g).\n', ...
            S, calib, fsc_ds(S));
else
    S = double(in.S);
end

L = S + Lspan;

[u, v, w, y] = solver_IBL_spectral(nx, ny, S, L, y_i, H, nu, Ue, We);

nxout = size(u, 2);
Y     = y(:);
Sgrid = S:(L - S)/(nx - 1):L;
Xgrid = Sgrid(1:nxout) - S;          % measured from the inlet (0 -> Lspan)

save(outfile, 'u', 'v', 'w', 'Y', 'Xgrid', '-v7');
end

% -------------------------------------------------------------------------
function ds = local_fsc_dstar(S, dUe1, Ue1, nu, We)
% Displacement thickness of the Falkner-Skan-Cooke inflow at virtual origin S,
% using the same β and FSC settings (eta_inf=30, eta_n=2000) as
% solver_IBL_spectral's inflow block. δ* is integrated from the wall to δ99
% (first u/Ue=0.99 crossing) to MATCH the Julia compute_bl_integrals
% definition used for the reported metrics — so the calibrated inlet δ*
% coincides with the comparison value.
m    = dUe1 * S / Ue1;
m    = min(m, 1);
beta = 2*m/(m + 1);
beta = max(beta, -0.19884);
[y, u, ~, ~] = FalknerSkanCookeReal(30, 2000, beta, nu, S, Ue1, We, true);
y = y(:); r = u(:) / Ue1;                    % FSC starts at wall (y=0, u=0)

% δ99: first crossing r >= 0.99 (linearly interpolated)
k = find(r(1:end-1) < 0.99 & r(2:end) >= 0.99, 1, 'first');
if isempty(k)
    ds = trapz(y, 1 - r);                    % never reaches 0.99 (fallback)
    return;
end
t   = (0.99 - r(k)) / (r(k+1) - r(k));
y99 = y(k) + t * (y(k+1) - y(k));

% Trapezoid of (1 - r) from 0 to δ99
ds = 0;
for i = 2:numel(y)
    if y(i) >= y99
        ri = r(i-1) + (y99 - y(i-1)) / (y(i) - y(i-1)) * (r(i) - r(i-1));
        ds = ds + 0.5 * ((1 - r(i-1)) + (1 - ri)) * (y99 - y(i-1));
        break;
    end
    ds = ds + 0.5 * ((1 - r(i-1)) + (1 - r(i))) * (y(i) - y(i-1));
end
end
