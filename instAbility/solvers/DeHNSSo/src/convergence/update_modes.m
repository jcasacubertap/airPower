function [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Unpack solution, compute amplitudes, activate modes, build NLT.
%
%   [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%
%   Reshapes the solution vector q into the primitive fields (recovering the
%   pressure perturbation from the equation of state in the compressible case),
%   computes the modal amplitudes, decides which new harmonics have grown large
%   enough to be advanced, and evaluates the nonlinear forcing term (skipped for
%   linear runs).
%

%% Inputs
% Name     size unit   explanation
% StabGrid      [-]    grid structure (ny, nx, D1, base flow, gamma, Ma)
% Bound         [-]    boundary/buffer structure (ibuff, bufcf)
% Mode          [-]    mode structure (RunJ, Aactive, nf)
% Opt           [-]    options (compressible, linear, allmodes, Sweep, thresholds)
% Stab          [-]    stability data (harmonic-balancing HB, A0, forcing f)
% StabRes       [-]    stability results (updated with u, v, w, p, rho, T, A)

%% Outputs
% Name    size unit   explanation
% StabRes      [-]    stability results (updated with u, v, w, p, rho, T, A)
% Mode         [-]    mode structure (RunJ, Aactive, nf)
% Stab         [-]    stability data (harmonic-balancing HB, A0, forcing f)

ny_      = StabGrid.ny;
nx_      = StabGrid.nx;
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + double(is_compr);

%% Split solution q into u, v, w and either p (incompressible) or rho, T (compressible).
% For compressible runs, p' is not stored directly and is recovered from the equation of state (EOS).
for j = Mode.RunJ
    qj = reshape(StabRes.q(j,:), nvar*ny_, nx_);
    StabRes.u(j,:,:) = qj(1:ny_,        :);
    StabRes.v(j,:,:) = qj(ny_+1:2*ny_,  :);
    StabRes.w(j,:,:) = qj(2*ny_+1:3*ny_,:);
    if is_compr
        StabRes.rho(j,:,:) = qj(3*ny_+1:4*ny_,:);
        StabRes.T(j,:,:)   = qj(4*ny_+1:5*ny_,:);
        % p' = (ρ'·T̄ + ρ̄·T') / (γ·Ma²)
        gM2 = StabGrid.gamma * StabGrid.Ma^2;
        StabRes.p(j,:,:) = (squeeze(StabRes.rho(j,:,:)) .* StabGrid.T ...
                          + StabGrid.rho .* squeeze(StabRes.T(j,:,:))) / gM2;
    else
        StabRes.p(j,:,:) = qj(3*ny_+1:4*ny_,:);
    end
end

%% Compute amplitudes
for j = Mode.RunJ
    StabRes.A(j,:) = max(abs(reshape(StabRes.u(j,:,:), ny_, nx_)), [], 1);
end

%% Determine active modes
if isfield(Opt,'allmodes') && strcmpi(Opt.allmodes,'on')
    % For the 2D cascade (oblique breakdown), activate every (m,n) mode. The
    % usual 1D harmonic ladder would miss off-axis modes such as (2,0), (0,2), (2,2).
    Mode.Aactive(:) = 1;
    Mode.RunJ = 1:Mode.nf;
else
% Outer product of per-mode amplitude maxima
Amax = max(StabRes.A, [], 2);
Amat = Amax * Amax';

[~, index] = max(StabRes.A(:,:), [], 2);
index = max(index);

% Determine which modes are allowed to enter
if Opt.Sweep == 1
    [nzA, ] = find(StabRes.A(:,index));
    nfm = max(nzA) + 1;
    nfl = min(nzA) - 1;
else
    % check_convergence increments k1 AFTER this function runs each iteration,
    % so k1 is still 0 on the first call; clamp so dalsave is never indexed at column 0.
    if Opt.AF >= 1 && max(abs(Opt.dalsave(:,max(Opt.k1,1)))) <= Opt.Conv*Opt.ConvF && ~isnan(Opt.dalsave(max(Mode.RunJ),max(Opt.k1,1)))
        [nzA, ] = find(StabRes.A(:,index));
        nfm = max(nzA) + 1;
        nfl = min(nzA) - 1;
    else
        if any(Stab.A0)
            [nzA, ] = find(Stab.A0(:));
            nfm = max(nzA) + 1;
            nfl = min(nzA) - 1;
        else
            [nzA, ] = [round(Mode.nf/2)-1, round(Mode.nf/2)+1];
            nfm = max(nzA) + 1;
            nfl = min(nzA) - 1;
        end
    end
end

% Predict amplitude of modes to be introduced
Amode = zeros(size(StabRes.A(:,1)));
for j = 1:Mode.nf
    if j > nfm || j < nfl
        Amode(j) = 0;
    else
        Amode(j) = sum(sum(Amat .* Stab.HB(:,:,j)));
    end
end

% Update active modes
Aactive_old = Mode.Aactive;
for j = 1:Mode.nf
    Mode.Aactive(j) = Aactive_old(j) || any((StabRes.A(j,:) + Amode(j,:)) >= Opt.TH);
end
NewModes = find(Mode.Aactive - Aactive_old > 0);
Mode.RunJ = sort([Mode.RunJ NewModes]);
end

%% Calculate nonlinear terms (skip for linear runs)
is_linear = isfield(Opt,'linear') && strcmpi(Opt.linear,'on');
if ~is_linear
    [Stab.f] = NLT_HNS(StabGrid, Mode.RunJ, StabRes, StabGrid.D1, Stab.HB, strcmpi(Opt.parfor,'on'));

    % Apply buffer attenuation (vectorized over active modes)
    ib = Bound.ibuff;
    bufcf_vec = repelem(Bound.bufcf(ib:nx_), nvar*ny_);
    idx = (ib-1)*nvar*ny_+1 : nx_*nvar*ny_;
    Stab.f(Mode.RunJ, idx) = Stab.f(Mode.RunJ, idx) .* bufcf_vec;
end

end
