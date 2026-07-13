function [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%UPDATE_MODES  Extract solution, compute amplitudes, activate new modes, compute NLT.

ny_ = StabGrid.ny;
nx_ = StabGrid.nx;

%% Extract u,v,w,p from solution q
for j = Mode.RunJ
    qj = reshape(StabRes.q(j,:), 4*ny_, nx_);
    StabRes.u(j,:,:) = qj(1:ny_,        :);
    StabRes.v(j,:,:) = qj(ny_+1:2*ny_,  :);
    StabRes.w(j,:,:) = qj(2*ny_+1:3*ny_,:);
    StabRes.p(j,:,:) = qj(3*ny_+1:4*ny_,:);
end

%% Compute amplitudes
for j = Mode.RunJ
    StabRes.A(j,:) = max(abs(reshape(StabRes.u(j,:,:), ny_, nx_)), [], 1);
end

%% Determine active modes
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
    if Opt.AF >= 1 && max(abs(Opt.dalsave(:,max(Opt.k1,1)))) <= Opt.Conv*Opt.ConvF && ~isnan(Opt.dalsave(max(Mode.RunJ),Opt.k1))
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
    Mode.Aactive(j) = any((StabRes.A(j,:) + Amode(j,:)) >= Opt.TH);
end
NewModes = find(Mode.Aactive - Aactive_old > 0);
Mode.RunJ = sort([Mode.RunJ NewModes]);

%% Calculate nonlinear terms
[Stab.f] = NLT_HNS(StabGrid, Mode.RunJ, StabRes, StabGrid.D1, Stab.HB);


% Apply buffer attenuation (vectorised over active modes)
ib = Bound.ibuff;
bufcf_vec = repelem(Bound.bufcf(ib:nx_), 4*ny_);
idx = (ib-1)*4*ny_+1 : nx_*4*ny_;
Stab.f(Mode.RunJ, idx) = Stab.f(Mode.RunJ, idx) .* bufcf_vec;

end
