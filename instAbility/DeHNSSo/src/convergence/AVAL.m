function [Opt,StabRes] = AVAL(Grid,Opt,Mode,StabRes)
%AVAL  Amplitude value check and inflow amplitude damping.
%
%   [Opt, StabRes] = AVAL(Grid, Opt, Mode, StabRes)
%
%   Checks if the maximum perturbation amplitude exceeds Opt.AMAX.
%   If so, reduces the inflow amplitude by half and resets the solution.
%   Tracks the number of damping events in Opt.AVAL.

    ny_   = Grid.ny;
    nx_   = Grid.nx;
    nvar  = size(StabRes.phi, 2) / ny_;
    n_vel = 3;                              % u, v, w slots
    for j = Mode.RunJ
        StabRes.phi(j,:,:) = reshape(StabRes.q(j,:), nvar*ny_, nx_);
    end

    [Opt.AVAL,~] = max(max(abs(squeeze(StabRes.phi(round(Mode.nf/2)+1,1:n_vel*ny_,:)))));

    if Opt.AVAL >= Opt.AMAX
        Opt.AF = Opt.AMAX/Opt.AVAL;
    end

    StabRes.q = Opt.AF*StabRes.q;
end
