function [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes)
%SAVE_SWEEP  Save intermediate results during convergence iterations.
%
%   [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes)
%
%   When Opt.Sweep is enabled, saves the current amplitude and wavenumber
%   data to disk at each convergence iteration for post-hoc analysis.

iu = sqrt(-1);
Opt.k2 = Opt.k2 + 1;
Opt.AF(Opt.k2) = Opt.AF;

% Compute du/dxi on the (possibly non-uniform) ξ grid
D_xi = StabGrid.D_xi;
for j = Mode.RunJ
    dux(j,:,:) = squeeze(StabRes.u(j,:,:)) * D_xi.';
end

yinteg = linspace(0, max(StabGrid.y), 4000);

for j = Mode.RunJ
    for i = 1:StabGrid.nx
        [~, yimax] = max(abs(interp1(StabGrid.y, StabRes.u(j,:,i), yinteg, 'spline')));
        dumax = interp1(StabGrid.y, dux(j,:,i), yinteg(yimax), 'spline');
        umax  = interp1(StabGrid.y, StabRes.u(j,:,i), yinteg(yimax), 'spline');

        StabRes.alphasweep(j, i, Opt.k2) = (dumax / umax) / iu;
        StabRes.usweep(j,:,i, Opt.k2) = StabRes.u(j,:,i) ./ abs(umax);
        StabRes.vsweep(j,:,i, Opt.k2) = StabRes.v(j,:,i) ./ abs(umax);
        StabRes.wsweep(j,:,i, Opt.k2) = StabRes.w(j,:,i) ./ abs(umax);
        StabRes.psweep(j,:,i, Opt.k2) = StabRes.p(j,:,i) ./ abs(umax);

        if j == round(Mode.nf/2)
            StabRes.Asweep(j, i, Opt.k2) = abs(umax);
        else
            StabRes.Asweep(j, i, Opt.k2) = 2 * abs(umax);
        end
    end
end

end
