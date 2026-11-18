function G2 = destreak_perturb_grid(G, DY)
%DESTREAK_PERTURB_GRID  Rebuild a StabGrid with the wall bent by DY.
%
%   G2 = destreak_perturb_grid(G, DY)
%
%   G   : flat StabGrid (from destreak_make_flatgrid / main_gridgen)
%   DY  : [ny x nx] wall-normal mesh displacement = h(x) .* phi(eta),
%         phi=1 at the wall, ->0 at the freestream (see destreak_relief_basis).
%
%   Installs the passive relief h(x) by displacing the body-fitted mesh Y and
%   letting build_stab_grid recompute the curvilinear metric tensors AND the
%   chain-ruled base-flow derivatives (dyU, ...). The base-flow VALUES U,V,W are
%   frozen in computational (xi,eta) coordinates — this is the frozen-base-flow
%   inner one-shot: no-slip stays exact (U at the wall node is untouched), and
%   the disturbance operator's sensitivity to geometry enters purely through the
%   metrics. How U itself responds to the wall is the deferred OpenFOAM outer loop.
%
%   DY = 0 reproduces G exactly.

here   = fileparts(mfilename('fullpath'));
ggsrc  = fullfile(here, '..', 'DeHNSSo', 'gridgen', 'src');
addpath(genpath(ggsrc));

BF.mesh_type = 'curvilinear';
BF.nx = G.nx;  BF.ny = G.ny;  BF.H = G.H;  BF.y_i = G.y_i;
% build_stab_grid does Grid.x = BF.xi1D + Re; undo the offset so Grid.x == G.x.
BF.xi1D = G.x(:).' - G.Re;
BF.eta1D = G.y(:);

% Flat swept-wing input is classified cartesian (Y constant per column), so a
% flat G carries only 1-D axes. Synthesize the 2-D body-fitted mesh: X=repmat(x),
% Y=repmat(y) gives identity metrics and reproduces the flat operator exactly,
% while giving us a Y surface to bend. A curvilinear G (already has 2-D X,Y) is
% used as-is.
if isfield(G,'X') && ~isempty(G.X) && ismatrix(G.X) && size(G.X,1)==G.ny && size(G.X,2)==G.nx
    BF.X = G.X;   Ybase = G.Y;
else
    BF.X  = repmat(G.x(:).', G.ny, 1);
    Ybase = repmat(G.y(:),   1, G.nx);
end
BF.Y = Ybase + DY;
BF.U = G.U;  BF.V = G.V;  BF.W = G.W;
BF.Re = G.Re;
for f = {'lref','Uref','nu'}
    if isfield(G, f{1}); BF.(f{1}) = G.(f{1}); end
end

gopt = struct('Re', G.Re, ...
              'FD_xi_order_1', 4, 'FD_xi_order_2', 2, ...
              'FD_eta_method', 'cheb', 'FD_eta_order', 4, 'plot', false);

% Silence the per-call gridgen banner/prints.
G2 = evalc_build(BF, gopt);
end

function G2 = evalc_build(BF, gopt)
    [~, G2] = evalc('build_stab_grid(BF, gopt)');
end
