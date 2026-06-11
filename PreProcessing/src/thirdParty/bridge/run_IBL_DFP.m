function run_IBL_DFP(infile, outfile)
% RUN_IBL_DFP  Batch driver for the spectral integral boundary-layer solver.
%
%   run_IBL_DFP(infile, outfile)
%
% Reads solver inputs from the MAT-file `infile` (written by the Julia
% orchestrator, runIBL.jl), runs solver_IBL_spectral, and writes the velocity
% fields to `outfile` for read-back in Julia. Intended to be invoked headless:
%   matlab -batch "run_IBL_DFP('in.mat','out.mat')"
%
% The spectral solver is used (rather than the Thomas-algorithm solver_IBL)
% because it seeds u, v AND w at the inlet from the Falkner-Skan-Cooke
% similarity solution — so the crossflow is imposed at x = 0 and matches the
% DFP at every station, instead of being grown from a uniform inflow.
%
% Expected fields in `infile`:
%   nx    number of streamwise grid points (= blockMesh streamwise cells)
%   ny    number of wall-normal Chebyshev nodes
%   S     streamwise starting point   [m]  (= xInlet)
%   L     streamwise domain end       [m]  (= xInlet + domainLength)
%   H     wall-normal domain height   [m]  (= domainHeight)
%   y_i   Chebyshev node median       [m]  (half the nodes lie below y_i)
%   nu    kinematic viscosity         [m^2/s]
%   Ue    edge velocity distribution, 1 x nx on the uniform grid S:dx:L
%   We    spanwise freestream velocity (scalar) [m/s]
%
% Output fields in `outfile`:
%   u,v,w   ny x nxout velocity components (row 1 = freestream, row ny = wall,
%           the native ordering of the Chebyshev grid)
%   Y       ny x 1 wall-normal coordinate [m] (Chebyshev nodes, H -> 0)
%   Xgrid   1 x nxout streamwise coordinate in the DFP/OpenFOAM frame
%           (x_OF = S_grid - S, i.e. 0 -> domainLength), truncated to nxout
%           if the marching solver stopped early at separation.

% Make the original solver and its dependencies importable, independent of the
% current working directory. The pristine vendored codes live in ../external
% (this driver is the project glue, in ../bridge).
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', 'external', 'IBL')));
addpath(genpath(fullfile(here, '..', 'external', 'SimilarFlow')));

in = load(infile);

nx  = double(in.nx);
ny  = double(in.ny);
S   = double(in.S);
L   = double(in.L);
H   = double(in.H);
y_i = double(in.y_i);
nu  = double(in.nu);
We  = double(in.We);
Ue  = double(in.Ue(:)).';        % force 1 x nx row vector

if numel(Ue) ~= nx
    error('run_IBL_DFP:UeSize', ...
        'Ue has %d entries but nx = %d.', numel(Ue), nx);
end

% Solve the non-similar incompressible boundary layer with the spectral
% (Chebyshev in y) solver. Inflow profiles (u, v, w) are the Falkner-Skan-Cooke
% similarity solution at the local Hartree parameter. The marching loop may
% stop early at separation, reducing the column count.
[u, v, w, y] = solver_IBL_spectral(nx, ny, S, L, y_i, H, nu, Ue, We);

nxout = size(u, 2);

% Wall-normal grid: the Chebyshev nodes returned by the solver (H -> 0).
Y = y(:);

% Streamwise grid in the DFP/OpenFOAM frame (origin at the inlet).
Sgrid = S:(L-S)/(nx-1):L;
Xgrid = Sgrid(1:nxout) - S;

save(outfile, 'u', 'v', 'w', 'Y', 'Xgrid', '-v7');
end
