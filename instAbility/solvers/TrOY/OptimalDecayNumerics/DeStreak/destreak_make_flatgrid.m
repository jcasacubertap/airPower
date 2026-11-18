function StabGrid = destreak_make_flatgrid(n_eta, n_xi, outfile)
%DESTREAK_MAKE_FLATGRID  Build the flat swept-wing StabGrid at chosen resolution.
%
%   StabGrid = destreak_make_flatgrid(n_eta, n_xi, outfile)
%
%   Mirrors gridgen/benchmark/sweptwing_flat.m but exposes the resolution
%   (n_eta wall-normal, n_xi streamwise) and writes to a caller-chosen file so
%   the DeStreak feasibility study never clobbers input/DeHNSSo_input.mat.
%
%   Uses ONLY the DeHNSSo tree inside OptimalDecayNumerics (never the
%   DeHNSSo_Antonis copy in the parent TrOY folder).
%
%   This is the FROZEN base flow for the de-streaking inner one-shot: a
%   flat-plate swept boundary layer. The wall geometry is flat (Y constant per
%   column); the DeStreak perturbation later bends Y to install h(x).

if nargin < 1 || isempty(n_eta);  n_eta = 60;   end
if nargin < 2 || isempty(n_xi);   n_xi  = 1200;  end
if nargin < 3 || isempty(outfile); outfile = fullfile(fileparts(mfilename('fullpath')), 'flat_grid.mat'); end

here    = fileparts(mfilename('fullpath'));
ggroot  = fullfile(here, '..', 'DeHNSSo', 'gridgen');
addpath(ggroot, genpath(fullfile(ggroot, 'src')));

%% Input (flat swept-wing base flow, body-fitted CFD scatter)
input.folder     = fullfile(here, '..', 'DeHNSSo', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_flat.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters (identical to the benchmark except resolution)
params.n_eta_new = n_eta;
params.n_xi_new  = n_xi;
params.Uref      = 12.41792766;
params.lref      = 2.3403348e-04;    % Blasius length at R=400
params.Re        = 199.518703662;
params.y_i       = [];
params.H         = [];
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = true;
params.plot      = false;

%% Grid-generation options — Chebyshev wall-normal (flat, smooth)
params.FD_xi_order_1  = 4;
params.FD_xi_order_2  = 2;
params.FD_eta_method  = 'cheb';
params.FD_eta_order   = 4;

%% Output
[outdir, outname, outext] = fileparts(outfile);
if isempty(outdir); outdir = here; end
output.folder   = outdir;
output.filename = [outname outext];

StabGrid = main_gridgen(input, params, output);
fprintf('[destreak] flat grid: %d eta x %d xi -> %s\n', StabGrid.ny, StabGrid.nx, outfile);
end
