function StabGrid = destreak_make_grid(csvname, n_eta, n_xi, eta_method, outfile, opts)
%DESTREAK_MAKE_GRID  Build a DeHNSSo StabGrid from a swept-wing base-flow CSV.
%
%   StabGrid = destreak_make_grid(csvname, n_eta, n_xi, eta_method, outfile, opts)
%
%   Generalizes destreak_make_flatgrid to any base-flow CSV in
%   DeHNSSo/baseflow/output/benchmark (e.g. 'bf_sweptwing_flat.csv' flat, or
%   'bf_sweptwing_hump.csv' = the DirectFlatPlate bump output). eta_method is
%   'cheb' or 'fd'. `opts` (optional struct) overrides:
%     .H, .y_i        force common Malik wall-normal extent/clustering so grids
%                     from different base flows are index-aligned (REQUIRED when
%                     differencing two base flows).
%     .mug/.sig/.ag   streamwise Gaussian clustering around a bump centre.
%
%   Uses ONLY the DeHNSSo tree inside OptimalDecayNumerics.

if nargin < 2 || isempty(n_eta);      n_eta = 60;    end
if nargin < 3 || isempty(n_xi);       n_xi  = 1200;  end
if nargin < 4 || isempty(eta_method); eta_method = 'cheb'; end
if nargin < 5 || isempty(outfile);    outfile = fullfile(fileparts(mfilename('fullpath')), 'grid.mat'); end
if nargin < 6; opts = struct(); end

here   = fileparts(mfilename('fullpath'));
ggroot = fullfile(here, '..', 'DeHNSSo', 'gridgen');
addpath(ggroot, genpath(fullfile(ggroot, 'src')));

if isfield(opts,'input_folder') && ~isempty(opts.input_folder)
    input.folder = opts.input_folder;
else
    input.folder = fullfile(here, '..', 'DeHNSSo', 'baseflow', 'output', 'benchmark');
end
input.filename   = csvname;
input.format     = 'csv';
input.structured = false;

params.n_eta_new = n_eta;
params.n_xi_new  = n_xi;
params.Uref      = 12.41792766;
params.lref      = 2.3403348e-04;
params.Re        = 199.518703662;
params.y_i       = [];   params.H = [];   params.xi_range = [];
params.xi_trim_inflow = [];  params.xi_trim_outflow = [];
params.rescale   = true;   params.plot = false;
params.FD_xi_order_1 = 4;  params.FD_xi_order_2 = 2;
params.FD_eta_method = eta_method;  params.FD_eta_order = 4;
if isfield(opts,'H')   && ~isempty(opts.H);   params.H   = opts.H;   end
if isfield(opts,'y_i') && ~isempty(opts.y_i); params.y_i = opts.y_i; end
if isfield(opts,'mug'); params.mug = opts.mug; params.sig = opts.sig; params.ag = opts.ag; end

[outdir, outname, outext] = fileparts(outfile);
if isempty(outdir); outdir = here; end
output.folder = outdir;  output.filename = [outname outext];

StabGrid = main_gridgen(input, params, output);
fprintf('[destreak] grid %s: %d eta x %d xi -> %s\n', csvname, StabGrid.ny, StabGrid.nx, outfile);
end
