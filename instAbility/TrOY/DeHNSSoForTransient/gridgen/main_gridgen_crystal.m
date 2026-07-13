function StabGrid = main_gridgen_crystal(input, params, output)
%MAIN_GRIDGEN  Build a DeHNSSo-ready StabGrid from a raw base-flow file.
%
%   StabGrid = main_gridgen(input, params, output)
%
%   Three input paths are supported, auto-selected by classify_input:
%     'cartesian'              — 1-D X, 1-D Y axes (e.g. bf_blasius.mat)
%     'curvilinear_structured' — 2-D X, 2-D Y body-fitted grid
%                                (e.g. bf_sweptwing_flat.mat, bf_sweptwing_hump.mat)
%     'unstructured'           — flat scatter cloud (e.g. m3j.csv)
%
%   input   struct with fields:
%     .folder     directory containing the BF file
%     .filename   BF filename
%     .format     'mat' | 'bin' | 'csv'
%     .structured (optional) — hint passed to classify_input
%
%   params  struct with the union of fields consumed by the chosen resampler
%           (n_eta_new / n_xi_new / y_i / H / xi_range / xi_trim_* / Re /
%            lref / Uref / nu / rescale / elliptic / plot / plot_opts ...)
%           plus grid-generation knobs (FD1_order, FD2_order, eta_method,
%           FD_eta_order). See the per-case callers in gridgen/benchmark/.
%
%   output  struct with fields:
%     .folder     directory where the output .mat is written
%     .filename   output filename (e.g. 'DeHNSSo_input_hump.mat')
%
%   Returns the StabGrid struct and writes it to disk.

%% 0. Banner
fprintf('\n');
fprintf(' ______     _   _  _   _  _____ _____       \n');
fprintf(' |  _  \\   | | | || \\ | |/  ___/  ___|      \n');
fprintf(' | | | |___| |_| ||  \\| |\\ `--.\\ `--.  ___  \n');
fprintf(' | | | / _ \\  _  || . ` | `--. \\`--. \\/ _ \\ \n');
fprintf(' | |/ /  __/ | | || |\\  |/\\__/ /\\__/ / (_) |\n');
fprintf(' |___/ \\___\\_| |_/\\_| \\_/\\____/\\____/ \\___/ \n');
fprintf('                                       v1.2  \n');
fprintf('  Delft Harmonic Navier-Stokes Solver\n');
fprintf('  =========================================\n');

%% 1. Load raw BF
BF = load_bf(input);

%% 2. Classify and dispatch to the appropriate resampler
class = classify_input(BF, input);

if isfield(params,'plot') && params.plot
    plot_raw_bf(BF, class);
end

switch class
    case 'cartesian'
        if isfield(params,'n_eta_new') && ~isempty(params.n_eta_new)
            BF_rs = resample_cartesian(BF, params);
        else
            % Pass-through: BF is already at the DeHNSSo resolution.
            BF_rs = BF;
            BF_rs.BC_mask = false;
        end

    case 'curvilinear_structured'
        BF_rs = resample_curvilinear(BF, params);

    case 'unstructured'
        BF_rs = resample_unstructured(BF, params);

    otherwise
        error('interpol: unknown input class ''%s''.', class);
end

%% 3. Guarantee Re/lref/Uref/nu are on BF_rs before build_stab_grid sees it
BF_rs = ensure_metadata(BF, BF_rs, params);

%% 4. Optional plot of the resampled BF
if isfield(params, 'plot') && params.plot && isfield(BF_rs, 'mesh_type') && ...
        strcmp(BF_rs.mesh_type, 'curvilinear')
    plot_opts = struct();
    if isfield(params, 'plot_opts');  plot_opts = params.plot_opts;  end
    plot_resampled_bf(BF_rs, plot_opts);
end

%% 5. Build the DeHNSSo StabGrid
fprintf('Generating StabGrid ...\n');
StabGrid = build_stab_grid(BF_rs, params_to_grid(params));
fprintf('StabGrid ready.\n');

%% 6. Write to disk
if ~exist(output.folder, 'dir');  mkdir(output.folder);  end
filepath = fullfile(output.folder, output.filename);
save(filepath, 'StabGrid');
fprintf('Saved to %s\n', filepath);

end


function Grid = params_to_grid(params)
% Extract the grid-generation subset of params into a struct that
% build_stab_grid understands.
Grid = struct();
passthrough = {'Re','plot','FD1_order','FD2_order','eta_method','FD_eta_order'};
for k = 1:numel(passthrough)
    f = passthrough{k};
    if isfield(params, f) && ~isempty(params.(f))
        Grid.(f) = params.(f);
    end
end
% Sensible defaults when the caller doesn't set them
if ~isfield(Grid, 'FD1_order');   Grid.FD1_order   = 4;      end
if ~isfield(Grid, 'FD2_order');   Grid.FD2_order   = 2;      end
if ~isfield(Grid, 'eta_method');  Grid.eta_method  = 'cheb'; end
if ~isfield(Grid, 'FD_eta_order');Grid.FD_eta_order= 4;      end
if ~isfield(Grid, 'plot');        Grid.plot        = false;  end
end
