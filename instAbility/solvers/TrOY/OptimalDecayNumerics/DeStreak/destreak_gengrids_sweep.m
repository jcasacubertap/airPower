function destreak_gengrids_sweep(ids, n_eta, n_xi)
%DESTREAK_GENGRIDS_SWEEP  Grid each swept base flow in bumps/ (clean gridgen session).
%
%   destreak_gengrids_sweep({'flat','p11','p15','p19','p23','d15'}, 60,1200)
%
%   Reads bumps/bf_<id>.csv (DirectFlatPlate output copied by run_bump_sweep.sh),
%   grids each onto a COMMON Malik wall-normal extent (H,y_i) so all grids are
%   index-aligned for differencing, saves gsweep_<id>_<neta>x<nxi>.mat.
%   Must run with ONLY the gridgen path (boundary() shadow — see destreak_gengrids).

if nargin<2||isempty(n_eta); n_eta=60;  end
if nargin<3||isempty(n_xi);  n_xi =1200; end
restoredefaultpath;
here = fileparts(mfilename('fullpath'));  addpath(here);
common = struct('H',84.53,'y_i',4.2,'input_folder',fullfile(here,'bumps'));
for k=1:numel(ids)
    id = ids{k};
    out = fullfile(here, sprintf('gsweep_%s_%dx%d.mat', id, n_eta, n_xi));
    csv = sprintf('bf_%s.csv', id);
    if ~isfile(fullfile(common.input_folder,csv))
        warning('missing %s — skipping', csv); continue;
    end
    destreak_make_grid(csv, n_eta,n_xi,'fd', out, common);
    fprintf('[gengrids_sweep] %s done\n', id);
end
end
