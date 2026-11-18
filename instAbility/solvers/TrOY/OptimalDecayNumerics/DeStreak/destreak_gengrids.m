function destreak_gengrids(n_eta, n_xi)
%DESTREAK_GENGRIDS  Generate flat + bump StabGrids in a CLEAN gridgen-only session.
%
%   Must run with ONLY the gridgen path (restoredefaultpath below): the solver's
%   src/init/boundary.m shadows MATLAB's built-in boundary() used by the
%   unstructured resampler, so grid generation and solving must not share a path.

if nargin<1||isempty(n_eta); n_eta=60;  end
if nargin<2||isempty(n_xi);  n_xi =1200; end
restoredefaultpath;
here = fileparts(mfilename('fullpath'));  addpath(here);
gf = fullfile(here,sprintf('gflat_%dx%d.mat',n_eta,n_xi));
gh = fullfile(here,sprintf('ghump_%dx%d.mat',n_eta,n_xi));
% Force COMMON Malik wall-normal extent so the two grids are index-aligned
% (identical eta nodes; arclength xi aligns to O(bump) which is localized
% upstream of the design window). Without this the resamplers auto-pick
% different H (84.5 vs 171) and window differencing is meaningless.
common = struct('H', 84.53, 'y_i', 4.2);
destreak_make_grid('bf_sweptwing_flat.csv', n_eta,n_xi,'fd',gf, common);
destreak_make_grid('bf_sweptwing_hump.csv', n_eta,n_xi,'fd',gh, common);
fprintf('[gengrids] done: %s , %s\n', gf, gh);
end
