%% Compile dgbsv_solve MEX.
%   Builds src/hns/dgbsv_solve.{c,mexmaci64,mexmaca64} linking against MATLAB's LAPACK.
%   Run once from any working directory.

here = fileparts(mfilename('fullpath'));
src = fullfile(here, 'dgbsv_solve.c');
fprintf('Compiling %s ...\n', src);
mex('-R2018a', src, '-lmwlapack', '-outdir', here);
fprintf('Done. Built into %s\n', here);
