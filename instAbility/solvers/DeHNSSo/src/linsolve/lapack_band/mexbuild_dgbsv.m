%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Compile the banded LU (LU factorization) MEX (MATLAB executable, compiled
%   C) files, linking against MATLAB's own LAPACK.
%   - dgbsv_solve        : factor and solve in one call (*gbsv)
%   - dgbsv_factor       : factor only (*gbtrf) -> a cacheable band matrix + pivots
%   - dgbsv_solve_cached : solve reusing a cached factor, no re-factoring (*gbtrs)
%   Run this once, from any working directory.

here = fileparts(mfilename('fullpath'));
srcs = {'dgbsv_solve.c', 'dgbsv_factor.c', 'dgbsv_solve_cached.c'};
for i = 1:numel(srcs)
    src = fullfile(here, srcs{i});
    fprintf('Compiling %s ...\n', src);
    mex('-R2018a', src, '-lmwlapack', '-outdir', here);
end
fprintf('Done. Built into %s\n', here);
