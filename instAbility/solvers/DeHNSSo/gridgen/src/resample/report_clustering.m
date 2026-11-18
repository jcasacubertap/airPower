function report_clustering(mug, sig, ag, coord)
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
% Print the streamwise clustering summary (all resamplers).
%
%   report_clustering(mug, sig, ag, coord)
%
%   Reports the requested clustering and, when coord is given, what the mesh
%   actually does: the finest dxi against the equidistant one, and the worst
%   ratio between neighbouring dxi (the stretching, 1.00 = none).

%% Inputs
% Name  size   unit   explanation
% mug   (1,nc) [-]    cluster centres in physical x
% sig   (1,nc) [-]    Gaussian width of each cluster
% ag    (1,nc) [-]    peak spacing factor of each cluster
% coord (1,n)  [-]    resulting stations; when given, the finest and worst spacing ratios are reported

%% Outputs
% None.

if isempty(mug) || all(ag(:) == 1)
    fprintf('xi distribution: equidistant (no refinement)\n');
    return
end

mug = mug(:).';
if isscalar(sig);  sig = sig * ones(size(mug));  end
if isscalar(ag);   ag  = ag  * ones(size(mug));  end

fprintf('xi distribution: Gaussian clustering, %d centre%s\n', ...
        numel(mug), repmat('s', 1, numel(mug) > 1));
for k = 1:numel(mug)
    fprintf('    centre %d: x_c = %-8.3f  sig = %-5.3f  a_g = %-5.3f  (dx ratio = %.2f)\n', ...
            k, mug(k), sig(k), ag(k), 1 / ag(k));
end
% Measured, from the stations actually produced: the finest spacing against the
% equidistant one, and the worst ratio between neighbouring cells.
if nargin >= 4 && ~isempty(coord) && numel(coord) > 2
    d  = diff(coord(:).');
    du = (coord(end) - coord(1)) / (numel(coord) - 1);
    r  = d(2:end) ./ d(1:end-1);
    fprintf('    measured: finest dxi = %.2f x equidistant (%.2fx finer) ; stretching = %.4f\n', ...
            min(d)/du, du/min(d), max(max(r), max(1./r)));
end
end
