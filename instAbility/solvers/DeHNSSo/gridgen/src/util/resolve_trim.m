function n = resolve_trim(params, fname, n_total)
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
% Convert params.<fname> to an integer station count.
%   []          -> 0
%   integer ≥ 1 -> that many stations
%   0 < x < 1   -> round(x * n_total)

%% Inputs
% Name    size         unit   explanation
% params  (1,1) struct [-]    parameter struct holding the trim specification
% fname   char         [-]    field of params to resolve
% n_total (1,1)        [-]    total number of stations, used to resolve a fractional request

%% Outputs
% Name size  unit   explanation
% n    (1,1) [-]    station count: 0 if empty, the value itself if an integer, round(x*n_total) if a fraction

n = 0;
if ~isfield(params, fname) || isempty(params.(fname));  return;  end
v = params.(fname);
if ~isscalar(v) || ~isnumeric(v) || v < 0
    error('%s must be a non-negative scalar.', fname);
end
if v < 1
    n = round(v * n_total);
else
    n = round(v);
end
end
