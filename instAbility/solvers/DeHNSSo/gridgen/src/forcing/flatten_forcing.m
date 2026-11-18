function specs = flatten_forcing(forcing)
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
% Normalise params.forcing to a cell of single-source specs.
%
%   specs = flatten_forcing(forcing)
%
%   Accepts a struct, a cell of structs, or an actuator struct with a vector xc
%   (and scalar or matching-length w) describing several actuators at once.
%   Returns a flat cell where each entry is one actuator or one distribution.

%% Inputs
% Name    size           unit   explanation
% forcing struct or cell [-]    forcing description: a struct, a cell of structs, or an actuator struct with a vector of centres

%% Outputs
% Name  size unit   explanation
% specs cell [-]    one normalized single-source spec per element

if ~iscell(forcing); forcing = {forcing}; end

specs = {};
for i = 1:numel(forcing)
    s = forcing{i};
    if isfield(s, 'xc') && numel(s.xc) > 1          % vector actuator -> split
        xc = s.xc(:).';
        w  = s.w;  if isscalar(w); w = w * ones(size(xc)); end
        assert(numel(w) == numel(xc), 'flatten_forcing: w must be scalar or match xc.');
        for k = 1:numel(xc)
            sk = s;  sk.xc = xc(k);  sk.w = w(k);
            specs{end+1} = sk;
        end
    else
        specs{end+1} = s;
    end
end
end
