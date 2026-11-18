function BF_new = ensure_metadata(BF_raw, BF_new, params)
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
% Guarantee Re/lref/Uref/nu are present on the resampled BF.
%
%   BF_new = ensure_metadata(BF_raw, BF_new, params)
%
%   Fills missing reference scales on BF_new, in order of precedence: already set
%   on BF_new, carried from BF_raw, given in params, or derived from the other
%   three (nu = Uref*lref/Re). Errors if Re cannot be established, since
%   build_stab_grid needs it.

%% Inputs
% Name   size         unit   explanation
% BF_raw (1,1) struct [-]    base flow as read from file, used as a fallback source of reference scales
% BF_new (1,1) struct [-]    base flow resampled onto the stability grid
% params (1,1) struct [-]    grid-generation parameters

%% Outputs
% Name   size         unit   explanation
% BF_new (1,1) struct [-]    base flow resampled onto the stability grid

fields = {'Re','lref','Uref','nu'};

% 1–3: inherit from BF_raw / params
for k = 1:numel(fields)
    f = fields{k};
    if isfield(BF_new, f) && ~isempty(BF_new.(f));  continue;  end
    if isfield(BF_raw, f) && ~isempty(BF_raw.(f))
        BF_new.(f) = BF_raw.(f);
    elseif isfield(params, f) && ~isempty(params.(f))
        BF_new.(f) = params.(f);
    end
end

% 4: derive the one missing scale from the other three, if possible
has = @(f) isfield(BF_new, f) && ~isempty(BF_new.(f));
if ~has('nu')   && has('Uref') && has('lref') && has('Re')
    BF_new.nu   = BF_new.Uref * BF_new.lref / BF_new.Re;
end
if ~has('Re')   && has('Uref') && has('lref') && has('nu')
    BF_new.Re   = BF_new.Uref * BF_new.lref / BF_new.nu;
end
if ~has('lref') && has('Uref') && has('Re')   && has('nu')
    BF_new.lref = BF_new.Re * BF_new.nu / BF_new.Uref;
end
if ~has('Uref') && has('lref') && has('Re')   && has('nu')
    BF_new.Uref = BF_new.Re * BF_new.nu / BF_new.lref;
end

if ~has('Re')
    error('ensure_metadata:missingRe', ...
        'Re is missing from BF and params. Set params.Re or include BF.Re.');
end
end
