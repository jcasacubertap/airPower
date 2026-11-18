function Grid = step_geometry(Grid)
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

%
%% Description
% Normalise the solid geometry carried on a stability grid.
%
%   Grid = step_geometry(Grid)
%
%   Grid.step is the single source of truth: its seg field holds one
%   [x_lo, x_hi, h] box per run of constant height, and seg_tol the tolerance on
%   each. The scalars StepType, StepX, StepH and the hasMask flag are derived
%   from it here, in one place, so they cannot drift apart. A grid that carries
%   only those scalars, as a caller may write by hand for a single step, is
%   accepted and the box is built from them. Calling this twice changes nothing.
%
%% Inputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    stability grid, with Grid.step or the legacy scalars
%
%% Outputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    same grid with step, StepType, StepX, StepH, hasMask consistent

have_step = isfield(Grid, 'step') && isstruct(Grid.step) && isfield(Grid.step, 'seg');
if ~have_step
    Grid.step = struct('type', 'FLAT', 'seg', zeros(0,3), 'seg_tol', zeros(0,2), ...
                       'x', 0, 'h', 0, 'x_tol', 0, 'h_tol', 0, 'nface', 0, 'snapped', false);
    % A hand-written single step: one box from the face to the end of the domain.
    if isfield(Grid, 'StepType') && ~isempty(Grid.StepType) && ~strcmpi(Grid.StepType, 'FLAT')
        assert(isfield(Grid,'StepX') && isfield(Grid,'StepH') && ~isempty(Grid.StepX) && ~isempty(Grid.StepH), ...
            'step_geometry: Grid.StepType = %s needs Grid.StepX and Grid.StepH.', Grid.StepType);
        x_end = max(Grid.xun(:));
        Grid.step.type    = Grid.StepType;
        Grid.step.seg     = [Grid.StepX, x_end, Grid.StepH];
        Grid.step.seg_tol = [0, 0];
        Grid.step.x       = Grid.StepX;
        Grid.step.h       = Grid.StepH;
        Grid.step.nface   = 1;
    end
end
if ~isfield(Grid.step, 'seg_tol') || size(Grid.step.seg_tol,1) ~= size(Grid.step.seg,1)
    Grid.step.seg_tol = zeros(size(Grid.step.seg,1), 2);     % written before seg_tol existed
end
if ~isfield(Grid.step, 'snapped');  Grid.step.snapped = false;  end

Grid.StepType = Grid.step.type;
Grid.StepX    = Grid.step.x;
Grid.StepH    = Grid.step.h;
Grid.hasMask  = ~strcmpi(Grid.step.type, 'FLAT') && ~isempty(Grid.step.seg);
end
