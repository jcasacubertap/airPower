function step = detect_wall_solid(BF, class)
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
% Read a solid on the wall off a Cartesian base flow.
%
%   step = detect_wall_solid(BF, class)
%
%   A solid shows up in the base flow as nodes with zero velocity above the wall.
%   Each run of columns with the same top node becomes one box [x_lo, x_hi, h],
%   so one step, several steps, a backward-facing step and a block are all lists
%   of boxes. The solid must rest on the wall: an overhang, a cavity or a body in
%   the flow raises an error. A face lies between the last fluid column and the
%   first solid one, so a box edge is put halfway between them.
%
%   Only Cartesian input is read. On the other paths a solid would be treated as
%   fluid, so the function warns if it sees one and returns a flat wall.

%% Inputs
% Name  size         unit   explanation
% BF    (1,1) struct [-]    base flow with X, Y, U, V, W
% class char         [-]    result of classify_input; 'cartesian' when omitted

%% Outputs
% Name size         unit   explanation
% step (1,1) struct [-]    type, boxes with their tolerances, first rising face

% Default: a flat wall. Every early return below hands this back.
step = struct('type', 'FLAT', 'seg', zeros(0,3), 'seg_tol', zeros(0,2), ...
              'x', 0, 'h', 0, 'x_tol', 0, 'h_tol', 0, 'nface', 0);
if nargin < 2 || isempty(class);  class = 'cartesian';  end
if ~strcmp(class, 'cartesian');  warn_if_solid(BF, class);  return;  end
if ~(isvector(BF.X) && isvector(BF.Y)); return; end
x = BF.X(:).';  y = BF.Y(:);  [ny, nx] = size(BF.U);
if numel(x) ~= nx || numel(y) ~= ny; return; end

%% Solid nodes
tol   = 1e-12 * max(abs(BF.U(:)));
solid = abs(BF.U) < tol & abs(BF.V) < tol & abs(BF.W) < tol;

[ys, iy] = sort(y, 'ascend');            % work bottom-up whatever the storage order
S        = solid(iy, :);
wall     = ys <= 1e-9 * max(ys);
S(wall, :) = false;                      % the wall row is zero anyway
if ~any(S(:)); return; end

%% Top node of the solid in each column
% The solid must fill the column from the first node above the wall to its top.
k0   = find(~wall, 1);
jtop = zeros(1, nx);
bad  = false(1, nx);
for i = 1:nx
    r = find(S(:, i));
    if isempty(r); continue; end
    jtop(i) = r(end);
    bad(i)  = ~isequal(r(:).', k0:r(end));
end
if any(bad)
    error('detect_wall_solid:overhang', ...
        'The solid does not rest on the wall in %d of %d columns (first at x = %.4g). Overhangs, cavities and detached bodies are not supported.', ...
        nnz(bad), nx, x(find(bad, 1)));
end

%% Boxes: runs of columns with the same top node
starts = find([true, diff(jtop) ~= 0]);
stops  = [starts(2:end) - 1, nx];
is_solid = jtop(starts) > 0;
starts = starts(is_solid);  stops = stops(is_solid);
seg    = zeros(numel(starts), 3);
segtol = zeros(numel(starts), 2);
for r = 1:numel(starts)
    a = starts(r);  b = stops(r);
    x_lo = x(1);   xt = 0;
    if a > 1;  x_lo = 0.5*(x(a-1) + x(a));  xt = 0.5*(x(a) - x(a-1));  end
    x_hi = x(nx);
    if b < nx; x_hi = 0.5*(x(b) + x(b+1));  end
    seg(r, :)    = [x_lo, x_hi, ys(jtop(a))];
    segtol(r, :) = [xt, 0.5*local_spacing(ys, jtop(a))];
end
step.seg     = seg;
step.seg_tol = segtol;
step.nface   = nnz(diff(jtop) ~= 0);

%% Classify
nseg      = size(seg, 1);
touch_in  = any(seg(:,1) <= x(1));
touch_out = any(seg(:,2) >= x(nx));
if touch_in && touch_out && nseg == 1
    error('detect_wall_solid:raised', ...
        'The solid spans the whole domain at h = %.4g. That is a raised wall, not a step.', seg(1,3));
end
if nseg > 1;         step.type = 'STEPS';
elseif touch_out;    step.type = 'FFS';
elseif touch_in;     step.type = 'BFS';
else;                step.type = 'BLOCK';
end

%% First rising face
inc = find(diff(jtop) > 0, 1);
if ~isempty(inc)
    jt         = jtop(inc + 1);
    step.x     = 0.5 * (x(inc) + x(inc+1));
    step.h     = ys(jt);
    step.x_tol = 0.5 * (x(inc+1) - x(inc));
    step.h_tol = 0.5 * local_spacing(ys, jt);
end

%% Warn when the surface row is not part of the solid
% Some solvers leave a wall-normal or spanwise velocity on the surface row. That
% row then counts as fluid and the height comes out one node low.
Uc = abs(BF.U(iy, :));  Vc = abs(BF.V(iy, :));  Wc = abs(BF.W(iy, :));
susp = false(1, nx);
for i = 1:nx
    j = jtop(i) + 1;
    if jtop(i) == 0 || j > ny; continue; end
    susp(i) = Uc(j,i) < tol && (Vc(j,i) >= tol || Wc(j,i) >= tol);
end
if any(susp)
    warning('detect_wall_solid:surfaceRow', ...
        'In %d of %d solid columns the node above the top has zero streamwise velocity but not zero wall-normal or spanwise velocity (first at x = %.4g). The height may be one node low.', ...
        nnz(susp), nnz(jtop > 0), x(find(susp, 1)));
end
end

function d = local_spacing(ys, j)
%LOCAL_SPACING  Smaller gap between node j and its two neighbours.
d = inf;
if j > 1;          d = ys(j) - ys(j-1);          end
if j < numel(ys);  d = min(d, ys(j+1) - ys(j));  end
if ~isfinite(d);   d = 0;                        end
end

function warn_if_solid(BF, class)
%WARN_IF_SOLID  Warn when a solid seems present on a path that cannot read it.
if ~all(isfield(BF, {'U','V','W'})); return; end
tol  = 1e-12 * max(abs(BF.U(:)));
zero = abs(BF.U) < tol & abs(BF.V) < tol & abs(BF.W) < tol;
switch class
    case 'curvilinear_structured'
        zero([1 end], :) = false;            % either boundary row may be the wall
        zero(:, [1 end]) = false;
        n = nnz(zero);
    case 'unstructured'
        % The wall points cannot be told apart here: allow one zero per station.
        tolX = max(1e-6 * (max(BF.X(:)) - min(BF.X(:))), 1e-12);
        n = nnz(zero) - numel(uniquetol(BF.X(:), tolX, 'DataScale', 1));
    otherwise
        return
end
if n > 0
    warning('main_gridgen:unreadSolid', ...
        '%d zero-velocity nodes above the wall look like a solid. Solids are read from Cartesian input only; here they are treated as fluid.', n);
end
end
