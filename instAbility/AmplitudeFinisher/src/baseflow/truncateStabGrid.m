function [StabGrid, info] = truncateStabGrid(StabGrid, xEndPhys, in)
% TRUNCATESTABGRID  Keep only the streamwise columns up to xEndPhys.
%
%   [StabGrid, info] = truncateStabGrid(StabGrid, xEndPhys, in)
%
% HNS cost scales with the number of streamwise (xi) stations. When only the
% first few PIV planes are matched, running the full domain is wasteful. This
% slices every nx-indexed StabGrid field to the columns whose physical
% arc-length (from the LE) is <= xEndPhys, so the solve runs on the shorter
% domain. All wall-normal (ny) fields, scalars, and differentiation matrices
% are left untouched; the xi spacing (dxi) is preserved.
%
% info.nxOld / info.nxNew / info.xEndPhys report what was kept.

origNx = StabGrid.nx;
lref   = StabGrid.lref;

% Physical arc-length of each column, from the LE, in the PIV frame.
xphys = StabGrid.x(:).';
if max(abs(xphys)) > 1; xphys = xphys * lref; end
xphys = xphys + in.match.xOffset;

nxKeep = find(xphys <= xEndPhys, 1, 'last');
if isempty(nxKeep) || nxKeep < 10
    error('truncateStabGrid:tooShort', ...
          'xEndPhys=%.4g m keeps only %d columns — check nPlanes/xOffset.', ...
          xEndPhys, max(nxKeep,0));
end
if nxKeep >= origNx
    info = struct('nxOld',origNx,'nxNew',origNx,'xEndPhys',xphys(end));
    return;   % nothing to trim
end

% Slice every field whose 2nd dim (matrix) or length (vector) equals origNx.
fn = fieldnames(StabGrid);
for i = 1:numel(fn)
    v = StabGrid.(fn{i});
    if ~isnumeric(v); continue; end
    if isrow(v) && numel(v) == origNx
        StabGrid.(fn{i}) = v(1:nxKeep);
    elseif iscolumn(v) && numel(v) == origNx
        StabGrid.(fn{i}) = v(1:nxKeep);
    elseif ismatrix(v) && size(v,2) == origNx && size(v,1) ~= origNx
        StabGrid.(fn{i}) = v(:, 1:nxKeep);
    end
end

StabGrid.nx = nxKeep;
if isfield(StabGrid,'L'); StabGrid.L = StabGrid.x(end); end

info = struct('nxOld',origNx,'nxNew',nxKeep,'xEndPhys',xphys(nxKeep));
end
