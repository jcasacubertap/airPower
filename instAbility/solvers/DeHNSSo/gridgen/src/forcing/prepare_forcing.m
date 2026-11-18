function f = prepare_forcing(spec)
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
% Wall-forcing spec -> refinement extent + profile.
%
%   f = prepare_forcing(spec)
%
%   Turns a description of what is at the wall into the centre / half-width used
%   to cluster the mesh and a unit-peak profile f.eval(x). spec is either a
%   parameterized actuator (fields xc, w, optional shape = 'cosine') or a
%   distribution (fields x, v), whose xc/w are taken from the support of v.
%
%   A tabulated profile supplies the SHAPE only: spec.v is stretched onto
%   [xc-w, xc+w], so xc and w place it exactly as they do a cosine and the data
%   may be in any units and any order,
%
%       struct('xc', Xc, 'w', w, 'v', h)          % samples assumed equispaced
%       struct('xc', Xc, 'w', w, 'x', xh, 'v', h) % xh gives their relative spacing
%
%   spec.x is only ever used for spacing, never for position. Given x and v with
%   no xc/w, the centre and half-width are instead taken from the support of v,
%   in the units of x.
%
%   Returns f.xc and f.w, the centre and half-width in physical x; f.eval, a
%   unit-peak complex profile that is zero outside [xc-w, xc+w]; and f.kind,
%   'cosine' or 'distribution'.

%% Inputs
% Name size         unit   explanation
% spec (1,1) struct [-]    description of what sits at the wall for one forcing source

%% Outputs
% Name size         unit   explanation
% f    (1,1) struct [-]    refinement extent and wall profile: centre, half-width and shape

% Distribution: spec.v is a shape sampled across the source. With xc and w it is
% stretched onto [xc-w, xc+w]; without them the support of v sets both, in the
% units of spec.x.
if isfield(spec, 'v') && ~isempty(spec.v)
    vd = spec.v(:).';
    assert(numel(vd) >= 2, 'prepare_forcing: spec.v needs at least two samples.');
    peak = max(abs(vd));
    assert(peak > 0, 'prepare_forcing: spec.v is all zero.');
    if isfield(spec, 'x') && ~isempty(spec.x)
        xs = spec.x(:).';
        assert(numel(xs) == numel(vd), ...
            'prepare_forcing: spec.x and spec.v must be equal-length vectors.');
        [xs, isort] = sort(xs);   vd = vd(isort);
        assert(all(diff(xs) > 0), 'prepare_forcing: spec.x has repeated points.');
    else
        xs = linspace(0, 1, numel(vd));                   % equispaced samples
    end

    has_place = isfield(spec,'xc') && isfield(spec,'w') && ...
                ~isempty(spec.xc) && ~isempty(spec.w);
    if has_place
        % Stretch the samples onto [xc-w, xc+w]: only their relative spacing is
        % used, so the tabulated units and origin never enter.
        f.xc = spec.xc;   f.w = spec.w;
        assert(f.w > 0, 'prepare_forcing: w must be positive.');
        t  = 2 * (xs - xs(1)) / (xs(end) - xs(1)) - 1;    % -> [-1, 1]
        xd = f.xc + f.w * t;
    else
        % No placement given: read it off the support, in the units of spec.x.
        assert(isfield(spec,'x') && ~isempty(spec.x), ...
            'prepare_forcing: spec.v needs either xc and w, or x to place it.');
        act  = find(abs(vd) > 1e-6 * peak);
        f.xc = 0.5 * (xs(act(end)) + xs(act(1)));
        f.w  = 0.5 * (xs(act(end)) - xs(act(1)));
        assert(f.w > 0, 'prepare_forcing: degenerate support.');
        xd   = xs;
    end
    f.eval = @(x) interp1(xd, vd / peak, x, 'pchip', 0);  % unit peak, phase kept
    f.kind = 'distribution';
    return
end

assert(isfield(spec,'xc') && isfield(spec,'w') && ~isempty(spec.xc) && ~isempty(spec.w), ...
    'prepare_forcing: actuator spec needs xc and w (or v for a tabulated shape).');
shape = 'cosine';
if isfield(spec,'shape') && ~isempty(spec.shape); shape = lower(spec.shape); end
assert(strcmp(shape,'cosine'), 'prepare_forcing: unknown shape ''%s'' (only ''cosine'').', shape);

% Actuator: raised cosine, unit peak at xc and zero at xc +- w, so w is the
% HALF-width and the footprint spans 2w.
xc = spec.xc;  w = spec.w;
f.xc   = xc;
f.w    = w;
f.eval = @(x) (abs(x - xc) <= w) .* 0.5 .* (1 + cos(pi * (x - xc) / w));
f.kind = shape;
end
