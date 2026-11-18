function [Stab, info] = add_wall_forcing(Stab, StabGrid, Opt, spec)
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
% Apply a wall forcing defined in gridgen (or overridden here).
%
%   [Stab, info] = add_wall_forcing(Stab, StabGrid, Opt, spec)
%
%   Grid-first wall forcing. gridgen already decides WHAT sits at the wall (an
%   actuator or a v-distribution) via params.forcing, refines the mesh around
%   it, and stores the unit-peak profile as StabGrid.forcing. This function
%   scales that profile by an amplitude and writes it into Stab.bcw, handling
%   the details: which velocity component, which mode index, sizing,
%   Stab.bcx, the conjugate mode for nonlinear runs, and adding to (not
%   overwriting) any forcing already stored there.
%
%   spec fields:
%     A      peak amplitude, physical 2|u| convention like Stab.A0 (required);
%            exact on the solver grid, max|Stab.bcw| == A
%     phi    phase [rad] (default 0)
%     comp   velocity component 1 u / 2 v / 3 w (default 2)
%     mode   [m n] mode the forcing drives (default [1 0], fundamental)
%   Where the profile comes from (checked in this order):
%     1. spec.x + spec.v        a distribution given here (override)
%     2. spec.xc + spec.w       a raised-cosine actuator built here (override)
%     3. StabGrid.forcing       the profile gridgen stored (default, recommended)
%   spec.which picks one if gridgen stored several (params.forcing a cell array).
%
%   The conjugate mode is added automatically for nonlinear runs (Opt.linear
%   'off'); for linear runs only the fundamental is set (reconstruction supplies
%   the conjugate later). Call this once per forcing to stack several of them.
%
%   Example (profile from gridgen):
%     Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 2e-4));

%% Inputs
% Name     size         unit   explanation
% Stab     (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Opt      (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% spec     (1,1) struct [-]    forcing specification: amplitude A, and optionally phase, velocity component and mode

%% Outputs
% Name size         unit   explanation
% Stab (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% info (1,1) struct [-]    what was applied: amplitude, phase, component, mode, the affected stations and the sampled peak before renormalisation

x  = StabGrid.xun(1, :);
nx = numel(x);

if ~isfield(spec, 'A') || isempty(spec.A)
    error('add_wall_forcing: spec.A (amplitude) is required.');
end
phi  = 0;     if isfield(spec,'phi')  && ~isempty(spec.phi);  phi  = spec.phi;  end
comp = 2;     if isfield(spec,'comp') && ~isempty(spec.comp); comp = spec.comp; end
mode = [1 0]; if isfield(spec,'mode') && ~isempty(spec.mode); mode = spec.mode; end

% ---- resolve the unit-peak profile on the grid ------------------------------
if (isfield(spec,'x') && isfield(spec,'v')) || (isfield(spec,'xc') && isfield(spec,'w'))
    f     = prepare_forcing(spec);           % caller override (distribution or actuator)
    vprof = f.eval(x);
    xc = f.xc;  w = f.w;  kind = f.kind;
elseif isfield(StabGrid, 'forcing') && ~isempty(StabGrid.forcing)
    F = StabGrid.forcing;
    if numel(F) > 1                              % grid stores several forcings
        assert(isfield(spec,'which') && ~isempty(spec.which), ...
            'add_wall_forcing: grid has %d forcings: pass spec.which to pick one.', numel(F));
        F = F(spec.which);
    end
    vprof = interp1(F.x, F.v, x, 'pchip', 0);   % F.x is already the grid, so this interp is a no-op
    xc = F.xc;  w = F.w;  kind = F.kind;
else
    error(['add_wall_forcing: no StabGrid.forcing and no override. Define ' ...
           'params.forcing in gridgen, or pass spec.x/spec.v (distribution) ' ...
           'or spec.xc/spec.w (actuator).']);
end

% ---- renormalise on the solver grid -------------------------------------------
pk = max(abs(vprof));
assert(pk > 0, 'add_wall_forcing:offGrid', ...
    'Forcing profile is zero on every grid station (xc=%.4g, w=%.4g).', xc, w);
vprof = vprof / pk;

% ---- resolution check: does the grid resolve the profile? -------------------
active = abs(vprof) > 1e-6;
if nnz(active) < 4
    warning('add_wall_forcing:underResolved', ...
        ['Wall forcing spans only %d grid stations. Refine the mesh in gridgen: ' ...
         'params.forcing = struct(...) sizes params.mug=%.4g, spread=%.4g.'], ...
        nnz(active), xc, w);
end

% ---- scale and write additively into Stab.bcw -------------------------------
v = spec.A * vprof * exp(1i*phi);

n_modes  = (2*Stab.N + 1) * (2*Stab.M + 1);
% u,v,w [,T]: the compressible operator carries a fourth wall channel, and init
% asserts the row count, so size it from Opt rather than assuming incompressible.
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
n_bc     = 3 + double(is_compr);
if ~isfield(Stab, 'bcw') || isempty(Stab.bcw) || ...
        size(Stab.bcw, 2) ~= nx || size(Stab.bcw, 1) ~= n_bc
    Stab.bcw = zeros(n_bc, nx, n_modes);
end
assert(comp <= n_bc, 'add_wall_forcing:badComp', ...
    'spec.comp = %d but only %d wall channels exist (Opt.compressible=%s).', ...
    comp, n_bc, string(is_compr));
Stab.bcx = x;

j = ModeToModeNumber(mode(1), mode(2), Stab.M, Stab.N);
Stab.bcw(comp, :, j) = Stab.bcw(comp, :, j) + reshape(v, 1, nx);

is_nonlinear = isfield(Opt,'linear') && strcmpi(Opt.linear,'off');
jc = [];
if is_nonlinear
    jc = ModeToModeNumber(-mode(1), -mode(2), Stab.M, Stab.N);
    Stab.bcw(comp, :, jc) = Stab.bcw(comp, :, jc) + reshape(conj(v), 1, nx);
end

if nargout > 1
    info = struct('xc', xc, 'w', w, 'kind', kind, 'active', active, ...
                  'idx', find(active), 'x', x, 'comp', comp, 'mode', mode, ...
                  'j', j, 'jc', jc, 'peak_sampled', pk);
end

fprintf(['Wall forcing (%s): comp %d, mode [%d %d], xc=%.4g, w=%.4g, A=%.3e, ' ...
         '%d stations%s\n'], kind, comp, mode(1), mode(2), xc, w, spec.A, ...
        nnz(active), ternary(is_nonlinear, ' (+conjugate)', ''));
end

function s = ternary(c, a, b)
if c; s = a; else; s = b; end
end
