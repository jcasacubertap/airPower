function [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
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
% Solve the harmonic Navier-Stokes (HNS) system for all active modes.
%
%   Every mode with the same (omega, beta) shares one operator, so the modes are
%   grouped by operator, each group is factorised once and its modes are solved
%   against that factor. Opt.lu_mode picks the back-end:
%     'full'           sparse LU, every factor kept
%     'auto'           sparse LU, factors kept up to Opt.lu_max_cache (default)
%     'seq'            sparse LU, one factor at a time, dropped after use
%     'none'           backslash on the assembled operator, nothing kept
%     'lapack_band'    banded LU (compiled MEX), factors kept as for 'auto'
%     'block_thomas'   block Thomas (Paredes 2024), compressible only, kept as for 'auto'
%   Opt.lu_cache = 'off' makes the kept back-ends factor and solve one-shot.
%   Opt.iterref = n adds n passes of iterative refinement wherever a factor is held.
%   Opt.parfor = 'on' solves the modes in parallel (sparse LU) or factorises the
%   block-Thomas groups in parallel.
%   Opt.bt_apply chooses how block_thomas applies the diagonal station block:
%     'lu' (default) stores its pivoted LU factors, 'inverse' stores G^-1.

%% Inputs
% Name     size              unit   explanation
% BF       (1,1) struct      [-]    base flow on its own grid
% Bound    (1,1) struct      [-]    boundary conditions and outflow-buffer coefficients
% StabGrid (1,1) struct      [-]    stability grid, coordinate metrics and base-flow fields on it
% Init     (1,1) struct      [-]    inflow initialization data
% Mode     (1,1) struct      [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt      (1,1) struct      [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab     (1,1) struct      [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes  (1,1) struct      [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

%% Outputs
% Name     size              unit   explanation
% StabRes  (1,1) struct      [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% Opt      (1,1) struct      [-]    solver options, carrying the factor cache
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

%% Options
if ~isfield(Opt, 'lu_mode');  Opt.lu_mode = 'auto'; end
if ~isfield(Opt, 'parfor');   Opt.parfor  = 'off';  end
if ~isfield(Opt, 'lu_ordering') || isempty(Opt.lu_ordering); Opt.lu_ordering = 'auto'; end
if ~isfield(Opt, 'bt_apply')    || isempty(Opt.bt_apply);    Opt.bt_apply    = 'lu';   end
if ~isfield(Opt, 'bt_np')       || isempty(Opt.bt_np);       Opt.bt_np       = 2;      end
if ~isfield(Opt, 'LU_cache');   Opt.LU_cache = struct();  end
if ~isfield(Opt, 'lu_order');   Opt.lu_order = {};        end   % kept factors, oldest first
use_scaling = isfield(Opt,'lu_scaling') && strcmpi(Opt.lu_scaling,'on');
niter = 0;
if isfield(Opt,'iterref') && ~isempty(Opt.iterref); niter = Opt.iterref; end
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible, 'on');
be       = lower(Opt.lu_mode);
phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));
nphys      = numel(phys_modes);
if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = max(8, nphys); end

switch be
    case {'full','auto','seq','none'}
    case 'lapack_band'
        if exist(['dgbsv_factor.' mexext],'file') ~= 3 || exist(['dgbsv_solve_cached.' mexext],'file') ~= 3
            error('solve_modes:no_mex', 'banded-LU MEX not built. Run src/linsolve/lapack_band/mexbuild_dgbsv.m first.');
        end
    case 'block_thomas'
        if ~is_compr
            error('solve_modes:bt_incompressible', ...
                ['lu_mode=''block_thomas'' is for the COMPRESSIBLE operator only. The ' ...
                 'incompressible saddle-point (zero-diagonal continuity rows) makes the ' ...
                 'unpivoted block factorization near-singular and the result SILENTLY ' ...
                 'WRONG (verified ~2x error). Use lu_mode=''lapack_band'', ''seq'', or ''auto''.']);
        end
    otherwise
        error('solve_modes:lu_mode', 'Unknown Opt.lu_mode ''%s''.', Opt.lu_mode);
end
% Factors are kept between solves by 'full' and 'auto' always, by the two MEX
% back-ends unless caching is switched off, and never by 'seq' and 'none'.
switch be
    case {'full','auto'};  keep = true;
    case {'seq','none'};   keep = false;
    otherwise;             keep = Opt.lu_max_cache > 0 && ~(isfield(Opt,'lu_cache') && strcmpi(Opt.lu_cache,'off'));
end
bt_ooc = strcmp(be,'block_thomas') && isfield(Opt,'bt_ooc') && strcmpi(Opt.bt_ooc,'on');
% Out-of-core factors go under Opt.bt_ooc_dir, tempdir by default. A large case
% wants a disk with room; MATLAB's own tempdir must stay on a disk it can run from.
ooc_root = tempdir;
if isfield(Opt,'bt_ooc_dir') && ~isempty(Opt.bt_ooc_dir);  ooc_root = Opt.bt_ooc_dir;  end
%% Phase 1: Build LHS matrices (once per run)
needs_beta = any(StabRes.betavec ~= 0);
if Opt.mmry == 1 && needs_beta && (isempty(MLA3) || isempty(MLA4))
    fprintf(' MLA3/MLA4 missing for 3D run, rebuilding LHS.\n');
    Opt.mmry = 0;
end
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible, 'on');
if Opt.mmry ~= 1
    Opt.mmry = 1;
    Opt.hr_in_MLA1 = false;
    Opt.hr_rows    = [];
    fprintf(' Constructing LHS matrices ...\n');
    tic
    if is_compr
        [MLAMFD, MLA1, MLA2, MLA3, MLA4, MLA5] = LHS_compressible(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
    else
        [MLAMFD, MLA1, MLA2, MLA3, MLA4]       = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
        MLA5 = sparse(size(MLA1,1), size(MLA1,2));
    end
    MFD_diff = MLAMFD - MLA1;
    clear MLAMFD
    fprintf(' LHS done in %.1f s\n', toc);
    if ~isempty(getenv('MFD_OPDUMP'))
        ML_mfd = MLA1 + MFD_diff;                       % MLAMFD (omega=beta=0)
        mfd_j  = find(Stab.Mvec==0 & Stab.Nvec==0, 1);
        R_mfd  = RHS(Bound, StabGrid, Mode, Opt, StabRes, mfd_j);
        f_mfd  = Stab.f(mfd_j,:);
        q_mfd  = StabRes.q(mfd_j,:);
        Xg = StabGrid.x(:).'; Yg = StabGrid.y(:); nyg = numel(Yg); nxg = numel(Xg);
        save('/tmp/mfd_op.mat','ML_mfd','R_mfd','f_mfd','q_mfd','Xg','Yg','nyg','nxg','mfd_j','-v7.3');
        fprintf(' [MFD_OPDUMP] saved MFD operator+forcing -> /tmp/mfd_op.mat\n');
    end
end

%% Phase 1b: Resonator wall rows (also applied to a restored LHS)
% Helmholtz-resonator impedance wall boundary condition (BC), plus an optional
% free-slip u' over the mouth. Applied to MLA1 only, so the (0,0) mean-flow
% distortion (MFD) mode keeps no-slip while the fluctuations see the resonator
% (MFD_diff = MLAMFD - MLA1 restores no-slip for the MFD mode).
want_hr = isfield(Bound,'hr') && ~isempty(Bound.hr);
if ~Opt.hr_in_MLA1 && (want_hr || ~isempty(Opt.hr_rows))
    MLAMFD = MLA1 + MFD_diff;                        % invariant under the wall rows
    if ~isempty(Opt.hr_rows)
        MLA1(Opt.hr_rows, :) = MLAMFD(Opt.hr_rows, :);
        Opt.hr_rows = [];
    end
    if want_hr
        [MLA1, Opt.hr_rows] = apply_resonator_wall(MLA1, Bound.hr, BF, StabGrid, is_compr);
    end
    MFD_diff = MLAMFD - MLA1;
    Opt.hr_in_MLA1 = want_hr;   % Picard calls solve_modes per iteration; write once
    clear MLAMFD
end

%% Forcing: zero it on the rows that carry a boundary condition
q_out  = zeros(nphys, size(StabRes.q, 2));
Stab_f = Stab.f;
bc_dir     = full(sum(MLA1 ~= 0, 2)) == 1;                  % single-entry rows
bc_dir_mfd = full(sum((MLA1 + MFD_diff) ~= 0, 2)) == 1;
for jj = 1:size(Stab_f, 1)
    if Stab.Mvec(jj) == 0 && Stab.Nvec(jj) == 0
        Stab_f(jj, bc_dir_mfd) = 0;
    else
        Stab_f(jj, bc_dir) = 0;
    end
end
% Ghost rows have several entries, so the count above misses them.
if any(Bound.ghost_idx)
    ny_g = StabGrid.ny;  nvar_g = size(Stab_f,2) / (ny_g*StabGrid.nx);
    for i_g = find(Bound.ghost_idx)
        rows = (i_g-1)*nvar_g*ny_g + Bound.ghost_idx(i_g) + (0:2)*ny_g;   % u, v, w
        Stab_f(:, rows) = 0;
    end
end

%% Group the modes by operator
groups = group_modes(phys_modes, StabRes.omegavec, StabRes.betavec, Stab.Mvec, Stab.Nvec);
ngrp   = numel(groups);
needed = {groups.key};
ctx = struct('be', be, 'use_symrcm', pick_ordering(Opt, is_compr), 'use_scaling', use_scaling, ...
             'b', round(size(MLA1,1) / StabGrid.nx), 'np', Opt.bt_np, 'amode', Opt.bt_apply, 'ooc', bt_ooc);

%% Phase 2: factors kept in Opt.LU_cache
if keep
    % Block Thomas can factorise the missing groups in parallel; the factor is the
    % expensive step. Needs a Threads pool so the operator blocks are shared.
    if strcmp(be,'block_thomas') && strcmpi(Opt.parfor,'on') && ~bt_ooc
        miss = find(~arrayfun(@(g) has_factor(Opt.LU_cache, g.key), groups));
        if ~isempty(miss)
            fprintf(' Block-Thomas: parallel factor of %d group(s) [np=%d] ...\n', numel(miss), Opt.bt_np);
            tf = tic;  Fnew = cell(numel(miss),1);  gm = groups(miss);  ctx_p = ctx;
            parfor ii = 1:numel(miss)
                Fnew{ii} = be_factor(ctx_p, assemble_ML(gm(ii).omega, gm(ii).beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff), '');
            end
            fprintf(' Block-Thomas parallel factor done (%.1f s)\n', toc(tf));
            for ii = 1:numel(miss);  Opt.LU_cache.(gm(ii).key) = Fnew{ii};  Fnew{ii} = [];  end
        end
    end
    for g = 1:ngrp
        G = groups(g);
        if has_factor(Opt.LU_cache, G.key);  continue;  end
        Opt = evict(Opt, needed, G.label, nphys);
        fprintf(' %s (%s) [%s, %d/%d kept] ...', factor_name(be), G.label, be, ...
                numel(fieldnames(Opt.LU_cache)) + 1, Opt.lu_max_cache);
        tic_f = tic;
        oocdir = '';
        if bt_ooc;  oocdir = [tempname(ooc_root) '_btcache'];  block_thomas('free', oocdir, 'reg');  end
        Opt.LU_cache.(G.key) = be_factor(ctx, assemble_ML(G.omega, G.beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff), oocdir);
        Opt.lu_order{end+1} = G.key;
        fprintf(' done (%.1f s)\n', toc(tic_f));
    end
    % Once a factor is held for every operator the run can ever need, not only
    % those active in this call, the blocks are not needed again unless a
    % refinement residual has to be assembled.
    all_keys = unique(arrayfun(@(j) ml_key(StabRes.omegavec(j), StabRes.betavec(j)), 1:Mode.nf, 'UniformOutput', false));
    if niter == 0 && all(isfield(Opt.LU_cache, all_keys)) && any(strcmp(be, {'full','auto'}))
        MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MLA5 = []; MFD_diff = [];
    end
end

%% Phase 3: solve
if keep && any(strcmp(be, {'full','auto','none'})) && strcmpi(Opt.parfor, 'on')
    % Sparse LU only: one mode per worker, with the factors sliced in.
    fprintf(' Solving %d modes in parallel ...\n', nphys);
    tic_s = tic;
    F_k = cell(nphys,1);  key_k = cell(nphys,1);
    for g = 1:ngrp;  for k = groups(g).k;  key_k{k} = groups(g).key;  end;  end
    for k = 1:nphys;  F_k{k} = Opt.LU_cache.(key_k{k});  end
    Stab_fk = Stab_f(phys_modes, :);
    parfor k = 1:nphys
        j   = phys_modes(k);
        R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
        f_j = Stab_fk(k,:);  f_j(R_j ~= 0) = 0;
        q_out(k,:) = lu_solve(F_k{k}, R_j + f_j.').';
    end
    fprintf(' Solve done in %.1f s\n', toc(tic_s));
else
    for g = 1:ngrp
        G = groups(g);
        tic_g = tic;
        if keep
            F = Opt.LU_cache.(G.key);  ML = [];
            fprintf(' Solving (%s) [%s, kept]', G.label, be);
        else
            ML = assemble_ML(G.omega, G.beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
            oocdir = '';
            if bt_ooc;  oocdir = fullfile(ooc_root, sprintf('btooc_g%d', g));  end
            F = be_factor_oneshot(ctx, ML, niter, oocdir);
            fprintf(' Factor+solve (%s) [%s, %d mode(s)]', G.label, be, numel(G.k));
        end
        if niter > 0 && isempty(ML)                             % the residual needs the operator
            ML = assemble_ML(G.omega, G.beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
        end
        if strcmp(be, 'block_thomas')
            rhs = zeros(size(MLA1,1), numel(G.k));
            for ki = 1:numel(G.k)
                rhs(:,ki) = mode_rhs(Bound, StabGrid, Mode, Opt, StabRes, phys_modes(G.k(ki)), Stab_f);
            end
            X = be_solve(ctx, F, rhs);
            for ki = 1:numel(G.k);  q_out(G.k(ki),:) = X(:,ki).';  end
        else
            for k = G.k
                rhs = mode_rhs(Bound, StabGrid, Mode, Opt, StabRes, phys_modes(k), Stab_f);
                if strcmp(be, 'lapack_band');  rhs = full(rhs);  end
                x = be_solve(ctx, F, rhs);
                for it = 1:niter
                    x = x + be_solve(ctx, F, rhs - ML*x);
                end
                q_out(k,:) = x.';
            end
        end
        fprintf(' ... %.1f s\n', toc(tic_g));
        clear F ML rhs X
    end
end

%% Collect results (including conjugate modes)
for k = 1:nphys
    j = phys_modes(k);
    StabRes.q(j,:)            = q_out(k,:);
    StabRes.q(Mode.nf+1-j,:)  = conj(q_out(k,:));
end

if ~isempty(getenv('MFD_QDUMP'))
    mfd_j = find(Stab.Mvec==0 & Stab.Nvec==0, 1);
    kmfd  = find(phys_modes==mfd_j, 1);
    if ~isempty(kmfd)
        q_mfd_conv = q_out(kmfd,:);                     % CONVERGED MFD field (operator-ordered)
        Xg = StabGrid.x(:).'; Yg = StabGrid.y(:); nyg = numel(Yg); nxg = numel(Xg);
        save('/tmp/mfd_qconv.mat','q_mfd_conv','Xg','Yg','nyg','nxg','mfd_j','-v7.3');  % overwritten each Picard iter -> last = converged
    end
end
end

%% Local helpers

function groups = group_modes(phys_modes, omegavec, betavec, Mvec, Nvec)
%GROUP_MODES  One entry per distinct (omega, beta), with the modes that share it.
groups = struct('key', {}, 'omega', {}, 'beta', {}, 'k', {}, 'label', {});
for k = 1:numel(phys_modes)
    j   = phys_modes(k);
    key = ml_key(omegavec(j), betavec(j));
    g   = find(strcmp({groups.key}, key), 1);
    if isempty(g)
        groups(end+1) = struct('key', key, 'omega', omegavec(j), 'beta', betavec(j), ...
                               'k', k, 'label', sprintf('%d,%d', Mvec(j), Nvec(j)));
    else
        groups(g).k(end+1) = k;
    end
end
end

function key = ml_key(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function rhs = mode_rhs(Bound, StabGrid, Mode, Opt, StabRes, j, Stab_f)
%MODE_RHS  Boundary forcing plus nonlinear forcing, the latter zeroed where a
% boundary value is imposed.
R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
f_j = Stab_f(j,:);
f_j(R_j ~= 0) = 0;
rhs = R_j + f_j.';
end

function tf = has_factor(cache, key)
tf = isfield(cache, key) && ~isempty(cache.(key));
end

function name = factor_name(be)
switch be
    case 'lapack_band';   name = 'Band-LU factor';
    case 'block_thomas';  name = 'Block-Thomas factor';
    case 'none';          name = 'Operator';
    otherwise;            name = 'LU';
end
end

function Opt = evict(Opt, needed, label, nphys)
%EVICT  Make room in the cache, dropping the oldest factor not needed this call.
% 'full' never drops. When every kept factor is needed the cache grows and warns.
if strcmpi(Opt.lu_mode, 'full');  return;  end
Opt.lu_order = Opt.lu_order(isfield(Opt.LU_cache, Opt.lu_order));   % drop stale keys
if numel(fieldnames(Opt.LU_cache)) < Opt.lu_max_cache;  return;  end
old = find(~ismember(Opt.lu_order, needed), 1);
if isempty(old)
    if ~isfield(Opt,'lu_grew_warned') || ~Opt.lu_grew_warned
        warning(['LU cache grown from lu_max_cache=%d to %d to hold all %d active modes; ' ...
                 'set lu_max_cache accordingly or use lu_mode=''seq''/''lapack_band'' for a bounded footprint.'], ...
                Opt.lu_max_cache, numel(fieldnames(Opt.LU_cache)) + 1, nphys);
        Opt.lu_grew_warned = true;
    end
    return
end
key = Opt.lu_order{old};
if isfield(Opt.LU_cache.(key), 'btF');  block_thomas('free', Opt.LU_cache.(key).btF);  end
Opt.LU_cache = rmfield(Opt.LU_cache, key);
Opt.lu_order(old) = [];
fprintf(' [dropping factor %s for %s]', key, label);
end

function F = be_factor(ctx, ML, oocdir)
%BE_FACTOR  Factor kept for reuse, in the form the back-end needs.
switch ctx.be
    case 'lapack_band'
        [ABf, IPIV, KLb, KUb] = dgbsv_factor(ML);
        F = struct('ABf', ABf, 'IPIV', IPIV, 'KL', KLb, 'KU', KUb);
    case 'block_thomas'
        F = struct('btF', block_thomas('factor', ML, ctx.b, ctx.np, oocdir, [], ctx.amode));
    case 'none'
        F = ML;
    otherwise
        F = do_lu(ML, ctx.use_symrcm, ctx.use_scaling);
end
end

function F = be_factor_oneshot(ctx, ML, niter, oocdir)
%BE_FACTOR_ONESHOT  Factor for a single group that is dropped after use. The
% banded and block-Thomas MEX routines have a combined factor-and-solve entry
% that is used when no refinement is asked for; it is marked by the field 'ML'.
switch ctx.be
    case 'lapack_band'
        if niter > 0;  F = be_factor(ctx, ML, '');
        else;          F = struct('ML', ML);  end
    case 'block_thomas'
        F = struct('ML', ML, 'oocdir', oocdir);
    otherwise
        F = be_factor(ctx, ML, '');
end
end

function x = be_solve(ctx, F, rhs)
%BE_SOLVE  Solve with a factor from be_factor or be_factor_oneshot.
switch ctx.be
    case 'lapack_band'
        if isfield(F, 'ML');  x = dgbsv_solve(F.ML, rhs);
        else;                 x = dgbsv_solve_cached(F.ABf, F.IPIV, F.KL, F.KU, rhs);  end
    case 'block_thomas'
        if isfield(F, 'ML');  x = block_thomas('solve', F.ML, ctx.b, rhs, ctx.np, F.oocdir, [], ctx.amode);
        else;                 x = block_thomas('apply', F.btF, rhs);  end
    otherwise
        x = lu_solve(F, rhs);
end
end
function x = lu_solve(lu_s, rhs)
%LU_SOLVE  Solve using cached LU.
if isstruct(lu_s) && isfield(lu_s, 'D')
    % Five-output lu: P*(D\ML)*Q = L*U  =>  x = Q*(U\(L\(P*(D\rhs))))
    if isfield(lu_s, 'r')
        y = lu_s.Q * (lu_s.U \ (lu_s.L \ (lu_s.P * (lu_s.D \ rhs(lu_s.r)))));
        x = zeros(size(rhs));
        x(lu_s.r) = y;
    else
        x = lu_s.Q * (lu_s.U \ (lu_s.L \ (lu_s.P * (lu_s.D \ rhs))));
    end
elseif isstruct(lu_s) && isfield(lu_s, 'q')
    % [L,U,p,q] = lu(ML,'vector'):  ML(p,q) = L*U  =>  x(q) = U\(L\rhs(p))
    x = zeros(size(rhs));
    x(lu_s.q) = lu_s.U \ (lu_s.L \ rhs(lu_s.p));
elseif isstruct(lu_s) && isfield(lu_s, 'dML')
    % symrcm + decomposition
    x = zeros(size(rhs));
    x(lu_s.r) = lu_s.dML \ rhs(lu_s.r);
elseif isstruct(lu_s) && isfield(lu_s, 'L')
    % symrcm + raw L,U,p
    rhs_r = rhs(lu_s.r);
    x = zeros(size(rhs));
    x(lu_s.r) = lu_s.U \ (lu_s.L \ rhs_r(lu_s.p));
else
    % decomposition object or raw sparse matrix
    x = lu_s \ rhs;
end
end

function s = do_lu(ML, use_symrcm, use_scaling)
%DO_LU  LU with a fill-reducing ordering.
if nargin < 2 || isempty(use_symrcm);  use_symrcm  = false; end
if nargin < 3 || isempty(use_scaling); use_scaling = false; end
piv = {};
if isreal(ML); piv = {1.0}; end
if use_scaling
    % Five-output lu: P*(D\ML)*Q = L*U, so the row scaling D is kept.
    if use_symrcm
        r = symrcm(ML);
        [L, U, P, Q, D] = lu(ML(r,r), piv{:});
        s = struct('L', L, 'U', U, 'P', P, 'Q', Q, 'D', D, 'r', r);
    else
        [L, U, P, Q, D] = lu(ML, piv{:});
        s = struct('L', L, 'U', U, 'P', P, 'Q', Q, 'D', D);
    end
elseif use_symrcm
    r = symrcm(ML);
    [L, U, p] = lu(ML(r,r), piv{:}, 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'r', r);
else
    [L, U, p, q] = lu(ML, piv{:}, 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'q', q);
end
end

function use_symrcm = pick_ordering(Opt, is_compr)
%PICK_ORDERING  Resolve Opt.lu_ordering to the flag do_lu expects.
switch lower(Opt.lu_ordering)
    case 'auto';   use_symrcm = is_compr;
    case 'symrcm'; use_symrcm = true;
    case 'colamd'; use_symrcm = false;
    otherwise
        error('solve_modes:lu_ordering', ...
            'Opt.lu_ordering must be ''auto'', ''symrcm'' or ''colamd'' (got ''%s'')', ...
            Opt.lu_ordering);
end
end
