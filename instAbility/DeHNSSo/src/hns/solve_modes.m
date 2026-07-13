function [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%SOLVE_MODES  Solve the HNS system for all active modes.
%
%   LU caching controlled by Opt.lu_mode:
%     'lapack_band' — LAPACK zgbsv/dgbsv via MEX, fastest single solve
%                     (~3x vs seq for compressible TS); needs ~12-20 GB.
%                     Build once via src/hns/mexbuild_dgbsv.m.
%     'full'        — cache all LU factorisations (fastest reuse, most memory)
%     'auto'        — LRU cache with soft limit Opt.lu_max_cache. Eviction skips
%                     keys still needed by the current iteration, so the cache
%                     grows to max(lu_max_cache, nphys) when necessary (one-time
%                     warning on growth).
%     'seq'         — symrcm + lu one mode-group at a time, discard after solve
%                     (medium memory, required for compressible 5-var on ~32 GB)
%     'none'        — use backslash each time, no caching (slowest, least memory)
%
%   Parallel solve controlled by Opt.parfor (set automatically in init):
%     true  — parfor over modes (requires Parallel Computing Toolbox)
%     false — standard for loop
%
%   GPU solve controlled by Opt.gpu (set automatically in init):
%     true  — transfer ML to GPU, solve with gpuArray \ (requires PCT + NVIDIA GPU)
%             Note: GPU does not support decomposition caching; ML is stored as
%             gpuArray and refactorised each solve. Beneficial for large ny (>>80).
%     false — CPU solve (default)

if ~isfield(Opt, 'lu_mode'); Opt.lu_mode = 'full'; end
if ~isfield(Opt, 'parfor');  Opt.parfor  = 'off';  end
if ~isfield(Opt, 'gpu');     Opt.gpu     = 'off';  end

phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));

%% Phase 1: Build LHS matrices (once per run)
needs_beta = any(StabRes.betavec ~= 0);
if Opt.mmry == 1 && needs_beta && (isempty(MLA3) || isempty(MLA4))
    fprintf(' MLA3/MLA4 missing for 3D run — rebuilding LHS.\n');
    Opt.mmry = 0;
end
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible, 'on');
if Opt.mmry ~= 1
    Opt.mmry = 1;
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
end

%% Phase 2 & 3: Factorise and solve
nphys  = numel(phys_modes);
q_out  = zeros(nphys, size(StabRes.q, 2));
Stab_f = Stab.f;
Mvec   = Stab.Mvec;
Nvec   = Stab.Nvec;

if strcmpi(Opt.lu_mode, 'seq')
    % Sequential mode: factorise one matrix, solve all modes that use it,
    % discard.  Only one LU in memory at a time — minimal memory footprint.
    unique_keys = {};
    mode_groups = {};
    for k = 1:nphys
        j = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key; %#ok<AGROW>
            mode_groups{end+1} = k;   %#ok<AGROW>
        else
            mode_groups{idx}(end+1) = k; %#ok<AGROW>
        end
    end

    for g = 1:numel(unique_keys)
        k_list = mode_groups{g};
        j0 = phys_modes(k_list(1));
        omega = StabRes.omegavec(j0);
        beta  = StabRes.betavec(j0);
        fprintf(' LU+solve (%d,%d) [seq, %d mode(s)] ...', Mvec(j0), Nvec(j0), numel(k_list));
        tic_g = tic;

        ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
        r  = symrcm(ML);
        [L, U, p] = lu(ML(r,r), 'vector');
        clear ML

        for ki = 1:numel(k_list)
            k = k_list(ki);
            j = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            rhs   = R_j + f_j.';
            rhs_r = rhs(r);
            x_r   = U \ (L \ rhs_r(p));
            x     = zeros(size(rhs));
            x(r)  = x_r;
            q_out(k,:) = x.';
        end
        clear L U p r
        fprintf(' done (%.1f s)\n', toc(tic_g));
    end

elseif strcmpi(Opt.lu_mode, 'lapack_band')
    % Banded LU via LAPACK zgbsv MEX (Apple Accelerate-backed).
    % ~3x faster than 'seq' on compressible 5-var, but uses ~20 GB for nx=800
    % at bsize=405 (LDAB ~ 3*bsize complex doubles).
    if exist(['dgbsv_solve.' mexext], 'file') ~= 3
        error('solve_modes:no_mex', ...
              'dgbsv_solve MEX not built. Run src/hns/mexbuild_dgbsv.m first.');
    end
    unique_keys = {};
    mode_groups = {};
    for k = 1:nphys
        j = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key; %#ok<AGROW>
            mode_groups{end+1} = k;   %#ok<AGROW>
        else
            mode_groups{idx}(end+1) = k; %#ok<AGROW>
        end
    end

    for g = 1:numel(unique_keys)
        k_list = mode_groups{g};
        j0 = phys_modes(k_list(1));
        omega = StabRes.omegavec(j0);
        beta  = StabRes.betavec(j0);
        fprintf(' Band-LU solve (%d,%d) [lapack_band, %d mode(s)] ...', ...
                Mvec(j0), Nvec(j0), numel(k_list));
        tic_g = tic;

        ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);

        for ki = 1:numel(k_list)
            k = k_list(ki);
            j = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            rhs = full(R_j + f_j.');
            q_out(k,:) = dgbsv_solve(ML, rhs).';
        end
        clear ML
        fprintf(' done (%.1f s)\n', toc(tic_g));
    end

else
    % Cache-based modes: 'full', 'auto', 'none'
    if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end

    % Keys needed by the solve loop this iteration — 'auto' eviction must not
    % remove any of these (otherwise the solve below reads a missing entry).
    needed_keys = cell(1, nphys);
    for k = 1:nphys
        jk = phys_modes(k);
        needed_keys{k} = ml_key(StabRes.omegavec(jk), StabRes.betavec(jk));
    end
    needed_set = unique(needed_keys);

    for j = phys_modes
        omega   = StabRes.omegavec(j);
        beta    = StabRes.betavec(j);
        ML_key  = ml_key(omega, beta);

        if ~isfield(Opt.LU_cache, ML_key)
            ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);

            if strcmpi(Opt.gpu, 'on')
                fprintf(' ML (%d,%d) [GPU] ...', Stab.Mvec(j), Stab.Nvec(j));
                Opt.LU_cache.(ML_key) = gpuArray(ML);
            else
                switch Opt.lu_mode
                    case 'full'
                        fprintf(' LU (%d,%d) [cache all] ...', Stab.Mvec(j), Stab.Nvec(j));
                        tic_lu = tic;
                        Opt.LU_cache.(ML_key) = do_lu(ML);
                        fprintf(' (%.1f s, L+U: %.0f MB)', toc(tic_lu), ...
                            (nnz(Opt.LU_cache.(ML_key).L)+nnz(Opt.LU_cache.(ML_key).U))*16/1e6);

                    case 'auto'
                        if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end
                        if ~isfield(Opt, 'lu_order');     Opt.lu_order = {}; Opt.lu_order_label = {}; end
                        n_cached = numel(fieldnames(Opt.LU_cache));
                        if n_cached >= Opt.lu_max_cache
                            % Evict the oldest cached key that's NOT needed this iteration.
                            evict_idx = [];
                            for ii = 1:numel(Opt.lu_order)
                                if ~any(strcmp(Opt.lu_order{ii}, needed_set))
                                    evict_idx = ii;
                                    break;
                                end
                            end
                            if ~isempty(evict_idx)
                                old_key  = Opt.lu_order{evict_idx};
                                old_mode = Opt.lu_order_label{evict_idx};
                                Opt.LU_cache = rmfield(Opt.LU_cache, old_key);
                                Opt.lu_order(evict_idx) = [];
                                Opt.lu_order_label(evict_idx) = [];
                                fprintf(' LU (%d,%d) [replacing %s, %d cached] ...', Stab.Mvec(j), Stab.Nvec(j), old_mode, Opt.lu_max_cache);
                            else
                                % Cache saturated with needed keys — must grow to avoid
                                % thrash. lu_max_cache is a soft limit in this case.
                                if ~isfield(Opt,'lu_grew_warned') || ~Opt.lu_grew_warned
                                    warning(['LU cache grown from lu_max_cache=%d to %d ' ...
                                             'to hold all %d active modes; set lu_max_cache ' ...
                                             'accordingly or use lu_mode=''seq''/''lapack_band'' ' ...
                                             'for a bounded footprint.'], ...
                                            Opt.lu_max_cache, n_cached+1, nphys);
                                    Opt.lu_grew_warned = true;
                                end
                                fprintf(' LU (%d,%d) [%d cached, grown] ...', Stab.Mvec(j), Stab.Nvec(j), n_cached+1);
                            end
                        else
                            fprintf(' LU (%d,%d) [%d/%d cached] ...', Stab.Mvec(j), Stab.Nvec(j), n_cached+1, Opt.lu_max_cache);
                        end
                        Opt.LU_cache.(ML_key) = decomposition(ML, 'lu');
                        Opt.lu_order{end+1}       = ML_key;
                        Opt.lu_order_label{end+1} = sprintf('(%d,%d)', Stab.Mvec(j), Stab.Nvec(j));

                    case 'none'
                        fprintf(' LU (%d,%d) [direct] ...', Stab.Mvec(j), Stab.Nvec(j));
                        Opt.LU_cache.(ML_key) = ML;
                end
            end
            fprintf(' done\n');
            clear ML
        end
    end

    % Free MLA matrices once all unique modes are cached
    n_unique = numel(fieldnames(Opt.LU_cache));
    n_needed = numel(unique(arrayfun(@(j) ml_key(StabRes.omegavec(j), StabRes.betavec(j)), ...
        1:Mode.nf, 'UniformOutput', false)));
    if n_unique >= n_needed
        MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MLA5 = []; MFD_diff = [];
    end

    % Solve
    LU_cell = cell(nphys, 1);
    for k = 1:nphys
        j = phys_modes(k);
        LU_cell{k} = Opt.LU_cache.(ml_key(StabRes.omegavec(j), StabRes.betavec(j)));
    end

    if strcmpi(Opt.parfor, 'on')
        fprintf(' Solving %d modes in parallel ...\n', nphys);
        tic_solve = tic;
        parfor k = 1:nphys
            j   = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';   %#ok<PFBNS>
        end
        fprintf(' Solve done in %.1f s\n', toc(tic_solve));
    elseif strcmpi(Opt.gpu, 'on')
        fprintf(' Solving %d modes on GPU ...\n', nphys);
        tic_solve = tic;
        for k = 1:nphys
            j   = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';
        end
        fprintf(' Solve done in %.1f s\n', toc(tic_solve));
    else
        for k = 1:nphys
            j = phys_modes(k);
            fprintf(' Solving (%d,%d)', Mvec(j), Nvec(j));
            tic_j = tic;
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';
            fprintf(' ... %.1fs\n', toc(tic_j));
        end
    end
end


%% Collect results (including conjugate modes)
for k = 1:nphys
    j = phys_modes(k);
    StabRes.q(j,:)              = q_out(k,:);
    StabRes.q(Mode.nf+1-j,:)   = conj(q_out(k,:));
end

end

%% Local helpers

function key = ml_key(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function x = lu_solve(lu_s, rhs)
%LU_SOLVE  Solve using cached LU.
if isstruct(lu_s) && isfield(lu_s, 'q')
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

function s = do_lu(ML)
%DO_LU  Factorise in a separate function scope.
[L, U, p, q] = lu(ML, 'vector');
s = struct('L', L, 'U', U, 'p', p, 'q', q);
end

function ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
% MLA5 = iβ·∂_x cross-viscous block (compressible only; zero sparse otherwise).
if omega == 0 && beta == 0
    ML = MLA1 + MFD_diff;
elseif omega == 0
    ML = MLA1 + MLA3*beta + MLA4*beta^2 + MLA5*beta;
elseif beta == 0
    ML = MLA1 + MLA2*omega;
else
    ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2 + MLA5*beta;
end
end
