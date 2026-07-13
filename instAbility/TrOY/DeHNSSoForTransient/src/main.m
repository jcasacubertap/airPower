function [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)
%MAIN  DeHNSSo main solver loop.
%
%   [StabRes, StabGrid]      = main(Stab, StabGrid, Opt)
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt)       % also return LHS
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)  % reuse LHS (Opt.rerun='on')
%
%   Solves the nonlinear Harmonic Navier-Stokes equations iteratively:
%   1. Initialise grid, operators, inflow, and boundary conditions
%   2. Convergence loop: solve modes, update amplitudes, check convergence
%   3. Post-process: normalize shape functions, compute wavenumbers
%
%   LHS struct fields: MLA1, MLA2, MFD_diff — base LHS matrices that depend
%   only on the base flow, not on omega/beta. Pass back in with Opt.rerun='on'
%   to skip matrix assembly (useful for omega/beta/amplitude sweeps).

%% Input validation
if ~exist('StabGrid','var') || ~isstruct(StabGrid); StabGrid = struct(); end
if ~exist('Opt','var')      || ~isstruct(Opt);      Opt = struct();      end
if ~exist('Stab','var')     || ~isstruct(Stab);     Stab = struct();     end
if ~exist('LHS','var');     LHS = struct(); end
if ~isfield(Opt,'rerun');   Opt.rerun = 'off'; end

fprintf('\n');
fprintf(' ______     _   _  _   _  _____ _____       \n');
fprintf(' |  _  \\   | | | || \\ | |/  ___/  ___|      \n');
fprintf(' | | | |___| |_| ||  \\| |\\ `--.\\ `--.  ___  \n');
fprintf(' | | | / _ \\  _  || . ` | `--. \\`--. \\/ _ \\ \n');
fprintf(' | |/ /  __/ | | || |\\  |/\\__/ /\\__/ / (_) |\n');
fprintf(' |___/ \\___\\_| |_/\\_| \\_/\\____/\\____/ \\___/ \n');
fprintf('                                       v1.2  \n');
fprintf('  Delft Harmonic Navier-Stokes Solver\n');
fprintf('  =========================================\n');

%% Extract base-flow struct from StabGrid
fprintf('Extracting base flow ... ');
BF.Re  = StabGrid.Re;
BF.U   = StabGrid.U;    BF.V   = StabGrid.V;    BF.W   = StabGrid.W;
BF.dxU = StabGrid.dxU;  BF.dxV = StabGrid.dxV;  BF.dxW = StabGrid.dxW;
BF.dyU = StabGrid.dyU;  BF.dyV = StabGrid.dyV;  BF.dyW = StabGrid.dyW;
fprintf('done (9 fields, Re = %.2f)\n', BF.Re);

%% 1. Initialisation
[StabGrid, Init, Mode, Opt, Stab, StabRes] = init(StabGrid, Opt, Stab);
[Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes, Opt);
[Bound] = boundary(StabGrid, Opt, Stab);

%% Print configuration (after init sets defaults)
print_config(StabGrid, Stab, Opt);

%% 2. Convergence loop
% Reuse LHS matrices from a previous run (Opt.rerun='on') or build fresh.
if strcmpi(Opt.rerun,'on') && isfield(LHS,'MLA1') && ~isempty(LHS.MLA1)
    MLA1     = LHS.MLA1;
    MLA2     = LHS.MLA2;
    if isfield(LHS, 'MLA3'); MLA3 = LHS.MLA3; else; MLA3 = []; end
    if isfield(LHS, 'MLA4'); MLA4 = LHS.MLA4; else; MLA4 = []; end
    MFD_diff = LHS.MFD_diff;
    Opt.mmry = 1;   % signal solve_modes to skip LHS construction
    % LU cache persists in Opt between runs (not in LHS, to avoid duplication).
    fprintf(' Reusing LHS matrices from previous run.\n');
else
    MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MFD_diff = [];
end

is_linear = isfield(Opt,'linear') && strcmpi(Opt.linear,'on');
is_newton = isfield(Opt,'solver') && strcmpi(Opt.solver, 'newton');

if is_newton
    %% Newton solve (JFNK — see solve_newton.m)
    fprintf('\n--- Newton solve (mode: %s) ---\n', Opt.solver);
    t_loop = tic;
    [StabRes, Mode, info] = solve_newton(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes);
    fprintf('--- Newton %s after %d iterations (%.1f s) ---\n\n', ...
            info.status, info.iters, toc(t_loop));
    MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MFD_diff = [];
elseif is_linear
    %% Linear solve (single iteration, no NLT, no amplitude ramping)
    fprintf('\n--- Linear solve ---\n');
    t_loop = tic;
    Opt.dal = zeros(Mode.nf, 1);
    Opt.AF  = 1;

    [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MFD_diff] = solve_modes( ...
        BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MFD_diff);
    [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);

    if strcmpi(Opt.plot,'on')
        plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt);
    end

    fprintf('--- Linear solve done (%.1f s) ---\n\n', toc(t_loop));
else
    %% Nonlinear convergence loop
    % AF-continuation Picard with damping, optional Anderson acceleration, and JFNK Newton fallback
    fprintf('\n--- Convergence loop ---\n');
    iter = 0;
    t_loop = tic;
    picard_stall_count = 0;
    AF_stall_check = -1;
    use_newton = false;
    Opt.AF_history = [];
    iters_at_AF1 = 0;          % Picard iters spent at AF=1
    max_A_prev   = 0;          % previous iter's max_A — for runaway detection
    stall_msg_printed = false; % only print stall-disabled message once

    % Anderson acceleration state (Walker-Ni 2011). Set Opt.anderson = m to
    % enable (history depth, typical 3-6); 0 disables.
    if ~isfield(Opt,'anderson'); Opt.anderson = 0; end
    AA_on     = Opt.anderson > 0;d
    AA_m      = Opt.anderson;
    AA_Xhist  = [];
    AA_Fhist  = [];
    AA_AF_ref = [];   % AF at which history was recorded — reset on AF change

    while max(abs(Opt.dal)) >= Opt.Conv || Opt.AF ~= 1
        iter = iter + 1;
        Opt.dal = zeros(Mode.nf, 1);
        t_iter = tic;

        if ~use_newton
            % === Damped Picard iteration ===
            q_old = StabRes.q;

            % Full Picard step: q_picard = ML⁻¹ * (rBC + f_NLT)
            [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MFD_diff] = solve_modes( ...
                BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MFD_diff);

            % AVAL rescales the first iterate down to AMAX based on the
            % linear-regime prediction. Skip it when LOAD has warm-started
            % the full domain — the warm-start is already at physical saturation
            % and scaling it down throws away the continuation benefit.
            loaded_full_domain = strcmpi(Stab.IC,'LOAD');
            if nnz(Stab.A0) ~= 1 && Opt.k1 == 0 && ~loaded_full_domain
                [Opt, StabRes] = AVAL(StabGrid, Opt, Mode, StabRes);
            elseif loaded_full_domain && Opt.k1 == 0
                % With warm-start, default AF=1 from init.m stands; no damping.
                Opt.AVAL = max(abs(StabRes.q(:)));
            end

            % --- Anderson acceleration (optional) ------
            % Least-squares fixed-point accelerator over last m iterates.
            % Falls back to plain Picard if the accelerated step is worse.
            if AA_on && Opt.k1 > 0
                q_picard  = StabRes.q;                  % g(q_old)
                f_k_vec   = q_picard(:) - q_old(:);     % residual of fixed-point

                % Reset history if AF changed since last iterate
                if ~isempty(AA_AF_ref) && AA_AF_ref ~= Opt.AF
                    AA_Xhist = []; AA_Fhist = [];
                end

                % Try Anderson step if we have at least 1 previous iterate
                if ~isempty(AA_Xhist)
                    dX = [AA_Xhist(:, 2:end), q_old(:)]   - [AA_Xhist(:, 1:end)];
                    dF = [AA_Fhist(:, 2:end), f_k_vec]    - [AA_Fhist(:, 1:end)];
                    % Tikhonov-regularized LS: (dFᴴdF + λI)γ = dFᴴ f_k.
                    % λ scales with ||dF||²_F so it's dimensionally consistent
                    % and degrades gracefully to plain Picard when dF is
                    % rank-deficient (parallel residuals near saturation).
                    lambda = 1e-10 * norm(dF, 'fro')^2;
                    gamma  = (dF'*dF + lambda * eye(size(dF,2))) \ (dF' * f_k_vec);
                    q_AA_vec = q_picard(:) - (dX + dF) * gamma;

                    % Safety: use AA only if result is finite AND not a
                    % bigger blowup than plain Picard. Otherwise keep q_picard.
                    max_AA = max(abs(q_AA_vec));
                    max_P  = max(abs(q_picard(:)));
                    if all(isfinite(q_AA_vec)) && max_AA < 1.5 * max_P
                        StabRes.q = reshape(q_AA_vec, size(q_picard));
                        fprintf('   [Anderson: m=%d, γ=(%s)]\n', ...
                            size(dX,2), sprintf('%.2f ', gamma));
                    end
                end

                % Push new entry, drop oldest if full
                AA_Xhist(:, end+1) = q_old(:);        %#ok<AGROW>
                AA_Fhist(:, end+1) = f_k_vec;         %#ok<AGROW>
                if size(AA_Xhist, 2) > AA_m
                    AA_Xhist(:, 1) = [];
                    AA_Fhist(:, 1) = [];
                end
                AA_AF_ref = Opt.AF;
            end
            % ----------------------------------------------------------------

            % Damp if amplitudes grow too fast (prevents blow-up at high AF).
            % Applied AFTER Anderson, as a final safety net.
            if Opt.k1 > 0
                q_picard = StabRes.q;
                dq = q_picard - q_old;
                max_A_old = max(max(abs(q_old(Mode.RunJ, :))));
                max_A_new = max(max(abs(q_picard(Mode.RunJ, :))));

                if max_A_new > 2 * max(max_A_old, 1e-6)
                    alpha = max_A_old / max(max_A_new, eps);
                    alpha = max(alpha, 0.1);
                    StabRes.q = q_old + alpha * dq;
                    fprintf('   [damped Picard: alpha=%.2f]\n', alpha);
                end
            end
        else
            % === Newton iteration: converge current AF, then ramp ===
            fprintf('  [Newton] Converging AF=%.3f ...\n', Opt.AF);
            [StabRes, Mode, ninfo] = solve_newton( ...
                BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes);
            fprintf('  [Newton] %s in %d iters (%.1f s)\n', ...
                    ninfo.status, ninfo.iters, toc(t_iter));

            % If Newton stalled (‖F‖ unchanged), don't re-invoke with the same
            % broken state — bail out with an explanation. This usually means
            % the warm-start state is too far from the physical solution
            % (grid under-discretized for this amplitude regime).
            if strcmpi(ninfo.status, 'stalled')
                fprintf(['\n  ** Newton cannot progress from this state ' ...
                         '(likely grid under-discretized for max_A=%.2f). ' ...
                         'Reduce Stab.A0_fund or refine the grid. **\n\n'], ...
                         max(StabRes.A(:)));
                break;
            end
        end

        % Update modes, amplitudes, NLT forcing
        [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);

        % After Newton: reset Aold so check_convergence ramps AF
        if use_newton
            StabRes.Aold = StabRes.A;
        end

        % Check convergence and ramp AF
        AF_before = Opt.AF;
        [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes);
        Opt.AF_history(end+1) = Opt.AF;

        % Status
        n_active = numel(Mode.RunJ(Mode.RunJ >= round(Mode.nf/2)));
        max_A = max(StabRes.A(:));
        max_dal = max(abs(Opt.dal));
        solver_tag = 'P';
        if use_newton; solver_tag = 'N'; end
        fprintf('  Iter %3d [%s] | AF=%.3f | %d modes | max(A)=%.3e | max(dal)=%.3e | %.1f s\n', ...
            iter, solver_tag, Opt.AF, n_active, max_A, max_dal, toc(t_iter));

        % --- Picard stall detection ---
        if ~use_newton
            % Reset counter when AF changes or number of active modes changes
            if ~isfield(Opt, 'n_active_prev'); Opt.n_active_prev = n_active; end
            modes_changed = (n_active ~= Opt.n_active_prev);
            Opt.n_active_prev = n_active;

            if modes_changed
                picard_stall_count = 0;
                picard_af_count = 0;
                AF_stall_check = Opt.AF;
                dal_prev = inf;
            elseif Opt.AF ~= AF_before
                % AF just ramped — reset
                picard_stall_count = 0;
                picard_af_count = 0;
                AF_stall_check = Opt.AF;
                dal_prev = inf;
            elseif max_dal > Opt.Conv
                if ~exist('dal_prev','var'); dal_prev = inf; end
                if ~exist('picard_af_count','var'); picard_af_count = 0; end
                picard_af_count = picard_af_count + 1;  % total iters at this AF
                if max_dal >= dal_prev
                    picard_stall_count = picard_stall_count + 1;
                else
                    picard_stall_count = 0;
                end
                dal_prev = max_dal;
            else
                picard_stall_count = 0;
                picard_af_count = 0;
                dal_prev = max_dal;
            end

            % Track time at AF=1 and detect runaway early, regardless of dal
            % (which gets masked to 0 when AVAL >= 2). These catch Picard
            % *before* it loses the physical saturation state.
            if Opt.AF >= 1 - eps
                iters_at_AF1 = iters_at_AF1 + 1;
            else
                iters_at_AF1 = 0;
            end
            runaway_jump = (iters_at_AF1 > 1) && (max_A > 1.5 * max_A_prev) && (max_A > 0.15);
            max_A_prev = max_A;

            % Trigger Newton if:
            %  - dal not decreasing for 5 iters (true stall/divergence), OR
            %  - too many iters at same AF (Picard too slow), OR
            %  - amplitudes blowing up (stall_count + max_A cross threshold), OR
            %  - Picard has been at AF=1 for 5 iters without converging (done
            %    its AF-ramp job, Newton should polish), OR
            %  - runaway signature: max_A jumps >50% in one iter at AF=1, OR
            %  - amplitudes unphysically large (max_A > 1 or NaN/Inf).
            if picard_stall_count >= 5 || picard_af_count >= 15 || ...
               (picard_stall_count >= 3 && max_A > 0.1) || ...
               (iters_at_AF1 >= 5 && max_dal > Opt.Conv) || ...
               runaway_jump || ...
               max_A > 1 || ~isfinite(max_A)
                picard_stall_count = 10;
            end

            if picard_stall_count >= 10
                if isfield(Opt,'newton_takeover') && strcmpi(Opt.newton_takeover,'off')
                    % User has disabled Newton takeover — keep running Picard.
                    % Reset the stall counter so we don't trip every iteration.
                    picard_stall_count = 0;
                    if ~stall_msg_printed
                        fprintf(['\n  ** Picard stalled at AF=%.3f ' ...
                                 '(Newton takeover disabled, continuing Picard) **\n\n'], Opt.AF);
                        stall_msg_printed = true;
                    end
                else
                    fprintf('\n  ** Picard stalled at AF=%.3f, switching to Newton permanently **\n\n', Opt.AF);
                    use_newton = true;
                    % Reduce AF growth for Newton (per-mode LU needs smaller steps at high AF)
                    Opt.AFg = 1.02;
                    % Free Picard's LU cache and LHS pieces — Newton builds its own
                    % preconditioner and doesn't reuse these, so holding both blows up RAM.
                    if isfield(Opt,'LU_cache')
                        Opt.LU_cache = struct();
                        Opt.lu_order = {}; Opt.lu_order_label = {};
                    end
                    MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MFD_diff = [];
                    Opt.mmry = 0;   % force LHS rebuild on the next call if needed
                end
            end
        end

        % Plotting
        if strcmpi(Opt.plot,'on')
            plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt);
        end
        if Opt.Sweep == 1 && max(abs(Opt.dal)) <= Opt.Sweep * Opt.Conv
            [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes);
        end
        if strcmpi(Opt.plot,'on')
            plot_convergence(StabGrid, Mode, Stab, Opt);
            % Newton convergence subplot
            if use_newton && exist('ninfo','var') && isfield(ninfo,'mode_F_history')
                % Accumulate Newton history across AF levels
                n_new = size(ninfo.mode_F_history, 2);
                if ~isfield(Opt, 'newton_F_history')
                    Opt.newton_F_history = ninfo.mode_F_history;
                    Opt.newton_AF_history = Opt.AF * ones(1, n_new);
                else
                    Opt.newton_F_history = [Opt.newton_F_history, ninfo.mode_F_history];
                    Opt.newton_AF_history = [Opt.newton_AF_history, Opt.AF * ones(1, n_new)];
                end
                plot_convergence(StabGrid, Mode, Stab, Opt);
            end
        end
    end

    % Report. Distinguish true convergence from the Newton-stall bail-out.
    converged = (max(abs(Opt.dal)) < Opt.Conv) && (Opt.AF == 1);
    if converged; label = 'Converged'; else; label = 'Stopped (did NOT converge)'; end
    try
        [~, rss] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
        fprintf('--- %s in %d iterations (%.1f s, peak %.1f GB) ---\n\n', ...
            label, iter, toc(t_loop), str2double(strtrim(rss))/1e6);
    catch
        fprintf('--- %s in %d iterations (%.1f s) ---\n\n', label, iter, toc(t_loop));
    end
end

%% 3. Save checkpoint (optional, before postprocess)
% Dumps StabRes with full mode set + .q intact, suitable for Stab.IC='LOAD'
% in a subsequent run (amplitude continuation, parameter sweep). Triggered by
% Opt.save_checkpoint = 'path/to/file.mat'. Only saves if solution is finite.
if isfield(Opt,'save_checkpoint') && ~isempty(Opt.save_checkpoint) ...
   && all(isfinite(StabRes.q(:)))
    ckpt_path = Opt.save_checkpoint;
    ckpt_dir  = fileparts(ckpt_path);
    if ~isempty(ckpt_dir) && ~exist(ckpt_dir,'dir'); mkdir(ckpt_dir); end
    save(ckpt_path, 'StabRes', 'StabGrid', 'Stab', 'Mode');
    fprintf('--- Checkpoint saved: %s ---\n', ckpt_path);
end

%% 4. Post-processing
half = round(Mode.nf/2);
% Skip postprocess if the solve diverged (StabRes contains NaN/Inf). Running
% interp1 / spline fits on NaN would crash. Just return the raw (broken)
% StabRes so the user can inspect.
if any(~isfinite(StabRes.q(:)))
    warning('main:postprocess_skipped', ...
        'Solve diverged (non-finite values in StabRes.q) — skipping postprocess.');
else
    StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes);
end

%% Pack LHS output (for reuse in subsequent runs)
LHS.MLA1     = MLA1;
LHS.MLA2     = MLA2;
LHS.MLA3     = MLA3;
LHS.MLA4     = MLA4;
LHS.MFD_diff = MFD_diff;
% Note: LU_cache is NOT copied into LHS to avoid doubling memory.
% For rerun mode, the cache persists in Opt (returned to workspace).

end
