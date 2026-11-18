function varargout = mumps_solve(op, varargin)
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
% MUMPS back-end for the Newton solver (factor / bsolve / destroy).
%
%   id      = mumps_solve('factor',  Jr, Opt, step, lam)
%   [X, id] = mumps_solve('bsolve',  id, Jr, RHS)
%   id      = mumps_solve('destroy', id)
%
%   Single entry point (as block_thomas); the ooc_* disk handling is internal.
%   JOB=4 analyze+factor, JOB=3 solve, JOB=-2 destroy. The factor is held across
%   steps, so id/Jr are threaded through the caller rather than captured in a
%   closure: a closure holding a destroyed id would segfault on the next solve.

%% Inputs
% Name     size unit   explanation
% op       char [-]    operation to perform: 'factor', 'bsolve' or 'destroy'
% varargin cell [-]    arguments forwarded to the selected operation

%% Outputs
% Name      size unit   explanation
% varargout cell [-]    values returned by the selected operation

switch lower(op)
    case 'factor';  [varargout{1:nargout}] = mumps_factor(varargin{:});
    case 'bsolve';  [varargout{1:nargout}] = mumps_bsolve(varargin{:});
    case 'destroy'; [varargout{1:nargout}] = mumps_destroy(varargin{:});
    otherwise
        error('mumps_solve:badOp', ...
              'Unknown op ''%s'': use factor, bsolve or destroy.', op);
end
end

% get_o and tern are one-liners kept here too, so this file stands alone.
function v = get_o(Opt, name, dflt)
if isfield(Opt, name) && ~isempty(Opt.(name)); v = Opt.(name); else; v = dflt; end
end
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function id = mumps_factor(Jr, Opt, step, lam)
%MUMPS_FACTOR  Factorise the real 2N x 2N Jacobian, retrying before going to disk.
%   step and lam only label the statistics record; they default to 0/NaN for
%   callers with no branch position, such as the natural-parameter seed solves.
% A NaN/Inf in Jr (state overshooting into nonphysical rho'/T') crashes dmumps
% hard. Check the stored nonzeros only: ~isfinite on the full sparse densifies it.
if nargin < 3 || isempty(step); step = 0;   end
if nargin < 4 || isempty(lam);  lam  = NaN; end
nbad = nnz(~isfinite(nonzeros(Jr)));
if nbad > 0
    error('mumps_factor:badJacobian', ...
          ['Jacobian has %d non-finite entries before factorisation: the state ' ...
           'is nonphysical (likely a cold-start overshoot). Lower Opt.newton_start ' ...
           'or newton_dsmax.'], nbad);
end
% Attempt schedule, columns = [ICNTL(14) relax, ICNTL(22) out-of-core]. In-core
% first, OOC last, since only steps near saturation need the disk. BLR and
% preconditioned Krylov did not survive measurement (project_m4_memory_wall).
sched = [100 0; 200 0; 100 1];
if isfield(Opt,'newton_mumps_relax') && ~isempty(Opt.newton_mumps_relax)
    r = Opt.newton_mumps_relax(:);
    sched = [[r, zeros(numel(r),1)]; [r(1), 1]];   % keep an OOC last resort
end
if isfield(Opt,'newton_mumps_ooc') && ~isempty(Opt.newton_mumps_ooc)
    sched(:,2) = Opt.newton_mumps_ooc;                    % force OOC on/off throughout
end
info1 = NaN;
fprintf('   [mumps_factor] Jr: %dx%d nnz=%.3e square=%d sympat=%d maxabs=%.3e\n', ...
        size(Jr,1), size(Jr,2), nnz(Jr), size(Jr,1)==size(Jr,2), ...
        nnz(Jr-Jr.')==0, full(max(abs(nonzeros(Jr))))+0);
for attempt = 1:size(sched,1)
    relax = sched(attempt,1); ooc = sched(attempt,2);
    crashed = false; info2 = NaN; t_fac = NaN;
    try
        id = initmumps; id.SYM = 0; id = dmumps(id);
        % ICNTL(7)=0 (AMD) by default, not METIS: the saturated M=4 Jacobian fails
        % under METIS with INFOG(1)=-13. Never size this solver from a seed probe,
        % which is far sparser. ICNTL(7)=4 (PORD) segfaults this build.
        vrb = get_o(Opt,'newton_mumps_verbose',0);   % 0=silent; >0 -> MUMPS diagnostics
        id.ICNTL(3) = tern(vrb>0,6,0); id.ICNTL(4) = vrb; id.ICNTL(16) = 1;
        id.ICNTL(7) = get_o(Opt, 'newton_mumps_ordering', 0);
        id.ICNTL(14) = relax;
        id.ICNTL(22) = ooc;
        % ICNTL(23): hard cap on working memory, MB per process. Needed because
        % ICNTL(22)=1 streams the factors but keeps the frontal stack in core, so
        % OOC alone does not bound the footprint. 0 = MUMPS sizes it itself.
        mxm = get_o(Opt, 'newton_mumps_maxmem', 0);
        if mxm > 0, id.ICNTL(23) = mxm; end
        % Block low-rank compression: ICNTL(35) picks the variant (1 throughout,
        % 2 full-rank solve, 3 compress after), CNTL(7) the drop tolerance. The
        % factor is then APPROXIMATE - validate against a BLR-off run first.
        blr = get_o(Opt, 'newton_mumps_blr', 0);
        if blr > 0
            id.ICNTL(35) = blr;
            id.CNTL(7)   = get_o(Opt, 'newton_mumps_blr_tol', 1e-9);
        end
        % JOB=-2 closes the OOC files but does not unlink them, and a killed process
        % leaves everything, so each factorisation gets its own MUMPS_OOC_PREFIX.
        if ooc
            ooc_pfx = sprintf('dehnsso_%d_%s', feature('getpid'), ...
                              char('a' + floor(26*rand(1,8))));
            setenv('MUMPS_OOC_PREFIX', ooc_pfx);
            if isempty(getenv('MUMPS_OOC_TMPDIR')); setenv('MUMPS_OOC_TMPDIR', tempdir); end
            ooc_check_space(getenv('MUMPS_OOC_TMPDIR'));
            ooc_track(ooc_pfx, 'push');
        end
        % Analysis only, when asked: JOB=1 returns the estimates without
        % building or writing anything, so the memory a case would need can be
        % checked in minutes before committing hours and hundreds of GB of
        % scratch to it. INFOG(16)/(17) are the in-core estimates and
        % INFOG(19)/(20) the out-of-core ones; the entry counts INFOG(9)/(29)
        % are only filled by a real factorisation and stay zero here.
        if isfield(Opt,'newton_analyze_only') && Opt.newton_analyze_only
            t_fac = tic; id.JOB = 1; id = dmumps(id, Jr); t_fac = toc(t_fac);
            fprintf(['   [mumps_analyze] STATS nnzJ=%d n=%d t_analyze=%.1f ' ...
                     'INFOG16=%g INFOG17=%g INFOG19=%g INFOG20=%g ' ...
                     'INFOG26=%g INFOG27=%g\n'], ...
                    nnz(Jr), size(Jr,1), t_fac, id.INFOG(16), id.INFOG(17), ...
                    id.INFOG(19), id.INFOG(20), id.INFOG(26), id.INFOG(27));
            % Not an error: throwing here would be caught by the retry loop
            % below, which reads any exception as a crash and starts over with
            % more relaxation. Close MUMPS and hand back an empty handle, which
            % the caller treats as "no factor".
            id.JOB = -2; dmumps(id); id = [];
            if ooc, ooc_track(ooc_pfx, 'pop'); end
            return
        end
        t_fac = tic; id.JOB = 4; id = dmumps(id, Jr); t_fac = toc(t_fac);
    catch ME
        crashed = true;
        fprintf('   [mumps_factor] exception: %s | %s\n', ME.identifier, ME.message);
    end
    if ~crashed && isstruct(id) && ~isempty(id.INFOG)
        info1 = id.INFOG(1); info2 = id.INFOG(2);
        if info1 >= 0
            if attempt > 1
                fprintf('   [mumps_factor] OK on attempt %d (relax=%d ooc=%d)\n', attempt, relax, ooc);
            end
            % What the factorisation actually cost, not what a failed attempt
            % asked for. MUMPS fills these after JOB=4:
            %   INFOG(9)/(10) real and integer entries in the retained factor
            %   INFOG(16)     in-core estimate from the analysis (what the
            %                 failed attempts above report)
            %   INFOG(21)     memory used during factorisation, MB
            %   INFOG(29)     factor entries after factorisation
            %   RINFOG(16)    disk written for the OOC factor, MB
            % Negative INFOG(9)/(10)/(29) counts millions of entries (MUMPS
            % convention) and is left as returned. Bytes per entry are build
            % properties: real double here, integer width from Opt (8 with
            % -DINTSIZE64, 4 in a default build), NaN when unset.
            rec = sprintf(['   [mumps_factor] STATS step=%d urms=%.6e nnzJ=%d ' ...
                     't_factor=%.1f relax=%d ooc=%d INFOG9=%g INFOG10=%g ' ...
                     'INFOG16=%g INFOG21=%g INFOG29=%g RINFOG16=%g ' ...
                     'scalar_bytes=%d int_bytes=%g\n'], ...
                    step, lam/sqrt(2), nnz(Jr), t_fac, relax, ooc, ...
                    id.INFOG(9), id.INFOG(10), id.INFOG(16), id.INFOG(21), ...
                    id.INFOG(29), id.RINFOG(16), 8, ...
                    get_o(Opt,'newton_mumps_int_bytes',NaN));
            fprintf('%s', rec);
            % A transcript can be lost or truncated; the records are the point of
            % the run, so they also go to a file when one is named.
            sfile = get_o(Opt,'newton_stats_file','');
            if ~isempty(sfile)
                fid = fopen(sfile,'a');
                if fid > 0; fprintf(fid,'%s',rec); fclose(fid); end
            end
            return
        end
        fprintf('   [mumps_factor] attempt %d relax=%d ooc=%d FAILED: INFOG(1)=%d INFOG(2)=%d\n', ...
                attempt, relax, ooc, info1, info2);
        try, id.JOB = -2; id = dmumps(id); catch, end
        ooc_purge();   % a failed attempt leaks its scratch otherwise
    else
        fprintf('   [mumps_factor] attempt %d relax=%d ooc=%d HARD CRASH\n', attempt, relax, ooc);
        try, if isstruct(id), id.JOB = -2; id = dmumps(id); end; catch, end
        ooc_purge();   % a failed attempt leaks its scratch otherwise
    end
    pause(20);
end
error('mumps_factor: failed after %d attempts (last INFOG(1)=%g)', size(sched,1), info1);
end

function [X, id] = mumps_bsolve(id, Jr, RHS, Opt)
%MUMPS_BSOLVE  Apply the retained factor to one right-hand side.
%   The corrector pays this once per iteration while the factorisation is paid
%   once per step, so the two together are the cost of a continuation step. The
%   time is recorded when Opt.newton_stats_file is set, in the same file as the
%   factorisation statistics, so a run reports both halves of that cost.
t_bs = tic; id.JOB = 3; id.RHS = RHS; id = dmumps(id, Jr); X = double(id.SOL);
t_bs = toc(t_bs);
if nargin >= 4
    sfile = get_o(Opt,'newton_stats_file','');
    if ~isempty(sfile)
        fid = fopen(sfile,'a');
        if fid > 0
            fprintf(fid,'   [mumps_bsolve] STATS t_solve=%.2f nrhs=%d\n', t_bs, size(RHS,2));
            fclose(fid);
        end
    end
end
end

function id = mumps_destroy(id)
%MUMPS_DESTROY  Release the instance and its out-of-core scratch.
%   JOB=-2 closes the OOC files but does not unlink them, so the scratch would
%   accumulate. Deletes only files carrying THIS instance's prefix.
id.JOB = -2; id = dmumps(id);
ooc_purge();
end

function ooc_check_space(d)
%OOC_CHECK_SPACE  Refuse to start an out-of-core factorisation on a nearly full disk.
%   The scratch defaults to tempdir, usually the system volume, and a large
%   factorisation writes tens of GB. Point MUMPS_OOC_TMPDIR at a roomy filesystem.
MIN_FREE_GB = 50;
try
    if ispc
        return                                   % no cheap portable check; skip
    end
    [st, out] = system(sprintf('df -Pk "%s" 2>/dev/null | tail -1', d));
    if st ~= 0 || isempty(out); return; end
    tok = regexp(strtrim(out), '\s+', 'split');
    if numel(tok) < 4; return; end
    free_gb = str2double(tok{4}) / 2^20;         % df -Pk reports 1K blocks
    if isnan(free_gb); return; end
    if free_gb < MIN_FREE_GB
        error('solve_newton:oocDiskSpace', ...
            ['Out-of-core scratch directory "%s" has only %.1f GB free (want >= %d GB). '  ...
             'A factorisation of this size writes tens of GB and filling this volume can ' ...
             'wedge the machine. Set MUMPS_OOC_TMPDIR to a filesystem with room, or free ' ...
             'space, before rerunning.'], d, free_gb, MIN_FREE_GB);
    end
    fprintf('   [mumps_factor] OOC scratch: %s (%.0f GB free)\n', d, free_gb);
catch ME
    if strcmp(ME.identifier, 'solve_newton:oocDiskSpace'); rethrow(ME); end
    % any probing failure is non-fatal: never block a run over a df quirk
end
end

function ooc_purge()
%OOC_PURGE  Delete the out-of-core scratch of the most recent OOC factorisation.
%   Pops this instance's prefix and unlinks only files carrying it, so a
%   concurrently live factorisation (natnewton holds its own) is untouched.
pfx = ooc_track('', 'pop');
if isempty(pfx); return; end
d = getenv('MUMPS_OOC_TMPDIR'); if isempty(d); d = tempdir; end
f = dir(fullfile(d, [pfx '*']));
for k = 1:numel(f)
    if ~f(k).isdir
        try
            delete(fullfile(d, f(k).name));
        catch
        end
    end
end
end

function out = ooc_track(pfx, op)
%OOC_TRACK  LIFO registry of OOC prefixes, so destroy removes the right set.
persistent stack
if isempty(stack); stack = {}; end
out = '';
switch op
    case 'push'; stack{end+1} = pfx;
    case 'pop'
        if ~isempty(stack); out = stack{end}; stack(end) = []; end
end
end
