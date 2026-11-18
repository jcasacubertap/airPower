function varargout = block_thomas(op, varargin)
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
% Dual-Schur block-pentadiagonal Thomas solver (Paredes et al. 2024).
%   Direct solver for the station-major block-pentadiagonal HNS operator ML,
%   block size b = Nq*Ny. The domain is split into np partitions separated by
%   two-station interfaces; each is factored by an unpivoted block-Thomas
%   recursion, and the partitions are coupled through a Schur complement on the
%   interfaces. The last partition is swept in reverse, a UL rather than an LU
%   decomposition, so both extreme partitions carry the same reduced cost.
%
%   The recursion is unpivoted, so the diagonal station blocks G_i must be
%   nonsingular. The stored factor is three dense blocks per station: Y_i, Z_i
%   and, depending on AMODE, either (G_i)^-1 or the pivoted LU factors of G_i.
%   Default np = 2: beyond two partitions the interior subdomains raise the
%   operation count and data movement and the scaling degrades.
%
%   Operations:
%     x  = block_thomas('solve',  ML, b, rhs, np, oocdir, par, amode)  factor+apply+free
%     F  = block_thomas('factor', ML, b, np, oocdir, par, amode)       reusable factor
%     x  = block_thomas('apply',  F, rhs)                              solve ML*x = rhs
%          block_thomas('free',   F)                                   delete F's disk files
%          block_thomas('free',   OOCDIR, 'reg')                       register a dir
%     mk = block_thomas('free',   [], 'mark')                          registry checkpoint
%          block_thomas('free',   mk, 'sweep')                         delete dirs after mk
%
%   With OOCDIR the factor is streamed to disk, so only a few blocks are ever in
%   RAM and resolutions too large to factor in core stay reachable. rhs may have
%   several columns. The registry lets main() clean up streamed files on exit.

%% Inputs
% Name     size unit   explanation
% op       char [-]    operation to perform: 'solve', 'factor', 'apply' or 'free'
% varargin cell [-]    arguments forwarded to the selected operation

%% Outputs
% Name      size unit   explanation
% varargout cell [-]    values returned by the selected operation

switch lower(op)
    case 'solve';  [varargout{1:nargout}] = bt_solve(varargin{:});
    case 'factor'; [varargout{1:nargout}] = bt_factor(varargin{:});
    case 'apply';  [varargout{1:nargout}] = bt_apply(varargin{:});
    case 'free';   [varargout{1:nargout}] = bt_free(varargin{:});
    otherwise
        error('block_thomas:unknownOp', ...
            'unknown operation ''%s'' (use solve/factor/apply/free)', op);
end
end

%% One-shot solve
function x = bt_solve(ML, b, rhs, np, oocdir, par, amode)
if nargin < 4 || isempty(np); np = 2; end
if nargin < 5; oocdir = ''; end
if nargin < 6 || isempty(par); par = false; end
if nargin < 7 || isempty(amode); amode = 'lu'; end
F = bt_factor(ML, b, np, oocdir, par, amode);
x = bt_apply(F, rhs);
bt_free(F);
end

%% Factor: build a right-hand-side (RHS) independent dual-Schur factor of ML (in-core, or streamed to OOCDIR)
function F = bt_factor(ML, b, np, oocdir, par, amode)
%BT_FACTOR  Build a right-hand-side independent dual-Schur factor of ML.
%   AMODE sets what is stored for each diagonal station block G:
%     'lu' (default)      its pivoted LU factors, applied by triangular solves
%     'inverse'           G^-1, as in Paredes et al. (2024)
%   'lu' keeps the pivoting that forming the inverse discards, which matters
%   where G is close to singular; see docs/DeHNSSo_guide.md for the measured
%   difference. Both forms stream, the LU record carrying b extra doubles per
%   station for the pivots.
if nargin < 3 || isempty(np); np = 2; end
if nargin < 4; oocdir = ''; end
if nargin < 5 || isempty(par); par = false; end
if nargin < 6 || isempty(amode); amode = 'lu'; end
use_lu = strcmpi(amode, 'lu');
ooc = ~isempty(oocdir);
% Out-of-core streams whichever form AMODE selects; the LU record carries the
% packed factors and the pivot vector, b doubles more per station than G^-1.
par = par && np > 1;
if ooc && ~exist(oocdir,'dir'); mkdir(oocdir); end
n = size(ML,1); N = round(n/b);

% stations whose stencil reaches past the pentadiagonal band (near a boundary) go into the interface
[r,c] = find(ML); bi = ceil(r/b); bj = ceil(c/b); far = abs(bi-bj) > 2;
bad = unique([bi(far); bj(far)]);
istart = 1; iend = N;
if ~isempty(bad)
    lo = bad(bad <= N/2); if ~isempty(lo); istart = max(lo) + 1; end
    hi = bad(bad >  N/2); if ~isempty(hi); iend  = min(hi) - 1; end
end

% partition the interior stations; the boundary blocks plus the seams between partitions form the interface
label = zeros(N,1);
Nin = iend - istart + 1; nint = Nin - 2*(np-1); chunk = floor(nint/np); pos = istart;
for p = 1:np
    csz = chunk + double(p <= mod(nint,np));
    label(pos:pos+csz-1) = p; pos = pos + csz;
    if p < np; label(pos:pos+1) = 0; pos = pos + 2; end
end
sd = @(st) reshape(((st(:)-1)*b + (1:b)).', [], 1);
ifDof = sd(find(label==0));

% Each partition builds its own Schur contribution Sloc instead of reducing into
% St, so the loop has no cross-iteration dependence and can run in parallel.
St = ML(ifDof,ifDof);
Cp = cell(np,1); W = cell(np,1); wcol = cell(np,1); pdof = cell(np,1);
part = cell(np,1); Dp = cell(np,1); wfile = cell(np,1); wsize = zeros(np,2);
Sloc = cell(np,1); Bc = cell(np,1); ffile = cell(np,1);
% Slice on the client: the slices are returned in the factor anyway, and it stops
% a process pool broadcasting the whole operator to every worker.
for p = 1:np
    % Last partition swept in reverse (UL): its interface coupling then sits in the
    % last two station-rows, as for the first. pdof undoes the reversal on scatter.
    st = find(label==p); if p == np && np > 1; st = flip(st); end
    pd = sd(st); pdof{p} = pd;
    Dp{p} = ML(pd,pd);
    Bp = ML(pd,ifDof); Cp{p} = ML(ifDof,pd);
    cols = find(any(Bp,1)); wcol{p} = cols;
    Bc{p} = full(Bp(:,cols));
    wsize(p,:) = [numel(pd), numel(cols)];
    if ooc
        ffile{p} = fullfile(oocdir, sprintf('fac_p%03d.bin', p));
        wfile{p} = fullfile(oocdir, sprintf('w_p%03d.bin', p));
    end
end
nif = numel(ifDof);
if par
parfor p = 1:np
    [Sloc{p}, W{p}, part{p}] = bt_partfactor(Dp{p}, Bc{p}, Cp{p}, b, nif, ooc, ffile{p}, wfile{p}, use_lu);
end
else
for p = 1:np
    [Sloc{p}, W{p}, part{p}] = bt_partfactor(Dp{p}, Bc{p}, Cp{p}, b, nif, ooc, ffile{p}, wfile{p}, use_lu);
end
end
Bc = [];
for p = 1:np
    St(:,wcol{p}) = St(:,wcol{p}) - Sloc{p};
end
Sloc = [];
[StL, StU, Stp] = lu(full(St), 'vector');

F = struct('b',b, 'n',n, 'np',np, 'ifDof',ifDof, 'pdof',{pdof}, 'wcol',{wcol}, ...
           'Cp',{Cp}, 'part',{part}, 'Dp',{Dp}, ...
           'StL',StL, 'StU',StU, 'Stp',Stp, 'ooc',ooc, 'oocdir',oocdir, 'par',par);
if ooc
    F.wfile = wfile; F.wsize = wsize;
else
    F.W = W;
end
end

function varargout = local_thomas_factor(D, b, Bcols, facfile, wfile, use_lu)
%LOCAL_THOMAS_FACTOR  Factor one partition and produce its interface response W.
%   In-core (3 args): returns [W, part], part = {L_i,U_i,pv_i,Y_i,Z_i}.
%   Out-of-core (5 args): streams the station factor blocks [L,U,Y,Z] to FACFILE
%   and W block-by-block to WFILE (pv kept in RAM), holding only a few blocks in
%   memory at a time; returns part = {pv, N, facfile}. Nothing the size of a
%   whole partition is ever materialized in memory.
stream = nargin >= 4 && ~isempty(facfile);
if nargin < 6 || isempty(use_lu); use_lu = false; end
n = size(D,1); N = round(n/b); zb = sparse(b,b); nc = size(Bcols,2);
if stream
    fid = fopen(facfile,'w');
    wfwd = [facfile '.wf']; fidWf = fopen(wfwd,'w');   % scratch file for the forward W sweep, deleted after use
else
    Gc = cell(N,1); Lc = cell(N,1); Uc = cell(N,1); Pc = cell(N,1);
    Yc = cell(N,1); Zc = cell(N,1);
    Wc = cell(N,1);
end
Ym1=zb; Ym2=zb; Zm1=zb; Zm2=zb; rm1=zeros(b,nc); rm2=zeros(b,nc);
blk = @(i,j) D((i-1)*b+(1:b),(j-1)*b+(1:b));
% Bcols is nonzero only in the two station-rows next to an interface, so W is
% identically zero until station i0 and those solves (and their I/O) are skipped.
nzst = any(reshape(any(Bcols ~= 0, 2), b, N), 1);
i0 = find(nzst, 1); if isempty(i0); i0 = N + 1; end
zW = zeros(b,nc);
for i = 1:N
    if i>=3, A = blk(i,i-2); else, A = zb; end
    if i>=2, Bb = blk(i,i-1); else, Bb = zb; end
    C = blk(i,i);
    if i<=N-1, Dd = blk(i,i+1); else, Dd = zb; end
    if i<=N-2, E  = blk(i,i+2); else, E  = zb; end
    K = Bb - A*Ym2;
    G = full(C - K*Ym1 - A*Zm2);
    if use_lu
        % Pivoted factors, applied by triangular solves.
        [Lg, Ug, pg] = lu(G, 'vector');
        gsolve = @(V) Ug \ (Lg \ V(pg,:));
        Yi = gsolve(full(Dd - K*Zm1));
        Zi = gsolve(full(E));
    else
        Gi = inv(G);                                        % stored inverse (Paredes et al. 2024)
        % Gi*E is a GEMM, so a sparse E (~1% dense) costs nnz(E)*b, not b^3.
        gsolve = @(V) Gi*V;
        Yi = Gi*(Dd - K*Zm1);
        Zi = Gi*E;
    end
    if i < i0
        Wi = zW;                                            % W is still identically zero
    else
        Wi = gsolve(full(Bcols((i-1)*b+(1:b),:) - K*rm1 - A*rm2));
    end
    if stream
        % L is unit lower triangular, so L and U pack into one block without loss.
        if use_lu; Gs = triu(Ug) + tril(Lg,-1); else; Gs = Gi; end
        blkrec = [real(Gs(:));imag(Gs(:)); ...
                  real(Yi(:));imag(Yi(:));real(Zi(:));imag(Zi(:))];
        if use_lu; blkrec = [blkrec; pg(:)]; end
        bt_fwrite(fid,blkrec,facfile);
        if i >= i0                                          % record (i-i0) holds station i
            bt_fwrite(fidWf,[real(Wi(:));imag(Wi(:))],wfwd);
        end
    else
        if use_lu
            Lc{i}=Lg; Uc{i}=Ug; Pc{i}=pg;
        else
            Gc{i}=Gi;
        end
        Yc{i}=Yi; Zc{i}=Zi; Wc{i}=Wi;
    end
    Ym2=Ym1; Ym1=Yi; Zm2=Zm1; Zm1=Zi; rm2=rm1; rm1=Wi;
end
if stream; fclose(fid); fclose(fidWf); end

% back-substitution for the interface response W = D \ Bcols
qp1=zeros(b,nc); qp2=zeros(b,nc);
if stream
    fid   = fopen(facfile,'r'); bb = b*b; rec = 6*bb + use_lu*b;
    fidWf = fopen(wfwd,'r');    recWf = 2*b*nc;
    fidW  = fopen(wfile,'w');   recW = 2*b*nc;   % W is written out contiguously, in reverse station order
    for i = N:-1:1                               % record (N-i) in the file holds station i
        fseek(fid,(i-1)*rec*8,'bof'); v = bt_fread(fid,rec,facfile);
        Yi = reshape(complex(v(2*bb+1:3*bb), v(3*bb+1:4*bb)), b, b);
        Zi = reshape(complex(v(4*bb+1:5*bb), v(5*bb+1:6*bb)), b, b);
        if i >= i0
            fseek(fidWf,(i-i0)*recWf*8,'bof'); vw = bt_fread(fidWf,recWf,wfwd);
            Wi = reshape(complex(vw(1:b*nc), vw(b*nc+1:end)), b, nc);
        else
            Wi = zW;
        end
        qi = Wi - Yi*qp1 - Zi*qp2;
        bt_fwrite(fidW,[real(qi(:));imag(qi(:))],wfile);
        qp2=qp1; qp1=qi;
    end
    fclose(fid); fclose(fidWf); fclose(fidW); delete(wfwd);
    part = struct('N',N, 'facfile',facfile, 'use_lu',use_lu);
    varargout = {part};
else
    W = zeros(n,nc);
    for i = N:-1:1
        qi = Wc{i} - Yc{i}*qp1 - Zc{i}*qp2;
        W((i-1)*b+(1:b),:) = qi; qp2=qp1; qp1=qi;
    end
    if use_lu
        part = struct('L',{Lc}, 'U',{Uc}, 'P',{Pc}, 'Y',{Yc}, 'Z',{Zc});
    else
        part = struct('G',{Gc}, 'Y',{Yc}, 'Z',{Zc});
    end
    varargout = {W, part};
end
end

%% Apply: solve ML*x = rhs from a precomputed factor F (rhs may have many columns)
function x = bt_apply(F, rhs)
b = F.b; n = F.n; np = F.np; nc = size(rhs,2);
ft = rhs(F.ifDof,:);
up = cell(np,1); fcon = cell(np,1);
par = isfield(F,'par') && F.par;
if par
parfor p = 1:np
    up{p} = local_thomas_apply(F.part{p}, F.Dp{p}, b, rhs(F.pdof{p},:), F.ooc);
    fcon{p} = F.Cp{p}*up{p};        % deferred, as in the factorization
end
else
for p = 1:np
    up{p} = local_thomas_apply(F.part{p}, F.Dp{p}, b, rhs(F.pdof{p},:), F.ooc);
    fcon{p} = F.Cp{p}*up{p};
end
end
for p = 1:np
    ft = ft - fcon{p};
end
fcon = [];
qs = F.StU \ (F.StL \ ft(F.Stp,:));
x = zeros(n, nc); x(F.ifDof,:) = qs;
xp = cell(np,1);
if par
parfor p = 1:np
    xp{p} = bt_backsub(F, p, b, up{p}, qs, nc);
end
else
for p = 1:np
    xp{p} = bt_backsub(F, p, b, up{p}, qs, nc);
end
end
for p = 1:np
    x(F.pdof{p},:) = xp{p};
end
end

function xp = bt_backsub(F, p, b, upp, qs, nc)
%BT_BACKSUB  Recover partition P's solution once the interface is known.
    if true
    if F.ooc
        ncols = F.wsize(p,2); Np = round(F.wsize(p,1)/b);
        qcols = qs(F.wcol{p},:); xp = zeros(Np*b, nc);
        fidW = fopen(F.wfile{p},'r'); recW = 2*b*ncols;
        for i = 1:Np
            fseek(fidW,(Np-i)*recW*8,'bof'); v = bt_fread(fidW,recW,F.wfile{p});
            qi = reshape(complex(v(1:b*ncols), v(b*ncols+1:end)), b, ncols);
            rows = (i-1)*b+(1:b);
            xp(rows,:) = upp(rows,:) - qi*qcols;
        end
        fclose(fidW);
    else
        xp = upp - F.W{p}*qs(F.wcol{p},:);
    end
    end
end

function up = local_thomas_apply(part, D, b, rhsp, ooc)
%LOCAL_THOMAS_APPLY  Compute the particular solution up = D \ rhsp on one partition.
%   K_i and A_i are recomputed from the sparse partition matrix D and a rolling
%   Y_{i-2} buffer, so only the G_i LU factors and Y_i, Z_i need to be kept
%   around (read from disk when OOC).
n = size(rhsp,1); N = round(n/b); zb = sparse(b,b); nc = size(rhsp,2);
rm1 = zeros(b,nc); rm2 = zeros(b,nc); rr = cell(N,1);
blk = @(i,j) D((i-1)*b+(1:b),(j-1)*b+(1:b));
Yb1 = zb; Yb2 = zb;                      % Y_{i-1}, Y_{i-2}
if ooc; fid = fopen(part.facfile,'r'); bb = b*b; rec = 6*bb + part.use_lu*b; end
for i = 1:N
    if i>=3, A = blk(i,i-2); Ym2 = Yb2; else, A = zb; Ym2 = zb; end
    if i>=2, Bb = blk(i,i-1); else, Bb = zb; end
    K = Bb - A*Ym2;
    if ooc
        fseek(fid,(i-1)*rec*8,'bof'); v = bt_fread(fid,2*bb,part.facfile);   % G^-1 or packed LU
        Gs = reshape(complex(v(1:bb), v(bb+1:2*bb)), b, b);
        fseek(fid,((i-1)*rec+2*bb)*8,'bof'); vy = bt_fread(fid,2*bb,part.facfile); % Y
        Yi = reshape(complex(vy(1:bb), vy(bb+1:2*bb)), b, b);
        if part.use_lu
            fseek(fid,((i-1)*rec+6*bb)*8,'bof'); vp = bt_fread(fid,b,part.facfile); % pivots
            Lg = tril(Gs,-1) + speye(b); Ug = triu(Gs); pg = round(vp(:)).';
        else
            Gi = Gs;
        end
    else
        Yi = part.Y{i};
        if isfield(part,'G'); Gi = part.G{i}; end
    end
    V = rhsp((i-1)*b+(1:b),:) - K*rm1 - A*rm2;
    if ooc
        if part.use_lu
            V = full(V);
            rr{i} = Ug \ (Lg \ V(pg,:));               % streamed pivoted factors
        else
            rr{i} = Gi*V;                               % streamed inverse
        end
    elseif isfield(part,'G')
        rr{i} = Gi*V;                                   % stored inverse
    else
        V = full(V);                                    % stored pivoted factors
        rr{i} = part.U{i} \ (part.L{i} \ V(part.P{i},:));
    end
    Yb2 = Yb1; Yb1 = Yi; rm2 = rm1; rm1 = rr{i};
end
up = zeros(n,nc); qp1 = zeros(b,nc); qp2 = zeros(b,nc);
for i = N:-1:1
    if ooc
        fseek(fid,((i-1)*rec+2*bb)*8,'bof'); v = bt_fread(fid,4*bb,part.facfile); % Y,Z
        Yi = reshape(complex(v(1:bb),      v(bb+1:2*bb)),   b, b);
        Zi = reshape(complex(v(2*bb+1:3*bb), v(3*bb+1:4*bb)), b, b);
    else
        Yi = part.Y{i}; Zi = part.Z{i};
    end
    qi = rr{i} - Yi*qp1 - Zi*qp2;
    up((i-1)*b+(1:b),:) = qi; qp2 = qp1; qp1 = qi;
end
if ooc; fclose(fid); end
end

%% Free / out-of-core temporary-file registry (see the header above for what each mode does)
function out = bt_free(F, mode)
persistent reg
if isempty(reg); reg = {}; end
out = [];
if nargin >= 2
    switch mode
        case 'reg'
            reg{end+1} = F;                          % F is an oocdir string
            return
        case 'mark'
            out = numel(reg);
            return
        case {'sweep','sweep_all'}
            mk = 0;
            if strcmp(mode,'sweep') && ~isempty(F); mk = F; end
            for i = numel(reg):-1:mk+1
                if ~isempty(reg{i}) && exist(reg{i},'dir'); rmdir(reg{i},'s'); end
                reg(i) = [];
            end
            return
    end
end

% default: free a factor struct's files
if isstruct(F) && isfield(F,'ooc') && F.ooc && isfield(F,'oocdir') ...
        && ~isempty(F.oocdir) && exist(F.oocdir,'dir')
    rmdir(F.oocdir,'s');
end
end

%% Guarded out-of-core I/O (turns a disk-full truncation into a clear error message)
function v = bt_fread(fid, n, fname)
v = fread(fid, n, 'double');
if numel(v) ~= n
    error('block_thomas:oocTruncated', ...
        ['out-of-core read from ''%s'' returned %d of %d values: the file is ' ...
         'truncated, most likely from a full disk during factorization. Free disk ' ...
         'space, lower the resolution, or use lu_mode=''lapack_band''.'], ...
        fname, numel(v), n);
end
end

function bt_fwrite(fid, data, fname)
n = numel(data);
if fwrite(fid, data, 'double') ~= n
    error('block_thomas:oocDiskFull', ...
        ['out-of-core write to ''%s'' was truncated: the disk is probably full. ' ...
         'The block-Thomas factor is streamed to disk and needs several GB per ' ...
         'partition at this resolution. Free disk space, lower the resolution, or ' ...
         'use lu_mode=''lapack_band''.'], fname);
end
end

function [Sloc, W, part] = bt_partfactor(Dp, Bcols, Cp, b, nif, ooc, facfile, wfile, use_lu)
%BT_PARTFACTOR  Factor one partition and return its interface contribution SLOC.
%   Receives only its own slice, so the loop that calls it carries no shared
%   state and may be a PARFOR on either a thread- or a process-based pool.
W = [];
if ooc
    part = local_thomas_factor(Dp, b, Bcols, facfile, wfile, use_lu);
    % Form Sloc = Cp*W by streaming W back from disk one block at a time, so the
    % partition never holds more than a single block of it and the out-of-core
    % footprint is unchanged by the deferral.
    ncols = size(Bcols,2); Np = round(size(Dp,1)/b);
    Sloc = zeros(nif, ncols);
    fidW = fopen(wfile,'r'); recW = 2*b*ncols;
    for i = 1:Np
        fseek(fidW,(Np-i)*recW*8,'bof'); v = bt_fread(fidW,recW,wfile);
        qi = reshape(complex(v(1:b*ncols), v(b*ncols+1:end)), b, ncols);
        rows = (i-1)*b+(1:b);
        Sloc = Sloc + Cp(:,rows)*qi;
    end
    fclose(fidW);
else
    [W, part] = local_thomas_factor(Dp, b, Bcols, [], [], use_lu);
    Sloc = Cp*W;
end
end
