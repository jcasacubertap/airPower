function [MLAMFD,MLA1,MLA2,MLA3,MLA4] = LHS(BF,Bound,Grid,Init,Stab,StabRes,Opt)
%LHS  Assemble the global sparse LHS matrices for the HNS system.
%
%   [MLAMFD, MLA1, MLA2, MLA3, MLA4] = LHS(BF, Bound, Grid, Init, Stab, StabRes, Opt)
%
%   Constructs the block-pentadiagonal LHS operator from local coefficient
%   matrices (Ai..Gi) at each streamwise station using 4th-order FD stencils.
%   Sparse triplets are accumulated in cell arrays (O(n) cost).
%
%   Parallel assembly controlled by Opt.parfor:
%     true  — parfor over stations i=2:nx (requires Parallel Computing Toolbox)
%     false — standard for loop with progress reporting
%
%   Outputs:
%     MLAMFD  [N x N] sparse — LHS for the mean-flow-distortion (omega=0, beta=0)
%     MLA1    [N x N] sparse — base LHS (common to all modes)
%     MLA2    [N x N] sparse — omega-dependent part  (ML = MLA1 + omega*MLA2)
%     MLA3    [N x N] sparse — beta-dependent part   (ML = MLA1 + ... + beta*MLA3)
%     MLA4    [N x N] sparse — beta^2-dependent part (ML = MLA1 + ... + beta^2*MLA4)
%   where N = 4 * ny * nx.
%   For 2D cases (N=0, all beta=0), MLA3 and MLA4 are empty sparse matrices.

nx    = Grid.nx;
ny    = Grid.ny;
bsize = 4*ny;

if ~isfield(Opt,'parfor'); Opt.parfor = 'off'; end

% Pre-allocate cell arrays for sparse triplets (avoids O(n^2) concatenation)
cell_i1MFD = cell(nx,1);  cell_j1MFD = cell(nx,1);  cell_s1MFD = cell(nx,1);
cell_i1    = cell(nx,1);  cell_j1    = cell(nx,1);  cell_s1    = cell(nx,1);
cell_i2    = cell(nx,1);  cell_j2    = cell(nx,1);  cell_s2    = cell(nx,1);
cell_i3    = cell(nx,1);  cell_j3    = cell(nx,1);  cell_s3    = cell(nx,1);
cell_i4    = cell(nx,1);  cell_j4    = cell(nx,1);  cell_s4    = cell(nx,1);

% Chebyshev/Malik grid: index 1 = freestream (y=H), index ny = wall (y=0).
bc_top     = [1, ny+1, 2*ny+1];
bc_top_mfd = [1, 2*ny+1];
bc_wall    = [ny, 2*ny, 3*ny];
bc_top_type = Opt.bc_top;  % cell array {'H_DR','H_DR','H_DR'} per component (u,v,w)

iu  = Init.iu;
D1  = Grid.D1;
D2  = Grid.D2;
I   = Stab.I;
Z   = Stab.Z;
Re  = BF.Re;

% Per-station Fornberg weights for ξ derivatives (1st and 2nd order).
% Reduces to classical 4th-order coefficients on uniform grids; on
% non-uniform grids it stays consistent. Boundary stencils are one-sided
% (5/4/3 points at i = 2 / nx-1 / nx).
x_xi = Grid.x(:);
W1 = zeros(nx, 5);
W2 = zeros(nx, 5);
if nx >= 5
    c = fornberg_weights(x_xi(2), x_xi(1:5), 2);
    W1(2, 1:5) = c(:, 2).';
    W2(2, 1:5) = c(:, 3).';
    for ii = 3 : nx-2
        c = fornberg_weights(x_xi(ii), x_xi(ii-2:ii+2), 2);
        W1(ii, 1:5) = c(:, 2).';
        W2(ii, 1:5) = c(:, 3).';
    end
    c = fornberg_weights(x_xi(nx-1), x_xi(nx-3:nx), 2);
    W1(nx-1, 1:4) = c(:, 2).';
    W2(nx-1, 1:4) = c(:, 3).';
    c = fornberg_weights(x_xi(nx), x_xi(nx-2:nx), 2);
    W1(nx, 1:3) = c(:, 2).';
    W2(nx, 1:3) = c(:, 3).';
end

% Use full (dense) arithmetic for large Chebyshev operators
D1_is_sparse = issparse(Grid.D1) && (nnz(Grid.D1) < 0.2 * ny^2);
use_full = (ny >= 60) && ~D1_is_sparse;
if use_full
    D1f = full(D1);  D2f = full(D2);  If = full(I);  Zf = full(Z);
else
    D1f = D1;  D2f = D2;  If = I;  Zf = Z;
end

%% Station i=1: inflow — identity block (no PDE, just stores IC)
[ii,jj,ss] = find(speye(bsize));
cell_i1MFD{1} = reshape(ii,[],1);  cell_j1MFD{1} = reshape(jj,[],1);  cell_s1MFD{1} = reshape(ss,[],1);
cell_i1{1}    = reshape(ii,[],1);  cell_j1{1}    = reshape(jj,[],1);  cell_s1{1}    = reshape(ss,[],1);
cell_i2{1}    = zeros(0,1);  cell_j2{1}    = zeros(0,1);  cell_s2{1}    = zeros(0,1);
cell_i3{1}    = zeros(0,1);  cell_j3{1}    = zeros(0,1);  cell_s3{1}    = zeros(0,1);
cell_i4{1}    = zeros(0,1);  cell_j4{1}    = zeros(0,1);  cell_s4{1}    = zeros(0,1);

%% Stations i=2:nx
if strcmpi(Opt.parfor, 'on')
    fprintf('   LHS: assembling stations 2:%d in parallel ...\n', nx);
    parfor i = 2:nx
        [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4] = ...
            lhs_station(i, nx, bsize, ny, BF, Grid, Bound, use_full, D1f, D2f, If, Zf, Re, iu, ...
                        W1(i,:), W2(i,:), bc_top, bc_top_mfd, bc_wall, bc_top_type);
        cell_i1MFD{i}=ci1MFD; cell_j1MFD{i}=cj1MFD; cell_s1MFD{i}=cs1MFD;
        cell_i1{i}=ci1;       cell_j1{i}=cj1;       cell_s1{i}=cs1;
        cell_i2{i}=ci2;       cell_j2{i}=cj2;       cell_s2{i}=cs2;
        cell_i3{i}=ci3;       cell_j3{i}=cj3;       cell_s3{i}=cs3;
        cell_i4{i}=ci4;       cell_j4{i}=cj4;       cell_s4{i}=cs4;
    end
    fprintf('   LHS: parallel assembly done\n');
else
    t_lhs = tic;
    lastPerc = 0;
    for i = 2:nx
        [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4] = ...
            lhs_station(i, nx, bsize, ny, BF, Grid, Bound, use_full, D1f, D2f, If, Zf, Re, iu, ...
                        W1(i,:), W2(i,:), bc_top, bc_top_mfd, bc_wall, bc_top_type);
        cell_i1MFD{i}=ci1MFD; cell_j1MFD{i}=cj1MFD; cell_s1MFD{i}=cs1MFD;
        cell_i1{i}=ci1;       cell_j1{i}=cj1;       cell_s1{i}=cs1;
        cell_i2{i}=ci2;       cell_j2{i}=cj2;       cell_s2{i}=cs2;
        cell_i3{i}=ci3;       cell_j3{i}=cj3;       cell_s3{i}=cs3;
        cell_i4{i}=ci4;       cell_j4{i}=cj4;       cell_s4{i}=cs4;

        perc = floor(i * 100 / nx);
        if perc >= lastPerc + 25 && perc < 100
            elapsed = toc(t_lhs);
            eta = elapsed / i * (nx - i);
            fprintf('   LHS: %d%% (%.0fs elapsed, ~%.0fs remaining)\n', perc, elapsed, eta);
            lastPerc = perc;
        end
    end
    fprintf('   LHS: 100%% (%.1fs total)\n', toc(t_lhs));
end

%% Assemble global sparse matrices from triplets
N = bsize * nx;

MLAMFD = sparse(cell2mat(cell_i1MFD), cell2mat(cell_j1MFD), cell2mat(cell_s1MFD), N, N);
clear cell_i1MFD cell_j1MFD cell_s1MFD

MLA1 = sparse(cell2mat(cell_i1), cell2mat(cell_j1), cell2mat(cell_s1), N, N);
clear cell_i1 cell_j1 cell_s1

if any(StabRes.omegavec)
    MLA2 = sparse(cell2mat(cell_i2), cell2mat(cell_j2), cell2mat(cell_s2), N, N);
else
    MLA2 = sparse(N, N);
end
clear cell_i2 cell_j2 cell_s2

if any(StabRes.betavec)
    MLA3 = sparse(cell2mat(cell_i3), cell2mat(cell_j3), cell2mat(cell_s3), N, N);
    MLA4 = sparse(cell2mat(cell_i4), cell2mat(cell_j4), cell2mat(cell_s4), N, N);
else
    MLA3 = sparse(N, N);
    MLA4 = sparse(N, N);
end
clear cell_i3 cell_j3 cell_s3 cell_i4 cell_j4 cell_s4

end


%% Local function: assemble triplets for one station i (i >= 2)
function [ci1MFD,cj1MFD,cs1MFD, ci1,cj1,cs1, ci2,cj2,cs2, ci3,cj3,cs3, ci4,cj4,cs4] = ...
    lhs_station(i, nx, bsize, ny, BF, Grid, Bound, use_full, D1f, D2f, If, Zf, Re, iu, ...
                w1, w2, bc_top, bc_top_mfd, bc_wall, bc_top_type)

% Base flow diagonal matrices
if use_full
    U   = diag(BF.U(:,i));    dxU = diag(BF.dxU(:,i));  dyU = diag(BF.dyU(:,i));
    V   = diag(BF.V(:,i));    dxV = diag(BF.dxV(:,i));  dyV = diag(BF.dyV(:,i));
    W   = diag(BF.W(:,i));    dxW = diag(BF.dxW(:,i));  dyW = diag(BF.dyW(:,i));
else
    U   = spdiags(BF.U(:,i),   0,ny,ny);  dxU = spdiags(BF.dxU(:,i),0,ny,ny);  dyU = spdiags(BF.dyU(:,i),0,ny,ny);
    V   = spdiags(BF.V(:,i),   0,ny,ny);  dxV = spdiags(BF.dxV(:,i),0,ny,ny);  dyV = spdiags(BF.dyV(:,i),0,ny,ny);
    W   = spdiags(BF.W(:,i),   0,ny,ny);  dxW = spdiags(BF.dxW(:,i),0,ny,ny);  dyW = spdiags(BF.dyW(:,i),0,ny,ny);
end

% Metric coefficients
xix  = Grid.xix(:,i);   etax = Grid.etax(:,i);
xiy  = Grid.xiy(:,i);   etay = Grid.etay(:,i);
xixx = Grid.xixx(:,i);  xiyy = Grid.xiyy(:,i);
etaxx= Grid.etaxx(:,i); etayy= Grid.etayy(:,i);

% Buffer
bufs_i  = Bound.bufs(:,i);
bufsp_i = Bound.bufsp(:,i);
bufc_i  = Bound.bufc(i);
bufp_i  = Bound.bufp(i);

% Common terms
Act = U.*etax*D1f + V.*etay*D1f ...
    - 1/Re*(etax.^2+etay.^2).*D2f - 1/Re*(etaxx+etayy).*D1f;
Ect = U.*xix + V.*xiy - 1/Re*(xixx+xiyy).*If ...
    - 1/Re*(2*etax.*xix+2*etay.*xiy).*D1f;

% Coefficient matrices
Ai = [(Act+dxU).*bufs_i    dyU.*bufs_i       Zf    etax.*D1f.*bufs_i
       dxV.*bufs_i         (Act+dyV).*bufs_i  Zf    etay.*D1f.*bufs_i
       dxW.*bufs_i          dyW.*bufs_i        Act.*bufs_i  Zf
      (etax.*D1f).*bufsp_i (etay.*D1f).*bufsp_i  Zf   Zf];

Bi = [-iu*If.*bufs_i  Zf             Zf             Zf
       Zf            -iu*If.*bufs_i  Zf             Zf
       Zf             Zf            -iu*If.*bufs_i  Zf
       Zf             Zf             Zf             Zf];

Ci = [iu*W.*bufs_i  Zf           Zf           Zf
       Zf           iu*W.*bufs_i  Zf           Zf
       Zf           Zf           iu*W.*bufs_i  iu*If.*bufs_i
       Zf           Zf           iu*If.*bufsp_i  Zf];

Di = [If/Re.*bufs_i  Zf           Zf           Zf
       Zf           If/Re.*bufs_i  Zf           Zf
       Zf           Zf           If/Re.*bufs_i  Zf
       Zf           Zf           Zf            Zf];

Ei = [Ect.*bufs_i   Zf          Zf          If.*xix.*bufs_i
       Zf           Ect.*bufs_i  Zf          If.*xiy.*bufs_i
       Zf           Zf          Ect.*bufs_i  Zf
      If.*xix.*bufsp_i  If.*xiy.*bufsp_i  Zf  Zf];

Fi = [-If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i  Zf  Zf  Zf
       Zf  -If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i  Zf  Zf
       Zf  Zf  -If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i  Zf
       Zf  Zf  Zf  Zf];

if use_full
    Pres = diag(1-bufsp_i) + diag(-1*(1-Bound.bufsp(2:end,i)), -1);
else
    Pres = If.*(1-bufsp_i) + spdiags(-1*(1-Bound.bufsp(2:end,i)), -1, ny, ny);
end

% Build Gi: identity rows for masked points
Gi_u = If.*(1-bufs_i) + If.*Bound.bufs(:,end)*(1-bufc_i);
Gi_v = If.*(1-bufs_i);
Gi_w = If.*(1-bufs_i) + If.*Bound.bufs(:,end)*(1-bufc_i);

Gi = [Gi_u  Zf   Zf   Zf
       Zf   Gi_v  Zf   Zf
       Zf   Zf   Gi_w  Zf
       Zf   Zf   Zf   Pres];

% Determine stencil: stenr = global column indices, cur_blk = 1-index of current station
sten = bsize*(i-1)+1 : bsize*i;   % global row (= global col of diagonal block)

% Per-station FD stencil weights (w1, w2) are precomputed in the parent
% function using Fornberg's algorithm on Grid.x, so xrefined / non-uniform
% streamwise grids are handled correctly. On a uniform grid these reduce
% to the classical 4th-order coefficients.
if i == 2
    stenr   = [sten-bsize, sten, sten+bsize, sten+2*bsize, sten+3*bsize];
    cur_blk = 2;
    blk1MFD = [ Ei*w1(1) + Fi*w2(1), ...
              (Ai + Ei*w1(2) + Fi*w2(2)) + Gi/bufc_i, ...
                Ei*w1(3) + Fi*w2(3), ...
                Ei*w1(4) + Fi*w2(4), ...
                Ei*w1(5) + Fi*w2(5)] * bufc_i;

elseif i == nx-1
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize];
    cur_blk = 3;
    blk1MFD = [ Ei*w1(1) + Fi*w2(1), ...
                Ei*w1(2) + Fi*w2(2), ...
              (Ai + Ei*w1(3) + Fi*w2(3)) + Gi/bufc_i, ...
                Ei*w1(4) + Fi*w2(4)] * bufc_i;

elseif i == nx
    stenr   = [sten-2*bsize, sten-bsize, sten];
    cur_blk = 3;
    blk1MFD = [ Ei*w1(1) + Fi*w2(1), ...
                Ei*w1(2) + Fi*w2(2), ...
              (Ai + Ei*w1(3) + Fi*w2(3)) + Gi/bufc_i] * bufc_i;

else   % interior (centred 5-point)
    stenr   = [sten-2*bsize, sten-bsize, sten, sten+bsize, sten+2*bsize];
    cur_blk = 3;
    blk1MFD = [ Ei*w1(1) + Fi*w2(1), ...
                Ei*w1(2) + Fi*w2(2), ...
              (Ai + Ei*w1(3) + Fi*w2(3)) + Gi/bufc_i, ...
                Ei*w1(4) + Fi*w2(4), ...
                Ei*w1(5) + Fi*w2(5)] * bufc_i;
end

ncols     = numel(stenr);          % = nblocks * bsize
blk1      = blk1MFD;              % ML1 identical to ML1MFD except for BCs
cur_off   = (cur_blk-1)*bsize;    % local column offset for the current station block

% ML2/ML3/ML4 only write to the diagonal (current station) block
blk2 = zeros(bsize, ncols);
blk3 = zeros(bsize, ncols);
blk4 = zeros(bsize, ncols);
blk2(:, cur_off+1:cur_off+bsize) = Bi*bufc_i;
blk3(:, cur_off+1:cur_off+bsize) = Ci*bufc_i;
blk4(:, cur_off+1:cur_off+bsize) = Di*bufc_i;

% Apply freestream BCs per component (u, v, w)
%   'H_DR' / 'IH_DR' → identity row (Dirichlet)
%   'H_NM'            → D1 row (Neumann: ∂/∂y = 0)
blk1MFD(bc_top_mfd,:) = 0;
blk1(bc_top,:)   = 0;
blk2(bc_top,:)   = 0;
blk3(bc_top,:)   = 0;
blk4(bc_top,:)   = 0;

for k = 1:3
    r  = (k-1)*ny + 1;
    is_mfd = ismember(r, bc_top_mfd);
    if strcmpi(bc_top_type{k}, 'H_NM')
        cs = cur_off + (k-1)*ny + 1;
        ce = cur_off + k*ny;
        blk1(r, cs:ce) = D1f(1,:);
        if is_mfd; blk1MFD(r, cs:ce) = D1f(1,:); end
    else  % H_DR or IH_DR → identity row
        blk1(r, cur_off+r) = 1;
        if is_mfd; blk1MFD(r, cur_off+r) = 1; end
    end
end

% Apply wall BCs
blk1MFD(bc_wall,:) = 0;
blk1(bc_wall,:)    = 0;
blk2(bc_wall,:)    = 0;
blk3(bc_wall,:)    = 0;
blk4(bc_wall,:)    = 0;
for r = bc_wall
    blk1MFD(r, cur_off+r) = 1;
    blk1(r,    cur_off+r) = 1;
end

% Extract sparse triplets; map local column indices to global via stenr.
% reshape(...,[],1) guarantees column vectors even when find returns empty 0x0.
row_off = bsize*(i-1);

[ii,jj,ss] = find(blk1MFD);  ci1MFD=reshape(ii+row_off,[],1); cj1MFD=reshape(stenr(jj),[],1); cs1MFD=reshape(ss,[],1);
[ii,jj,ss] = find(blk1);     ci1   =reshape(ii+row_off,[],1); cj1   =reshape(stenr(jj),[],1); cs1   =reshape(ss,[],1);
[ii,jj,ss] = find(blk2);     ci2   =reshape(ii+row_off,[],1); cj2   =reshape(stenr(jj),[],1); cs2   =reshape(ss,[],1);
[ii,jj,ss] = find(blk3);     ci3   =reshape(ii+row_off,[],1); cj3   =reshape(stenr(jj),[],1); cs3   =reshape(ss,[],1);
[ii,jj,ss] = find(blk4);     ci4   =reshape(ii+row_off,[],1); cj4   =reshape(stenr(jj),[],1); cs4   =reshape(ss,[],1);

end
