function x = block_penta_solve(ML, rhs, bsize, nx)
%BLOCK_PENTA_SOLVE  Solve a block-pentadiagonal system Ax = b.
%
%   x = block_penta_solve(ML, rhs, bsize, nx)
%
%   Solves ML*x = rhs where ML is block-pentadiagonal with blocks of size
%   bsize x bsize and nx block-rows. Uses block LU elimination with
%   pre-extracted dense blocks for fast BLAS operations.
%
%   Memory: O(bsize^2 * nx) — no global LU factorisation.
%
%   The stencil at each station is determined from the sparse structure.
%   Interior stations have 5 blocks (4th-order FD), boundaries may differ.

N = bsize * nx;

%% Pre-extract all blocks into a 3D array
% For each block-row i, store up to 5 blocks in D(:,:,i,1:5)
% and their block-column indices in jmap(i,1:5). nblk(i) = number of blocks.
max_bw = 5;  % max blocks per row
D    = zeros(bsize, bsize, nx, max_bw);
jmap = zeros(nx, max_bw);  % block-column index for each stored block
nblk = zeros(nx, 1);       % number of blocks per row
diag_loc = zeros(nx, 1);   % position of diagonal block within jmap(i,:)
f    = zeros(bsize, nx);   % RHS per block-row

for i = 1:nx
    rows = (i-1)*bsize+1 : i*bsize;
    f(:,i) = rhs(rows);

    % Find block-column range from sparse structure
    [~, jj] = find(ML(rows,:));
    if isempty(jj); continue; end
    j1 = ceil(min(jj) / bsize);
    j2 = ceil(max(jj) / bsize);
    jidx = j1:j2;
    nb = numel(jidx);
    nblk(i) = nb;
    jmap(i,1:nb) = jidx;
    diag_loc(i) = find(jidx == i);

    for kb = 1:nb
        j = jidx(kb);
        cols = (j-1)*bsize+1 : j*bsize;
        D(:,:,i,kb) = full(ML(rows, cols));
    end
end

%% Forward elimination
for i = 1:nx-1
    dl = diag_loc(i);
    if dl == 0; continue; end

    % LU factorize diagonal block in-place
    [Ld, Ud, pd] = lu(D(:,:,i,dl), 'vector');

    % Eliminate column i from rows below
    for ii = i+1 : min(i+4, nx)
        % Find where column i appears in row ii
        col_pos = 0;
        for kb = 1:nblk(ii)
            if jmap(ii,kb) == i; col_pos = kb; break; end
        end
        if col_pos == 0; continue; end

        % Multiplier: L_block = D(ii,i) * inv(D(i,i))
        Dii_i = D(:,:,ii,col_pos);
        L_block = (Ud \ (Ld \ Dii_i(pd,:).')).';

        % Update all blocks in row ii
        for kb_i = 1:nblk(i)
            j = jmap(i,kb_i);
            % Find matching column in row ii
            for kb_ii = 1:nblk(ii)
                if jmap(ii,kb_ii) == j
                    D(:,:,ii,kb_ii) = D(:,:,ii,kb_ii) - L_block * D(:,:,i,kb_i);
                    break
                end
            end
        end

        % Update RHS
        f(:,ii) = f(:,ii) - L_block * f(:,i);

        % Zero eliminated block
        D(:,:,ii,col_pos) = 0;
    end
end

%% Back substitution
x = zeros(N, 1);
for i = nx:-1:1
    dl = diag_loc(i);
    if dl == 0; continue; end

    r = f(:,i);

    % Subtract known contributions from blocks right of diagonal
    for kb = 1:nblk(i)
        j = jmap(i,kb);
        if j == i; continue; end
        cols = (j-1)*bsize+1 : j*bsize;
        r = r - D(:,:,i,kb) * x(cols);
    end

    % Solve diagonal block
    rows = (i-1)*bsize+1 : i*bsize;
    x(rows) = D(:,:,i,dl) \ r;
end

end
