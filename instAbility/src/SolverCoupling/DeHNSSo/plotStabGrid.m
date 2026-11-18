%% plotStabGrid — visual assessment of a DeHNSSo StabGrid after `mesh`
%
% Runs after the gridgen case script, reading StabGrid from the base workspace
% (or from gridgen's saved input/*.mat as a fallback). Read-only with respect to
% the DeHNSSo repository.
%
% Expected in the workspace (all optional; panels degrade gracefully):
%   StabGrid  gridgen output. X, Y are [ny x nx] with eta(1) = H (freestream) and
%             eta(end) = 0 (wall), so the wall is the LAST row.
%   figfile   where to save the figure (no extension). Supplied by run.jl.
%   stabCaller path of the solver caller; beta_0 and N are read from it here
%             (beta_0 is usually an expression in StabGrid.lref, evaluated with
%             the StabGrid above). Supplied by run.jl. A ready-made struct
%             `stab` (beta_0, N) is used instead when present.
%
% Panels
%   1 full-domain grid           gross body-fitting failures          (top left)
%   2 dxi vs x                   clustering in the wrong place,       (top right)
%                                trim artefacts
%   3 stretching dxi(i+1)/dxi(i) smoothness for the 4th-order xi      (middle right,
%                                operators                            under panel 2)
%   4 eta distribution mid-domain the wall-normal resolution question (middle left)
%   5 resolution summary         whether nx meets the streamwise      (bottom left)
%                                criteria

if ~exist('StabGrid', 'var')
    error('plotStabGrid: StabGrid not in the workspace; run the gridgen case first.');
end
if ~exist('stab', 'var')
    stab = struct();
    if exist('stabCaller', 'var') && isfile(stabCaller)
        stab = readCallerStab(stabCaller, StabGrid);
    end
end

gv = @(s, f, d) (isfield(s, f) && ~isempty(s.(f))) * 1;   % presence test helper
getf = @(s, f, d) ternary_local(gv(s, f, d), s, f, d);

% A Cartesian (flat-wall) StabGrid carries only the 1-D x, y; the 2-D X, Y are
% written for curvilinear grids alone, so rebuild them for a Cartesian one.
if isfield(StabGrid, 'X') && isfield(StabGrid, 'Y')
    X = StabGrid.X;  Y = StabGrid.Y;
else
    [X, Y] = meshgrid(StabGrid.x(:).', StabGrid.y(:));
end
[ny, nx] = size(X);
lref = getf(StabGrid, 'lref', NaN);

% Wall row and its arc length (grid units).
xw = X(end, :).';   yw = Y(end, :).';
sw = [0; cumsum(hypot(diff(xw), diff(yw)))];
L  = sw(end);

% Station where the wall-normal distribution is shown: mid-domain.
[~, jc] = min(abs(sw - 0.5 * L));

figure('Name', 'StabGrid assessment', 'Position', [60 60 1250 860]);
tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 1. Full-domain grid -----------------------------------------------------
nexttile; hold on
kx = max(1, round(nx / 60));   ke = max(1, round(ny / 20));
plot(X(:, 1:kx:nx), Y(:, 1:kx:nx), 'Color', [.75 .75 .75], 'LineWidth', .25);
plot(X(1:ke:ny, :).', Y(1:ke:ny, :).', 'Color', [.75 .75 .75], 'LineWidth', .25);
plot(xw, yw, 'k', 'LineWidth', 1.2);
plot(xw(jc), yw(jc), 'r.', 'MarkerSize', 14);
daspect([1 1 1]); box on
xlabel('$x/\ell$', 'Interpreter', 'latex');
ylabel('$y/\ell$', 'Interpreter', 'latex');
title(sprintf('grid  (%d \\times %d, every %d\\times%d line)', nx, ny, kx, ke));
xlim([min(xw) max(xw)]);

%% 2. dxi vs x -------------------------------------------------------------
nexttile; hold on
dxi = diff(sw);   xm = 0.5 * (xw(1:end-1) + xw(2:end));
du  = L / (nx - 1);
plot(xm, dxi / du, 'k', 'LineWidth', 1);
yline(1, 'r--');
% With no clustering the spread is round-off (1e-5), and autoscale would magnify
% it into spurious-looking structure. Pin the axis unless it is really varying.
if (max(dxi) - min(dxi)) / du < 1e-3
    ylim([0.9 1.1]);
    text(0.02, 0.92, 'uniform to <0.1%', 'Units', 'normalized', 'FontSize', 8);
end
grid on; box on
xlabel('$x/\ell$', 'Interpreter', 'latex');
ylabel('$\Delta\xi/\Delta\xi_{\rm unif}$', 'Interpreter', 'latex');
title(sprintf('streamwise spacing   (min %.3f of uniform)', min(dxi) / du));
xlim([min(xw) max(xw)]);

%% 3. Stretching -----------------------------------------------------------
nexttile(4)                     % right column, under the streamwise spacing
r = dxi(2:end) ./ dxi(1:end-1);
r = max(r, 1 ./ r);                       % worst direction at each interface
plot(xm(2:end), r, 'k', 'LineWidth', 1); hold on
yline(1.01, 'r--');
grid on; box on
xlabel('$x/\ell$', 'Interpreter', 'latex');
ylabel('$\max(\Delta\xi_{i+1}/\Delta\xi_i,\,{\rm inv})$', 'Interpreter', 'latex');
title(sprintf('stretching   (worst %.5f; 1.01 target)', max(r)));
xlim([min(xw) max(xw)]);

%% 4. eta distribution mid-domain -------------------------------
nexttile(3); hold on             % left column, middle row
d = hypot(X(:, jc) - X(end, jc), Y(:, jc) - Y(end, jc));   % height above wall
d = flipud(d);                                             % wall -> freestream
U = flipud(StabGrid.U(:, jc));
Ue = U(end);
d99 = NaN;
if isfinite(Ue) && Ue ~= 0
    k = find(U >= 0.99 * Ue, 1, 'first');
    if ~isempty(k) && k > 1
        t = (0.99 * Ue - U(k-1)) / (U(k) - U(k-1));
        d99 = d(k-1) + t * (d(k) - d(k-1));
    end
end
% Distinct colours AND line styles: delta_99 and y_i are often close, so the two
% reference lines must stay tellable apart even where they overlap.
h = plot(1:ny, max(d, 1e-6), 'k.-', 'MarkerSize', 8);
hs = h;  lbl = {'grid points'};
if isfinite(d99)
    hs(end+1) = yline(d99, '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.2);
    lbl{end+1} = sprintf('\\delta_{99} = %.2f  (%d pts below)', d99, sum(d <= d99));
end
yi = getf(StabGrid, 'y_i', NaN);
if isfinite(yi)
    hs(end+1) = yline(yi, '-.', 'Color', [0.95 0.55 0.00], 'LineWidth', 1.2);
    lbl{end+1} = sprintf('y_i = %.2f  (%d pts below)', yi, sum(d <= yi));
end
legend(hs, lbl, 'Location', 'northwest');
grid on; box on
xlabel('$\eta$ index (from wall)', 'Interpreter', 'latex');
ylabel('wall distance $/\ell$', 'Interpreter', 'latex');
title('wall-normal points, mid-domain');

%% 5. Resolution summary ---------------------------------------------------
nexttile(5); axis off
txt = {};
txt{end+1} = sprintf('nx = %d   ny = %d', nx, ny);
txt{end+1} = sprintf('L  = %.1f x/l  =  %.2f mm', L, L * lref * 1e3);
txt{end+1} = sprintf('dxi = %.4f x/l  =  %.5f mm', du, du * lref * 1e3);
b0 = getf(stab, 'beta_0', NaN);
if isfinite(b0) && b0 > 0
    lz = 2 * pi / b0;                       % spanwise wavelength, grid units
    txt{end+1} = '';
    txt{end+1} = sprintf('beta_0 = %.4f  ->  lambda_z = %.3f mm', b0, lz * lref * 1e3);
    txt{end+1} = 'lambda_x taken as lambda_z (proxy; alpha_r known only after the solve)';
    txt{end+1} = '   N    pts/lambda    accum. 4th-order phase err';
    Nmax = max(6, getf(stab, 'N', 6));
    for n = 1:Nmax
        ln = lz / n;   pp = ln / du;   Nw = L / ln;
        acc = 360 * Nw * (2 * pi / pp)^4 / 30;    % degrees
        flag = '';
        if pp < 10;  flag = '   <-- under 10 pts';
        elseif acc > 30; flag = '   <-- phase error high'; end
        txt{end+1} = sprintf('  %2d    %8.1f    %12.1f deg%s', n, pp, acc, flag);
    end
end
if isfield(StabGrid, 'cluster') && ~isempty(StabGrid.cluster)
    txt{end+1} = '';
    txt{end+1} = 'clustering ACTIVE (StabGrid.cluster present)';
else
    txt{end+1} = '';
    txt{end+1} = 'clustering: none (uniform xi)';
end
text(0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'Courier', 'FontSize', 9, 'Interpreter', 'none');
title('resolution summary');

title(tl, sprintf('StabGrid assessment   (Re = %.4g)', getf(StabGrid, 'Re', NaN)));

% `matlab -batch` has no display, so the figure must be written to disk or the
% user would never see it. Saved as png.
if exist('figfile', 'var') && ~isempty(figfile)
    [fd, ~, ~] = fileparts(figfile);
    if ~isempty(fd) && ~isfolder(fd); mkdir(fd); end
    try
        exportgraphics(gcf, [figfile '.png'], 'Resolution', 150);
    catch
        print(gcf, [figfile '.png'], '-dpng', '-r150');
    end
    fprintf('\n[airPower] grid assessment figure -> %s.png\n', figfile);
end

% Stab.beta_0 and Stab.N from a solver caller: the right-hand side of their last
% assignment, evaluated with StabGrid in scope (beta_0 is written in terms of
% StabGrid.lref). A field that cannot be read or evaluated is left out.
function stab = readCallerStab(caller, StabGrid) %#ok<INUSD>
    stab = struct();
    lines = splitlines(fileread(caller));
    for f = {'beta_0', 'N'}
        rhs = '';
        for k = 1:numel(lines)
            ln  = regexprep(lines{k}, '%.*$', '');
            tok = regexp(ln, ['Stab\.' f{1} '\s*=\s*([^;]+);'], 'tokens', 'once');
            if ~isempty(tok); rhs = tok{1}; end
        end
        if isempty(rhs); continue; end
        try
            stab.(f{1}) = eval(rhs);
        catch
        end
    end
end

function v = ternary_local(present, s, f, d)
    if present; v = s.(f); else; v = d; end
end
