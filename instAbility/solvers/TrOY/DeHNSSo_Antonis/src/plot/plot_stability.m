function plot_stability(StabGrid, StabRes, Opt)
%PLOT_STABILITY  Plot N-factor, growth rates, wavenumber, and eigenfunctions.
%
%   plot_stability(StabGrid, StabRes, Opt)

x  = StabGrid.xun;
nx = StabGrid.nx;
nf = size(StabRes.A, 1);

% R = sqrt(Re_x) = sqrt(x_norm * R0), with R0 = sqrt(x_norm(1))
R0 = sqrt(x(1));
R  = sqrt(x) * R0;

% Exclude buffer region from plots
if nargin >= 3 && isfield(Opt, 'xb')
    i_buf = round(nx * Opt.xb / 100);
else
    i_buf = round(nx * 0.85);  % default 85%
end
xr = 1:i_buf;  % plot range (indices before buffer)

% Mode labels
omega_vec = StabRes.omegavec;
beta_vec  = StabRes.betavec;
omega_f = min(omega_vec(omega_vec > 0));
beta_f  = min(beta_vec(beta_vec > 0));

% Skip only MFD; flag modes with invalid alpha for N-factor only
has_valid_alpha = true(1, nf);
for j = 1:nf
    if any(~isfinite(StabRes.alpha(j, xr))); has_valid_alpha(j) = false; end
end

% Skip first few inflow points (ILST transient)
xr(1:min(10, numel(xr))) = [];

colors = lines(nf);

%% 1. N-factor and growth rate
figure('Name', 'N-factor & growth rate', 'Position', [100 100 900 500]);

subplot(2,1,1)
hold on
for j = 1:nf
    [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f);
    if (m == 0 && n == 0) || ~has_valid_alpha(j); continue; end
    alpha_i = -imag(StabRes.alpha(j,:));
    N_fac = cumtrapz(x, alpha_i);
    plot(R(xr), N_fac(xr), 'LineWidth', 1.5, 'Color', colors(j,:), ...
        'DisplayName', sprintf('(%d,%d)', m, n));
end
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N$', 'Interpreter', 'latex', 'FontSize', 12)
title('N-factor', 'FontSize', 13)
if nf > 2; legend('Location', 'best'); end
grid on; box on

subplot(2,1,2)
hold on
for j = 1:nf
    [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f);
    if m == 0 && n == 0; continue; end
    alpha_i = -imag(StabRes.alpha(j,:));
    plot(R(xr), alpha_i(xr), 'LineWidth', 1.2, 'Color', colors(j,:), ...
        'DisplayName', sprintf('(%d,%d)', m, n));
end
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$-\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
title('Spatial growth rate', 'FontSize', 13)
if nf > 2; legend('Location', 'best'); end
grid on; box on

%% 2. Wavenumber and phase speed
figure('Name', 'Wavenumber & phase speed', 'Position', [100 600 900 500]);

subplot(2,1,1)
hold on
for j = 1:nf
    [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f);
    if m == 0 && n == 0; continue; end
    alpha_r = real(StabRes.alpha(j,:));
    plot(R(xr), alpha_r(xr), 'LineWidth', 1.2, 'Color', colors(j,:), ...
        'DisplayName', sprintf('(%d,%d)', m, n));
end
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
title('Streamwise wavenumber', 'FontSize', 13)
if nf > 2; legend('Location', 'best'); end
grid on; box on

subplot(2,1,2)
hold on
for j = 1:nf
    [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f);
    if m == 0 && n == 0; continue; end
    omega = omega_vec(j);
    alpha_r = real(StabRes.alpha(j,:));
    c_ph = omega ./ alpha_r;
    c_ph(abs(alpha_r) < 1e-12) = NaN;
    plot(R(xr), c_ph(xr), 'LineWidth', 1.2, 'Color', colors(j,:), ...
        'DisplayName', sprintf('(%d,%d)', m, n));
end
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$c = \omega / \alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
title('Phase speed', 'FontSize', 13)
if nf > 2; legend('Location', 'best'); end
grid on; box on

%% 3. Eigenfunctions at selected stations
x_stations = x(1) + [0.2 0.4 0.6] * (x(i_buf) - x(1));
j_fund = find_fundamental(omega_vec, beta_vec);

if ~isempty(j_fund)
    y = StabGrid.y(:);

    figure('Name', 'Eigenfunctions', 'Position', [100 100 1100 400]);
    vars   = {'u', 'v', 'w', 'p'};
    labels = {'$|\hat{u}|$', '$|\hat{v}|$', '$|\hat{w}|$', '$|\hat{p}|$'};
    for iv = 1:4
        subplot(1, 4, iv)
        hold on
        for k = 1:numel(x_stations)
            [~, ix] = min(abs(x - x_stations(k)));
            phi = StabRes.(vars{iv})(j_fund, :, ix);
            plot(abs(phi(:)), y, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('$x=%.0f$', x(ix)));
        end
        xlabel(labels{iv}, 'Interpreter', 'latex', 'FontSize', 12)
        if iv == 1; ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 12); end
        title(labels{iv}, 'Interpreter', 'latex', 'FontSize', 13)
        legend('Interpreter', 'latex', 'Location', 'best')
        grid on; box on
    end
    sgtitle(sprintf('Eigenfunctions — mode (1,0), omega = %.4g', omega_vec(j_fund)), ...
        'FontSize', 13)
end

%% Print N-factors
i_nfac = i_buf;   % report N-factors at buffer start
fprintf('\n=== N-factors ===\n');
for j = 1:nf
    [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f);
    if (m == 0 && n == 0) || ~has_valid_alpha(j); continue; end
    alpha_i = -imag(StabRes.alpha(j,:));
    N_x = cumtrapz(x, alpha_i);
    N_end = N_x(i_nfac);
    [N_max, i_max] = max(N_x(1:i_nfac));
    fprintf('  Mode (%d,%d): N(x=%.1f) = %.2f,  N_max = %.2f at x = %.1f\n', ...
        m, n, x(i_nfac), N_end, N_max, x(i_max));
end
fprintf('\n');

end

%% Local helpers

function [m, n] = mode_mn(j, omega_vec, beta_vec, omega_f, beta_f)
if ~isempty(omega_f); m = round(omega_vec(j)/omega_f); else; m = 0; end
if ~isempty(beta_f);  n = round(beta_vec(j)/beta_f);   else; n = 0; end
end

function j = find_fundamental(omega_vec, beta_vec)
j = [];
for k = 1:numel(omega_vec)
    if omega_vec(k) ~= 0 || beta_vec(k) ~= 0
        j = k; return;
    end
end
end
