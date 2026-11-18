%% DESTREAK_PLOT  Feasibility summary figure (frozen-base-flow inner one-shot).
clear; clc;
here = fileparts(mfilename('fullpath'));
L = load(fullfile(here,'feas_60x1200.mat'));  R = L.R;
Sw = load(fullfile(here,'basis_sweep.mat'));
G = R.G;  stn = R.stn;  K = numel(R.c_g);

set(groot,'defaultTextInterpreter','latex', ...
          'defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');

% achieved roller field
vA = R.v0;  uA = R.u0;
for k=1:K; vA = vA + R.c_g(k)*R.Vfull(:,:,k); uA = uA + R.c_g(k)*R.Sfull(:,:,k); end

x  = G.x(:).';  xn = (x - x(1))/(x(end)-x(1));    % chordwise (window arclength)
xw = xn(stn);

fig = figure('Color','w','Position',[100 100 980 760]);

%% (a) singular values of the reachable roller subspace
subplot(2,2,1);
sv = R.sv/max(R.sv);
bar(1:numel(sv), sv, 0.6, 'FaceColor',[0.2 0.4 0.7]); hold on
xlabel('mode index $k$'); ylabel('$\sigma_k/\sigma_1$');
title('\textbf{(a)} reachable roller subspace $V$');
xlim([0.4 numel(sv)+0.6]); ylim([0 1.05]); box on

%% (b) drain reach R^2 vs basis size K (mechanism vs basis limited)
subplot(2,2,2);
mk = {'-o','-s','-^'};
for it=1:numel(Sw.tapers)
    plot(Sw.Ylist, Sw.R2tab(it,:), mk{it}, 'LineWidth',1.3, ...
         'MarkerFaceColor','w'); hold on
end
xlabel('number of relief modes $K$'); ylabel('drain reach $R^2$');
title('\textbf{(b)} reach saturates $\Rightarrow$ mechanism-limited');
legend(arrayfun(@(s) sprintf('$s_0=%.2f$',s), Sw.tapers,'uni',0), ...
       'Location','northwest');
xlim([0 26]); ylim([0 0.5]); box on

%% (c) transfer T(x) over the window: flat / ceiling / achieved
subplot(2,2,3);
plot(xw, R.T_flat, '-',  'LineWidth',1.5,'Color',[0.85 0.33 0.10]); hold on
plot(xw, R.T_ach,  '--', 'LineWidth',1.5,'Color',[0.2 0.4 0.7]);
plot(xw, R.T_ceil, '-.', 'LineWidth',1.5,'Color',[0.1 0.6 0.3]);
yline(0,'k:');
xlabel('chordwise $\xi/L$'); ylabel('$T=-I_2$');
title('\textbf{(c)} transfer: $T<0$ drains, $T>0$ sustains');
legend({'flat (self-roller)','achieved $c^\ast$','sealed-Poisson ceiling'}, ...
       'Location','east');
xlim([min(xw) max(xw)]); box on

%% (d) roller profiles at a representative window station
subplot(2,2,4);
idx = round(numel(stn)/2);  i = stn(idx);
Y = G.Y(:,i);
plot(abs(R.v0(:,i)),       Y, '-',  'LineWidth',1.5,'Color',[0.85 0.33 0.10]); hold on
plot(abs(vA(:,i)),         Y, '--', 'LineWidth',1.5,'Color',[0.2 0.4 0.7]);
plot(abs(R.Vopt_hom(:,idx)),Y,'-.', 'LineWidth',1.5,'Color',[0.1 0.6 0.3]);
xlabel('$|\hat v|$'); ylabel('$y/\ell$');
title(sprintf('\\textbf{(d)} roller at $\\xi/L=%.2f$', xw(idx)));
legend({'$|\hat v_{\rm hom}|$','$|\hat v(c^\ast)|$','$|\hat v_{\rm opt}|$'}, ...
       'Location','northeast');
ylim([0 8]); box on

sgtitle(sprintf(['Passive de-streaking, frozen base flow (60$\\times$1200): ' ...
    'drain reach $R^2=%.3f$, progress $=%.1f\\%%$'], R.R2_E, 100*R.progress), ...
    'Interpreter','latex');

exportgraphics(fig, fullfile(here,'destreak_feasibility.png'), 'Resolution',180);
fprintf('[plot] saved destreak_feasibility.png\n');
