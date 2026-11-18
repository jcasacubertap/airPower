%% DESTREAK_PLOT_BFRESP  Base-flow response vs frozen-metric channel (the verdict).
clear; clc;
here = fileparts(mfilename('fullpath'));
R = getfield(load(fullfile(here,'bfresp_60x1200.mat')),'R');
G = R.Gf;  stn = R.stn;
set(groot,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');

x = G.x(:).'; xn=(x-x(1))/(x(end)-x(1)); xw=xn(stn);
Vopt = R.r + R.v0(:,stn);                         % ceiling roller over window

fig = figure('Color','w','Position',[100 100 1000 340]);

%% (a) single-mode drain reach: frozen vs full response
subplot(1,3,1);
b = bar([1 2],[R.cos_frozen^2, R.cos_full^2],0.55); b.FaceColor='flat';
b.CData(1,:)=[0.75 0.4 0.2]; b.CData(2,:)=[0.15 0.5 0.35];
set(gca,'XTickLabel',{'frozen','full resp.'});
text(1,R.cos_frozen^2+0.03,sprintf('%.3f',R.cos_frozen^2),'HorizontalAlignment','center','Interpreter','latex');
text(2,R.cos_full^2+0.03,sprintf('%.2f',R.cos_full^2),'HorizontalAlignment','center','Interpreter','latex');
ylabel('single-bump drain reach $\cos^2(\Delta\hat v, r)$');
title('\textbf{(a)} base-flow response drains');
ylim([0 0.9]); box on

%% (b) roller profiles at mid-window: flat / frozen / full / optimal
subplot(1,3,2);
idx=round(numel(stn)/2); i=stn(idx); Y=G.Y(:,i);
plot(abs(R.v0(:,i)),            Y,'-','LineWidth',1.5,'Color',[0.85 0.33 0.10]); hold on
plot(abs(R.v0(:,i)+R.dV_frozen(:,idx)),Y,':','LineWidth',1.5,'Color',[0.75 0.4 0.2]);
plot(abs(R.v0(:,i)+R.dV_full(:,idx)),  Y,'--','LineWidth',1.6,'Color',[0.15 0.5 0.35]);
plot(abs(Vopt(:,idx)),          Y,'-.','LineWidth',1.5,'Color',[0.1 0.2 0.7]);
xlabel('$|\hat v|$'); ylabel('$y/\ell$');
title(sprintf('\\textbf{(b)} roller at $\\xi/L=%.2f$',xw(idx)));
legend({'flat','+bump (frozen)','+bump (full resp.)','optimal (ceiling)'}, ...
       'Location','northeast'); ylim([0 8]); box on

%% (c) transfer over the window
subplot(1,3,3);
plot(xw,R.Tf,'-','LineWidth',1.4,'Color',[0.85 0.33 0.10]); hold on
plot(xw,R.Th,'--','LineWidth',1.4,'Color',[0.15 0.5 0.35]);
plot(xw,R.Tc,'-.','LineWidth',1.4,'Color',[0.1 0.2 0.7]); yline(0,'k:');
xlabel('chordwise $\xi/L$'); ylabel('$T=-I_2$');
title('\textbf{(c)} transfer (protrusion amplifies)');
legend({'flat','+bump (full)','ceiling'},'Location','east');
xlim([min(xw) max(xw)]); box on

sgtitle(sprintf(['Passive de-streaking: geometry alone $\\cos^2=%.3f$, ' ...
  'base-flow response $\\cos^2=%.2f$ (60$\\times$1200)'], ...
  R.cos_frozen^2, R.cos_full^2),'Interpreter','latex');
exportgraphics(fig, fullfile(here,'destreak_bfresp.png'),'Resolution',180);
fprintf('[plot] saved destreak_bfresp.png\n');
