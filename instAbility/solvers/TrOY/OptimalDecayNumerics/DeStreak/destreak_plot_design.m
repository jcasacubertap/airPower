function destreak_plot_design(n_eta, n_xi)
%DESTREAK_PLOT_DESIGN  Summary figure for the corrected base-oriented design.
if nargin<1||isempty(n_eta); n_eta=60;  end
if nargin<2||isempty(n_xi);  n_xi =1200; end
here=fileparts(mfilename('fullpath'));
R=getfield(load(fullfile(here,sprintf('design_base_%dx%d.mat',n_eta,n_xi))),'R');
set(groot,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');
K=numel(R.basis);
fig=figure('Color','w','Position',[100 100 1050 340]);

%% (a) per-bump: does each bump drain (dI2<0) and align with g?
subplot(1,3,1);
yyaxis left;  b=bar(1:K, R.dI2, 0.5); b.FaceColor=[0.2 0.4 0.7];
ylabel('$\Delta I_2$ per bump  ($<0$ drains)');
yyaxis right; plot(1:K, R.cos1,'o-','LineWidth',1.4,'MarkerFaceColor','w','Color',[0.85 0.33 0.1]);
ylabel('$\cos(\Delta v_n'', g)$');
set(gca,'XTick',1:K,'XTickLabel',R.basis); xtickangle(0);
title('\textbf{(a)} per-bump base-flow response'); box on

%% (b) transfer I2(x): flat / achieved / ceiling
subplot(1,3,2);
plot(R.xw,R.I2f_x,'-','LineWidth',1.5,'Color',[0.85 0.33 0.1]); hold on
plot(R.xw,R.I2ach_x,'--','LineWidth',1.5,'Color',[0.2 0.4 0.7]);
plot(R.xw,R.I2ceil_x,'-.','LineWidth',1.5,'Color',[0.1 0.6 0.3]); yline(0,'k:');
xlabel('chordwise $\xi/L$'); ylabel('$I_2$  (base-oriented)');
title('\textbf{(b)} transfer: $I_2<0$ drains');
legend({'flat','designed $c^\ast$','ceiling $\|\Pi g\|$'},'Location','best');
xlim([min(R.xw) max(R.xw)]); box on

%% (c) design coefficients + scores
subplot(1,3,3);
b=bar(1:K, R.c_g, 0.5); b.FaceColor=[0.4 0.4 0.45];
set(gca,'XTick',1:K,'XTickLabel',R.basis);
ylabel('$c^\ast_i$ (relief weight, signed)');
title('\textbf{(c)} designed relief'); box on
txt=sprintf(['reach $R^2=%.2f$\nprogress $=%.1f\\%%$\n' ...
  '$\\sum I_2$: flat $%.2f$\n$\\to$ $c^\\ast$ $%.2f$\nceiling $%.2f$\nsign-lin $%s$'], ...
  R.R2_g, 100*R.prog, R.sumI2f, R.I2_ach, R.I2_ceil, R.lin);
xl=xlim; yl=ylim;
text(xl(1)+0.05*diff(xl), yl(2)-0.03*diff(yl), txt, 'Interpreter','latex', ...
     'VerticalAlignment','top','FontSize',9);

sgtitle(sprintf('Corrected base-oriented de-streaking design (%d bumps, %d$\\times$%d)', ...
  K, n_eta, n_xi),'Interpreter','latex');
exportgraphics(fig, fullfile(here,'destreak_design.png'),'Resolution',180);
fprintf('[plot] saved destreak_design.png\n');
end
