%% DESTREAK_SHAPES_PREVIEW  Design forward-facing-step-like wall reliefs.
%  Replicates the DFP Julia shape functions (wall_bump_esn / wall_bump_sigmoidal)
%  to pick asymmetric protrusions: steep (rounded) leading edge + gentle trailing
%  edge, to induce near-wall rollers (reverse lift-up) without separating.
clear; clc;
set(groot,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');
A=1e-3; xc=0.16; R=50; yTol=1e-5;
x=linspace(0.10,0.24,2000);

esn = @(eps) arrayfun(@(xx) wbe(xx,A,xc,eps,R,yTol), x);
sig = @(xs,xp,xe,p,q) arrayfun(@(xx) wbs(xx,A,xs,xp,xe,p,q), x);

fig=figure('Color','w','Position',[100 100 960 380]);

subplot(1,2,1);  % ESN skew family
plot(x*1e3, esn(0.0)*1e3,'k-','LineWidth',1.6); hold on
plot(x*1e3, esn(-0.4)*1e3,'-','LineWidth',1.6,'Color',[0.85 0.33 0.1]);
plot(x*1e3, esn(-0.7)*1e3,'-','LineWidth',1.6,'Color',[0.2 0.4 0.7]);
plot(x*1e3, esn(+0.5)*1e3,'--','LineWidth',1.2,'Color',[0.4 0.6 0.3]);
xline(xc*1e3,'k:'); xlabel('$x$ [mm]'); ylabel('$h$ [mm]');
title('\textbf{(a)} ESN skew ($\epsilon<0$ = forward-facing)');
legend({'$\epsilon=0$ (sym)','$\epsilon=-0.4$','$\epsilon=-0.7$','$\epsilon=+0.5$ (bwd)'}, ...
       'Location','northeast'); xlim([120 210]); box on
text(0.02,0.98,'flow $\rightarrow$','Units','normalized','VerticalAlignment','top');

subplot(1,2,2);  % sigmoidal forward-facing step
plot(x*1e3, sig(0.150,0.160,0.200,4,3)*1e3,'-','LineWidth',1.6,'Color',[0.85 0.33 0.1]); hold on
plot(x*1e3, sig(0.150,0.158,0.210,6,3)*1e3,'-','LineWidth',1.6,'Color',[0.2 0.4 0.7]);
plot(x*1e3, sig(0.150,0.170,0.190,3,3)*1e3,'k--','LineWidth',1.2);
xlabel('$x$ [mm]'); ylabel('$h$ [mm]');
title('\textbf{(b)} sigmoidal: steep rise + gentle fall');
legend({'rise10/fall40, $p4q3$','rise8/fall52, $p6q3$','sym rise20/fall20 $p3q3$'}, ...
       'Location','northeast'); xlim([120 210]); box on

exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')),'destreak_shapes.png'),'Resolution',170);
fprintf('[shapes] saved destreak_shapes.png\n');

%% --- local replicas of the Julia shape functions ---
function y=ss(t,n); y=t.^n./(t.^n+(1-t).^n); end
function y=wbs(x,A,xStart,xPeak,xEnd,p,q)
if x<=xStart||x>=xEnd, y=0; return; end
if x<=xPeak, t=(x-xStart)/(xPeak-xStart); y=A*ss(t,p);
else, s=(x-xPeak)/(xEnd-xPeak); y=A*(1-ss(s,q)); end
end
function y=wbe(x,A,xCenter,epsilon,R,yTol)
sl=1+epsilon; sr=1-epsilon; tol=yTol/abs(A);
if tol>=1, y=0; return; end
xtl=-sl*sqrt(-2*log(tol)); xtr=sr*sqrt(-2*log(tol)); besn=xtr-xtl;
scale=(R*abs(A))/besn; bw=0.05*besn; xi=(x-xCenter)/scale;
if xi>xtl && xi<xtr
    s=(xi<0)*sl+(xi>=0)*sr; y=A*exp(-xi^2/(2*s^2));
elseif xi>=xtl-bw && xi<=xtl
    g=A*exp(-xi^2/(2*sl^2)); t=(xi-(xtl-bw))/bw; y=g*ss(t,3);
elseif xi>=xtr && xi<=xtr+bw
    g=A*exp(-xi^2/(2*sr^2)); t=(xi-xtr)/bw; y=g*(1-ss(t,3));
else, y=0; end
end
