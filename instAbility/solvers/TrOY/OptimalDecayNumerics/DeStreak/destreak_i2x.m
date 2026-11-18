function destreak_i2x(ids, n_eta, n_xi)
%DESTREAK_I2X  Streamwise I2(x) profiles for the shape family (local drain vs sum).
if nargin<2||isempty(n_eta); n_eta=60; end
if nargin<3||isempty(n_xi); n_xi=1200; end
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';
gn=@(id) fullfile(here,sprintf('gsweep_%s_%dx%d.mat',id,n_eta,n_xi));
Gflat=getfield(load(gn('flat')),'StabGrid');
Gfc=destreak_perturb_grid(Gflat, zeros(Gflat.ny,Gflat.nx));
beta0=2*pi*Gflat.lref/7.5e-3;
[u0,v0,w0]=destreak_run_case(Gfc,beta0,ropt);
x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1));
rng = xn>=0.33 & xn<=0.78; stn=find(rng); xw=xn(stn);   % step region + downstream
[I2f,~]=destreak_i2base(Gfc,u0,v0,w0,stn);

set(groot,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');
fig=figure('Color','w','Position',[100 100 900 400]); cols=lines(numel(ids));
for j=1:numel(ids)
    id=ids{j};
    if strcmp(id,'flat'); I2=I2f;
    else; G=getfield(load(gn(id)),'StabGrid'); [u,v,w]=destreak_run_case(G,beta0,ropt);
        p=sum(conj(u(:)).*u0(:)); p=p/abs(p); u=u*p;v=v*p;w=w*p;
        [I2,~]=destreak_i2base(G,u,v,w,stn); end
    plot(xw, I2, '-','LineWidth',1.6,'Color',cols(j,:)); hold on
end
yline(0,'k:'); xline(0.44,'k--');   % ~step location
xlabel('chordwise $\xi/L$'); ylabel('$I_2(x)$  (base-oriented; $<0$ drains)');
title('Streamwise $I_2$: step vs symmetric vs flat'); legend(ids,'Location','best');
xlim([min(xw) max(xw)]); box on
exportgraphics(fig, fullfile(here,'destreak_i2x.png'),'Resolution',170);
fprintf('[i2x] saved. per-shape min I2(x) and its location:\n');
end
