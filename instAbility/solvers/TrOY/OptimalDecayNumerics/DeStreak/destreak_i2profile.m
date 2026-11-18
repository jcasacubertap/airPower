function destreak_i2profile(ids, n_eta, n_xi)
%DESTREAK_I2PROFILE  Wall-normal-resolved I2: isolate the near-wall roller from
%   the (over-predicted) separated region on top. Cumulative I2 from the wall out.
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
x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1)); win=xn>=0.30&xn<=0.72; stn=find(win); ny=Gfc.ny;
Yrep=mean(Gfc.Y(:,stn),2);            % representative wall-normal coord (η rows)

set(groot,'defaultTextInterpreter','latex','defaultLegendInterpreter','latex', ...
          'defaultAxesTickLabelInterpreter','latex');
fig=figure('Color','w','Position',[100 100 900 380]); cols=lines(numel(ids));
yb=3.0;                                % near-wall band threshold [lref]
fprintf('%-5s  I2_total   I2_nearwall(y<%.1f)  I2_outer\n','id',yb);
for j=1:numel(ids)
    id=ids{j};
    if strcmp(id,'flat'); G=Gfc; u=u0;v=v0;w=w0;
    else; G=getfield(load(gn(id)),'StabGrid'); [u,v,w]=destreak_run_case(G,beta0,ropt);
        p=sum(conj(u(:)).*u0(:)); p=p/abs(p); u=u*p;v=v*p;w=w*p; end
    [~,o]=destreak_i2base(G,u,v,w,stn);
    % per (row,station) I2 contribution; sum over window stations -> profile(row)
    prof=zeros(ny,1);
    for k=1:numel(stn)
        wq=loctrap(G.Y(:,stn(k)));
        c=-wq.*real(conj(o.gx(:,k)).*o.vnx(:,k)+conj(o.gy(:,k)).*o.vny(:,k));
        prof=prof+c;
    end
    % cumulative from wall (row ny) outward to freestream (row 1)
    cumw=flipud(cumsum(flipud(prof)));      % cumw(row)=sum from row..ny
    I2tot=sum(prof); nw=Yrep<yb; I2nw=sum(prof(nw)); I2out=I2tot-I2nw;
    fprintf('%-5s  %+8.3f   %+8.3f            %+8.3f\n', id, I2tot, I2nw, I2out);
    plot(cumw, Yrep, '-', 'LineWidth',1.6,'Color',cols(j,:)); hold on
end
yline(yb,'k:'); xline(0,'k:');
xlabel('cumulative $I_2$ from wall $\to y$'); ylabel('$y/\ell$ (wall at 0)');
title('Wall-normal-resolved $I_2$ (left of 0 = draining)');
legend(ids,'Location','southeast'); ylim([0 12]); box on
exportgraphics(fig, fullfile(here,'destreak_i2profile.png'),'Resolution',170);
fprintf('[i2profile] saved destreak_i2profile.png\n');
end
function w=loctrap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
