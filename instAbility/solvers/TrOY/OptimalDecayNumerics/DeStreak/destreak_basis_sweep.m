%% DESTREAK_BASIS_SWEEP  Is the ~10%% drain reach mechanism-limited or basis-limited?
%  Build a RICH relief basis (many narrow bumps) on the cheap grid, then measure
%  how the reachable fraction R^2 of the required drain roller r=(P u_hom - v_hom)
%  grows with basis size K and with the wall-normal taper s0. Saturation => the
%  frozen-base-flow geometric channel is intrinsically weak (mechanism-limited).

clear; clc;
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));

S = load(fullfile(here,'flat_grid_40x300.mat'));  G = S.StabGrid;
beta0 = 2*pi*G.lref/7.5e-3;
Gc = destreak_perturb_grid(G, zeros(G.ny,G.nx));
ropt.lu_mode = 'none';
[u0,v0] = destreak_run_case(Gc, beta0, ropt);

% window + rich set of narrow bump centres
x = Gc.x(:).'; xn=(x-x(1))/(x(end)-x(1));
win = xn>=0.30 & xn<=0.72;  stn=find(win);  ny=Gc.ny;
Krich = 24;
xc_n = linspace(0.33,0.69,Krich);
sig  = 0.6*(xc_n(2)-xc_n(1))*(x(end)-x(1));    % narrow, ~half the spacing
xc   = x(1)+xc_n*(x(end)-x(1));

% quadrature weight (energy) over window
Wc = zeros(ny,numel(stn));
for idx=1:numel(stn); Wc(:,idx)=abs(local_trap(Gc.Y(:,stn(idx)))); end
sw = sqrt(Wc(:));

% required correction r (taper-independent)
[~,Vopt] = destreak_Popt(Gc,beta0,stn,u0);
r = Vopt - v0(:,stn);  rr=[real(sw.*r(:));imag(sw.*r(:))];

tapers = [0.08 0.18 0.35];   % near-wall / medium / broad
Ylist  = [4 8 12 16 20 24];
R2tab  = zeros(numel(tapers), numel(Ylist));

for it=1:numel(tapers)
    s0 = tapers(it);
    Yc = Gc.Y(:,round(Gc.nx/2)); s=(Yc-Yc(ny))/(Yc(1)-Yc(ny));
    phi = exp(-(s/s0).^2); phi(end)=1;
    epsA = 5e-4*(Yc(1)-Yc(ny));
    Vmat = zeros(ny*numel(stn), Krich);
    for k=1:Krich
        g = exp(-((x-xc(k))/sig).^2).'; g=g/max(g);
        DY = epsA*(phi*g.');
        Gk = destreak_perturb_grid(Gc, DY);
        [uk,vk] = destreak_run_case(Gk,beta0,ropt);
        phk=sum(conj(uk(:)).*u0(:)); vk=vk*phk/abs(phk);
        Vk=(vk-v0)/epsA; Vmat(:,k)=reshape(Vk(:,stn),[],1);
    end
    for iK=1:numel(Ylist)
        sel = round(linspace(1,Krich,Ylist(iK)));
        Vs  = Vmat(:,sel);
        Ar  = [real(sw.*Vs); imag(sw.*Vs)];
        c   = Ar\rr;  R2tab(it,iK)=1-norm(Ar*c-rr)^2/norm(rr)^2;
    end
    fprintf('[sweep] s0=%.2f : R2(K=%s) = %s\n', s0, mat2str(Ylist), mat2str(R2tab(it,:),3));
end

save(fullfile(here,'basis_sweep.mat'),'R2tab','tapers','Ylist');
fprintf('[sweep] DONE\n');

function w = local_trap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
