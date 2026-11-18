function destreak_rescore(n_eta, n_xi)
%DESTREAK_RESCORE  Re-score the existing flat/hump base flows with the CORRECTED
%   base-flow-oriented I2 (Casacuberta et al. 2024), replacing the parallel
%   reduction. No new CFD — reuses gflat_*/ghump_* grids on disk.

if nargin<1||isempty(n_eta); n_eta=60;  end
if nargin<2||isempty(n_xi);  n_xi =1200; end
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';
gn=@(id) fullfile(here,sprintf('g%s_%dx%d.mat',id,n_eta,n_xi));

Gflat=getfield(load(gn('flat')),'StabGrid');
Gfc  =destreak_perturb_grid(Gflat, zeros(Gflat.ny,Gflat.nx));
Gh   =getfield(load(gn('hump')),'StabGrid');
beta0=2*pi*Gflat.lref/7.5e-3;

[u0,v0,w0]=destreak_run_case(Gfc,beta0,ropt);
[uh,vh,wh]=destreak_run_case(Gh, beta0,ropt);
p=sum(conj(uh(:)).*u0(:)); p=p/abs(p); uh=uh*p; vh=vh*p; wh=wh*p;

x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1));
win=xn>=0.30 & xn<=0.72; stn=find(win); ns=numel(stn);

% --- base-flow tilt (why parallel fails): mean flow angle from chordwise ---
ang = atan2d(hypot(Gfc.V(:,stn),Gfc.W(:,stn)), Gfc.U(:,stn));
angW= atan2d(Gfc.W(:,stn), Gfc.U(:,stn));      % crossflow (w/u) angle
fprintf('[rescore] base-flow tilt in window: median |v,w|/u angle=%.1f deg, crossflow atan(W/U)=%.1f deg\n', ...
        median(ang(:)), median(angW(:)));

% --- corrected base-oriented I2 ---
[I2f,of]=destreak_i2base(Gfc,u0,v0,w0,stn);
[I2h,oh]=destreak_i2base(Gh, uh,vh,wh,stn);
chk=max(abs([I2f-of.I2_check, I2h-oh.I2_check]));
fprintf('[rescore] I2 consistency <vt,Lvn> vs <g,vn>: max diff=%.2e\n', chk);

% --- parallel (old) transfer for comparison ---
Tf=zeros(1,ns); Th=zeros(1,ns);
for k=1:ns; i=stn(k); wq=loctrap(Gfc.Y(:,i));
  Tf(k)=-sum(wq.*Gfc.dyU(:,i).*real(conj(u0(:,i)).*v0(:,i)));
  Th(k)=-sum(wq.*Gh.dyU(:,i) .*real(conj(uh(:,i)).*vh(:,i))); end

% --- corrected drain alignment: induced roller change vs full g ---
W=zeros(Gfc.ny,ns); dxv=gradient(x);
for k=1:ns; W(:,k)=abs(loctrap(Gfc.Y(:,stn(k))))*abs(dxv(stn(k))); end
sw=sqrt([W(:);W(:)]);
gvec =[of.gx(:); of.gy(:)];
dvn  =[(oh.vnx(:)-of.vnx(:)); (oh.vny(:)-of.vny(:))];   % induced cross-stream change
ip=@(a,b) real(sum(sw.^2.*conj(a).*b)); nrm=@(a)sqrt(ip(a,a));
cosC = ip(gvec,dvn)/(nrm(gvec)*nrm(dvn));

% parallel alignment (old cos, for reference): dv_parallel vs r_parallel
% r_par ~ (P u0 - v0); reuse scalar Popt
[~,Vopt]=destreak_Popt(Gfc,beta0,stn,u0);
rpar=Vopt-v0(:,stn); dvpar=vh(:,stn)-v0(:,stn);
ipp=@(a,b) real(sum(W(:).*conj(a(:)).*b(:))); nrmp=@(a)sqrt(ipp(a,a));
cosP=ipp(dvpar,rpar)/(nrmp(dvpar)*nrmp(rpar));

fprintf('\n=============== CORRECTED (base-oriented) RE-SCORE ===============\n');
fprintf(' sum I2_flat  (corrected) = %+.4e\n', sum(I2f));
fprintf(' sum I2_hump  (corrected) = %+.4e   (I2<0 drains)\n', sum(I2h));
fprintf(' sum T_flat   (parallel)  = %+.4e\n', sum(Tf));
fprintf(' sum T_hump   (parallel)  = %+.4e\n', sum(Th));
fprintf(' drain alignment  cos(dvn'', g)  CORRECTED = %+.3f\n', cosC);
fprintf(' drain alignment  cos(dv, r)     parallel  = %+.3f\n', cosP);
fprintf('=================================================================\n');

save(fullfile(here,sprintf('rescore_%dx%d.mat',n_eta,n_xi)), ...
     'I2f','I2h','Tf','Th','cosC','cosP','of','oh','stn','-v7.3');
end

function w=loctrap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
