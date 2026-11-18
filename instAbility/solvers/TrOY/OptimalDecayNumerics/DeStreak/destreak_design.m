function R = destreak_design(basis, n_eta, n_xi)
%DESTREAK_DESIGN  Base-flow-response inverse design of a passive de-streaking relief.
%
%   R = destreak_design({'p11','p15','p19','p23'}, 60,1200)
%
%   Uses the DirectFlatPlate base-flow sweep (gsweep_*.mat) to build the induced-
%   roller operator WITH the base-flow response: column i = (v_bump_i - v_flat)/1,
%   the beta0 roller change from bump i (U responds to the geometry). Then forms
%   the sealed-Poisson ceiling P, solves the real-c least squares for the relief
%   h = sum c_i bump_i, and reports achieved net transfer + reach vs the ceiling.
%   Also checks sign-linearity against the depression run 'd15' if present.

if nargin<2||isempty(n_eta); n_eta=60;  end
if nargin<3||isempty(n_xi);  n_xi =1200; end
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';
gname=@(id) fullfile(here,sprintf('gsweep_%s_%dx%d.mat',id,n_eta,n_xi));

%% flat baseline
Gflat = getfield(load(gname('flat')),'StabGrid');
Gfc   = destreak_perturb_grid(Gflat, zeros(Gflat.ny,Gflat.nx));
beta0 = 2*pi*Gflat.lref/7.5e-3;
[u0,v0] = destreak_run_case(Gfc, beta0, ropt);

x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1));
win = xn>=0.30 & xn<=0.72; stn=find(win); ny=Gflat.ny;
[applyP,Vopt] = destreak_Popt(Gfc,beta0,stn,u0);
r = Vopt - v0(:,stn);                              % required roller correction

%% base-flow-response columns for each basis bump
K=numel(basis); dV=zeros(ny,Gflat.nx,K); dU=zeros(ny,Gflat.nx,K);
for k=1:K
    Gk = getfield(load(gname(basis{k})),'StabGrid');
    [uk,vk]=destreak_run_case(Gk,beta0,ropt);
    p=sum(conj(uk(:)).*u0(:)); p=p/abs(p); uk=uk*p; vk=vk*p;
    dV(:,:,k)=vk-v0; dU(:,:,k)=uk-u0;
end

%% weights
W=zeros(ny,numel(stn));
for idx=1:numel(stn); W(:,idx)=abs(loc_trap(Gfc.Y(:,stn(idx)))); end
dxv=gradient(x); for idx=1:numel(stn); W(:,idx)=W(:,idx)*abs(dxv(stn(idx))); end
wqv=W(:);
dvec=reshape(-Gfc.dyU(:,stn).*u0(:,stn),[],1); gw=abs(dvec)/max(abs(dvec));
W_E=wqv; W_g=wqv.*gw;

%% assemble A c ~= r,  real c
flat=@(M) reshape(M,[],1);
Amat=zeros(ny*numel(stn),K);
for k=1:K; Amat(:,k)=flat(dV(:,stn,k) - applyP(dU(:,:,k))); end
rr=flat(r);
stack=@(M)[real(M);imag(M)];
solve_c=@(Wm)(sqrt([Wm;Wm]).*stack(Amat))\(sqrt([Wm;Wm]).*stack(rr));
c_E=solve_c(W_E); c_g=solve_c(W_g);
R2=@(c,Wm)1-sum([Wm;Wm].*stack(Amat*c-rr).^2)/sum([Wm;Wm].*stack(rr).^2);
R2_E=R2(c_E,W_E); R2_g=R2(c_g,W_g);

%% achieved fields + transfer
uA=u0; vA=v0; for k=1:K; uA=uA+c_g(k)*dU(:,:,k); vA=vA+c_g(k)*dV(:,:,k); end
Tf=xfer(Gfc,u0,v0,stn);
Tc=xfer(Gfc,u0,embedw(Vopt,stn,Gfc),stn);
Ta=xfer(Gfc,uA,vA,stn);
prog=(sum(Tf)-sum(Ta))/(sum(Tf)-sum(Tc));

%% per-bump single-mode drain cosine + sign-linearity check
ip=@(a,b) real(sum(W(:).*conj(a(:)).*b(:))); nrm=@(a)sqrt(ip(a,a));
cos1=zeros(1,K);
for k=1:K; d=dV(:,stn,k); cos1(k)=ip(d,r)/(nrm(d)*nrm(r)); end
lin='n/a';
if any(strcmp(basis,'p15')) && isfile(gname('d15'))
    Gd=getfield(load(gname('d15')),'StabGrid');
    [ud,vd]=destreak_run_case(Gd,beta0,ropt); p=sum(conj(ud(:)).*u0(:)); vd=vd*p/abs(p);
    dVd=vd-v0; kp=find(strcmp(basis,'p15'),1);
    a=dVd(:,stn); b=-dV(:,stn,kp); lin=sprintf('%.3f',ip(a,b)/(nrm(a)*nrm(b)));
end

%% report
fprintf('\n================= BASE-FLOW-RESPONSE DESIGN =================\n');
fprintf(' basis bumps: %s\n', strjoin(basis,', '));
fprintf(' per-bump cos(dV_i, r): %s\n', mat2str(cos1,3));
fprintf(' depression sign-linearity cos(dV_d15,-dV_p15) = %s  (1=linear)\n', lin);
fprintf(' reach R^2 (energy) = %.3f   (g-weighted) = %.3f\n', R2_E, R2_g);
fprintf(' c* (g-weighted) = %s\n', mat2str(c_g(:)',3));
fprintf(' sum T_flat = %+.4e [+sustain]   T_ceil = %+.4e [-drain]\n', sum(Tf), sum(Tc));
fprintf(' sum T_ach  = %+.4e   progress flat->ceiling = %.3f\n', sum(Ta), prog);
fprintf('============================================================\n');

R=struct('Gfc',Gfc,'beta0',beta0,'stn',stn,'u0',u0,'v0',v0,'Vopt',Vopt,'r',r, ...
  'dV',dV,'dU',dU,'basis',{basis},'Amat',Amat,'W_E',W_E,'W_g',W_g, ...
  'c_E',c_E,'c_g',c_g,'R2_E',R2_E,'R2_g',R2_g,'cos1',cos1,'lin',lin, ...
  'Tf',Tf,'Tc',Tc,'Ta',Ta,'prog',prog);
save(fullfile(here,sprintf('design_%dx%d.mat',n_eta,n_xi)),'R','-v7.3');
end

function w=loc_trap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
function T=xfer(G,u,v,stn)
ns=numel(stn); T=zeros(1,ns);
for k=1:ns; i=stn(k); T(k)=-sum(loc_trap(G.Y(:,i)).*(G.dyU(:,i).*real(conj(u(:,i)).*v(:,i)))); end
end
function F=embedw(Vw,stn,G); F=zeros(G.ny,G.nx); F(:,stn)=Vw; end
