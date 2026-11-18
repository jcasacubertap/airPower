function R = destreak_design_base(basis, n_eta, n_xi)
%DESTREAK_DESIGN_BASE  Corrected (base-flow-oriented) inverse design of a passive
%   de-streaking relief. Uses the DirectFlatPlate base-flow sweep and the
%   base-oriented I2 / representer g / general ceiling P (Casacuberta et al. 2024).
%
%   R = destreak_design_base({'p11','p15','p19','p23'}, 60,1200)
%
%   Skeleton (brief): A = V - P S, c* = (A'A)^{-1} A' r, but every operator is
%   base-oriented: V = induced cross-stream roller change Δvn' (base-flow response
%   from each bump), S = tangential change Δvt', P = general optimal roller. The
%   objective I2 is scored base-oriented; achieved I2(c) is the linear
%   superposition of per-bump ΔI2 (validate later with a CFD run of the relief).

if nargin<2||isempty(n_eta); n_eta=60;  end
if nargin<3||isempty(n_xi);  n_xi =1200; end
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';
gn=@(id) fullfile(here,sprintf('gsweep_%s_%dx%d.mat',id,n_eta,n_xi));

Gflat=getfield(load(gn('flat')),'StabGrid');
Gfc  =destreak_perturb_grid(Gflat, zeros(Gflat.ny,Gflat.nx));
beta0=2*pi*Gflat.lref/7.5e-3;
[u0,v0,w0]=destreak_run_case(Gfc,beta0,ropt);
x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1));
win=xn>=0.30 & xn<=0.72; stn=find(win); ns=numel(stn); ny=Gflat.ny;

[I2f,of]=destreak_i2base(Gfc,u0,v0,w0,stn);  sumI2f=sum(I2f);
applyP=destreak_Popt_base(Gfc,beta0,stn);
vnopt_f=applyP(struct('vtx',of.vtx,'vty',of.vty,'vtz',of.vtz));   % ceiling roller

%% per-bump base-flow-response columns
K=numel(basis);
dvnx=zeros(ny,ns,K); dvny=zeros(ny,ns,K); Pd=cell(1,K); dI2=zeros(1,K); cos1=zeros(1,K);
dI2x=zeros(ns,K);
% energy weight
W1=zeros(ny,ns); dxv=gradient(x);
for k=1:ns; W1(:,k)=abs(loctrap(Gfc.Y(:,stn(k))))*abs(dxv(stn(k))); end
w1=W1(:);
ipg=@(a,b) real(sum([w1;w1].*conj(a).*b)); nrmg=@(a)sqrt(ipg(a,a));
gvec=[of.gx(:); of.gy(:)];
for k=1:K
    Gk=getfield(load(gn(basis{k})),'StabGrid');
    [uk,vk,wk]=destreak_run_case(Gk,beta0,ropt);
    p=sum(conj(uk(:)).*u0(:)); p=p/abs(p); uk=uk*p; vk=vk*p; wk=wk*p;
    [I2k,ok]=destreak_i2base(Gk,uk,vk,wk,stn); dI2(k)=sum(I2k)-sumI2f; dI2x(:,k)=(I2k-I2f).';
    dvnx(:,:,k)=ok.vnx-of.vnx; dvny(:,:,k)=ok.vny-of.vny;
    Pd{k}=applyP(struct('vtx',ok.vtx-of.vtx,'vty',ok.vty-of.vty,'vtz',ok.vtz-of.vtz));
    dvn=[reshape(dvnx(:,:,k),[],1); reshape(dvny(:,:,k),[],1)];
    cos1(k)=ipg(dvn,gvec)/(nrmg(dvn)*nrmg(gvec));
end

%% assemble A c ~= r  (stacked x,y window comps), real c
flat=@(M) reshape(M,[],1);
r=[flat(vnopt_f.x - of.vnx); flat(vnopt_f.y - of.vny)];
A=zeros(2*ny*ns,K);
for k=1:K; A(:,k)=[flat(dvnx(:,:,k)-Pd{k}.x); flat(dvny(:,:,k)-Pd{k}.y)]; end
Wd=[w1;w1];
gm=sqrt(abs(of.gx).^2+abs(of.gy).^2); gm=gm(:)/max(gm(:)); Wg=[w1.*gm; w1.*gm];
stack=@(M)[real(M);imag(M)];
solve_c=@(Wm)(sqrt([Wm;Wm]).*stack(A))\(sqrt([Wm;Wm]).*stack(r));
c_E=solve_c(Wd); c_g=solve_c(Wg);
R2=@(c,Wm)1-sum([Wm;Wm].*stack(A*c-r).^2)/sum([Wm;Wm].*stack(r).^2);
R2_E=R2(c_E,Wd); R2_g=R2(c_g,Wg);

%% achieved I2 (linear superposition of per-bump changes) + ceiling
I2_ach = sumI2f + sum(c_g(:)'.*dI2);
I2ceil_x = -ceil_transfer(of,vnopt_f,Gfc,stn);        % per-station ceiling
I2_ceil = sum(I2ceil_x);
I2ach_x = I2f + (dI2x*c_g(:)).';                       % per-station achieved (linear)
prog=(sumI2f - I2_ach)/(sumI2f - I2_ceil);

%% sign-linearity check with depression
lin='n/a';
if any(strcmp(basis,'p15')) && isfile(gn('d15'))
    Gd=getfield(load(gn('d15')),'StabGrid');
    [ud,vd,wd]=destreak_run_case(Gd,beta0,ropt); p=sum(conj(ud(:)).*u0(:)); ud=ud*p/abs(p); vd=vd*p/abs(p); wd=wd*p/abs(p);
    [~,od]=destreak_i2base(Gd,ud,vd,wd,stn); kp=find(strcmp(basis,'p15'),1);
    a=[flat(od.vnx-of.vnx); flat(od.vny-of.vny)];
    b=-[flat(dvnx(:,:,kp)); flat(dvny(:,:,kp))];
    lin=sprintf('%.3f', ipg(a,b)/(nrmg(a)*nrmg(b)));
end

fprintf('\n=========== CORRECTED BASE-ORIENTED DESIGN ===========\n');
fprintf(' basis: %s\n', strjoin(basis,', '));
fprintf(' per-bump cos(dvn_i, g): %s\n', mat2str(cos1,3));
fprintf(' per-bump dI2 (each drains if <0): %s\n', mat2str(dI2,3));
fprintf(' depression sign-linearity = %s\n', lin);
fprintf(' reach R^2 (E)=%.3f  (g)=%.3f\n', R2_E, R2_g);
fprintf(' c* (g-weighted) = %s\n', mat2str(c_g(:)',3));
fprintf(' sum I2_flat = %+.4e   I2_ceil = %+.4e   I2_ach = %+.4e\n', sumI2f, I2_ceil, I2_ach);
fprintf(' progress flat->ceiling = %.3f\n', prog);
fprintf('======================================================\n');

R=struct('Gfc',Gfc,'beta0',beta0,'stn',stn,'basis',{basis},'of',of, ...
  'vnopt_f',vnopt_f,'A',A,'r',r,'c_E',c_E,'c_g',c_g,'R2_E',R2_E,'R2_g',R2_g, ...
  'cos1',cos1,'dI2',dI2,'lin',lin,'sumI2f',sumI2f,'I2_ceil',I2_ceil,'I2_ach',I2_ach,'prog',prog, ...
  'I2f_x',I2f,'I2ach_x',I2ach_x,'I2ceil_x',I2ceil_x,'xw',xn(stn));
save(fullfile(here,sprintf('design_base_%dx%d.mat',n_eta,n_xi)),'R','-v7.3');
end

function T=ceil_transfer(of,vnopt,G,stn)
% <g, vn_opt> per window station (trapezoid in y)
ns=numel(stn); T=zeros(1,ns);
for k=1:ns
    wq=loctrap(G.Y(:,stn(k)));
    T(k)=sum(wq.*real(conj(of.gx(:,k)).*vnopt.x(:,k)+conj(of.gy(:,k)).*vnopt.y(:,k)));
end
end
function w=loctrap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
