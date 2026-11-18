function R = destreak_score_shapes(ids, n_eta, n_xi)
%DESTREAK_SCORE_SHAPES  Corrected-I2 + separation + near-wall-roller scoring of a
%   shape family (symmetric vs forward-facing-step protrusions).
%
%   R = destreak_score_shapes({'flat','p16','e05','e07'}, 60,1200)
%
%   For each shape: corrected base-oriented sum I2 and dI2 vs flat; a separation
%   flag (min near-wall base-flow U over the window; <0 = reversed/separated); and
%   the induced near-wall roller's alignment with the drain g and its near-wall
%   energy fraction. Tests whether forward-facing asymmetry drains more (dI2<0)
%   while staying attached, and whether the productive roller is near-wall.

if nargin<2||isempty(n_eta); n_eta=60;  end
if nargin<3||isempty(n_xi);  n_xi =1200; end
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';
gn=@(id) fullfile(here,sprintf('gsweep_%s_%dx%d.mat',id,n_eta,n_xi));

% flat reference
Gflat=getfield(load(gn('flat')),'StabGrid');
Gfc  =destreak_perturb_grid(Gflat, zeros(Gflat.ny,Gflat.nx));
beta0=2*pi*Gflat.lref/7.5e-3;
[u0,v0,w0]=destreak_run_case(Gfc,beta0,ropt);
x=Gflat.x(:).'; xn=(x-x(1))/(x(end)-x(1)); win=xn>=0.30&xn<=0.72; stn=find(win); ny=Gfc.ny;
[I2f,of]=destreak_i2base(Gfc,u0,v0,w0,stn); sumI2f=sum(I2f);

% near-wall band (rows near the wall = last rows) and quad weight
nwb = (ny-round(ny/4)) : ny;         % near-wall quarter
W=zeros(ny,numel(stn)); for k=1:numel(stn); W(:,k)=abs(loctrap(Gfc.Y(:,stn(k)))); end
Wd=[W(:);W(:)];                                   % weights for stacked (x,y) fields
ip=@(a,b)real(sum(Wd.*conj(a(:)).*b(:))); nrm=@(a)sqrt(ip(a,a));
gvec=[of.gx(:);of.gy(:)];

rows={};
for j=1:numel(ids)
    id=ids{j};
    if strcmp(id,'flat'); G=Gfc; u=u0;v=v0;w=w0; oo=of; sI2=sumI2f;
    else
        G=getfield(load(gn(id)),'StabGrid');
        [u,v,w]=destreak_run_case(G,beta0,ropt);
        p=sum(conj(u(:)).*u0(:)); p=p/abs(p); u=u*p; v=v*p; w=w*p;
        [I2i,oo]=destreak_i2base(G,u,v,w,stn); sI2=sum(I2i);
    end
    % separation: min near-wall U over the window (base flow rides with mesh -> U
    % values frozen, so check the ACTUAL bump grid's interpolated U near wall)
    Uwin=G.U(nwb,stn); sepU=min(Uwin(:));
    % induced near-wall roller vs flat
    if strcmp(id,'flat'); dvn=[zeros(ny,numel(stn));zeros(ny,numel(stn))]; cosg=NaN; nwf=NaN;
    else
        dvnx=oo.vnx-of.vnx; dvny=oo.vny-of.vny;
        dvn=[dvnx(:);dvny(:)]; cosg=ip(dvn,gvec)/(nrm(dvn)*nrm(gvec));
        enw=sum(W(nwb,:).*(abs(dvnx(nwb,:)).^2+abs(dvny(nwb,:)).^2),'all');
        etot=sum(W.*(abs(dvnx).^2+abs(dvny).^2),'all');
        nwf=enw/etot;                  % near-wall energy fraction of induced roller
    end
    rows(end+1,:)={id, sI2, sI2-sumI2f, sepU, cosg, nwf}; %#ok<AGROW>
    fprintf('%-5s  sumI2=%+8.3f  dI2=%+8.3f  minU_nw=%+7.4f%s  cos(dvn,g)=%+.3f  nwFrac=%.2f\n', ...
        id, sI2, sI2-sumI2f, sepU, ternary(sepU<0,' [SEPARATED]','           '), cosg, nwf);
end
R=struct('ids',{ids},'table',{rows},'sumI2f',sumI2f,'stn',stn);
save(fullfile(here,sprintf('shapes_score_%dx%d.mat',n_eta,n_xi)),'R','-v7.3');
end

function s=ternary(c,a,b); if c, s=a; else, s=b; end; end
function w=loctrap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
