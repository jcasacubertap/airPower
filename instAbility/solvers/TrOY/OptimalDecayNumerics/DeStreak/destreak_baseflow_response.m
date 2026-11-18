function R = destreak_baseflow_response(n_eta, n_xi)
%DESTREAK_BASEFLOW_RESPONSE  Does the base-flow response drain, where geometry alone cannot?
%
%   Compares, for the SAME physical bump (the DirectFlatPlate ESN bump, already
%   on disk as bf_sweptwing_hump.csv), the induced roller from two channels:
%     - frozen-metric : bend the flat mesh by the bump's wall shape, base flow
%                       frozen (my V channel) -> Delta v_frozen
%     - full response : run the streak on the ACTUAL bumped base flow (U responds)
%                       -> Delta v_full = v_hump - v_flat
%   and measures each against the drain direction r = (P u_flat - v_flat).
%   If the full-response roller is much larger / better drain-aligned than the
%   frozen one, the de-streaking mechanism lives in the base-flow response.
%
%   No OpenFOAM run needed: reuses the committed flat and bump base-flow CSVs.

if nargin<1||isempty(n_eta); n_eta=60;  end
if nargin<2||isempty(n_xi);  n_xi =1200; end
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));
ropt.lu_mode='none';

%% 1. Load flat + bump grids (pre-built by destreak_gengrids in a clean session)
gf = fullfile(here,sprintf('gflat_%dx%d.mat',n_eta,n_xi));
gh = fullfile(here,sprintf('ghump_%dx%d.mat',n_eta,n_xi));
assert(isfile(gf) && isfile(gh), ...
    'Grids missing — run  destreak_gengrids(%d,%d)  first (separate MATLAB session).', n_eta, n_xi);
Gf = getfield(load(gf),'StabGrid');
Gh = getfield(load(gh),'StabGrid');
% lift flat to curvilinear (2-D X,Y) for the frozen-metric channel
Gfc = destreak_perturb_grid(Gf, zeros(Gf.ny,Gf.nx));
beta0 = 2*pi*Gf.lref/7.5e-3;

%% 2. Streak on each base flow
[u0,v0] = destreak_run_case(Gfc, beta0, ropt);      % flat
[uh,vh] = destreak_run_case(Gh,  beta0, ropt);      % bumped (full response)
ph = sum(conj(uh(:)).*u0(:)); ph=ph/abs(ph); uh=uh*ph; vh=vh*ph;

%% 3. Window + drain references (on flat)
x=Gf.x(:).'; xn=(x-x(1))/(x(end)-x(1));
win = xn>=0.30 & xn<=0.72;  stn=find(win);  ny=Gf.ny;
[~,Vopt] = destreak_Popt(Gfc,beta0,stn,u0);
r = Vopt - v0(:,stn);                                % required roller correction

%% 4. Frozen-metric channel: bend flat mesh by the bump's wall displacement
DYh = Gh.Y - Gfc.Y;                                  % same (xi,eta) index grid
Gfroz = destreak_perturb_grid(Gfc, DYh);
[uf2,vf2] = destreak_run_case(Gfroz, beta0, ropt);
pf=sum(conj(uf2(:)).*u0(:)); pf=pf/abs(pf); uf2=uf2*pf; vf2=vf2*pf;

dV_full   = vh(:,stn)  - v0(:,stn);
dV_frozen = vf2(:,stn) - v0(:,stn);

%% 4b. Alignment sanity: grids must share wall-normal extent, bump must sit
%  outside the window (so window differencing is valid).
fprintf('[bfresp] H_flat=%.3f H_hump=%.3f  max|DYh| all=%.3e in-window=%.3e\n', ...
        Gfc.H, Gh.H, max(abs(DYh(:))), max(abs(DYh(:,stn)),[],'all'));

%% 5. Energy-weighted inner products over the window
W=zeros(ny,numel(stn));
for idx=1:numel(stn); W(:,idx)=abs(loc_trap(Gfc.Y(:,stn(idx)))); end
ip=@(a,b) real(sum(W(:).*conj(a(:)).*b(:)));  nrm=@(a) sqrt(ip(a,a));
cosr=@(a) ip(a,r)/(nrm(a)*nrm(r));

%% 6. Transfers (each on its own base flow) over the window
Tf = xfer(Gfc, u0,  v0,  stn);      % flat streak
Th = xfer(Gh,  uh,  vh,  stn);      % bumped streak (full response)
Tz = xfer(Gfroz,uf2, vf2, stn);     % frozen-metric bump
Tc = xfer(Gfc, u0, embedw(Vopt,stn,Gf), stn);   % flat sealed-Poisson ceiling

prog = @(T) (sum(Tf)-sum(T))/(sum(Tf)-sum(Tc));

%% 7. Report
fprintf('\n============== BASE-FLOW RESPONSE vs FROZEN METRIC ==============\n');
fprintf(' bump wall height max|DY| = %.3e (lref units), window %d stations\n', max(abs(DYh(:))), numel(stn));
fprintf(' |dV_frozen|/|v_flat|  = %.3f   cos(dV_frozen, r) = %+.3f\n', nrm(dV_frozen)/nrm(v0(:,stn)), cosr(dV_frozen));
fprintf(' |dV_full|  /|v_flat|  = %.3f   cos(dV_full,   r) = %+.3f\n', nrm(dV_full)/nrm(v0(:,stn)),   cosr(dV_full));
fprintf(' single-DOF drain reach  frozen cos^2 = %.3f   full cos^2 = %.3f\n', cosr(dV_frozen)^2, cosr(dV_full)^2);
fprintf(' sum T_flat  = %+.4e  [+ sustains]\n', sum(Tf));
fprintf(' sum T_ceil  = %+.4e  [- drains]\n',   sum(Tc));
fprintf(' sum T_frozen= %+.4e   progress=%.3f\n', sum(Tz), prog(Tz));
fprintf(' sum T_hump  = %+.4e   progress=%.3f\n', sum(Th), prog(Th));
fprintf('================================================================\n');

R=struct('Gf',Gfc,'Gh',Gh,'beta0',beta0,'stn',stn,'u0',u0,'v0',v0, ...
  'uh',uh,'vh',vh,'uf2',uf2,'vf2',vf2,'r',r,'dV_full',dV_full,'dV_frozen',dV_frozen, ...
  'DYh',DYh,'Tf',Tf,'Th',Th,'Tz',Tz,'Tc',Tc, ...
  'cos_frozen',cosr(dV_frozen),'cos_full',cosr(dV_full), ...
  'prog_frozen',prog(Tz),'prog_hump',prog(Th));
save(fullfile(here,sprintf('bfresp_%dx%d.mat',n_eta,n_xi)),'R','-v7.3');
end

function w=loc_trap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
function T=xfer(G,u,v,stn)
ns=numel(stn); T=zeros(1,ns);
for k=1:ns; i=stn(k);
  T(k)=-sum(loc_trap(G.Y(:,i)).*(G.dyU(:,i).*real(conj(u(:,i)).*v(:,i)))); end
end
function F=embedw(Vw,stn,G)
F=zeros(G.ny,G.nx); F(:,stn)=Vw;
end
