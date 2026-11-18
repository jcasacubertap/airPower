function applyP = destreak_Popt_base(G, beta0, stn)
%DESTREAK_POPT_BASE  General (non-parallel) optimal de-streaking roller = ceiling P.
%
%   applyP = destreak_Popt_base(G, beta0, stn)
%
%   Builds the base-flow-oriented optimal cross-stream roller operator. Given a
%   tangential perturbation vt' = (vtx,vty,vtz) it forms the representer
%       g_x = vt'.dx vB ,   g_y = vt'.dy vB
%   and returns the energy-optimal solenoidal roller that maximises the drain
%   <g, vn'> at fixed energy (local cross-plane continuity dy vn'_y + i beta vn'_z = 0):
%       vn'_x = g_x                              (streamwise channel, unconstrained)
%       vn'_y = i beta psi ,  psi'' - beta^2 psi = i beta g_y ,  psi=0 at walls (sealed)
%   In the parallel limit (g_x->0, g_y->U' u') this collapses to solve_psihat.
%
%   applyP(VT) takes VT with fields vtx,vty,vtz (each [ny x ns], window-ordered
%   to match stn) and returns the optimal roller as a struct with .x and .y
%   ([ny x ns]).

ny=G.ny; ns=numel(stn); I=eye(ny);
decs=cell(1,ns);
dxU=cell(1,ns); dyU=cell(1,ns); dxV=cell(1,ns); dyV=cell(1,ns); dxW=cell(1,ns); dyW=cell(1,ns);
for k=1:ns
    i=stn(k);
    Dy=diag(G.etay(:,i))*G.D1;  L=Dy*Dy - beta0^2*I;
    L(1,:)=0; L(1,1)=1;  L(ny,:)=0; L(ny,ny)=1;      % sealed psi=0
    decs{k}=decomposition(L,'lu');
    dxU{k}=G.dxU(:,i); dyU{k}=G.dyU(:,i);
    dxV{k}=G.dxV(:,i); dyV{k}=G.dyV(:,i);
    dxW{k}=G.dxW(:,i); dyW{k}=G.dyW(:,i);
end
applyP=@(VT) impl(VT,decs,dxU,dyU,dxV,dyV,dxW,dyW,beta0,stn,ny);
end

function out=impl(VT,decs,dxU,dyU,dxV,dyV,dxW,dyW,beta0,stn,ny)
ns=numel(stn); vnox=zeros(ny,ns); vnoy=zeros(ny,ns);
for k=1:ns
    vtx=VT.vtx(:,k); vty=VT.vty(:,k); vtz=VT.vtz(:,k);
    gx=vtx.*dxU{k}+vty.*dxV{k}+vtz.*dxW{k};
    gy=vtx.*dyU{k}+vty.*dyV{k}+vtz.*dyW{k};
    rhs=1i*beta0*gy; rhs(1)=0; rhs(ny)=0;
    psi=decs{k}\rhs;
    vnox(:,k)=gx;
    vnoy(:,k)=1i*beta0*psi;
end
out=struct('x',vnox,'y',vnoy);
end
