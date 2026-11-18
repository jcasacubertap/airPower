function [I2, out] = destreak_i2base(G, u, v, w, stn)
%DESTREAK_I2BASE  Base-flow-oriented lift-up transfer I2 (Casacuberta et al. 2024).
%
%   [I2, out] = destreak_i2base(G, u, v, w, stn)
%
%   Correct (non-parallel) transfer over a 3-D / swept base flow. At each window
%   station it forms the local base-flow tangent t_hat = vB/|vB|, splits the full
%   perturbation q'=[u',v',w'] into tangential vt'=(q'.t_hat)t_hat and cross-stream
%   remainder vn'=q'-vt', and evaluates
%        I2 = -<vt', L vn'>,   L vn' = (vn'.grad) vB
%   with the FULL base-flow gradient tensor grad vB = [dx vB, dy vB] (dz vB = 0,
%   base flow spanwise-invariant). The representer is g = (grad vB)^T vt', with
%   components g_x = vt'.dx vB, g_y = vt'.dy vB (g_z = 0); then I2 = -<g, vn'>.
%   Both forms are computed and returned for a consistency check.
%
%   This REPLACES the parallel reduction -int U' Re(conj u' v') dy over the feature.
%
%   Returns I2 [1 x ns] (per window station) and out with g_x,g_y, vn' components,
%   and I2_check (the <vt',Lvn'> form).

ny = G.ny;  ns = numel(stn);
I2 = zeros(1,ns);  I2b = zeros(1,ns);
gx=zeros(ny,ns); gy=zeros(ny,ns); vnx=zeros(ny,ns); vny=zeros(ny,ns); vnz=zeros(ny,ns);
vtx=zeros(ny,ns); vty=zeros(ny,ns); vtz=zeros(ny,ns);
for k = 1:ns
    i = stn(k);
    U=G.U(:,i); V=G.V(:,i); W=G.W(:,i);
    S=sqrt(U.^2+V.^2+W.^2);  S(S==0)=eps;
    tx=U./S; ty=V./S; tz=W./S;                 % base-flow tangent t_hat
    up=u(:,i); vp=v(:,i); wp=w(:,i);           % full perturbation q'
    tau = up.*tx + vp.*ty + wp.*tz;            % q'.t_hat
    vtx_=tau.*tx; vty_=tau.*ty; vtz_=tau.*tz;  % vt' (tangential)
    vnx_=up-vtx_; vny_=vp-vty_; vnz_=wp-vtz_;  % vn' (cross-stream remainder)
    dxU=G.dxU(:,i); dyU=G.dyU(:,i);
    dxV=G.dxV(:,i); dyV=G.dyV(:,i);
    dxW=G.dxW(:,i); dyW=G.dyW(:,i);
    % representer g = (grad vB)^T vt'
    gx_ = vtx_.*dxU + vty_.*dxV + vtz_.*dxW;
    gy_ = vtx_.*dyU + vty_.*dyV + vtz_.*dyW;
    % L vn' = vn'_x dx vB + vn'_y dy vB   (3 components)
    Lx = vnx_.*dxU + vny_.*dyU;
    Ly = vnx_.*dxV + vny_.*dyV;
    Lz = vnx_.*dxW + vny_.*dyW;
    wq = loc_trap(G.Y(:,i));
    I2b(k) = -sum(wq.*real(conj(vtx_).*Lx + conj(vty_).*Ly + conj(vtz_).*Lz));
    I2(k)  = -sum(wq.*real(conj(gx_).*vnx_ + conj(gy_).*vny_));
    gx(:,k)=gx_; gy(:,k)=gy_; vnx(:,k)=vnx_; vny(:,k)=vny_; vnz(:,k)=vnz_;
    vtx(:,k)=vtx_; vty(:,k)=vty_; vtz(:,k)=vtz_;
end
out = struct('gx',gx,'gy',gy,'vnx',vnx,'vny',vny,'vnz',vnz, ...
             'vtx',vtx,'vty',vty,'vtz',vtz,'I2_check',I2b);
end

function w = loc_trap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
