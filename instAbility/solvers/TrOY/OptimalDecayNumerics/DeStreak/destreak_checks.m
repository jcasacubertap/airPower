%% DESTREAK_CHECKS  Verify the geometry-sensitivity operator on the cheap grid.
%  (1) linearity of V wrt eps (is the FD sensitivity in the linear regime?)
%  (2) alignment: v_hom vs v_opt, and required correction r vs range(V).

clear; clc;
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));

S = load(fullfile(here,'flat_grid_40x300.mat'));  G = S.StabGrid;
beta0 = 2*pi*G.lref/7.5e-3;
Gc = destreak_perturb_grid(G, zeros(G.ny,G.nx));
ropt.lu_mode = 'none';
[u0,v0] = destreak_run_case(Gc, beta0, ropt);

B = destreak_relief_basis(Gc, struct('K',6));
stn = find(B.win);
phaseal = @(uk,ref) uk * (sum(conj(uk(:)).*ref(:))/abs(sum(conj(uk(:)).*ref(:))));

%% (1) Linearity of mode 1: V at eps vs 2*eps
DY = B.DY{1};
G1 = destreak_perturb_grid(Gc, DY);    [u1,v1] = destreak_run_case(G1,beta0,ropt);
G2 = destreak_perturb_grid(Gc, 2*DY);  [u2,v2] = destreak_run_case(G2,beta0,ropt);
ph1 = sum(conj(u1(:)).*u0(:)); v1 = v1*ph1/abs(ph1);
ph2 = sum(conj(u2(:)).*u0(:)); v2 = v2*ph2/abs(ph2);
V1 = (v1-v0)/B.eps;        % 1st-order FD
V2 = (v2-v0)/(2*B.eps);    % half-slope (should match V1 if linear)
lin_err = norm(V2(:,stn)-V1(:,stn),'fro')/norm(V1(:,stn),'fro');
fprintf('[chk] linearity mode1: ||V(2eps)-V(eps)||/||V(eps)|| = %.3e  (eps=%.3e)\n', lin_err, B.eps);

%% (2) Alignment v_hom vs v_opt (energy-weighted over window)
[~, Vopt] = destreak_Popt(Gc, beta0, stn, u0);
wcol = @(i) abs(local_trap(Gc.Y(:,i)));
W = zeros(Gc.ny, numel(stn));
for idx=1:numel(stn); W(:,idx)=wcol(stn(idx)); end
ip = @(a,b) real(sum(W(:).*conj(a(:)).*b(:)));
nrm = @(a) sqrt(ip(a,a));
vh = v0(:,stn);  vo = Vopt;
cos_ho = ip(vh,vo)/(nrm(vh)*nrm(vo));
fprintf('[chk] cos(v_hom, v_opt) energy-weighted = %+.3f  (|v_hom|=%.2e |v_opt|=%.2e)\n', ...
        cos_ho, nrm(vh), nrm(vo));

%% (3) required correction r vs range(V) — real-c reachable
r = vo - vh;
K = 6;  Vmat = zeros(Gc.ny*numel(stn), K);
for k=1:K
    Gk = destreak_perturb_grid(Gc, B.DY{k});
    [uk,vk] = destreak_run_case(Gk,beta0,ropt);
    phk = sum(conj(uk(:)).*u0(:)); vk = vk*phk/abs(phk);
    Vk = (vk-v0)/B.eps;  Vmat(:,k) = reshape(Vk(:,stn),[],1);
end
sw = sqrt(W(:));
Ar = [real(sw.*Vmat); imag(sw.*Vmat)];
rr = [real(sw.*r(:)); imag(sw.*r(:))];
c  = Ar\rr;  R2 = 1 - norm(Ar*c-rr)^2/norm(rr)^2;
fprintf('[chk] real-c reach of r by range(V): R^2 = %.3f\n', R2);
fprintf('[chk] DONE\n');

function w = local_trap(Y)
Ys=flipud(Y(:)); d=diff(Ys); n=numel(Ys); w=zeros(n,1);
w(1)=d(1)/2; w(end)=d(end)/2; w(2:end-1)=(d(1:end-1)+d(2:end))/2; w=flipud(w);
end
