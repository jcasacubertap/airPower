function destreak_orient(id)
% Diagnose flow direction, wall-shape asymmetry, and separation location.
if nargin<1; id='e05'; end
here=fileparts(mfilename('fullpath'));
G=getfield(load(fullfile(here,sprintf('gsweep_%s_60x1200.mat',id))),'StabGrid');
ny=G.ny; x=G.x(:).'; xn=(x-x(1))/(x(end)-x(1));
fprintf('[%s] freestream U (row1 mean) = %+.3f  (flow +x if >0)\n', id, mean(G.U(1,:)));
wallY=G.Y(ny,:); [ypk,ipk]=max(wallY);
fprintf('[%s] wall peak at xi/L=%.3f  h_peak=%.4f h_min=%.4f\n', id, xn(ipk), ypk, min(wallY));
dY=gradient(wallY,x);
up=max(1,ipk-80):ipk; dn=ipk:min(numel(x),ipk+80);
fprintf('[%s] steepest slope upstream=%+.3f  downstream=%+.3f  (|up|>|dn| => forward-facing)\n', ...
        id, max(abs(dY(up))), max(abs(dY(dn))));
nwb=(ny-round(ny/4)):ny; mU=min(G.U(nwb,:),[],1);
isep=find(mU<-0.02);
if ~isempty(isep)
  fprintf('[%s] separation (U<-0.02) at xi/L=[%.3f .. %.3f]; peak at %.3f  (sep downstream => bump lee)\n', ...
          id, xn(isep(1)), xn(isep(end)), xn(ipk));
else
  fprintf('[%s] no significant near-wall reversal\n', id);
end
end
