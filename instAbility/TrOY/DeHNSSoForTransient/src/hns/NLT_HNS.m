function [ f ] = NLT_HNS(Grid,RunJ,StabRes,D1,HB)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:
% DOI:

%% Find dimensions
[nf,ny,nx] = size(StabRes.u);
iu = 1i;

% Reset f
f = zeros(nf, 4*ny*nx);

% Metric coefficients
etay = Grid.etay;
etax = Grid.etax;
xiy  = Grid.xiy;
xix  = Grid.xix;

% ξ-derivatives via sparse non-uniform FD matrix (built once in init.m).
% On a uniform grid this reduces to classical 4th-order coefficients; on
% xrefined / curvilinear grids it correctly accounts for variable Δξ.
D_xi = Grid.D_xi;

%% Calculate derivative fields (batched over modes)

% xi-derivatives: apply the sparse D_xi matrix along rows (streamwise)
dudxi = zeros(nf, ny, nx);
dvdxi = zeros(nf, ny, nx);
dwdxi = zeros(nf, ny, nx);
for j = RunJ
    dudxi(j,:,:) = squeeze(StabRes.u(j,:,:)) * D_xi.';
    dvdxi(j,:,:) = squeeze(StabRes.v(j,:,:)) * D_xi.';
    dwdxi(j,:,:) = squeeze(StabRes.w(j,:,:)) * D_xi.';
end

% eta-derivatives: apply D1 to all stations at once per mode
dudeta = zeros(nf, ny, nx);
dvdeta = zeros(nf, ny, nx);
dwdeta = zeros(nf, ny, nx);
for j = RunJ
    dudeta(j,:,:) = D1 * squeeze(StabRes.u(j,:,:));
    dvdeta(j,:,:) = D1 * squeeze(StabRes.v(j,:,:));
    dwdeta(j,:,:) = D1 * squeeze(StabRes.w(j,:,:));
end

%% Calculate forcing (local ny x nx accumulation per mode, avoids 3D indexing)

for j = RunJ(RunJ >= round(nf/2))

    % Find interacting mode pairs from Harmonic Balancing
    [gharray, jkarray] = find(HB(:,:,j));

    % Keep only active mode pairs
    active = ismember(gharray, RunJ) & ismember(jkarray, RunJ);
    gharray = gharray(active);
    jkarray = jkarray(active);

    % Local accumulation arrays — avoids squeeze into 3D array each interaction
    xm_j = zeros(ny, nx);
    ym_j = zeros(ny, nx);
    zm_j = zeros(ny, nx);

    for j2 = 1:length(gharray)
        gh = gharray(j2);
        jk = jkarray(j2);
        beta = StabRes.betavec(jk);

        % Extract fields for this mode pair (squeeze once)
        u_gh = squeeze(StabRes.u(gh,:,:));
        v_gh = squeeze(StabRes.v(gh,:,:));
        w_gh = squeeze(StabRes.w(gh,:,:));

        dudxi_jk  = squeeze(dudxi(jk,:,:));
        dvdxi_jk  = squeeze(dvdxi(jk,:,:));
        dwdxi_jk  = squeeze(dwdxi(jk,:,:));
        dudeta_jk = squeeze(dudeta(jk,:,:));
        dvdeta_jk = squeeze(dvdeta(jk,:,:));
        dwdeta_jk = squeeze(dwdeta(jk,:,:));
        u_jk = squeeze(StabRes.u(jk,:,:));
        v_jk = squeeze(StabRes.v(jk,:,:));
        w_jk = squeeze(StabRes.w(jk,:,:));

        % Accumulate momentum forcing into local arrays
        xm_j = xm_j ...
            - u_gh.*dudxi_jk.*xix - u_gh.*dudeta_jk.*etax ...
            - v_gh.*dudxi_jk.*xiy - v_gh.*dudeta_jk.*etay ...
            - iu*beta .* w_gh.*u_jk;

        ym_j = ym_j ...
            - u_gh.*dvdxi_jk.*xix - u_gh.*dvdeta_jk.*etax ...
            - v_gh.*dvdxi_jk.*xiy - v_gh.*dvdeta_jk.*etay ...
            - iu*beta .* w_gh.*v_jk;

        zm_j = zm_j ...
            - u_gh.*dwdxi_jk.*xix - u_gh.*dwdeta_jk.*etax ...
            - v_gh.*dwdxi_jk.*xiy - v_gh.*dwdeta_jk.*etay ...
            - iu*beta .* w_gh.*w_jk;
    end

    % Reshape momentum into f (vectorised, skip i=1)
    block = [xm_j(:,2:nx); ym_j(:,2:nx); zm_j(:,2:nx); zeros(ny, nx-1)];
    f(j, 4*ny+1:end) = block(:).';

end

end
