function plot_amplitudes(StabGrid, StabRes, ref_data, scholten_csv, plot_metric)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
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

%% Description
% Final amplitude development plot.
%
%   plot_amplitudes(StabGrid, StabRes)
%   plot_amplitudes(StabGrid, StabRes, ref_data)
%   plot_amplitudes(StabGrid, StabRes, ref_data, scholten_csv)
%   plot_amplitudes(StabGrid, StabRes, ref_data, scholten_csv, plot_metric)
%
%   ref_data (optional) struct with:
%     .file     path to .mat file with reference StabRes/StabGrid
%     .label    legend prefix (e.g. 'NPSE', 'LPSE', 'DNS')
%
%   scholten_csv (optional): set to true (or a path string) to overlay the
%   digitized Scholten Fig. 1a data. Pass [] or omit to skip.
%
%   plot_metric (optional, default 'umax'):
%     'umax'     peak |u'|        (StabRes.A directly)
%     'urms'     RMS u'           (StabRes.A / sqrt(2))
%     'energy'   sqrt(E); E is kinetic (incompressible) or Chu/Mack (compressible)
%
%   Mode labels use StabRes.omegavec/betavec to determine (m,n) indices.

%% Inputs
% Name         size         unit   explanation
% StabGrid     (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% StabRes      (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% ref_data     (1,1) struct [-]    optional reference data to overlay; needs a file field pointing at a .mat
% scholten_csv char         [-]    optional path to a Scholten CSV reference curve
% plot_metric  char         [-]    quantity to plot: 'umax' (default) or 'energy', the latter requiring StabRes.E

%% Outputs
% None.

if nargin < 5 || isempty(plot_metric); plot_metric = 'umax'; end
switch lower(plot_metric)
    case 'umax'
        A_get  = @(SR, j) SR.A(j, :);
        sqrt2  = 1;     % no /sqrt(2) for ref overlay
        ylab   = '$\max|u''|$';
    case 'urms'
        A_get  = @(SR, j) SR.A(j, :) / sqrt(2);
        sqrt2  = sqrt(2);
        ylab   = '$u''_{\mathrm{rms}}$';
    case 'energy'
        A_get  = @(SR, j) sqrt(max(SR.E(j, :), 0));
        sqrt2  = NaN;   % energy: ref overlay only if ref has .E (rare); skipped otherwise
        ylab   = '$\sqrt{E_{\mathrm{kin}}}$';
    otherwise
        error('plot_amplitudes:unknown_metric', 'Unknown plot_metric: %s', plot_metric);
end

if strcmpi(plot_metric, 'energy') && ~isfield(StabRes, 'E')
    error('plot_amplitudes:no_E', 'plot_metric=''energy'' requires StabRes.E (run postprocess).');
end

figure(4); clf; set(gcf, 'Position', [100 100 900 500]); hold on;

% R = sqrt(Re_x) = sqrt(x_norm) * R0, with R0 = sqrt(x_norm(1))
R0    = sqrt(StabGrid.xun(1));
R_vec = sqrt(StabGrid.xun) * R0;

marker = {'o','+','*','s','^','v','d','x','<','>','p','h'};
xmark  = linspace(R_vec(1), R_vec(end), 15);
h_list = [];
legs   = {};

% Derive fundamental frequency/wavenumber for integer (m,n) labels
omega_vec = StabRes.omegavec;
beta_vec  = StabRes.betavec;
omega_f = min(omega_vec(omega_vec > 0));  % fundamental frequency
beta_f  = min(beta_vec(beta_vec > 0));    % fundamental wavenumber

A_all = [];  % collect all plotted amplitudes for auto ylim

% HNS (harmonic Navier-Stokes) modes, plotted in black using the selected metric via A_get
mk_count = 0;
for j = 1:size(StabRes.A, 1)
    if ~any(StabRes.A(j,:)); continue; end
    Aj = A_get(StabRes, j);
    A_all = [A_all, Aj(Aj > 0)];

    if ~isempty(omega_f); m = round(omega_vec(j)/omega_f); else; m = 0; end
    if ~isempty(beta_f);  n = round(beta_vec(j)/beta_f);   else; n = 0; end
    lbl = mode_label('HNS', m, n);

    % Solver results: lines only (no markers). Distinct linestyle per mode.
    linestyles = {'-','--',':','-.'};
    if m == 0 && n == 0
        h = semilogy(R_vec, Aj, 'k--', 'LineWidth', 1.5);
    else
        mk_count = mk_count + 1;
        ls = linestyles{mod(mk_count-1, numel(linestyles))+1};
        h = semilogy(R_vec, Aj, ['k' ls], 'LineWidth', 1.5);
    end
    h_list(end+1) = h;
    legs{end+1}   = lbl;
end

% Reference data (red, optional)
if nargin >= 3 && ~isempty(ref_data) && isfile(ref_data.file)
    S = load(ref_data.file);
    lbl_prefix = 'REF';
    if isfield(ref_data, 'label'); lbl_prefix = ref_data.label; end

    [ref_res, ref_grid] = find_stab_structs(S);

    if ~isempty(ref_res) && ~isempty(ref_grid)
        % Reference (m,n) labels from its own omegavec/betavec if available
        if isfield(ref_res,'omegavec') && isfield(ref_res,'betavec')
            ref_omega_vec = ref_res.omegavec;
            ref_beta_vec = ref_res.betavec;
            ref_omega_f = min(ref_omega_vec(ref_omega_vec > 0));
            ref_beta_f = min(ref_beta_vec(ref_beta_vec > 0));
        else
            ref_omega_vec = []; ref_beta_vec = []; ref_omega_f = []; ref_beta_f = [];
        end

        mk_count_r = 0;
        for j = 1:size(ref_res.A, 1)
            if ~any(ref_res.A(j,:)); continue; end
            % Apply same metric to ref. For 'energy' fall back to peak if ref
            % lacks .E (heritage refs don't store energy).
            if strcmpi(plot_metric, 'energy') && isfield(ref_res, 'E')
                Aj = sqrt(max(ref_res.E(j, :), 0));
            elseif strcmpi(plot_metric, 'energy')
                Aj = ref_res.A(j, :);   % fall back
            else
                Aj = ref_res.A(j, :) / sqrt2;
            end
            A_all = [A_all, Aj(Aj > 0)];

            if ~isempty(ref_omega_f) && numel(ref_omega_vec) >= j
                m = round(ref_omega_vec(j)/ref_omega_f);
            else
                m = j - 1;
            end
            if ~isempty(ref_beta_f) && numel(ref_beta_vec) >= j
                n = round(ref_beta_vec(j)/ref_beta_f);
            else
                n = 0;
            end
            lbl = mode_label(lbl_prefix, m, n);

            R0_ref    = sqrt(ref_grid.xun(1));
            R_vec_ref = sqrt(ref_grid.xun) * R0_ref;
            if m == 0 && n == 0
                h = semilogy(R_vec_ref, Aj, 'r--', 'LineWidth', 1.5);
            else
                mk_count_r = mk_count_r + 1;
                semilogy(R_vec_ref, Aj, 'r-', 'LineWidth', 1.5);
                ym = interp1(R_vec_ref, Aj, xmark);
                mk = marker{mod(mk_count_r-1, numel(marker))+1};
                h = semilogy(xmark, ym, mk, 'Color', 'r', 'MarkerSize', 5);
            end
            h_list(end+1) = h;
            legs{end+1}   = lbl;
        end
    else
        warning('Could not find StabRes/StabGrid in %s', ref_data.file);
    end
end

% Scholten reference (from combined CSV)
if nargin >= 4 && ~isempty(scholten_csv)
    if islogical(scholten_csv) || isnumeric(scholten_csv)
        if scholten_csv
            % Default location: <project root>/input/literature/...
            % plot_amplitudes lives in <project root>/src/plot/, so go up 2.
            project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            scholten_path = fullfile(project_root, 'input', 'literature', ...
                                     'Scholten_0.002_combined.csv');
        else
            scholten_path = '';
        end
    else
        scholten_path = scholten_csv;
    end
    if ~isempty(scholten_path)
        if ~isfile(scholten_path)
            warning('Scholten CSV not found: %s', scholten_path);
        else
            T = readmatrix(scholten_path);   % cols: R_00 urms_00 R_10 urms_10 ...
            n_pairs = floor(size(T,2) / 2);
            tags = {'00','10','20','30','40','50'};
            % Scholten data is u_rms; for 'umax' scale non-MFD (mean-flow distortion)
            % modes by sqrt(2). MFD, urms, and energy metrics are left as-is.
            if strcmpi(plot_metric, 'umax')
                for kk = 1:n_pairs
                    if strcmp(tags{kk}, '00'); continue; end   % MFD: leave as-is
                    T(:, 2*kk) = T(:, 2*kk) * sqrt(2);
                end
            end
            mk_count_s = 0;
            for k = 1:n_pairs
                R_k = T(:, 2*k-1);
                u_k = T(:, 2*k);
                ok  = isfinite(R_k) & isfinite(u_k) & u_k > 0;
                if nnz(ok) < 2; continue; end
                R_k = R_k(ok);  u_k = u_k(ok);
                [R_k, idx] = sort(R_k);
                u_k = u_k(idx);
                [R_k, iu] = unique(R_k, 'stable');
                u_k = u_k(iu);

                tag = tags{k};
                m = str2double(tag(1));
                n = str2double(tag(2));
                lbl = mode_label('Scholten', m, n);

                % Scholten reference: markers only (no connecting line)
                mk_count_s = mk_count_s + 1;
                mk = marker{mod(mk_count_s-1, numel(marker))+1};
                h  = semilogy(R_k, u_k, mk, 'Color', 'b', ...
                              'MarkerSize', 6, 'LineStyle', 'none');
                h_list(end+1) = h;
                legs{end+1}   = lbl;
                A_all = [A_all, u_k(:)'];
            end
        end
    end
end

set(gca, 'YScale', 'log');
xlim(R_vec([1 end]))
if ~isempty(A_all)
    ymin = 10^floor(log10(min(A_all)));
    ymax = 10^ceil(log10(max(A_all)));
    ylim([max(ymin, 1e-12), min(ymax, 1)]);
end
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel(ylab, 'Interpreter', 'latex', 'FontSize', 12)
title('Amplitude development', 'FontSize', 13)
if ~isempty(h_list)
    legend(h_list, legs, 'Location', 'EastOutside', 'FontSize', 8)
end
grid on; box on;

end

%% Local helpers

function lbl = mode_label(prefix, m, n)
if n == 0
    lbl = sprintf('%s (%d,0)', prefix, m);
else
    lbl = sprintf('%s (%d,%d)', prefix, m, n);
end
end

function [res, grid] = find_stab_structs(S)
%FIND_STAB_STRUCTS  Find StabRes and StabGrid in a loaded .mat struct.
    res = []; grid = [];
    fnames = fieldnames(S);
    for k = 1:numel(fnames)
        v = S.(fnames{k});
        if isstruct(v) && isfield(v, 'A')          % results: amplitudes (alpha optional)
            res = v;
        end
        if isstruct(v) && (isfield(v, 'xun') || isfield(v, 'x'))   % grid: streamwise coord
            grid = v;
            if ~isfield(grid, 'xun'); grid.xun = grid.x; end       % minimal refs store .x only
        end
    end
end
