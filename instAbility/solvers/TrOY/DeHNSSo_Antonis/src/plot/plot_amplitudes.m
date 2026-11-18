function plot_amplitudes(StabGrid, StabRes, ref_data, scholten_csv)
%PLOT_AMPLITUDES  Final amplitude development plot.
%
%   plot_amplitudes(StabGrid, StabRes)
%   plot_amplitudes(StabGrid, StabRes, ref_data)
%   plot_amplitudes(StabGrid, StabRes, ref_data, scholten_csv)
%
%   ref_data (optional) struct with:
%     .file   — path to .mat file with reference StabRes/StabGrid
%     .label  — legend prefix (e.g. 'NPSE', 'LPSE', 'DNS')
%
%   scholten_csv (optional) — set to true (or a path string) to overlay the
%   digitised Scholten Fig. 1a data from a combined CSV with columns
%       R_00, urms_00, R_10, urms_10, R_20, urms_20, R_30, urms_30
%   If `true`, defaults to input/literature/Scholten_0.002_combined.csv.
%   Pass [] or omit to skip.
%
%   Mode labels use StabRes.omegavec/betavec to determine (m,n) indices,
%   so 3D cases (N>0) are labelled correctly.

figure(4); clf; set(gcf, 'Position', [100 100 900 500]); hold on;

% R = sqrt(Re_x) = sqrt(x_norm * R0), with R0 = sqrt(x_norm(1))
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

% HNS modes (black) — plot StabRes.A directly (peak |u'|, German eq. 4.7b)
mk_count = 0;
for j = 1:size(StabRes.A, 1)
    if ~any(StabRes.A(j,:)); continue; end
    Aj = StabRes.A(j,:);
    A_all = [A_all, Aj(Aj > 0)]; %#ok<AGROW>

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
    h_list(end+1) = h; %#ok<AGROW>
    legs{end+1}   = lbl; %#ok<AGROW>
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
            ref_ov = ref_res.omegavec;
            ref_bv = ref_res.betavec;
            ref_of = min(ref_ov(ref_ov > 0));
            ref_bf = min(ref_bv(ref_bv > 0));
        else
            ref_ov = []; ref_bv = []; ref_of = []; ref_bf = [];
        end

        mk_count_r = 0;
        for j = 1:size(ref_res.A, 1)
            if ~any(ref_res.A(j,:)); continue; end
            Aj = ref_res.A(j,:);
            A_all = [A_all, Aj(Aj > 0)]; %#ok<AGROW>

            if ~isempty(ref_of) && numel(ref_ov) >= j
                m = round(ref_ov(j)/ref_of);
            else
                m = j - 1;
            end
            if ~isempty(ref_bf) && numel(ref_bv) >= j
                n = round(ref_bv(j)/ref_bf);
            else
                n = 0;
            end
            lbl = mode_label(lbl_prefix, m, n);

            if isfield(ref_grid, 'xun')
                xun_ref = ref_grid.xun;
            else
                xun_ref = ref_grid.x;
            end
            R0_ref    = sqrt(xun_ref(1));
            R_vec_ref = sqrt(xun_ref) * R0_ref;
            if m == 0 && n == 0
                h = semilogy(R_vec_ref, Aj, 'r--', 'LineWidth', 1.5);
            else
                mk_count_r = mk_count_r + 1;
                semilogy(R_vec_ref, Aj, 'r-', 'LineWidth', 1.5);
                ym = interp1(R_vec_ref, Aj, xmark);
                mk = marker{mod(mk_count_r-1, numel(marker))+1};
                h = semilogy(xmark, ym, mk, 'Color', 'r', 'MarkerSize', 5);
            end
            h_list(end+1) = h; %#ok<AGROW>
            legs{end+1}   = lbl; %#ok<AGROW>
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
            % Convert Scholten u_rms -> peak amplitude on the German convention:
            %   non-MFD modes (eq. 4.7b/c): peak = sqrt(2) * u_rms
            %   MFD (0,0) (eq. 4.7a):       peak = u_rms        (factor 1, no doubling)
            for kk = 1:n_pairs
                tag = tags{kk};
                if strcmp(tag, '00')
                    continue   % MFD: leave as-is
                end
                T(:, 2*kk) = T(:, 2*kk) * sqrt(2);
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
                h_list(end+1) = h; %#ok<AGROW>
                legs{end+1}   = lbl; %#ok<AGROW>
                A_all = [A_all, u_k(:)']; %#ok<AGROW>
            end
        end
    end
end

set(gca, 'YScale', 'log');
%xlim(R_vec([1 end]))
if ~isempty(A_all)
    ymin = 10^floor(log10(min(A_all)));
    ymax = 10^ceil(log10(max(A_all)));
    ylim([max(ymin, 1e-12), min(ymax, 1)]);
end
xlabel('$R = \sqrt{Re_x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\max|u''|$', 'Interpreter', 'latex', 'FontSize', 12)
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
%   Accepts both modern schema (StabRes.A/alpha, StabGrid.xun) and legacy
%   NPSE exports that only carry StabRes.A and StabGrid.x.
    res = []; grid = [];
    fnames = fieldnames(S);
    for k = 1:numel(fnames)
        v = S.(fnames{k});
        if isstruct(v) && isfield(v, 'A')
            res = v;
        end
        if isstruct(v) && (isfield(v, 'xun') || isfield(v, 'x'))
            grid = v;
        end
    end
end
