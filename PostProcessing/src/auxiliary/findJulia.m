function [jl, tried] = findJulia(isWin)
% findJulia  Locate a julia executable, without relying on the PATH MATLAB sees.
%
%   jl          = findJulia()        % this platform
%   [jl, tried] = findJulia(isWin)   % force the Windows / POSIX branch (for tests)
%
% Returns '' if nothing was found; `tried` describes the places looked at, for
% the caller's error message.
%
% Why this exists: a MATLAB launched from Finder, the Dock or the Windows Start
% menu inherits a minimal environment, so a perfectly good julia is often absent
% from the PATH that `system()` can see. The ladder is
%   1. $AIRPOWER_JULIA        explicit override
%   2. the shell's own lookup (`where` on Windows, `command -v` elsewhere)
%   3. the standard install locations for the platform, newest version first
%
% Windows differs in every step: the probe is `where`, the binary is julia.exe,
% HOME is usually unset (USERPROFILE holds it), and versioned installers drop
% Julia-<version> folders that need a numeric sort to pick the newest.

    if nargin < 1 || isempty(isWin); isWin = ispc; end

    jl = '';
    tried = {};

    % --- 1. explicit override ---
    ov = getenv('AIRPOWER_JULIA');
    if ~isempty(ov) && isfile(ov); jl = ov; return; end
    tried{end+1} = '$AIRPOWER_JULIA';

    % --- 2. ask the shell ---
    if isWin; probe = 'where julia'; else; probe = 'command -v julia'; end
    tried{end+1} = probe;
    [st, out] = system(probe);
    if st == 0 && ~isempty(strtrim(out))
        lines = strsplit(strtrim(out), {sprintf('\n'), sprintf('\r')});
        lines = lines(~cellfun(@isempty, lines));       % `where` can return several
        if ~isempty(lines) && isfile(strtrim(lines{1}))
            jl = strtrim(lines{1});  return;
        end
    end

    % --- 3. standard install locations ---
    [cands, where] = juliaCandidates(isWin);
    tried = [tried, where];
    for k = 1:numel(cands)
        if isfile(cands{k}); jl = cands{k}; return; end
    end
end

% --- candidate binaries, in preference order, and a human-readable list of the
%     directories they came from (used verbatim in the caller's error message) ---
function [cands, where] = juliaCandidates(isWin)

    home = getenv('HOME');
    if isempty(home); home = getenv('USERPROFILE'); end     % Windows

    if ~isWin
        exe   = 'julia';
        where = {'~/.juliaup/bin', '~/.local/bin', '/opt/homebrew/bin', ...
                 '/usr/local/bin', '/usr/bin'};
        cands = { fullfile(home, '.juliaup', 'bin', exe), ...
                  fullfile(home, '.local',   'bin', exe), ...
                  '/opt/homebrew/bin/julia', '/usr/local/bin/julia', '/usr/bin/julia' };
        return;
    end

    exe   = 'julia.exe';
    pf    = getenv('PROGRAMFILES');   if isempty(pf);  pf  = 'C:\Program Files';       end
    pf86  = getenv('PROGRAMFILES(X86)');
    if isempty(pf86); pf86 = 'C:\Program Files (x86)'; end
    lad   = getenv('LOCALAPPDATA');   if isempty(lad); lad = fullfile(home, 'AppData', 'Local'); end

    where = {'%USERPROFILE%\.juliaup\bin', '%LOCALAPPDATA%\Microsoft\WindowsApps', ...
             '%LOCALAPPDATA%\Programs\Julia-*', ...
             'C:\Program Files\Julia-*', 'C:\Program Files (x86)\Julia-*'};

    cands = { fullfile(home, '.juliaup', 'bin', exe), ...        % juliaup (recommended)
              fullfile(lad, 'Microsoft', 'WindowsApps', exe) };  % Store shim

    % Versioned installers: Julia-1.11.2, Julia-1.9.4, ... Sort NUMERICALLY —
    % alphabetically 'Julia-1.9' beats 'Julia-1.11', which would pick the oldest.
    for r = { fullfile(lad, 'Programs'), pf, pf86 }
        d = dir(fullfile(r{1}, 'Julia-*'));
        d = d([d.isdir]);
        if isempty(d); continue; end
        v = zeros(numel(d), 3);
        for k = 1:numel(d)
            t = sscanf(erase(d(k).name, 'Julia-'), '%d.%d.%d');
            v(k, 1:numel(t)) = t(:).';
        end
        [~, ord] = sortrows(v, [-1 -2 -3]);                      % newest first
        for k = ord(:).'
            cands{end+1} = fullfile(r{1}, d(k).name, 'bin', exe); %#ok<AGROW>
        end
    end
end
