function [gen, caseId] = resolveGenCase(inputsJlPath)
% RESOLVEGENCASE  Read the validation Gen/Case selector from airPower inputs.jl.
%
%   [gen, caseId] = resolveGenCase(inputsJlPath)
%
% airPower's inputs.jl (Julia) defines, inside its VAL block:
%       Gen  = 0,   # validation generation
%       Case = 0,   # validation case
% These pick Validation/Gen{gen}/Experimental/Case{caseId}/. We parse them
% directly so this MATLAB tool stays downstream of the single Julia source of
% truth rather than duplicating the numbers.

txt = fileread(inputsJlPath);

gen    = local_scan_int(txt, 'Gen');
caseId = local_scan_int(txt, 'Case');
end

% -------------------------------------------------------------------------
function v = local_scan_int(txt, key)
% Match e.g. "Gen     = 0," / "Case = 12", ignoring trailing comment/comma.
tok = regexp(txt, ['(?m)^\s*' key '\s*=\s*(-?\d+)'], 'tokens', 'once');
if isempty(tok)
    error('resolveGenCase:notFound', ...
          'Could not find `%s = <int>` in %s', key, inputsJlPath);
end
v = str2double(tok{1});
end
