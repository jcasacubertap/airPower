function v = get_opt(opts, name, default)
%GET_OPT  Fetch opts.(name) if set and non-empty, otherwise return default.
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
