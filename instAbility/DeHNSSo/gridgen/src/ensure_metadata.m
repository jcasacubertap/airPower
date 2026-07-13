function BF_new = ensure_metadata(BF_raw, BF_new, params)
%ENSURE_METADATA  Guarantee Re/lref/Uref/nu are present on the resampled BF.
%
%   BF_new = ensure_metadata(BF_raw, BF_new, params)
%
%   Fills missing reference scales on BF_new using this precedence:
%       1. value already set on BF_new
%       2. value carried over from BF_raw
%       3. value provided explicitly in params
%       4. derived from the other three (nu = Uref*lref/Re, etc.)
%
%   Errors if Re cannot be established — downstream build_stab_grid needs it.

fields = {'Re','lref','Uref','nu'};

% 1–3: inherit from BF_raw / params
for k = 1:numel(fields)
    f = fields{k};
    if isfield(BF_new, f) && ~isempty(BF_new.(f));  continue;  end
    if isfield(BF_raw, f) && ~isempty(BF_raw.(f))
        BF_new.(f) = BF_raw.(f);
    elseif isfield(params, f) && ~isempty(params.(f))
        BF_new.(f) = params.(f);
    end
end

% 4: derive the one missing scale from the other three, if possible
has = @(f) isfield(BF_new, f) && ~isempty(BF_new.(f));
if ~has('nu')   && has('Uref') && has('lref') && has('Re')
    BF_new.nu   = BF_new.Uref * BF_new.lref / BF_new.Re;
end
if ~has('Re')   && has('Uref') && has('lref') && has('nu')
    BF_new.Re   = BF_new.Uref * BF_new.lref / BF_new.nu;
end
if ~has('lref') && has('Uref') && has('Re')   && has('nu')
    BF_new.lref = BF_new.Re * BF_new.nu / BF_new.Uref;
end
if ~has('Uref') && has('lref') && has('Re')   && has('nu')
    BF_new.Uref = BF_new.Re * BF_new.nu / BF_new.lref;
end

if ~has('Re')
    error('ensure_metadata:missingRe', ...
        'Re is missing from BF and params. Set params.Re or include BF.Re.');
end
end
