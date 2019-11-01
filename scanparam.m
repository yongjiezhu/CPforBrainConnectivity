function varargout = scanparam2(defoptions,options)

allfields = fieldnames(defoptions);
opts = defoptions;
for k = 1:numel(allfields)
    if isfield(options,allfields{k})
        if ~iscell(options.(allfields{k})) && numel(options.(allfields{k})) < numel(defoptions.(allfields{k}))
            opts.(allfields{k}) = repmat(options.(allfields{k}),...
                1,numel(defoptions.(allfields{k})));
        else
            opts.(allfields{k}) = options.(allfields{k});
        end
    end
end
if nargout > 1
    varargout = struct2cell(opts);
else
    varargout = {opts};
end