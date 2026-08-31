function useProjection = runtimeRateProjectionPolicy(cfg,fallback)
%RUNTIMERATEPROJECTIONPOLICY Resolve projection from the runtime body case.

if nargin < 2 || isempty(fallback)
    fallback = false;
end
validateattributes(fallback,{'logical','numeric'}, ...
    {'scalar','real','finite'},mfilename,'fallback');
useProjection = logical(fallback);

bodyCase = "";
if isstruct(cfg) && isfield(cfg,'sim') && isstruct(cfg.sim)
    if isfield(cfg.sim,'bodyCase') && ~isempty(cfg.sim.bodyCase)
        bodyCase = string(cfg.sim.bodyCase);
    elseif isfield(cfg.sim,'body_case') && ~isempty(cfg.sim.body_case)
        bodyCase = string(cfg.sim.body_case);
    end
end
bodyCase = lower(erase(erase(strtrim(bodyCase),"_"),"-"));

if bodyCase == "wingonly"
    useProjection = false;
elseif any(bodyCase == ["coupledfull","fullycoupled"])
    useProjection = true;
end
end
