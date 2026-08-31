function value = evaluateInternalCoupledNonlinearTerms(state,parConst,idx,transport)
%EVALUATEINTERNALCOUPLEDNONLINEARTERMS Evaluate nonlinear terms in V8 coordinates.

arguments
    state double
    parConst (1,1) struct
    idx (1,1) struct
    transport (1,1) struct = struct()
end

if ~isfield(transport,'enabled') || ~logical(transport.enabled)
    value=AeroFlex.sim.nonlinear_terms(state,parConst,idx);
    return
end
requireTransport(transport,numel(state));
physical=transport.internalToPhysical*state(:);
value=transport.physicalToInternal* ...
    AeroFlex.sim.nonlinear_terms(physical,parConst,idx);
end

function requireTransport(transport,count)
required={'physicalToInternal','internalToPhysical'};
if ~all(isfield(transport,required)) || ...
        ~isequal(size(transport.physicalToInternal),[count count]) || ...
        ~isequal(size(transport.internalToPhysical),[count count])
    error('AeroFlex:sched:InternalCoupledTransport', ...
        'The V8 internal coupled-coordinate transport is invalid.');
end
end
