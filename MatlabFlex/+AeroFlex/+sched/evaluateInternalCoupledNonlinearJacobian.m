function [Nq,Nu] = evaluateInternalCoupledNonlinearJacobian( ...
        state,idx,parConst,transport)
%EVALUATEINTERNALCOUPLEDNONLINEARJACOBIAN Evaluate V8 internal derivatives.

arguments
    state double
    idx (1,1) struct
    parConst (1,1) struct
    transport (1,1) struct = struct()
end

if ~isfield(transport,'enabled') || ~logical(transport.enabled)
    [Nq,Nu]=AeroFlex.sim.nonlinearJacobian(state,idx,parConst);
    return
end
count=numel(state);
if ~isfield(transport,'physicalToInternal') || ...
        ~isfield(transport,'internalToPhysical') || ...
        ~isequal(size(transport.physicalToInternal),[count count]) || ...
        ~isequal(size(transport.internalToPhysical),[count count])
    error('AeroFlex:sched:InternalCoupledTransport', ...
        'The V8 internal coupled-coordinate transport is invalid.');
end
physical=transport.internalToPhysical*state(:);
[NqPhysical,NuPhysical]=AeroFlex.sim.nonlinearJacobian(physical,idx,parConst);
Nq=transport.physicalToInternal*NqPhysical*transport.internalToPhysical;
Nu=transport.physicalToInternal*NuPhysical;
end
