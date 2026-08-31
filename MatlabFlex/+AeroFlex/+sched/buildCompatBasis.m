
function W = buildCompatBasis(P, allowIdentityAero)
%BUILDCOMPATBASIS Return stored bases for coordinate alignment.

if nargin < 2
    allowIdentityAero = false;
end

if isfield(P,'compat') && isfield(P.compat,'W')
    W = P.compat.W;
else
    W = struct();
end

if ~isfield(W,'q1') || isempty(W.q1)
    if isfield(P,'beam') && isfield(P.beam,'phi1')
        W.q1 = P.beam.phi1;
    else
        error('buildCompatBasis:MissingQ1Basis', ...
            'Missing q1 basis. Rebuild the library with point.compat.W.q1.');
    end
end

if ~isfield(W,'qxi') || isempty(W.qxi)
    if isfield(P,'base') && isfield(P.base,'phi_xi_modes')
        W.qxi = P.base.phi_xi_modes;
    else
        error('buildCompatBasis:MissingQxiBasis', ...
            'Missing qxi basis. Rebuild the library with point.compat.W.qxi.');
    end
end

if ~isfield(W,'qGam') || isempty(W.qGam)
    if allowIdentityAero
        W.qGam = eye(P.dims.Na);
    else
        error('buildCompatBasis:MissingQGamBasis', ...
            'Missing qGam basis. Rebuild the library with point.compat.W.qGam.');
    end
end

W.q1   = full(W.q1);
W.qxi  = full(W.qxi);
W.qGam = full(W.qGam);
end
% function W = buildCompatBasis(P, allowIdentityAero)
% %BUILDCOMPATBASIS Return bases used for coordinate alignment.
% 
% if nargin < 2
%     allowIdentityAero = false;
% end
% 
% W = struct();
% 
% W.q1 = AeroFlex.sched.getFirstNested(P, { ...
%     {'compat','W','q1'}, ...
%     {'beam','phi1'}, ...
%     {'beam','red','phi1'}, ...
%     {'beam','red','Phi1'} });
% 
% W.qxi = AeroFlex.sched.getFirstNested(P, { ...
%     {'compat','W','qxi'}, ...
%     {'base','phi_xi_modes'}, ...
%     {'base','phiXi_sA'} });
% 
% W.qGam = AeroFlex.sched.getFirstNested(P, { ...
%     {'compat','W','qGam'}, ...
%     {'aero','Z'}, ...
%     {'aero','basis','Z'}, ...
%     {'aero','DataMatrix','Z'}, ...
%     {'aero','DataMatrix','KrylovBasis'} });
% 
% if isempty(W.q1)
%     error('buildCompatBasis:MissingQ1Basis', ...
%         'Missing q1 basis. Store it in point.compat.W.q1 or point.beam.phi1.');
% end
% 
% if isempty(W.qxi)
%     error('buildCompatBasis:MissingQxiBasis', ...
%         'Missing qxi basis. Store it in point.compat.W.qxi or point.base.phi_xi_modes.');
% end
% 
% if isempty(W.qGam)
%     if allowIdentityAero
%         W.qGam = eye(P.dims.Na);
%     else
%         error('buildCompatBasis:MissingAeroBasis', ...
%             'Missing aerodynamic basis. Store it in point.compat.W.qGam.');
%     end
% end
% 
% W.q1 = full(W.q1);
% W.qxi = full(W.qxi);
% W.qGam = full(W.qGam);
% end