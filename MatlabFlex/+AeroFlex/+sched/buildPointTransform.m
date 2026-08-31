function tr = buildPointTransform(P, Pref, varargin)
%BUILDPOINTTRANSFORM Local-to-reference reduced coordinate map.

p = inputParser;
p.addParameter('allow_identity_aero', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

W    = AeroFlex.sched.buildCompatBasis(P, opt.allow_identity_aero);
Wref = AeroFlex.sched.buildCompatBasis(Pref, opt.allow_identity_aero);

tr = struct();

tr.q1.T      = AeroFlex.sched.constraintAlignedOrthogonalMap( ...
    W.q1, Wref.q1, P.beam.Pz, Pref.beam.Pz);
tr.q1.Tinv   = tr.q1.T.';
tr.q1.method = 'constraint_aligned_orthogonal_identity_metric';

% q1/q2 are paired intrinsic structural coordinates.
tr.q2.T    = tr.q1.T;
tr.q2.Tinv = tr.q1.Tinv;
tr.q2.method = 'paired_with_q1';

tr.qxi.T      = AeroFlex.sched.procrustesLocalToRef( ...
    W.qxi, Wref.qxi, 'orthogonal');
tr.qxi.Tinv   = tr.qxi.T.';
tr.qxi.method = 'orthogonal_procrustes';

% tr.qGam.T    = AeroFlex.sched.procrustesLocalToRef(W.qGam, Wref.qGam);
% tr.qGam.Tinv = pinv(tr.qGam.T);
% -------------------------------------------------------------------------
% Aerodynamic memory coordinates.
% -------------------------------------------------------------------------
if size(W.qGam,1) ~= size(Wref.qGam,1)
    error(['buildPointTransform:qGamAmbientMismatch\n', ...
           'qGam Krylov bases are not embedded in the same full aerodynamic state space.\n', ...
           'local mu=[%.3f %.3f], local size=%d x %d\n', ...
           'ref   mu=[%.3f %.3f], ref   size=%d x %d\n\n', ...
           'Rebuild the SHARPy ROM library with a consistent aerodynamic full-state layout.\n', ...
           'Recommended first fix: use mimo_rational_arnoldi with single_side=''observability''.'], ...
           P.mu(1), P.mu(2), size(W.qGam,1), size(W.qGam,2), ...
           Pref.mu(1), Pref.mu(2), size(Wref.qGam,1), size(Wref.qGam,2));
end

tr.qGam.T    = AeroFlex.sched.procrustesLocalToRef( ...
    W.qGam, Wref.qGam, 'least_squares');
tr.qGam.Tinv = tr.qGam.T \ eye(size(tr.qGam.T));
tr.qGam.method = 'basis_least_squares';


tr.chi.T    = eye(3);
tr.chi.Tinv = eye(3);

idx = P.idx;
nx = size(P.L,1);

T = eye(nx);
Ti = eye(nx);

T(idx.q1,idx.q1)       = tr.q1.T;
Ti(idx.q1,idx.q1)      = tr.q1.Tinv;

T(idx.q2,idx.q2)       = tr.q2.T;
Ti(idx.q2,idx.q2)      = tr.q2.Tinv;

T(idx.qxi,idx.qxi)     = tr.qxi.T;
Ti(idx.qxi,idx.qxi)    = tr.qxi.Tinv;

T(idx.qGam,idx.qGam)   = tr.qGam.T;
Ti(idx.qGam,idx.qGam)  = tr.qGam.Tinv;

if isfield(idx,'chi') && ~isempty(idx.chi)
    T(idx.chi,idx.chi)  = tr.chi.T;
    Ti(idx.chi,idx.chi) = tr.chi.Tinv;
end

tr.T = T;
tr.Tinv = Ti;
tr.err = struct();

tr.err.q1   = norm(Wref.q1*tr.q1.T - W.q1, 'fro') / max(1,norm(W.q1,'fro'));
tr.err.qxi  = norm(Wref.qxi*tr.qxi.T - W.qxi, 'fro') / max(1,norm(W.qxi,'fro'));
tr.err.qGam = norm(Wref.qGam*tr.qGam.T - W.qGam, 'fro') / max(1,norm(W.qGam,'fro'));
end
