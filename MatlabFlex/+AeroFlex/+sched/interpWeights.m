function [w, ids, info] = interpWeights(ROMlib, mu, cfgLibrary)
%INTERPWEIGHTS Linear interpolation weights for U-alpha library points.
%
% Supports scattered 2-D libraries using Delaunay triangulation.  If the
% requested point coincides with a training point, the returned weight is one
% at that point.

if nargin < 3 || isempty(cfgLibrary)
    cfgLibrary = struct('noExtrapolate',true,'interpTol',1e-10);
end
if ~isfield(cfgLibrary,'noExtrapolate'), cfgLibrary.noExtrapolate = true; end
if ~isfield(cfgLibrary,'interpTol'), cfgLibrary.interpTol = 1e-10; end

mu = double(mu(:).');
P = double(ROMlib.mu);
assert(size(P,2) == numel(mu), 'mu dimension mismatch.');

% Exact point match.
d = vecnorm(P - mu,2,2);
[dmin, imin] = min(d);
if dmin < cfgLibrary.interpTol
    w = 1;
    ids = imin;
    info = struct('mode','exact','simplex',[],'barycentric',1,'distance',dmin);
    return
end

if size(P,2) == 1
    [Ps, order] = sort(P(:,1));
    mu1 = mu(1);
    if cfgLibrary.noExtrapolate && (mu1 < Ps(1) || mu1 > Ps(end))
        error('Requested mu=%.6g outside library range [%.6g, %.6g].', mu1, Ps(1), Ps(end));
    end
    if mu1 <= Ps(1)
        ids = order(1); w = 1;
    elseif mu1 >= Ps(end)
        ids = order(end); w = 1;
    else
        j = find(Ps <= mu1,1,'last');
        if j == numel(Ps), j = j-1; end
        ids = order([j j+1]);
        lam = (mu1 - Ps(j))/(Ps(j+1)-Ps(j));
        w = [1-lam; lam];
    end
    info = struct('mode','linear1d','simplex',[],'barycentric',w,'distance',dmin);
    return
end

if size(P,2) ~= 2
    error('interpWeights currently supports 1-D or 2-D libraries. Use [U_inf, alpha_deg].');
end

DT = delaunayTriangulation(P(:,1),P(:,2));
tid = pointLocation(DT, mu(1), mu(2));

if isnan(tid)
    if cfgLibrary.noExtrapolate
        error('Requested mu=[%.6g %.6g] is outside the ROM library convex hull.', mu(1), mu(2));
    end
    [~, ids] = mink(d, min(3,numel(d)));
    invd = 1./max(d(ids),eps);
    w = invd/sum(invd);
    info = struct('mode','nearest_extrap','simplex',[],'barycentric',w,'distance',dmin);
    return
end

ids = DT.ConnectivityList(tid,:).';
w = cartesianToBarycentric(DT, tid, mu);
w = w(:);

% Remove tiny negative numerical artifacts.
w(abs(w) < 10*eps) = 0;
w = w/sum(w);

info = struct('mode','barycentric2d','simplex',tid,'barycentric',w,'distance',dmin);
end
