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
if any(~isfinite(mu))
    error('interpWeights:NonFiniteMu', ...
        'Requested schedule point contains NaN/Inf: [%s]', num2str(mu));
end
P = double(ROMlib.mu);
assert(size(P,2) == numel(mu), 'mu dimension mismatch.');

% Exact point match.
d = vecnorm(P - mu,2,2);
[dmin, imin] = min(d);
if isfield(cfgLibrary,'caseView') && ~isempty(cfgLibrary.caseView)
    [w,ids,info] = caseViewWeights(ROMlib,P,mu,d,dmin,imin, ...
        cfgLibrary.caseView,cfgLibrary.interpTol);
    return
end
if dmin < cfgLibrary.interpTol
    w = 1;
    ids = imin;
    info = struct('mode','exact','simplex',[],'barycentric',1,'distance',dmin);
    return
end

if isfield(cfgLibrary,'interpMode') && strcmpi(string(cfgLibrary.interpMode),"nearest")
    [~, imin] = min(d);
    w = 1;
    ids = imin;
    info = struct('mode','nearest','simplex',[],'barycentric',1,'distance',dmin);
    return
end

if size(P,2) == 1
    [Ps, order] = sort(P(:,1));
    mu1 = mu(1);
    if cfgLibrary.noExtrapolate && (mu1 < Ps(1) || mu1 > Ps(end))
        error('interpWeights:Extrapolation', ...
            'Requested mu=%.6g outside library range [%.6g, %.6g].', ...
            mu1, Ps(1), Ps(end));
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
        error('interpWeights:Extrapolation', ...
            'Requested mu=[%.6g %.6g] is outside the ROM library convex hull.', ...
            mu(1), mu(2));
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

% Exact-zero vertices do not belong to the active interpolation stencil.
% Retaining them changes certificate membership and can change the selected
% runtime time step even though they contribute nothing numerically. Preserve
% the original Delaunay simplex separately for provenance.
simplexIds = ids;
simplexWeights = w;
active = (w ~= 0);
ids = ids(active);
w = w(active);
w = w/sum(w);

info = struct('mode','barycentric2d','simplex',tid, ...
    'simplexIds',simplexIds,'simplexBarycentric',simplexWeights, ...
    'activeIds',ids,'barycentric',w,'distance',dmin);
end

function [w,ids,info] = caseViewWeights(ROMlib,P,mu,d,dmin,imin,view,tol)
required = {'testId','query','heldOutSourceIds','permittedSourceIds', ...
    'expectedInterpolationMode','caseManifestSha256'};
for k = 1:numel(required)
    if ~isfield(view,required{k})
        error('interpWeights:CaseViewMissingField', ...
            'Case view is missing required field %s.',required{k});
    end
end
if norm(double(view.query(:).')-mu,inf)>tol
    error('interpWeights:CaseViewQueryMismatch', ...
        'The runtime query does not match case view %s.',string(view.testId));
end
if ~isfield(ROMlib,'points') || ~isfield(ROMlib.points,'name')
    error('interpWeights:CaseViewSourceIdentity', ...
        'Case-view interpolation requires named library points.');
end
names = string({ROMlib.points.name});
permitted = string(view.permittedSourceIds(:)).';
heldOut = string(view.heldOutSourceIds(:)).';
if ~isequal(names,permitted)
    error('interpWeights:CaseViewSourceSetMismatch', ...
        'Library sources do not match the ordered permitted set for %s.', ...
        string(view.testId));
end
if any(ismember(names,heldOut))
    error('interpWeights:HeldoutSourcePresent', ...
        'A held-out source is present before weight evaluation for %s.', ...
        string(view.testId));
end

mode = string(view.expectedInterpolationMode);
switch mode
    case "exact"
        if dmin>=tol
            error('interpWeights:CaseViewExactMiss', ...
                'Exact-node case view %s did not resolve an exact source.', ...
                string(view.testId));
        end
        ids = imin;
        w = 1;
        modeOut = 'exact';
    case "case_barycentric2d"
        if size(P,1)~=3 || size(P,2)~=2
            error('interpWeights:CaseViewGeometry', ...
                'Barycentric case view %s requires three 2-D sources.', ...
                string(view.testId));
        end
        w = [P.';ones(1,3)]\[mu.';1];
        if any(w < -100*eps) || abs(sum(w)-1)>100*eps
            error('interpWeights:CaseViewGeometry', ...
                'Query lies outside the governed simplex for %s.', ...
                string(view.testId));
        end
        w(abs(w)<max(10*eps,tol)) = 0;
        ids = find(w~=0);
        w = w(ids);
        w = w/sum(w);
        modeOut = 'case_barycentric2d';
    case "case_guarded_idw4"
        if size(P,1)~=4 || size(P,2)~=2 || any(d<=tol)
            error('interpWeights:CaseViewGeometry', ...
                'Guarded-IDW case view %s requires four nonexact 2-D sources.', ...
                string(view.testId));
        end
        ids = (1:4).';
        w = 1./d.^2;
        w = w/sum(w);
        modeOut = 'case_guarded_idw4';
    case "case_linear2d"
        if size(P,1)~=2 || size(P,2)~=2
            error('interpWeights:CaseViewGeometry', ...
                'Linear case view %s requires two 2-D sources.', ...
                string(view.testId));
        end
        direction = P(2,:)-P(1,:);
        fraction = dot(mu-P(1,:),direction)/dot(direction,direction);
        projection = P(1,:)+fraction*direction;
        if norm(projection-mu,inf)>tol || fraction < -tol || fraction > 1+tol
            error('interpWeights:CaseViewGeometry', ...
                'Query is not on the governed segment for %s.', ...
                string(view.testId));
        end
        ids = [1;2];
        w = [1-fraction;fraction];
        w(abs(w)<max(10*eps,tol)) = 0;
        active = w~=0;
        ids = ids(active);
        w = w(active)/sum(w(active));
        modeOut = 'case_linear2d';
    otherwise
        error('interpWeights:UnsupportedCaseViewMode', ...
            'Unsupported case-view mode %s.',mode);
end
info = struct('mode',modeOut,'caseViewId',char(string(view.testId)), ...
    'caseViewHash',char(string(view.caseManifestSha256)), ...
    'activeIds',ids,'barycentric',w,'distance',dmin, ...
    'heldOutSourceIds',{cellstr(heldOut)});
end
