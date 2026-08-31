
function T = procrustesLocalToRef(Wloc, Wref, method)
%PROCRUSTESLOCALTOREF Local-to-reference coordinate map.

if nargin < 3 || isempty(method)
    method = 'orthogonal';
end

if isempty(Wloc) || isempty(Wref)
    error('procrustesLocalToRef:EmptyBasis', 'Basis matrix is empty.');
end

Wloc = full(double(Wloc));
Wref = full(double(Wref));

if size(Wloc,1) ~= size(Wref,1)
    error('procrustesLocalToRef:BadRows', ...
        'Basis row counts differ: %d vs %d.', size(Wloc,1), size(Wref,1));
end

if size(Wloc,2) ~= size(Wref,2)
    error('procrustesLocalToRef:BadCols', ...
        'Basis column counts differ: %d vs %d.', size(Wloc,2), size(Wref,2));
end

if isequal(Wloc,Wref)
    T = eye(size(Wloc,2));
    return
end

switch lower(string(method))
    case "orthogonal"
        [U,~,V] = svd(Wref.'*Wloc,'econ');
        T = U*V.';
    case "least_squares"
        T = Wref \ Wloc;
    otherwise
        error('procrustesLocalToRef:UnknownMethod', ...
            'Unknown alignment method "%s".', method);
end

if any(~isfinite(T(:)))
    error('procrustesLocalToRef:NonFinite', 'Transform contains NaN/Inf.');
end

if rcond(T) < 1e-10
    warning('procrustesLocalToRef:IllConditioned', ...
        'Coordinate map is ill-conditioned. rcond(T)=%.3e.', rcond(T));
end
end
