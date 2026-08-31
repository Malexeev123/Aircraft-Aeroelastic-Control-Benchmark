function T = constraintAlignedOrthogonalMap(Wloc,Wref,Ploc,Pref)
%CONSTRAINTALIGNEDORTHOGONALMAP Align q1 while preserving root constraints.

n = size(Ploc,1);
if ~isequal(size(Ploc),[n n]) || ~isequal(size(Pref),[n n])
    error('constraintAlignedOrthogonalMap:ProjectorSize', ...
        'Local and reference projectors must be square and equal-sized.');
end

if isequal(Wloc,Wref) && isequal(Ploc,Pref)
    T = eye(n);
    return
end

[Uloc,Sloc] = svd((Ploc+Ploc.')/2);
[Uref,Sref] = svd((Pref+Pref.')/2);
rLoc = sum(diag(Sloc)>0.5);
rRef = sum(diag(Sref)>0.5);
if rLoc ~= rRef
    error('constraintAlignedOrthogonalMap:ProjectorRank', ...
        'Local/reference projector ranks differ: %d versus %d.',rLoc,rRef);
end

Qloc = {Uloc(:,1:rLoc),Uloc(:,rLoc+1:end)};
Qref = {Uref(:,1:rRef),Uref(:,rRef+1:end)};
O = cell(1,2);
for k = 1:2
    A = Wref*Qref{k};
    B = Wloc*Qloc{k};
    [Us,~,Vs] = svd(A.'*B,'econ');
    O{k} = Us*Vs.';
end

T = [Qref{1} Qref{2}]*blkdiag(O{1},O{2})*[Qloc{1} Qloc{2}].';
if norm(T.'*T-eye(n),'fro') > 1e-12
    error('constraintAlignedOrthogonalMap:Orthogonality', ...
        'Constraint-aligned q1 map is not orthogonal.');
end
end
