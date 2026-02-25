function e = sub_to_elem(fem, chain, iSub)
% map sub-element index to finite element index by connectivity walk
% Assumes chain follows the mesh order.
% This simplistic version treats each consecutively connected pair as belonging
% to the element whose (n1,n3,n2) includes both endpoints; adjust if needed.
nA = chain(iSub);
nB = chain(iSub+1);
for e = 1:fem.num_elem
    c = fem.connectivities(e,:);
    if any(c==nA) && any(c==nB)
        return;
    end
end
% fallback
e = max(1, min(iSub, fem.num_elem));
end
