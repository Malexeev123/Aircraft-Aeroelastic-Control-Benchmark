function [S1, S2] = build_spatial_derivative_operators(fem, keepDofs, chain, elemLen, Eall)
% Build first- and second-derivative operators consistent with your psi2 stencil:
%   x1' ≈ -(x1_{i+1}-x1_i)/ds + 0.5*E_i^T (x1_{i+1}+x1_i)
% Then S2 is a centered difference of S1 along sub-elements.

ndof = length(keepDofs);
nSub = length(chain)-1;

rowBase = @(i) (i-1)*6 + (1:6);
S1 = sparse(6*nSub, ndof);
for i = 1:nSub
    r  = rowBase(i);
    ds = elemLen(i);
    E6 = Eall(r,:); % 6x6 at sub i

    nA = chain(i);
    nB = chain(i+1);
    ia = localIdx(keepDofs, (nA-1)*6 + (1:6));
    ib = localIdx(keepDofs, (nB-1)*6 + (1:6));

    if ~isempty(ia)
        S1(r, ia) = S1(r, ia) + (+eye(6)/ds) + 0.5*E6.';  % + because of the leading minus
    end
    if ~isempty(ib)
        S1(r, ib) = S1(r, ib) + (-eye(6)/ds) + 0.5*E6.';  % - for forward difference
    end
end

% Second-derivative: difference of first-derivative across neighbouring subs
S2 = sparse(6*nSub, ndof);
for i = 1:nSub
    ds = elemLen(i);
    if i>1
        S2(rowBase(i), :) = S2(rowBase(i), :) + ( S1(rowBase(i), :) - S1(rowBase(i-1), :) ) / ds;
    else
        % one-sided at the first sub
        S2(rowBase(i), :) = S2(rowBase(i), :) + ( S1(rowBase(i), :) ) / ds;
    end
end
end

function I = localIdx(keepDofs, dofRange)
I = [];
for d = dofRange(:).'
    k = find(keepDofs==d, 1, 'first');
    if ~isempty(k), I(end+1) = k; end
end
end
