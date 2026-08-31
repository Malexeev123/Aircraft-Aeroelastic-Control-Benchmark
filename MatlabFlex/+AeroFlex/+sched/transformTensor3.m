% function Tstar = transformTensor3(T, Tout, TaInv, TbInv)
% %TRANSFORMTENSOR3 Transform a bilinear tensor between reduced coordinates.
% %
% % Local bilinear map:
% %   y_i = sum_j sum_k T(i,j,k) a_j b_k
% %
% % Coordinate transforms:
% %   y_star = Tout*y,  a = TaInv*a_star,  b = TbInv*b_star
% %
% % Transformed tensor:
% %   Tstar(r,p,q) = sum_i,j,k Tout(r,i)*T(i,j,k)*TaInv(j,p)*TbInv(k,q)
% %
% % This routine is needed if structural coordinates q1, q2 or qxi are
% % transformed before interpolation.  If only qGam/aerodynamic states are
% % transformed, Gamma tensors can be interpolated directly.
% 
% [nout, na, nb] = size(T);
% assert(size(Tout,2)==nout, 'Tout size mismatch.');
% assert(size(TaInv,1)==na, 'TaInv size mismatch.');
% assert(size(TbInv,1)==nb, 'TbInv size mismatch.');
% 
% Tstar = zeros(size(Tout,1), size(TaInv,2), size(TbInv,2));
% for k = 1:nb
%     % Apply output and first-input transforms page by page.
%     page = Tout * T(:,:,k) * TaInv;
%     for q = 1:size(TbInv,2)
%         Tstar(:,:,q) = Tstar(:,:,q) + page * TbInv(k,q);
%     end
% end
% end
function Gout = transformTensor3(G, Tout, S1, S2)
%TRANSFORMTENSOR3 Apply Gref = Tout*Glocal*(S1,S2).

nout = size(Tout,1);
n1   = size(S1,2);
n2   = size(S2,2);

if size(G,1) ~= size(Tout,2)
    error('transformTensor3:BadOutputSize', ...
        'Tensor output size %d does not match transform size %d.', ...
        size(G,1), size(Tout,2));
end

if size(G,2) ~= size(S1,1)
    error('transformTensor3:BadInput1Size', ...
        'Tensor input-1 size %d does not match transform size %d.', ...
        size(G,2), size(S1,1));
end

if size(G,3) ~= size(S2,1)
    error('transformTensor3:BadInput2Size', ...
        'Tensor input-2 size %d does not match transform size %d.', ...
        size(G,3), size(S2,1));
end

Gout = zeros(nout,n1,n2,'like',G);

for k2 = 1:n2
    A = zeros(nout,n1,'like',G);
    for j2 = 1:size(S2,1)
        A = A + (Tout * G(:,:,j2) * S1) * S2(j2,k2);
    end
    Gout(:,:,k2) = A;
end
end