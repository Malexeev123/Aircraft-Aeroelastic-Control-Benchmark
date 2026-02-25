function UT = halfTangentialOperator( xibar )
% eq. (3.49) or eq. (3.50) from the reference:
%  U_t(xi_bar) = partial_x [ U(x)*xi_bar ] but the text states a closed form:
%  U_t( xi_bar ) = 1/2 [ -xi_bar_v^T;  xi_bar_0 I3 + tilde(xi_bar_v ) ]
%
% So we define:
  xi0 = xibar(1);
  xiv = xibar(2:4);
  UT = zeros(4,3);
  UT(1,1:3) = - (xiv.');
  skewv = [0 -xiv(3) xiv(2); xiv(3) 0 -xiv(1); -xiv(2) xiv(1) 0];
  UT(2:4,1:3) = xi0*eye(3) + skewv;
  UT = 0.5*UT;
end