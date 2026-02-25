function T_AB = T_phi_quat(xi)
% eq. (32) style from the reference:
%  T_phi(xi^{(a)}, xi^{(b)})= -2 (xi_v^(a))^T xi_v^(b)* I_3 + 2 xi_v^(a) (xi_v^(b))^T + ...
%                             + 2 xi_0^(a) ~xi_v^(b).
% We'll do a version for when A=B => so "xiB" = "xiA"? 
% Or we keep the general approach with different xiA, xiB. 
% This is the wrong eq I think
% Let xiA= [xi0_A; xiV_A], xiB= [xi0_B; xiV_B].

  xi0 = xi(1);
  xv  = xi(2:4);

  dotVV = xv.' * xv;  % scalar
  crossMat = AeroFlex.core.skew(xv);
  T_AB = (xi0^2 - dotVV)*eye(3) + 2*(xv*xv.') + 2*xi0*crossMat; 
end