function M_CG = shiftMomentToCG(F_B, M_A_B, r_A_CG_B)
%SHIFTMOMENTTOCG Move a force/moment resultant from point A to CG.
%
% r_A_CG_B is the vector from CG to the point A where M_A_B is reported.
M_CG = M_A_B(:) + cross(r_A_CG_B(:), F_B(:));
end
