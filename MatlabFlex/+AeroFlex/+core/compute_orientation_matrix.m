function [X,chi_bar] = compute_orientation_matrix(phi_A, xi_bar_array)
   quat_A = xi_bar_array(1,:);
   
    eul_A = quat2eul(quat_A);
    roll_A = eul_A(1); 
    pitch_A = eul_A(2); 
    yaw_A = eul_A(3); 
    eul_mat_A = [1, sin(roll_A)*tan(pitch_A), -cos(roll_A)*tan(pitch_A);
                0, cos(roll_A), sin(roll_A);
                0, -sin(roll_A)/cos(pitch_A), cos(roll_A)/cos(pitch_A)];
    mat_input_A = eul_mat_A*[zeros(3,3), eye(3)];
    X = mat_input_A * phi_A;
    chi_bar = [roll_A;pitch_A;yaw_A];
end