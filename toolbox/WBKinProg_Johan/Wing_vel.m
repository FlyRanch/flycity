function [ U_wing_L, U_wing_R , U_wing_rot_L, U_wing_rot_R, U_trans_L, U_trans_R, U_body_rot_L, U_body_rot_R] = Wing_vel(wing_loc_L, wing_loc_R, omega_b, omega_L, omega_R, J_loc_L, J_loc_R, qL, qR, U_body)

    DCM_L = quat2matNEW(qL);

    DCM_R = quat2matNEW(qR);
    
    N = size(wing_loc_L,2);
      
    
    for i = 1:N

        U_trans_L(:,i) = DCM_L'*U_body;

        U_trans_R(:,i) = DCM_R'*U_body;

        U_body_rot_L(:,i) = DCM_L'*cross(omega_b,J_loc_L+DCM_L*wing_loc_L(:,i));

        U_body_rot_R(:,i) = DCM_R'*cross(omega_b,J_loc_R+DCM_R*wing_loc_R(:,i));

%         U_body_rot_L(:,i) = DCM_L*cross(omega_b,J_loc_L+DCM_L'*wing_loc_L(:,i));
% 
%         U_body_rot_R(:,i) = DCM_R*cross(omega_b,J_loc_R+DCM_R'*wing_loc_R(:,i));


        U_wing_rot_L(:,i) = cross(omega_L,wing_loc_L(:,i));

        U_wing_rot_R(:,i) = cross(omega_R,wing_loc_R(:,i));

        U_wing_L(:,i) = U_trans_L(:,i) + U_body_rot_L(:,i) + U_wing_rot_L(:,i);

        U_wing_R(:,i) = U_trans_R(:,i) + U_body_rot_R(:,i) + U_wing_rot_R(:,i);
        
    end

end

