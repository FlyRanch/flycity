function [ Aero_f ] = Aerodynamic_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st, down, up )


    % Function which computes aerodynamic force on left and right wing as a
    % function of wing and body kinematics for time t_1 until time t_2.
    
    Aero_f = {};
    
    % Import body parameters and wing kinematics:
    
    wing_l = Body.wing_l;       % mm    
    JL = Body.joint_L;          % mm    
    JR = Body.joint_R;          % mm    
    cg_body = Body.cg_body;     % mm    
    y_sect = Body.y_sect;       % mm    
    chords = Body.chords;       % mm 
    
    w_L = wing_kin.w_L;         % rad/s
    w_R = wing_kin.w_R;         % rad/s
    w_dot_L = wing_kin.w_dot_L; % rad/s^2
    w_dot_R = wing_kin.w_dot_R; % rad/s^2
    R_L = wing_kin.Rot_L;
    R_R = wing_kin.Rot_R;
    
    
    % Calculate density:
    
    pres = 101.325e3;           % Pressure in Pascal
    Temp = 297.039;             % Average temperature in Kelvin
    R = 287.058;                % Gas constant for air
    Rho = 1e-9*pres/(R*Temp);   % Density in kg/mm^3
    
    
    nr_sect = length(y_sect);
    
    
    % Sections on left and right wing:
    
    y_sect_L = -y_sect;
    y_sect_R = y_sect;
    
    
    % Position wing hinges w.r.t. cg_body:
    
    arm_JL = JL-cg_body;
    arm_JR = JR-cg_body;
    
    
    % Body velocity and angular velocity

    beta = (55/180)*pi;
    
    R_beta = [cos(-beta) 0 -sin(-beta); ...
              0 1 0; ...
              sin(-beta) 0 cos(-beta)]; 
    
    u_b = zeros(3,n);
    w_b = zeros(3,n);
          
    for i = 1:n
    
        u_b(:,i) = R_beta'*v_st(:,i);    
        w_b(:,i) = R_beta'*w_st(:,i);
    
    end
    
    
    % Calculate the velocity at the different sections of the wings:
    
    u_L = zeros(3,nr_sect,n);
    u_R = zeros(3,nr_sect,n);
    

   
    
    for j = 1:n
    
        for i = 1:nr_sect

            u_L(:,i,j) = (R_L(:,:,j)*(u_b(:,j) + cross(w_b(:,j),arm_JL + R_L(:,:,j)'*y_sect_L(:,i)))+cross(w_L(:,j),y_sect_L(:,i)));
            u_R(:,i,j) = (R_R(:,:,j)*(u_b(:,j) + cross(w_b(:,j),arm_JR + R_R(:,:,j)'*y_sect_R(:,i)))+cross(w_R(:,j),y_sect_R(:,i)));

        end
    
    end
    
%     UL_plot = zeros(3,n);
%     UR_plot = zeros(3,n);
%     
%     for j = 1:n
%         
%         UL_plot(:,j) = u_L(:,nr_sect,j);
%         UR_plot(:,j) = u_R(:,nr_sect,j);
%         
%     end
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(1:n,UL_plot(1,:),1:n,UR_plot(1,:))
%     subplot(3,1,2); plot(1:n,UL_plot(2,:),1:n,UR_plot(2,:))
%     subplot(3,1,3); plot(1:n,UL_plot(3,:),1:n,UR_plot(3,:))
%     hold off
%     pause
    
    % Calculate angle of attack at the sections:
    
    alfa_L = zeros(nr_sect,n);
    alfa_R = zeros(nr_sect,n);
    
    for j = 1:n
    
        for i = 1:nr_sect

           alfa_L(i,j) = real(atan2(u_L(3,i,j),u_L(1,i,j)));
           alfa_R(i,j) = real(atan2(u_R(3,i,j),u_R(1,i,j)));

        end
    
    end
    
    
    % Calculate the local lift and drag forces:
    
    deltaR = y_sect_R(2,2)-y_sect_R(2,1);
    
    Lift_L = zeros(nr_sect,n);
    Drag_L = zeros(nr_sect,n);
    Lift_R = zeros(nr_sect,n);
    Drag_R = zeros(nr_sect,n);
    F_loc_L = zeros(3,nr_sect,n);
    M_loc_L = zeros(3,nr_sect,n);
    F_loc_R = zeros(3,nr_sect,n);
    M_loc_R  = zeros(3,nr_sect,n);
    
%     for j = 1:n
%     
%         for i = 1:nr_sect
% 
%                 alfa = radtodeg(alfa_L(i,j));
%                 Cl = 0.225+1.58*sind(2.13*alfa-7.2);
%                 Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
%                 Lift_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cl*deltaR;
%                 Drag_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cd*deltaR;
%                 F_loc_L(:,i,j) = [(Lift_L(i,j)*sind(alfa)-Drag_L(i,j)*cosd(alfa)); 0; -(Lift_L(i,j)*cosd(alfa)+Drag_L(i,j)*sind(alfa))];
%                 M_loc_L(:,i,j) = cross(y_sect_L(:,i),F_loc_L(:,i,j));
% 
%             clear Cl Cd alfa
% 
%         end
% 
%         for i = 1:nr_sect
% 
%                 alfa = radtodeg(alfa_R(i,j));
%                 Cl = 0.225+1.58*sind(2.13*alfa-7.2);
%                 Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
%                 Lift_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cl*deltaR;
%                 Drag_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cd*deltaR;
%                 F_loc_R(:,i,j) = [(Lift_R(i,j)*sind(alfa)-Drag_R(i,j)*cosd(alfa)); 0; -(Lift_R(i,j)*cosd(alfa)+Drag_R(i,j)*sind(alfa))];
%                 M_loc_R(:,i,j) = cross(y_sect_R(:,i),F_loc_R(:,i,j));
% 
%             clear Cl Cd alfa
% 
%         end    
% 
%     end

    for j = 1:n
    
        for i = 1:nr_sect

            if up(j) == 1

                alfa = -radtodeg(alfa_L(i,j));
                Cl = 0.225+1.58*sind(2.13*alfa-7.2);
                Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
                Lift_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cl*deltaR; % [ kg/mm^3 * mm * mm^2/s^2 * mm ]
                Drag_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cd*deltaR;
                F_loc_L(:,i,j) = [(Lift_L(i,j)*sind(alfa)-Drag_L(i,j)*cosd(alfa)); 0; (Lift_L(i,j)*cosd(alfa)+Drag_L(i,j)*sind(alfa))];
                M_loc_L(:,i,j) = cross(y_sect_L(:,i),F_loc_L(:,i,j));

            elseif down(j) == 1

                alfa = radtodeg(alfa_L(i,j));
                Cl = 0.225+1.58*sind(2.13*alfa-7.2);
                Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
                Lift_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cl*deltaR;
                Drag_L(i,j) = 0.5*Rho*chords(i)*(u_L(1,i,j)^2+u_L(3,i,j)^2)*Cd*deltaR;
                F_loc_L(:,i,j) = [(Lift_L(i,j)*sind(alfa)-Drag_L(i,j)*cosd(alfa)); 0; -(Lift_L(i,j)*cosd(alfa)+Drag_L(i,j)*sind(alfa))];
                M_loc_L(:,i,j) = cross(y_sect_L(:,i),F_loc_L(:,i,j));

            end

            clear Cl Cd alfa

        end

        for i = 1:nr_sect

            if up(j) == 1

                alfa = -radtodeg(alfa_R(i,j));
                Cl = 0.225+1.58*sind(2.13*alfa-7.2);
                Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
                Lift_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cl*deltaR;
                Drag_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cd*deltaR;
                F_loc_R(:,i,j) = [(Lift_R(i,j)*sind(alfa)-Drag_R(i,j)*cosd(alfa)); 0; (Lift_R(i,j)*cosd(alfa)+Drag_R(i,j)*sind(alfa))];
                M_loc_R(:,i,j) = cross(y_sect_R(:,i),F_loc_R(:,i,j));

            elseif down(j) == 1

                alfa = radtodeg(alfa_R(i,j));
                Cl = 0.225+1.58*sind(2.13*alfa-7.2);
                Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
                Lift_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cl*deltaR;
                Drag_R(i,j) = 0.5*Rho*chords(i)*(u_R(1,i,j)^2+u_R(3,i,j)^2)*Cd*deltaR;
                F_loc_R(:,i,j) = [(Lift_R(i,j)*sind(alfa)-Drag_R(i,j)*cosd(alfa)); 0; -(Lift_R(i,j)*cosd(alfa)+Drag_R(i,j)*sind(alfa))];
                M_loc_R(:,i,j) = cross(y_sect_R(:,i),F_loc_R(:,i,j));

            end

            clear Cl Cd alfa

        end    

    end
    
    
    c_pres_L = [ 0; -0.75*wing_l; 0];
    c_pres_R = [ 0; 0.75*wing_l; 0];
    
    % Calculate the combined lift and drag forces at the center of
    % pressure:
    
    F_L = zeros(3,n);
    F_R = zeros(3,n);
    
    for i = 1:n
    
        F_L(:,i) = sum(F_loc_L(:,:,i),2); % Force in Newton
        F_R(:,i) = sum(F_loc_R(:,:,i),2);% Force in Newton
        
    end
    

    F_L_st = zeros(3,n);
    F_R_st = zeros(3,n);
    
    for i = 1:n
    
        F_L_st(:,i) = R_beta*R_L(:,:,i)'*F_L(:,i);
        F_R_st(:,i) = R_beta*R_R(:,:,i)'*F_R(:,i);
        
    end
    
    
    M_L_st = zeros(3,n);
    M_R_st = zeros(3,n);    
    
    arm_cp_L = zeros(3,n);
    arm_cp_R = zeros(3,n);
    
    for i = 1:n
        
        arm_cp_L(:,i) = R_beta*(arm_JL+R_L(:,:,i)*c_pres_L);
        arm_cp_R(:,i) = R_beta*(arm_JR+R_R(:,:,i)*c_pres_R);
        
        M_L_st(:,i) = cross(arm_cp_L(:,i),F_L_st(:,i));
        M_R_st(:,i) = cross(arm_cp_R(:,i),F_R_st(:,i));
        
    end
    
    Aero_f.u_L = u_L;
    Aero_f.u_R = u_R;
    Aero_f.c_pres_L = c_pres_L;
    Aero_f.c_pres_R = c_pres_R;
    Aero_f.F_L = F_L;
    Aero_f.F_R = F_R;
    Aero_f.F_L_st = F_L_st;
    Aero_f.F_R_st = F_R_st;
    Aero_f.arm_cp_L = arm_cp_L;
    Aero_f.arm_cp_R = arm_cp_R;
    Aero_f.M_L_st = M_L_st;
    Aero_f.M_R_st = M_R_st;
        
end

