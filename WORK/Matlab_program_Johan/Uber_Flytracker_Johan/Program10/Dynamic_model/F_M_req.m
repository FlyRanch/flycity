function [ F_M_r ] = F_M_req( settings, pathDB, all_maneuvers )


    F_M_r = {};
   
    roll_maneuvers = all_maneuvers.roll_maneuvers;
    pitch_maneuvers = all_maneuvers.pitch_maneuvers;
    yaw_maneuvers = all_maneuvers.yaw_maneuvers;
    ax_maneuvers = all_maneuvers.ax_maneuvers;
    ay_maneuvers = all_maneuvers.ay_maneuvers;
    az_maneuvers = all_maneuvers.az_maneuvers;
    pure_roll = all_maneuvers.pure_roll;
    pure_pitch = all_maneuvers.pure_pitch;
    pure_yaw = all_maneuvers.pure_yaw;
    pure_ax = all_maneuvers.pure_ax;
    pure_ay = all_maneuvers.pure_ay;
    pure_az = all_maneuvers.pure_az;
    
    roll_man_name = fieldnames(roll_maneuvers);
    pitch_man_name = fieldnames(pitch_maneuvers);
    yaw_man_name = fieldnames(yaw_maneuvers);
    ax_man_name = fieldnames(ax_maneuvers);
    ay_man_name = fieldnames(ay_maneuvers);
    az_man_name = fieldnames(az_maneuvers);
    pure_roll_name = fieldnames(pure_roll);
    pure_pitch_name = fieldnames(pure_pitch);
    pure_yaw_name = fieldnames(pure_yaw);
    pure_ax_name = fieldnames(pure_ax);
    pure_ay_name = fieldnames(pure_ay);
    pure_az_name = fieldnames(pure_az);
    
    nr_roll_man = size(roll_man_name,1);
    nr_pitch_man = size(pitch_man_name,1);
    nr_yaw_man = size(yaw_man_name,1);
    nr_ax_man = size(ax_man_name,1);
    nr_ay_man = size(ay_man_name,1);
    nr_az_man = size(az_man_name,1);
    nr_pure_roll = size(pure_roll_name,1);
    nr_pure_pitch = size(pure_pitch_name,1);
    nr_pure_yaw = size(pure_yaw_name,1);
    nr_pure_ax = size(pure_ax_name,1);
    nr_pure_ay = size(pure_ay_name,1);
    nr_pure_az = size(pure_az_name,1);
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    %%
    
    roll_man = {};

    'roll maneuvers'
    
    for i = 1:nr_roll_man
        
        i
        
        temp_man = roll_maneuvers.(char(roll_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);

           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            roll_man.(char(roll_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;            
            
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        roll_man.(char(roll_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        roll_man.(char(roll_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        roll_man.(char(roll_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        roll_man.(char(roll_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        roll_man.(char(roll_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        roll_man.(char(roll_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;        
                
    end
        
    %%
    pitch_man = {};
    
    'pitch maneuvers'
    
    for i = 1:nr_pitch_man
        
        i
        
        temp_man = pitch_maneuvers.(char(pitch_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);

           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pitch_man.(char(pitch_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean; 
            
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pitch_man.(char(pitch_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end


        
    
    %%
    
    yaw_man = {};
    
    'yaw maneuvers'
    
    for i = 1:nr_yaw_man
        
        i
        
        temp_man = yaw_maneuvers.(char(yaw_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            yaw_man.(char(yaw_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean; 
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        yaw_man.(char(yaw_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end


        
    
    
    %%
    
    ax_man = {};
    
    'a_x maneuvers'
    
    for i = 1:nr_ax_man
        
        i
        
        temp_man = ax_maneuvers.(char(ax_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            ax_man.(char(ax_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean; 
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        ax_man.(char(ax_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        ax_man.(char(ax_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        ax_man.(char(ax_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        ax_man.(char(ax_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        ax_man.(char(ax_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        ax_man.(char(ax_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
        
       
                
    end


    
    
    
    %%
    
    ay_man = {};
    
    'a_y maneuvers'
    
    for i = 1:nr_ay_man
        
        i
        
        temp_man = ay_maneuvers.(char(ay_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            ay_man.(char(ay_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean; 
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        ay_man.(char(ay_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        ay_man.(char(ay_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        ay_man.(char(ay_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        ay_man.(char(ay_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        ay_man.(char(ay_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        ay_man.(char(ay_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end

 
    
    
    
    %%
    
    az_man = {};
    
    'a_z maneuvers'
    
    for i = 1:nr_az_man
        
        i
        
        temp_man = az_maneuvers.(char(az_man_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            az_man.(char(az_man_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean; 
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        az_man.(char(az_man_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        az_man.(char(az_man_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        az_man.(char(az_man_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        az_man.(char(az_man_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        az_man.(char(az_man_name(i))).F_M_avg.F_I      = Req_f.F_I;
        az_man.(char(az_man_name(i))).F_M_avg.M_I      = Req_f.M_I;
        az_man.(char(az_man_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        az_man.(char(az_man_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        az_man.(char(az_man_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        az_man.(char(az_man_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        az_man.(char(az_man_name(i))).F_M_avg.F_g      = Req_f.F_g;

        az_man.(char(az_man_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        az_man.(char(az_man_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        az_man.(char(az_man_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        az_man.(char(az_man_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        az_man.(char(az_man_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        az_man.(char(az_man_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        az_man.(char(az_man_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        az_man.(char(az_man_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        az_man.(char(az_man_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        az_man.(char(az_man_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end
    
    
    
    %%
    
    
    %%
    
    pure_roll_man = {};
    
    'pure roll'
    
    for i = 1:nr_pure_roll
        
        i
        
        temp_man = pure_roll.(char(pure_roll_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_roll_man.(char(pure_roll_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_roll_man.(char(pure_roll_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
        
        
        
        
                
    end


    
    
    
    %%
    
    
    pure_pitch_man = {};
    
    'pure pitch'
    
    for i = 1:nr_pure_pitch
        
        i
        
        temp_man = pure_pitch.(char(pure_pitch_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_pitch_man.(char(pure_pitch_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_pitch_man.(char(pure_pitch_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end


    
    
    
    %%
    
    
    pure_yaw_man = {};
    
    'pure yaw'
    
    for i = 1:nr_pure_yaw
        
        i
        
        temp_man = pure_yaw.(char(pure_yaw_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_yaw_man.(char(pure_yaw_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_yaw_man.(char(pure_yaw_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end



    
    
    %%
    
    
    pure_ax_man = {};
    
    'pure a_x'
    
    for i = 1:nr_pure_ax
        
        i
        
        temp_man = pure_ax.(char(pure_ax_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_ax_man.(char(pure_ax_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_ax_man.(char(pure_ax_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
        
        
        

                
    end


    
    
    
    %%
    
    
    pure_ay_man = {};
    
    'pure a_y'
    
    for i = 1:nr_pure_ay
        
        i
        
        temp_man = pure_ay.(char(pure_ay_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_ay_man.(char(pure_ay_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_ay_man.(char(pure_ay_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
        
        

                
    end
    
    
    
    %%
    
    
    pure_az_man = {};
    
    'pure a_z'
    
    for i = 1:nr_pure_az
        
        i
        
        temp_man = pure_az.(char(pure_az_name(i)));
       
        seq_nr = temp_man.seq_nr;       
        
        f_vect = temp_man.f;
        
        nr_man_wb = length(f_vect);
        
        
        for j = 1:nr_man_wb
            
            nr_points_vect = round((1./f_vect(j))./dt)+1;
            
            t_vect = [0 (nr_points_vect-1)*dt];

            down_up_vect = temp_man.down_up(j);
            
            a_sim = {};
            
            a_sim.theta_L = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_L1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_L2(:,j)];
            a_sim.eta_L = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_L1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_L2(:,j)];
            a_sim.phi_L = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_L1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_L2(:,j)];

            a_sim.theta_R = [temp_man.a_avg.theta_LR1+temp_man.a_dev.theta_R1(:,j); temp_man.a_avg.theta_LR2+temp_man.a_dev.theta_R2(:,j)];
            a_sim.eta_R = [temp_man.a_avg.eta_LR1+temp_man.a_dev.eta_R1(:,j); temp_man.a_avg.eta_LR2+temp_man.a_dev.eta_R2(:,j)];
            a_sim.phi_R = [temp_man.a_avg.phi_LR1+temp_man.a_dev.phi_R1(:,j); temp_man.a_avg.phi_LR2+temp_man.a_dev.phi_R2(:,j)];

            q_b_t = temp_man.q_body_wb(:,j);
            v_st_t = temp_man.V_strk(:,j);
            a_st_t = temp_man.A_strk(:,j);
            w_st_t = temp_man.Omega_strk(:,j);
            w_dot_st_t = temp_man.Omega_dot_strk(:,j);

            q_b = zeros(4,nr_points_vect);
            v_st = zeros(3,nr_points_vect);
            a_st = zeros(3,nr_points_vect);
            w_st = zeros(3,nr_points_vect);
            w_dot_st = zeros(3,nr_points_vect);
            
            q_b(1,:) = ones(1,nr_points_vect).*q_b_t(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_t(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_t(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_t(4);
            v_st(1,:) = ones(1,nr_points_vect).*v_st_t(1);
            v_st(2,:) = ones(1,nr_points_vect).*v_st_t(2);
            v_st(3,:) = ones(1,nr_points_vect).*v_st_t(3);
            a_st(1,:) = ones(1,nr_points_vect).*a_st_t(1);
            a_st(2,:) = ones(1,nr_points_vect).*a_st_t(2);
            a_st(3,:) = ones(1,nr_points_vect).*a_st_t(3);
            w_st(1,:) = ones(1,nr_points_vect).*w_st_t(1);
            w_st(2,:) = ones(1,nr_points_vect).*w_st_t(2);
            w_st(3,:) = ones(1,nr_points_vect).*w_st_t(3);
            w_dot_st(1,:) = ones(1,nr_points_vect).*w_dot_st_t(1);
            w_dot_st(2,:) = ones(1,nr_points_vect).*w_dot_st_t(2);
            w_dot_st(3,:) = ones(1,nr_points_vect).*w_dot_st_t(3);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
            

            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).t_man = Req_f.time;

            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_acc  = Req_f.F_I_acc;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_vel  = Req_f.F_I_vel;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_acc  = Req_f.M_I_acc;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_vel  = Req_f.M_I_vel;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I      = Req_f.F_I;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I      = Req_f.M_I;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).c_pres_L = Req_f.c_pres_L;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).c_pres_R = Req_f.c_pres_R;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_wL   = Req_f.F_a_wL;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_wR   = Req_f.F_a_wR;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_L    = Req_f.F_a_L;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_R    = Req_f.F_a_R;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_L    = Req_f.M_a_L;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_R    = Req_f.M_a_R;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_g      = Req_f.F_g;

            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_mean      = Req_f.F_I_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean  = Req_f.F_I_acc_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean  = Req_f.F_I_vel_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean    = Req_f.F_a_L_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean    = Req_f.F_a_R_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_g_mean      = Req_f.F_g_mean;

            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_mean      = Req_f.M_I_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean  = Req_f.M_I_acc_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean  = Req_f.M_I_vel_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean    = Req_f.M_a_L_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean    = Req_f.M_a_R_mean;
            
            
            
           phi_b = pi;
           theta_b = (55/180)*pi;

           q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                      cos(phi_b/2)*sin(theta_b/2); ...
                     -sin(phi_b/2)*sin(theta_b/2); ...
                      cos(phi_b/2)*cos(theta_b/2)];

            q_b(1,:) = ones(1,nr_points_vect).*q_b_c(1);
            q_b(2,:) = ones(1,nr_points_vect).*q_b_c(2);
            q_b(3,:) = ones(1,nr_points_vect).*q_b_c(3);
            q_b(4,:) = ones(1,nr_points_vect).*q_b_c(4);
            
        
            [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect(j), dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );            
            
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_mean_ss     = Req_f.F_I_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_acc_mean_ss  = Req_f.F_I_acc_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_I_vel_mean_ss  = Req_f.F_I_vel_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_L_mean_ss    = Req_f.F_a_L_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_a_R_mean_ss    = Req_f.F_a_R_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).F_g_mean_ss      = Req_f.F_g_mean;

            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_mean_ss      = Req_f.M_I_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_acc_mean_ss  = Req_f.M_I_acc_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_I_vel_mean_ss  = Req_f.M_I_vel_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_L_mean_ss    = Req_f.M_a_L_mean;
            pure_az_man.(char(pure_az_name(i))).(char(['wb_' int2str(j)])).M_a_R_mean_ss    = Req_f.M_a_R_mean;
            
        end
        
        % Compute hovering flight forces and moments:
        
        f_avg = temp_man.f_avg;
        
        nr_points_avg = round((1./f_avg)./dt)+1;
           
        t_avg = [0 (nr_points_avg-1)*dt];

        down_up_avg = temp_man.down_up_avg;
            
        a_sim = {};
            
        a_sim.theta_L = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_L = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_L = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];

        a_sim.theta_R = [temp_man.a_avg.theta_LR1; temp_man.a_avg.theta_LR2];
        a_sim.eta_R = [temp_man.a_avg.eta_LR1; temp_man.a_avg.eta_LR2];
        a_sim.phi_R = [temp_man.a_avg.phi_LR1; temp_man.a_avg.phi_LR2];      
        
        q_b = zeros(4,nr_points_avg);
        v_st = zeros(3,nr_points_avg);
        a_st = zeros(3,nr_points_avg);
        w_st = zeros(3,nr_points_avg);
        w_dot_st = zeros(3,nr_points_avg);
        
       phi_b = pi;
       theta_b = (55/180)*pi;
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];

        q_b(1,:) = ones(1,nr_points_avg).*q_b_c(1);
        q_b(2,:) = ones(1,nr_points_avg).*q_b_c(2);
        q_b(3,:) = ones(1,nr_points_avg).*q_b_c(3);
        q_b(4,:) = ones(1,nr_points_avg).*q_b_c(4);
              
        [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_avg, t_avg, down_up_avg, f_avg, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
        
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I_acc  = Req_f.F_I_acc;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I_vel  = Req_f.F_I_vel;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I_acc  = Req_f.M_I_acc;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I_vel  = Req_f.M_I_vel;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I      = Req_f.F_I;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I      = Req_f.M_I;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.c_pres_L = Req_f.c_pres_L;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.c_pres_R = Req_f.c_pres_R;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_wL   = Req_f.F_a_wL;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_wR   = Req_f.F_a_wR;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_L    = Req_f.F_a_L;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_R    = Req_f.F_a_R;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_a_L    = Req_f.M_a_L;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_a_R    = Req_f.M_a_R;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_g      = Req_f.F_g;

        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I_mean      = Req_f.F_I_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I_acc_mean  = Req_f.F_I_acc_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_I_vel_mean  = Req_f.F_I_vel_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_L_mean    = Req_f.F_a_L_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_a_R_mean    = Req_f.F_a_R_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.F_g_mean      = Req_f.F_g_mean;

        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I_mean      = Req_f.M_I_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I_acc_mean  = Req_f.M_I_acc_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_I_vel_mean  = Req_f.M_I_vel_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_a_L_mean    = Req_f.M_a_L_mean;
        pure_az_man.(char(pure_az_name(i))).F_M_avg.M_a_R_mean    = Req_f.M_a_R_mean;
                
    end 
    
    
    %%
    
        F_M_r.roll_maneuver = roll_man;
        F_M_r.pitch_maneuver = pitch_man;
        F_M_r.yaw_maneuver = yaw_man;
        F_M_r.ax_maneuver = ax_man;
        F_M_r.ay_maneuver = ay_man;
        F_M_r.az_maneuver = az_man;
        F_M_r.pure_roll_maneuver = pure_roll_man;
        F_M_r.pure_pitch_maneuver = pure_pitch_man;
        F_M_r.pure_yaw_maneuver = pure_yaw_man;
        F_M_r.pure_ax_maneuver = pure_ax_man;
        F_M_r.pure_ay_maneuver = pure_ay_man;
        F_M_r.pure_az_maneuver = pure_az_man;
        
        

end

