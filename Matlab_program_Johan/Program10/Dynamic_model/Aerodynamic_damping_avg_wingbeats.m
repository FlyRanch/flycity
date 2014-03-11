function [ F_exp, M_exp ] = Aerodynamic_damping_avg_wingbeats( settings, pathDB, a_avg, body_kin, avg_name, n_points, dt2, seq_nr )


    experiment = {};
    
    nr_of_exp = length(body_kin.phi_b);
    
    F_exp = zeros(18,nr_of_exp);
       
    M_exp = zeros(15,nr_of_exp);    
        
    a_theta_L = [ a_avg.theta_LR1; a_avg.theta_LR2 ];
    a_eta_L = [ a_avg.eta_LR1; a_avg.eta_LR2 ];
    a_phi_L = [ a_avg.phi_LR1; a_avg.phi_LR2 ];

    a_theta_R = [ a_avg.theta_LR1; a_avg.theta_LR2 ];
    a_eta_R = [ a_avg.eta_LR1; a_avg.eta_LR2 ];
    a_phi_R = [ a_avg.phi_LR1; a_avg.phi_LR2 ];
           
    a_sim.theta_L = a_theta_L;
    a_sim.eta_L = a_eta_L;
    a_sim.phi_L = a_phi_L;
           
    a_sim.theta_R = a_theta_R;
    a_sim.eta_R = a_eta_R;
    a_sim.phi_R = a_phi_R;
    
    f_vect = a_avg.f;
       
    down_up = a_avg.down_up;
    
    t_vect = a_avg.t_vect;
    
    for j = 1:nr_of_exp

            
       phi_b = body_kin.phi_b(j);
       theta_b = body_kin.theta_b(j);
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];
              
       v_st_c = body_kin.v_st_c(:,j);
       a_st_c = body_kin.a_st_c(:,j);
       w_st_c = body_kin.w_st_c(:,j);
       w_dot_st_c = body_kin.w_dot_st_c(:,j);
           
       q_b = zeros(4,n_points);
       v_st = zeros(3,n_points);
       a_st = zeros(3,n_points);
       w_st = zeros(3,n_points);
       w_dot_st = zeros(3,n_points);
           
       for k = 1:n_points
               
         q_b(:,k) = q_b_c;
         v_st(:,k) = v_st_c;
         a_st(:,k) = a_st_c;
         w_st(:,k) = w_st_c;
         w_dot_st(:,k) = w_dot_st_c;
               
       end

       Req_f = Required_forces( settings, pathDB, a_sim, n_points, t_vect, down_up, f_vect, dt2, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );

           
       t_man = Req_f.time;
       
       F_I_acc  = Req_f.F_I_acc;
       F_I_vel  = Req_f.F_I_vel;
       M_I_acc  = Req_f.M_I_acc;
       M_I_vel  = Req_f.M_I_vel;
       F_I      = Req_f.F_I;
       M_I      = Req_f.M_I;
       c_pres_L = Req_f.c_pres_L;
       c_pres_R = Req_f.c_pres_R;
       F_a_wL   = Req_f.F_a_wL;
       F_a_wR   = Req_f.F_a_wR;
       F_a_L    = Req_f.F_a_L;
       F_a_R    = Req_f.F_a_R;
       M_a_L    = Req_f.M_a_L;
       M_a_R    = Req_f.M_a_R;
       F_g      = Req_f.F_g;
       
       
        F_I_mean = Req_f.F_I_mean;
        F_I_acc_mean = Req_f.F_I_acc_mean;
        F_I_vel_mean = Req_f.F_I_vel_mean;
        F_a_L_mean = Req_f.F_a_L_mean;
        F_a_R_mean = Req_f.F_a_R_mean;
        F_g_mean = Req_f.F_g_mean;
        
        M_I_mean = Req_f.M_I_mean;
        M_I_acc_mean = Req_f.M_I_acc_mean;
        M_I_vel_mean = Req_f.M_I_vel_mean;
        M_a_L_mean = Req_f.M_a_L_mean;
        M_a_R_mean = Req_f.M_a_R_mean;

        F_exp(1:3,j) = F_I_mean;
        F_exp(4:6,j) = F_I_acc_mean;
        F_exp(7:9,j) = F_I_vel_mean;
        F_exp(10:12,j) = F_a_L_mean;
        F_exp(13:15,j) = F_a_R_mean;
        F_exp(16:18,j) = F_g_mean;
            
        M_exp(1:3,j) = M_I_mean;
        M_exp(4:6,j) = M_I_acc_mean;
        M_exp(7:9,j) = M_I_vel_mean;
        M_exp(10:12,j) = M_a_L_mean;
        M_exp(13:15,j) = M_a_R_mean;

        
    end
       

    
end

