function Aerodynamic_damping( settings, pathDB )


    % Test the effects of aerodynamic damping on several different
    % wingbeats:
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    dt2 = dt;
    
    nr_of_seq = size(pathDB.x,2);   
    
    
    
    % Set up different experiments:----------------------------------------
    
    body_kin = {};
    
    a_grid = [-30; -25; 20; 15; 10; 5; 0; 5; 10; 15; 20; 25; 30];
    
    w_dot_grid = [-30000; -20000; -15000; -10000; -5000; 0; 5000; 10000; 15000; 20000; 30000];
    
    v_grid = [-1500; -1250; -1000; -750; -500; -250; 0; 250; 500; 750; 1000; 1250; 1500];
    
    w_grid = [-150; -125; -100; -75; -50; -25; 0; 25; 50; 75; 100; 125; 150];
    
    n_v_points = length(v_grid);
    n_w_points = length(w_grid);
    n_a_points = length(a_grid);
    n_w_dot_points = length(w_dot_grid);
    
    v_exps = 3;
    w_exps = 3;
    a_exps = 3;
    w_dot_exps = 3;
    
    body_kin.phi_b = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    body_kin.theta_b = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    body_kin.v_st_c = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    body_kin.a_st_c = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    body_kin.w_st_c = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    body_kin.w_dot_st_c = zeros(1,1+v_exps*n_v_points+w_exps*n_w_points+a_exps*n_a_points+w_dot_exps*n_w_dot_points);
    
    % All experiments take place at a body pitch angle of 55 degrees:
    
    body_kin.phi_b(:) = pi;
    body_kin.theta_b(:) = -(55/180)*pi;
    
    % Experiment 1:
    
    exp_range1 = 2:(n_v_points+1);
    body_kin.v_st_c(1,exp_range1) = v_grid;
    
    % Experiment 2:
    
    exp_range2 = (exp_range1(end)+1):(exp_range1(end)+n_v_points);
    body_kin.v_st_c(2,exp_range2) = v_grid;
    
    % Experiment 3:
    
    exp_range3 = (exp_range2(end)+1):(exp_range2(end)+n_v_points);
    body_kin.v_st_c(3,exp_range3) = v_grid;
    
    % Experiment 4:
    
    exp_range4 = (exp_range3(end)+1):(exp_range3(end)+n_w_points);
    body_kin.w_st_c(1,exp_range4) = w_grid;
    
    % Experiment 5:
    
    exp_range5 = (exp_range4(end)+1):(exp_range4(end)+n_w_points);
    body_kin.w_st_c(2,exp_range5) = w_grid;
    
    % Experiment 6:
    
    exp_range6 = (exp_range5(end)+1):(exp_range5(end)+n_w_points);
    body_kin.w_st_c(3,exp_range6) = w_grid;
    
    % Experiment 7:
    
    exp_range7 = (exp_range6(end)+1):(exp_range6(end)+n_a_points);
    body_kin.a_st_c(1,exp_range7) = a_grid;
    
    % Experiment 8:
    
    exp_range8 = (exp_range7(end)+1):(exp_range7(end)+n_a_points);
    body_kin.a_st_c(2,exp_range8) = a_grid;
    
    % Experiment 9:
    
    exp_range9 = (exp_range8(end)+1):(exp_range8(end)+n_a_points);
    body_kin.a_st_c(3,exp_range9) = a_grid;
    
    % Experiment 10:
    
    exp_range10 = (exp_range9(end)+1):(exp_range9(end)+n_w_dot_points);
    body_kin.w_dot_st_c(1,exp_range10) = w_dot_grid;
    
    % Experiment 11:
    
    exp_range11 = (exp_range10(end)+1):(exp_range10(end)+n_w_dot_points);
    body_kin.w_dot_st_c(2,exp_range11) = w_dot_grid;
    
    % Experiment 12:
    
    exp_range12 = (exp_range11(end)+1):(exp_range11(end)+n_w_dot_points);
    body_kin.w_dot_st_c(3,exp_range12) = w_dot_grid;
    
    
    
    body_kin.exp_range1 = exp_range1;
    body_kin.exp_range2 = exp_range2;
    body_kin.exp_range3 = exp_range3;
    body_kin.exp_range4 = exp_range4;
    body_kin.exp_range5 = exp_range5;
    body_kin.exp_range6 = exp_range6;
    body_kin.exp_range7 = exp_range7;
    body_kin.exp_range8 = exp_range8;
    body_kin.exp_range9 = exp_range9;
    body_kin.exp_range10 = exp_range10;
    body_kin.exp_range11 = exp_range11;
    body_kin.exp_range12 = exp_range12;

    %----------------------------------------------------------------------
    
    
    
    % Test aerodynamic damping and inertial damping effects on average
    % wingbeats:-----------------------------------------------------------
    
    a_avg = pathDB.a_avg;
    
    avg_name = fieldnames(a_avg);
    
    experiment = {};
    
    for i = 1:nr_of_seq
        
        i
        
        temp_avg = a_avg.(char(avg_name(i)));
        
        if isstruct(temp_avg) == 1
        
            a_sim = {};
            
            a_sim.theta_LR1 = temp_avg.theta_LR1;
            a_sim.theta_LR2 = temp_avg.theta_LR2;
            a_sim.eta_LR1 = temp_avg.eta_LR1;
            a_sim.eta_LR2 = temp_avg.eta_LR2;
            a_sim.phi_LR1 = temp_avg.phi_LR1;
            a_sim.phi_LR2 = temp_avg.phi_LR2;
            
            a_sim.f = pathDB.f_avg.(char(['f_avg_' int2str(i)]));
            
            a_sim.down_up = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i)]));
            
            n_points = round((1/a_sim.f)/dt2)+1;
            
            a_sim.t_vect = [0 n_points*dt2];
            
            [ F_exp, M_exp ] = Aerodynamic_damping_avg_wingbeats( settings, pathDB, a_sim, body_kin, avg_name(i), n_points, dt2, i );
            
            experiment.(char(['exp_' int2str(i)])).F_exp = F_exp ;
            experiment.(char(['exp_' int2str(i)])).M_exp = M_exp ;
            experiment.(char(['exp_' int2str(i)])).w_length = pathDB.wing_l(i) ;
            experiment.(char(['exp_' int2str(i)])).freq = a_sim.f ;
            
        end

        
    end
    
    %----------------------------------------------------------------------
    
    
    exp_names = fieldnames(experiment)
    
    nr_exp = length(exp_names)

    
    C_Fx_vx = zeros(nr_exp,n_v_points);
    C_Fy_vx = zeros(nr_exp,n_v_points);
    C_Fz_vx = zeros(nr_exp,n_v_points);
    
    C_Fx_vy = zeros(nr_exp,n_v_points);
    C_Fy_vy = zeros(nr_exp,n_v_points);
    C_Fz_vy = zeros(nr_exp,n_v_points);
    
    C_Fx_vz = zeros(nr_exp,n_v_points);
    C_Fy_vz = zeros(nr_exp,n_v_points);
    C_Fz_vz = zeros(nr_exp,n_v_points);
    
    C_Mx_wx = zeros(nr_exp,n_w_points);
    C_My_wx = zeros(nr_exp,n_w_points);
    C_Mz_wx = zeros(nr_exp,n_w_points);
    
    C_Mx_wy = zeros(nr_exp,n_w_points);
    C_My_wy = zeros(nr_exp,n_w_points);
    C_Mz_wy = zeros(nr_exp,n_w_points);
    
    C_Mx_wz = zeros(nr_exp,n_w_points);
    C_My_wz = zeros(nr_exp,n_w_points);
    C_Mz_wz = zeros(nr_exp,n_w_points);    
    
    C_Fx_ax = zeros(nr_exp,n_a_points);
    C_Fy_ax = zeros(nr_exp,n_a_points);
    C_Fz_ax = zeros(nr_exp,n_a_points);
    
    C_Fx_ay = zeros(nr_exp,n_a_points);
    C_Fy_ay = zeros(nr_exp,n_a_points);
    C_Fz_ay = zeros(nr_exp,n_a_points);
    
    C_Fx_az = zeros(nr_exp,n_a_points);
    C_Fy_az = zeros(nr_exp,n_a_points);
    C_Fz_az = zeros(nr_exp,n_a_points);
    
    C_Mx_w_dot_x = zeros(nr_exp,n_w_dot_points);
    C_My_w_dot_x = zeros(nr_exp,n_w_dot_points);
    C_Mz_w_dot_x = zeros(nr_exp,n_w_dot_points);
    
    C_Mx_w_dot_y = zeros(nr_exp,n_w_dot_points);
    C_My_w_dot_y = zeros(nr_exp,n_w_dot_points);
    C_Mz_w_dot_y = zeros(nr_exp,n_w_dot_points);
    
    C_Mx_w_dot_z = zeros(nr_exp,n_w_dot_points);
    C_My_w_dot_z = zeros(nr_exp,n_w_dot_points);
    C_Mz_w_dot_z = zeros(nr_exp,n_w_dot_points);

    for j = 1:nr_exp
        
        exp_mat = experiment.(char(exp_names(j)));
        loc_wl = exp_mat.w_length;
        loc_f = exp_mat.freq;
        
        C_Fx_vx(j,:) = (exp_mat.F_exp(10,exp_range1)+exp_mat.F_exp(13,exp_range1)-exp_mat.F_exp(10,1)-exp_mat.F_exp(13,1))/(loc_f^2*loc_wl^4);
        C_Fy_vx(j,:) = (exp_mat.F_exp(11,exp_range1)+exp_mat.F_exp(14,exp_range1)-exp_mat.F_exp(11,1)-exp_mat.F_exp(14,1))/(loc_f^2*loc_wl^4);
        C_Fz_vx(j,:) = (exp_mat.F_exp(12,exp_range1)+exp_mat.F_exp(15,exp_range1)-exp_mat.F_exp(12,1)-exp_mat.F_exp(15,1))/(loc_f^2*loc_wl^4);

        C_Fx_vy(j,:) = (exp_mat.F_exp(10,exp_range2)+exp_mat.F_exp(13,exp_range2)-exp_mat.F_exp(10,1)-exp_mat.F_exp(13,1))/(loc_f^2*loc_wl^4);
        C_Fy_vy(j,:) = (exp_mat.F_exp(11,exp_range2)+exp_mat.F_exp(14,exp_range2)-exp_mat.F_exp(11,1)-exp_mat.F_exp(14,1))/(loc_f^2*loc_wl^4);
        C_Fz_vy(j,:) = (exp_mat.F_exp(12,exp_range2)+exp_mat.F_exp(15,exp_range2)-exp_mat.F_exp(12,1)-exp_mat.F_exp(15,1))/(loc_f^2*loc_wl^4);
        
        C_Fx_vz(j,:) = (exp_mat.F_exp(10,exp_range3)+exp_mat.F_exp(13,exp_range3)-exp_mat.F_exp(10,1)-exp_mat.F_exp(13,1))/(loc_f^2*loc_wl^4);
        C_Fy_vz(j,:) = (exp_mat.F_exp(11,exp_range3)+exp_mat.F_exp(14,exp_range3)-exp_mat.F_exp(11,1)-exp_mat.F_exp(14,1))/(loc_f^2*loc_wl^4);
        C_Fz_vz(j,:) = (exp_mat.F_exp(12,exp_range3)+exp_mat.F_exp(15,exp_range3)-exp_mat.F_exp(12,1)-exp_mat.F_exp(15,1))/(loc_f^2*loc_wl^4);

        C_Mx_wx(j,:) = (exp_mat.M_exp(10,exp_range4)+exp_mat.M_exp(13,exp_range4)-exp_mat.M_exp(10,1)-exp_mat.M_exp(13,1))/(loc_f^2*loc_wl^5);
        C_My_wx(j,:) = (exp_mat.M_exp(11,exp_range4)+exp_mat.M_exp(14,exp_range4)-exp_mat.M_exp(11,1)-exp_mat.M_exp(14,1))/(loc_f^2*loc_wl^5);
        C_Mz_wx(j,:) = (exp_mat.M_exp(12,exp_range4)+exp_mat.M_exp(15,exp_range4)-exp_mat.M_exp(12,1)-exp_mat.M_exp(15,1))/(loc_f^2*loc_wl^5);

        C_Mx_wy(j,:) = (exp_mat.M_exp(10,exp_range5)+exp_mat.M_exp(13,exp_range5)-exp_mat.M_exp(10,1)-exp_mat.M_exp(13,1))/(loc_f^2*loc_wl^5);
        C_My_wy(j,:) = (exp_mat.M_exp(11,exp_range5)+exp_mat.M_exp(14,exp_range5)-exp_mat.M_exp(11,1)-exp_mat.M_exp(14,1))/(loc_f^2*loc_wl^5);
        C_Mz_wy(j,:) = (exp_mat.M_exp(12,exp_range5)+exp_mat.M_exp(15,exp_range5)-exp_mat.M_exp(12,1)-exp_mat.M_exp(15,1))/(loc_f^2*loc_wl^5);

        C_Mx_wz(j,:) = (exp_mat.M_exp(10,exp_range6)+exp_mat.M_exp(13,exp_range6)-exp_mat.M_exp(10,1)-exp_mat.M_exp(13,1))/(loc_f^2*loc_wl^5);
        C_My_wz(j,:) = (exp_mat.M_exp(11,exp_range6)+exp_mat.M_exp(14,exp_range6)-exp_mat.M_exp(11,1)-exp_mat.M_exp(14,1))/(loc_f^2*loc_wl^5);
        C_Mz_wz(j,:) = (exp_mat.M_exp(12,exp_range6)+exp_mat.M_exp(15,exp_range6)-exp_mat.M_exp(12,1)-exp_mat.M_exp(15,1))/(loc_f^2*loc_wl^5);

        C_Fx_ax(j,:) = (exp_mat.F_exp(1,exp_range7)-exp_mat.F_exp(1,1))/(loc_f^2*loc_wl^4);
        C_Fy_ax(j,:) = (exp_mat.F_exp(2,exp_range7)-exp_mat.F_exp(2,1))/(loc_f^2*loc_wl^4);
        C_Fz_ax(j,:) = (exp_mat.F_exp(3,exp_range7)-exp_mat.F_exp(3,1))/(loc_f^2*loc_wl^4);

        C_Fx_ay(j,:) = (exp_mat.F_exp(1,exp_range8)-exp_mat.F_exp(1,1))/(loc_f^2*loc_wl^4);
        C_Fy_ay(j,:) = (exp_mat.F_exp(2,exp_range8)-exp_mat.F_exp(2,1))/(loc_f^2*loc_wl^4);
        C_Fz_ay(j,:) = (exp_mat.F_exp(3,exp_range8)-exp_mat.F_exp(3,1))/(loc_f^2*loc_wl^4);

        C_Fx_az(j,:) = (exp_mat.F_exp(1,exp_range9)-exp_mat.F_exp(1,1))/(loc_f^2*loc_wl^4);
        C_Fy_az(j,:) = (exp_mat.F_exp(2,exp_range9)-exp_mat.F_exp(2,1))/(loc_f^2*loc_wl^4);
        C_Fz_az(j,:) = (exp_mat.F_exp(3,exp_range9)-exp_mat.F_exp(3,1))/(loc_f^2*loc_wl^4);

        C_Mx_w_dot_x(j,:) = (exp_mat.M_exp(1,exp_range10)-exp_mat.M_exp(1,1))/(loc_f^2*loc_wl^5);
        C_My_w_dot_x(j,:) = (exp_mat.M_exp(2,exp_range10)-exp_mat.M_exp(2,1))/(loc_f^2*loc_wl^5);
        C_Mz_w_dot_x(j,:) = (exp_mat.M_exp(3,exp_range10)-exp_mat.M_exp(3,1))/(loc_f^2*loc_wl^5);

        C_Mx_w_dot_y(j,:) = (exp_mat.M_exp(1,exp_range11)-exp_mat.M_exp(1,1))/(loc_f^2*loc_wl^5);
        C_My_w_dot_y(j,:) = (exp_mat.M_exp(2,exp_range11)-exp_mat.M_exp(2,1))/(loc_f^2*loc_wl^5);
        C_Mz_w_dot_y(j,:) = (exp_mat.M_exp(3,exp_range11)-exp_mat.M_exp(3,1))/(loc_f^2*loc_wl^5);

        C_Mx_w_dot_z(j,:) = (exp_mat.M_exp(1,exp_range12)-exp_mat.M_exp(1,1))/(loc_f^2*loc_wl^5);
        C_My_w_dot_z(j,:) = (exp_mat.M_exp(2,exp_range12)-exp_mat.M_exp(2,1))/(loc_f^2*loc_wl^5);
        C_Mz_w_dot_z(j,:) = (exp_mat.M_exp(3,exp_range12)-exp_mat.M_exp(3,1))/(loc_f^2*loc_wl^5);
        
    end
    
    
    
    C_Fx_vx_p = polyfit(body_kin.v_st_c(1,exp_range1),mean(C_Fx_vx),3);
    C_Fy_vx_p = polyfit(body_kin.v_st_c(1,exp_range1),mean(C_Fy_vx),3);
    C_Fz_vx_p = polyfit(body_kin.v_st_c(1,exp_range1),mean(C_Fz_vx),3);
    
    C_Fx_vy_p = polyfit(body_kin.v_st_c(2,exp_range2),mean(C_Fx_vy),3);
    C_Fy_vy_p = polyfit(body_kin.v_st_c(2,exp_range2),mean(C_Fy_vy),3);
    C_Fz_vy_p = polyfit(body_kin.v_st_c(2,exp_range2),mean(C_Fz_vy),3);
    
    C_Fx_vz_p = polyfit(body_kin.v_st_c(3,exp_range3),mean(C_Fx_vz),3);
    C_Fy_vz_p = polyfit(body_kin.v_st_c(3,exp_range3),mean(C_Fy_vz),3);
    C_Fz_vz_p = polyfit(body_kin.v_st_c(3,exp_range3),mean(C_Fz_vz),3);
    
    C_Mx_wx_p = polyfit(body_kin.w_st_c(1,exp_range4),mean(C_Mx_wx),3);
    C_My_wx_p = polyfit(body_kin.w_st_c(1,exp_range4),mean(C_My_wx),3);
    C_Mz_wx_p = polyfit(body_kin.w_st_c(1,exp_range4),mean(C_Mz_wx),3);
   
    C_Mx_wy_p = polyfit(body_kin.w_st_c(2,exp_range5),mean(C_Mx_wy),3);
    C_My_wy_p = polyfit(body_kin.w_st_c(2,exp_range5),mean(C_My_wy),3);
    C_Mz_wy_p = polyfit(body_kin.w_st_c(2,exp_range5),mean(C_Mz_wy),3);
    
    C_Mx_wz_p = polyfit(body_kin.w_st_c(3,exp_range6),mean(C_Mx_wz),3);
    C_My_wz_p = polyfit(body_kin.w_st_c(3,exp_range6),mean(C_My_wz),3);
    C_Mz_wz_p = polyfit(body_kin.w_st_c(3,exp_range6),mean(C_Mz_wz),3);
    
    C_Fx_ax_p = polyfit(body_kin.a_st_c(1,exp_range7),mean(C_Fx_ax),1);
    C_Fy_ax_p = polyfit(body_kin.a_st_c(1,exp_range7),mean(C_Fy_ax),1);
    C_Fz_ax_p = polyfit(body_kin.a_st_c(1,exp_range7),mean(C_Fz_ax),1);
    
    C_Fx_ay_p = polyfit(body_kin.a_st_c(2,exp_range8),mean(C_Fx_ay),1);
    C_Fy_ay_p = polyfit(body_kin.a_st_c(2,exp_range8),mean(C_Fy_ay),1);
    C_Fz_ay_p = polyfit(body_kin.a_st_c(2,exp_range8),mean(C_Fz_ay),1);
    
    C_Fx_az_p = polyfit(body_kin.a_st_c(3,exp_range9),mean(C_Fx_az),1);
    C_Fy_az_p = polyfit(body_kin.a_st_c(3,exp_range9),mean(C_Fy_az),1);
    C_Fz_az_p = polyfit(body_kin.a_st_c(3,exp_range9),mean(C_Fz_az),1);
    
    C_Mx_w_dot_x_p = polyfit(body_kin.w_dot_st_c(1,exp_range10),mean(C_Mx_w_dot_x),1);
    C_My_w_dot_x_p = polyfit(body_kin.w_dot_st_c(1,exp_range10),mean(C_My_w_dot_x),1);
    C_Mz_w_dot_x_p = polyfit(body_kin.w_dot_st_c(1,exp_range10),mean(C_Mz_w_dot_x),1);
   
    C_Mx_w_dot_y_p = polyfit(body_kin.w_dot_st_c(2,exp_range11),mean(C_Mx_w_dot_y),1);
    C_My_w_dot_y_p = polyfit(body_kin.w_dot_st_c(2,exp_range11),mean(C_My_w_dot_y),1);
    C_Mz_w_dot_y_p = polyfit(body_kin.w_dot_st_c(2,exp_range11),mean(C_Mz_w_dot_y),1);
    
    C_Mx_w_dot_z_p = polyfit(body_kin.w_dot_st_c(3,exp_range12),mean(C_Mx_w_dot_z),1);
    C_My_w_dot_z_p = polyfit(body_kin.w_dot_st_c(3,exp_range12),mean(C_My_w_dot_z),1);
    C_Mz_w_dot_z_p = polyfit(body_kin.w_dot_st_c(3,exp_range12),mean(C_Mz_w_dot_z),1);    
    
    'C_Fx_vx_p'
    
    C_Fx_vx_p(1)
    C_Fx_vx_p(2)
    C_Fx_vx_p(3)
    C_Fx_vx_p(4)
    
    'C_Fy_vx_p'
    
    C_Fy_vx_p(1)
    C_Fy_vx_p(2)
    C_Fy_vx_p(3)
    C_Fy_vx_p(4)
    
    'C_Fz_vx_p'
   
    C_Fz_vx_p(1)
    C_Fz_vx_p(2)
    C_Fz_vx_p(3)
    C_Fz_vx_p(4)
    
    
    
    
    'C_Fx_vy_p'
    
    C_Fx_vy_p(1)
    C_Fx_vy_p(2)
    C_Fx_vy_p(3)
    C_Fx_vy_p(4)
    
    'C_Fy_vy_p'
    
    C_Fy_vy_p(1)
    C_Fy_vy_p(2)
    C_Fy_vy_p(3)
    C_Fy_vy_p(4)
    
    'C_Fz_vy_p'
   
    C_Fz_vy_p(1)
    C_Fz_vy_p(2)
    C_Fz_vy_p(3)
    C_Fz_vy_p(4)
    
    
    'C_Fx_vz_p'
    
    C_Fx_vz_p(1)
    C_Fx_vz_p(2)
    C_Fx_vz_p(3)
    C_Fx_vz_p(4)
    
    'C_Fy_vz_p'
    
    C_Fy_vz_p(1)
    C_Fy_vz_p(2)
    C_Fy_vz_p(3)
    C_Fy_vz_p(4)
    
    'C_Fz_vz_p'
   
    C_Fz_vz_p(1)
    C_Fz_vz_p(2)
    C_Fz_vz_p(3)
    C_Fz_vz_p(4)
    

    
    
    
    'C_Mx_wx_p'
    
    C_Mx_wx_p(1)
    C_Mx_wx_p(2)
    C_Mx_wx_p(3)
    C_Mx_wx_p(4)
    
    'C_My_wx_p'
    
    C_My_wx_p(1)
    C_My_wx_p(2)
    C_My_wx_p(3)
    C_My_wx_p(4)
    
    'C_Mz_wx_p'
   
    C_Mz_wx_p(1)
    C_Mz_wx_p(2)
    C_Mz_wx_p(3)
    C_Mz_wx_p(4)
    
    
    
    
    'C_Mx_wy_p'
    
    C_Mx_wy_p(1)
    C_Mx_wy_p(2)
    C_Mx_wy_p(3)
    C_Mx_wy_p(4)
    
    'C_My_wy_p'
    
    C_My_wy_p(1)
    C_My_wy_p(2)
    C_My_wy_p(3)
    C_My_wy_p(4)
    
    'C_Mz_wy_p'
   
    C_Mz_wy_p(1)
    C_Mz_wy_p(2)
    C_Mz_wy_p(3)
    C_Mz_wy_p(4)
    
    
    'C_Mx_wz_p'
    
    C_Mx_wz_p(1)
    C_Mx_wz_p(2)
    C_Mx_wz_p(3)
    C_Mx_wz_p(4)
    
    'C_My_wz_p'
    
    C_My_wz_p(1)
    C_My_wz_p(2)
    C_My_wz_p(3)
    C_My_wz_p(4)
    
    'C_Mz_wz_p'
   
    C_Mz_wz_p(1)
    C_Mz_wz_p(2)
    C_Mz_wz_p(3)
    C_Mz_wz_p(4)
    
    
    
    
    'C_Fx_ax_p'
    
    C_Fx_ax_p(1)
    C_Fx_ax_p(2)
    
    'C_Fy_ax_p'
    
    C_Fy_ax_p(1)
    C_Fy_ax_p(2)
    
    'C_Fz_ax_p'
    
    C_Fz_ax_p(1)
    C_Fz_ax_p(2)
    
    
    
    'C_Fx_ay_p'
    
    C_Fx_ay_p(1)
    C_Fx_ay_p(2)
    
    'C_Fy_ay_p'
    
    C_Fy_ay_p(1)
    C_Fy_ay_p(2)
    
    'C_Fz_ay_p'
    
    C_Fz_ay_p(1)
    C_Fz_ay_p(2)
    
    
    
    'C_Fx_az_p'
    
    C_Fx_az_p(1)
    C_Fx_az_p(2)
    
    'C_Fy_az_p'
    
    C_Fy_az_p(1)
    C_Fy_az_p(2)
    
    'C_Fz_az_p'
    
    C_Fz_az_p(1)
    C_Fz_az_p(2)
    
    
    
    
    
    
    'C_Mx_w_dot_x_p'
    
    C_Mx_w_dot_x_p(1)
    C_Mx_w_dot_x_p(2)
    
    'C_My_w_dot_x_p'
    
    C_My_w_dot_x_p(1)
    C_My_w_dot_x_p(2)
    
    'C_Mz_w_dot_x_p'
    
    C_Mz_w_dot_x_p(1)
    C_Mz_w_dot_x_p(2)
    
    
    
    'C_Mx_w_dot_y_p'
    
    C_Mx_w_dot_y_p(1)
    C_Mx_w_dot_y_p(2)
    
    'C_My_w_dot_y_p'
    
    C_My_w_dot_y_p(1)
    C_My_w_dot_y_p(2)
    
    'C_Mz_w_dot_y_p'
    
    C_Mz_w_dot_y_p(1)
    C_Mz_w_dot_y_p(2)
    
    
    
    'C_Mx_w_dot_z_p'
    
    C_Mx_w_dot_z_p(1)
    C_Mx_w_dot_z_p(2)
    
    'C_My_w_dot_z_p'
    
    C_My_w_dot_z_p(1)
    C_My_w_dot_z_p(2)
    
    'C_Mz_w_dot_z_p'
    
    C_Mz_w_dot_z_p(1)
    C_Mz_w_dot_z_p(2)

    
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.v_st_c(1,exp_range1),C_Fx_vx(j,:),'b')
    plot(body_kin.v_st_c(1,exp_range1),C_Fy_vx(j,:),'r')
    plot(body_kin.v_st_c(1,exp_range1),C_Fz_vx(j,:),'g')
    end
    plot(body_kin.v_st_c(1,exp_range1),polyval(C_Fx_vx_p,body_kin.v_st_c(1,exp_range1)),'k')
    plot(body_kin.v_st_c(1,exp_range1),polyval(C_Fy_vx_p,body_kin.v_st_c(1,exp_range1)),'k')
    plot(body_kin.v_st_c(1,exp_range1),polyval(C_Fz_vx_p,body_kin.v_st_c(1,exp_range1)),'k')
    hold off
    title('Aerodynamic damping forces for v_x')
    xlabel('v_x [mm/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.v_st_c(2,exp_range2),C_Fx_vy(j,:),'b')
    plot(body_kin.v_st_c(2,exp_range2),C_Fy_vy(j,:),'r')
    plot(body_kin.v_st_c(2,exp_range2),C_Fz_vy(j,:),'g')
    end
    plot(body_kin.v_st_c(2,exp_range2),polyval(C_Fx_vy_p,body_kin.v_st_c(2,exp_range2)),'k')
    plot(body_kin.v_st_c(2,exp_range2),polyval(C_Fy_vy_p,body_kin.v_st_c(2,exp_range2)),'k')
    plot(body_kin.v_st_c(2,exp_range2),polyval(C_Fz_vy_p,body_kin.v_st_c(2,exp_range2)),'k')
    hold off
    title('Aerodynamic damping forces for v_y')
    xlabel('v_y [mm/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.v_st_c(3,exp_range3),C_Fx_vz(j,:),'b')
    plot(body_kin.v_st_c(3,exp_range3),C_Fy_vz(j,:),'r')
    plot(body_kin.v_st_c(3,exp_range3),C_Fz_vz(j,:),'g')
    end
    plot(body_kin.v_st_c(3,exp_range3),polyval(C_Fx_vz_p,body_kin.v_st_c(3,exp_range3)),'k')
    plot(body_kin.v_st_c(3,exp_range3),polyval(C_Fy_vz_p,body_kin.v_st_c(3,exp_range3)),'k')
    plot(body_kin.v_st_c(3,exp_range3),polyval(C_Fz_vz_p,body_kin.v_st_c(3,exp_range3)),'k')
    hold off
    title('Aerodynamic damping forces for v_z')
    xlabel('v_z [mm/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_st_c(1,exp_range4),C_Mx_wx(j,:),'b')
    plot(body_kin.w_st_c(1,exp_range4),C_My_wx(j,:),'r')
    plot(body_kin.w_st_c(1,exp_range4),C_Mz_wx(j,:),'g')
    end
    plot(body_kin.w_st_c(1,exp_range4),polyval(C_Mx_wx_p,body_kin.w_st_c(1,exp_range4)),'k')
    plot(body_kin.w_st_c(1,exp_range4),polyval(C_My_wx_p,body_kin.w_st_c(1,exp_range4)),'k')
    plot(body_kin.w_st_c(1,exp_range4),polyval(C_Mz_wx_p,body_kin.w_st_c(1,exp_range4)),'k')
    hold off
    title('Aerodynamic damping moment for w_x')
    xlabel('w_x [rad/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_st_c(2,exp_range5),C_Mx_wy(j,:),'b')
    plot(body_kin.w_st_c(2,exp_range5),C_My_wy(j,:),'r')
    plot(body_kin.w_st_c(2,exp_range5),C_Mz_wy(j,:),'g')
    end
    plot(body_kin.w_st_c(2,exp_range5),polyval(C_Mx_wy_p,body_kin.w_st_c(2,exp_range5)),'k')
    plot(body_kin.w_st_c(2,exp_range5),polyval(C_My_wy_p,body_kin.w_st_c(2,exp_range5)),'k')
    plot(body_kin.w_st_c(2,exp_range5),polyval(C_Mz_wy_p,body_kin.w_st_c(2,exp_range5)),'k')
    hold off
    title('Aerodynamic damping moment for w_y')
    xlabel('w_y [rad/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_st_c(3,exp_range6),C_Mx_wz(j,:),'b')
    plot(body_kin.w_st_c(3,exp_range6),C_My_wz(j,:),'r')
    plot(body_kin.w_st_c(3,exp_range6),C_Mz_wz(j,:),'g')
    end
    plot(body_kin.w_st_c(3,exp_range6),polyval(C_Mx_wz_p,body_kin.w_st_c(3,exp_range6)),'k')
    plot(body_kin.w_st_c(3,exp_range6),polyval(C_My_wz_p,body_kin.w_st_c(3,exp_range6)),'k')
    plot(body_kin.w_st_c(3,exp_range6),polyval(C_Mz_wz_p,body_kin.w_st_c(3,exp_range6)),'k')
    hold off
    title('Aerodynamic damping moment for w_z')
    xlabel('w_z [rad/s]')
    ylabel('C_{ad}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.a_st_c(1,exp_range7),C_Fx_ax(j,:),'b')
    plot(body_kin.a_st_c(1,exp_range7),C_Fy_ax(j,:),'r')
    plot(body_kin.a_st_c(1,exp_range7),C_Fz_ax(j,:),'g')
    end
    plot(body_kin.a_st_c(1,exp_range7),polyval(C_Fx_ax_p,body_kin.a_st_c(1,exp_range7)),'k')
    plot(body_kin.a_st_c(1,exp_range7),polyval(C_Fy_ax_p,body_kin.a_st_c(1,exp_range7)),'k')
    plot(body_kin.a_st_c(1,exp_range7),polyval(C_Fz_ax_p,body_kin.a_st_c(1,exp_range7)),'k')
    hold off
    title('Inertial forces for a_x')
    xlabel('a_x [mm/s^2]')
    ylabel('C_{I}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.a_st_c(2,exp_range8),C_Fx_ay(j,:),'b')
    plot(body_kin.a_st_c(2,exp_range8),C_Fy_ay(j,:),'r')
    plot(body_kin.a_st_c(2,exp_range8),C_Fz_ay(j,:),'g')
    end
    plot(body_kin.a_st_c(2,exp_range8),polyval(C_Fx_ay_p,body_kin.a_st_c(2,exp_range8)),'k')
    plot(body_kin.a_st_c(2,exp_range8),polyval(C_Fy_ay_p,body_kin.a_st_c(2,exp_range8)),'k')
    plot(body_kin.a_st_c(2,exp_range8),polyval(C_Fz_ay_p,body_kin.a_st_c(2,exp_range8)),'k')
    hold off
    title('Inertial forces for a_y')
    xlabel('a_y [mm/s^2]')
    ylabel('C_{I}')    
    
   
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.a_st_c(3,exp_range9),C_Fx_az(j,:),'b')
    plot(body_kin.a_st_c(3,exp_range9),C_Fy_az(j,:),'r')
    plot(body_kin.a_st_c(3,exp_range9),C_Fz_az(j,:),'g')
    end
    plot(body_kin.a_st_c(3,exp_range9),polyval(C_Fx_az_p,body_kin.a_st_c(3,exp_range9)),'k')
    plot(body_kin.a_st_c(3,exp_range9),polyval(C_Fy_az_p,body_kin.a_st_c(3,exp_range9)),'k')
    plot(body_kin.a_st_c(3,exp_range9),polyval(C_Fz_az_p,body_kin.a_st_c(3,exp_range9)),'k')
    hold off
    title('Inertial forces for a_z')
    xlabel('a_z [mm/s^2]')
    ylabel('C_{I}') 
    
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_dot_st_c(1,exp_range10),C_Mx_w_dot_x(j,:),'b')
    plot(body_kin.w_dot_st_c(1,exp_range10),C_My_w_dot_x(j,:),'r')
    plot(body_kin.w_dot_st_c(1,exp_range10),C_Mz_w_dot_x(j,:),'g')
    end
    plot(body_kin.w_dot_st_c(1,exp_range10),polyval(C_Mx_w_dot_x_p,body_kin.w_dot_st_c(1,exp_range10)),'k')
    plot(body_kin.w_dot_st_c(1,exp_range10),polyval(C_My_w_dot_x_p,body_kin.w_dot_st_c(1,exp_range10)),'k')
    plot(body_kin.w_dot_st_c(1,exp_range10),polyval(C_Mz_w_dot_x_p,body_kin.w_dot_st_c(1,exp_range10)),'k')
    hold off
    title('Inertial moment for w dot x')
    xlabel('w dot x [rad/s^2]')
    ylabel('C_{I}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_dot_st_c(2,exp_range11),C_Mx_w_dot_y(j,:),'b')
    plot(body_kin.w_dot_st_c(2,exp_range11),C_My_w_dot_y(j,:),'r')
    plot(body_kin.w_dot_st_c(2,exp_range11),C_Mz_w_dot_y(j,:),'g')
    end
    plot(body_kin.w_dot_st_c(2,exp_range11),polyval(C_Mx_w_dot_y_p,body_kin.w_dot_st_c(2,exp_range11)),'k')
    plot(body_kin.w_dot_st_c(2,exp_range11),polyval(C_My_w_dot_y_p,body_kin.w_dot_st_c(2,exp_range11)),'k')
    plot(body_kin.w_dot_st_c(2,exp_range11),polyval(C_Mz_w_dot_y_p,body_kin.w_dot_st_c(2,exp_range11)),'k')
    hold off
    title('Inertial moment for w dot y')
    xlabel('w dot y [rad/s^2]')
    ylabel('C_{I}')
    
    figure()
    hold on
    for j = 1:nr_exp
    plot(body_kin.w_dot_st_c(3,exp_range12),C_Mx_w_dot_z(j,:),'b')
    plot(body_kin.w_dot_st_c(3,exp_range12),C_My_w_dot_z(j,:),'r')
    plot(body_kin.w_dot_st_c(3,exp_range12),C_Mz_w_dot_z(j,:),'g')
    end
    plot(body_kin.w_dot_st_c(3,exp_range12),polyval(C_Mx_w_dot_z_p,body_kin.w_dot_st_c(3,exp_range12)),'k')
    plot(body_kin.w_dot_st_c(3,exp_range12),polyval(C_My_w_dot_z_p,body_kin.w_dot_st_c(3,exp_range12)),'k')
    plot(body_kin.w_dot_st_c(3,exp_range12),polyval(C_Mz_w_dot_z_p,body_kin.w_dot_st_c(3,exp_range12)),'k')
    hold off
    title('Inertial moment for w dot z')
    xlabel('w dot z [rad/s^2]')
    ylabel('C_{I}')

%     figure()
%     hold on
%     subplot(3,1,1); plot(body_kin.v_st_c(2,exp_range2),
%     subplot(3,1,2); plot(body_kin.v_st_c(2,exp_range2),
%     subplot(3,1,3); plot(body_kin.v_st_c(2,exp_range2),
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(body_kin.v_st_c(3,exp_range3),
%     subplot(3,1,2); plot(body_kin.v_st_c(3,exp_range3),
%     subplot(3,1,3); plot(body_kin.v_st_c(3,exp_range3),
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(body_kin.v_st_c(4,exp_range4),
%     subplot(3,1,2); plot(body_kin.v_st_c(4,exp_range4),
%     subplot(3,1,3); plot(body_kin.v_st_c(4,exp_range4),
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(body_kin.v_st_c(5,exp_range5),
%     subplot(3,1,2); plot(body_kin.v_st_c(5,exp_range5),
%     subplot(3,1,3); plot(body_kin.v_st_c(5,exp_range5),
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(body_kin.v_st_c(6,exp_range6),
%     subplot(3,1,2); plot(body_kin.v_st_c(6,exp_range6),
%     subplot(3,1,3); plot(body_kin.v_st_c(6,exp_range6),
%     hold off

    
%     aero_damp = {};
%         
%     all_maneuvers = pathDB.all_maneuvers;
%     
%     roll_maneuvers = all_maneuvers.roll_maneuvers;
%     pitch_maneuvers = all_maneuvers.pitch_maneuvers;
%     yaw_maneuvers = all_maneuvers.yaw_maneuvers;
%     ax_maneuvers = all_maneuvers.ax_maneuvers;
%     ay_maneuvers = all_maneuvers.ay_maneuvers;
%     az_maneuvers = all_maneuvers.az_maneuvers;
%     pure_roll = all_maneuvers.pure_roll;
%     pure_pitch = all_maneuvers.pure_pitch;
%     pure_yaw = all_maneuvers.pure_yaw;
%     pure_ax = all_maneuvers.pure_ax;
%     pure_ay = all_maneuvers.pure_ay;
%     pure_az = all_maneuvers.pure_az;
%     
%     roll_man_name = fieldnames(roll_maneuvers);
%     pitch_man_name = fieldnames(pitch_maneuvers);
%     yaw_man_name = fieldnames(yaw_maneuvers);
%     ax_man_name = fieldnames(ax_maneuvers);
%     ay_man_name = fieldnames(ay_maneuvers);
%     az_man_name = fieldnames(az_maneuvers);
%     pure_roll_name = fieldnames(pure_roll);
%     pure_pitch_name = fieldnames(pure_pitch);
%     pure_yaw_name = fieldnames(pure_yaw);
%     pure_ax_name = fieldnames(pure_ax);
%     pure_ay_name = fieldnames(pure_ay);
%     pure_az_name = fieldnames(pure_az);
%     
%     nr_roll_man = size(roll_man_name,1);
%     nr_pitch_man = size(pitch_man_name,1);
%     nr_yaw_man = size(yaw_man_name,1);
%     nr_ax_man = size(ax_man_name,1);
%     nr_ay_man = size(ay_man_name,1);
%     nr_az_man = size(az_man_name,1);
%     nr_pure_roll = size(pure_roll_name,1);
%     nr_pure_pitch = size(pure_pitch_name,1);
%     nr_pure_yaw = size(pure_yaw_name,1);
%     nr_pure_ax = size(pure_ax_name,1);
%     nr_pure_ay = size(pure_ay_name,1);
%     nr_pure_az = size(pure_az_name,1);
%     
%     for i = 1:nr_roll_man
%         
%        temp_man = roll_maneuvers.(char(roll_man_name(i)));
%        
%        seq_nr = temp_man.seq_nr;
%         
% %        temp_theta_LR1_avg = temp_man.a_avg.theta_LR1;
% %        temp_theta_LR2_avg = temp_man.a_avg.theta_LR2;
% %        temp_eta_LR1_avg = temp_man.a_avg.eta_LR1;
% %        temp_eta_LR2_avg = temp_man.a_avg.eta_LR2;
% %        temp_phi_LR1_avg = temp_man.a_avg.phi_LR1;
% %        temp_phi_LR2_avg = temp_man.a_avg.phi_LR2;
% %        
% %        temp_theta_L1_dev = temp_man.a_dev.theta_L1;
% %        temp_theta_L2_dev = temp_man.a_dev.theta_L2;
% %        temp_eta_L1_dev = temp_man.a_dev.eta_L1;
% %        temp_eta_L2_dev = temp_man.a_dev.eta_L2;
% %        temp_phi_L1_dev = temp_man.a_dev.phi_L1;
% %        temp_phi_L2_dev = temp_man.a_dev.phi_L2;
% %         
% %        temp_theta_R1_dev = temp_man.a_dev.theta_R1;
% %        temp_theta_R2_dev = temp_man.a_dev.theta_R2;
% %        temp_eta_R1_dev = temp_man.a_dev.eta_R1;
% %        temp_eta_R2_dev = temp_man.a_dev.eta_R2;
% %        temp_phi_R1_dev = temp_man.a_dev.phi_R1;
% %        temp_phi_R2_dev = temp_man.a_dev.phi_R2;
%        
%        f_vect = temp_man.f;
%        
%        down_up = temp_man.down_up;
%        
%        nr_wb_man = length(f_vect);
%        
% %        for j = 1:nr_wb_man
% %            
% %            a_theta_L = [ temp_theta_LR1_avg+temp_theta_L1_dev(:,j); temp_theta_LR2_avg+temp_theta_L2_dev(:,j)];
% %            a_eta_L = [ temp_eta_LR1_avg+temp_eta_L1_dev(:,j); temp_eta_LR2_avg+temp_eta_L2_dev(:,j)];
% %            a_phi_L = [ temp_phi_LR1_avg+temp_phi_L1_dev(:,j); temp_phi_LR2_avg+temp_phi_L2_dev(:,j)];
% % 
% %            a_theta_R = [ temp_theta_LR1_avg+temp_theta_R1_dev(:,j); temp_theta_LR2_avg+temp_theta_R2_dev(:,j)];
% %            a_eta_R = [ temp_eta_LR1_avg+temp_eta_R1_dev(:,j); temp_eta_LR2_avg+temp_eta_R2_dev(:,j)];
% %            a_phi_R = [ temp_phi_LR1_avg+temp_phi_R1_dev(:,j); temp_phi_LR2_avg+temp_phi_R2_dev(:,j)];
% %            
% %            a_sim.theta_L(:,j) = a_theta_L;
% %            a_sim.eta_L(:,j) = a_eta_L;
% %            a_sim.phi_L(:,j) = a_phi_L;
% %            
% %            a_sim.theta_R(:,j) = a_theta_R;
% %            a_sim.eta_R(:,j) = a_eta_R;
% %            a_sim.phi_R(:,j) = a_phi_R;
% %            
% %        
% %        end
% %        
% %        t_vect = zeros(nr_wb_man+1,1);
% %           
% %        nr_points_vect = zeros(nr_wb_man,1);
% %     
% %        for k = 1:nr_wb_man
% % 
% %             nr_of_points = round((1/f_vect(k))/dt2)+1;
% %             nr_points_vect(k) = nr_of_points;
% %             t_vect(k+1) = t_vect(k)+(nr_of_points-1)*dt2;
% %             % t_vect(k+1) = k;
% %        end
% %             
% %        n_tot = sum(nr_points_vect)-nr_wb_man+1;
% %             
% %        phi_b = pi;
% %        theta_b = -(55/180)*pi;
% %        
% %        q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
% %                   cos(phi_b/2)*sin(theta_b/2); ...
% %                  -sin(phi_b/2)*sin(theta_b/2); ...
% %                   cos(phi_b/2)*cos(theta_b/2)];
% %               
% %        v_st_c = [0; 0; 0];
% %        a_st_c = [0; 0; 0];
% %        w_st_c = [0; 0; 0];
% %        w_dot_st_c = [0; 0; 0];
% %            
% %        q_b = zeros(4,n_tot);
% %        v_st = zeros(3,n_tot);
% %        a_st = zeros(3,n_tot);
% %        w_st = zeros(3,n_tot);
% %        w_dot_st = zeros(3,n_tot);
% %            
% %        for k = 1:n_tot
% %                
% %          q_b(:,k) = q_b_c;
% %          v_st(:,k) = v_st_c;
% %          a_st(:,k) = a_st_c;
% %          w_st(:,k) = w_st_c;
% %          w_dot_st(:,k) = w_dot_st_c;
% %                
% %        end
%        
%        body_kin.phi_b = 
%        body_kin.theta_b = 
%        body_kin.v_st_c = 
%        body_kin.a_st_c = 
%        body_kin.w_st_c = 
%        body_kin.w_dot_st_c = 
%        
%        
%        [ experiment ] = Aerodynamic_damping_maneuvers( nr_man, maneuvers, man_name, body_kin )
% 
% 
%        Req_f = Required_forces( settings,pathDB, a_sim, nr_points_vect, t_vect, down_up, f_vect, dt2, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
%        
%        Req_f
%            
%        t_man = Req_f.time;
%        
%        F_I_acc  = Req_f.F_I_acc;
%        F_I_vel  = Req_f.F_I_vel;
%        M_I_acc  = Req_f.M_I_acc;
%        M_I_vel  = Req_f.M_I_vel;
%        F_I      = Req_f.F_I;
%        M_I      = Req_f.M_I;
%        c_pres_L = Req_f.c_pres_L;
%        c_pres_R = Req_f.c_pres_R;
%        F_a_wL   = Req_f.F_a_wL;
%        F_a_wR   = Req_f.F_a_wR;
%        F_a_L    = Req_f.F_a_L;
%        F_a_R    = Req_f.F_a_R;
%        M_a_L    = Req_f.M_a_L;
%        M_a_R    = Req_f.M_a_R;
%        F_g      = Req_f.F_g;
% 
%        
%        arm_L = [0; -1; 0];
%        arm_R = [0; 1; 0];
%        
%        arm_L1 = [0.025; -1; 0];
%        arm_R1 = [0.025; 1; 0];
%        
%        arm_L2 = [-0.075; -1; 0];
%        arm_R2 = [-0.075; 1; 0];
%        
%        figure()
%        hold on
%        for k = 1:nr_points_vect(1)
%            Lwt = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L;
%            Rwt = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R;
%            Lwt1 = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L1;
%            Rwt1 = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R1;
%            Lwt2 = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L2;
%            Rwt2 = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R2;
%            scale = 10;
%            F_L = scale*Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*Req_f.F_a_wL(:,k);
%            F_R = scale*Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*Req_f.F_a_wR(:,k);
%            plot3([Lwt1(1) Lwt2(1)],[Lwt1(2) Lwt2(2)],[Lwt1(3) Lwt2(3)],'k')
%            plot3([Rwt1(1) Rwt2(1)],[Rwt1(2) Rwt2(2)],[Rwt1(3) Rwt2(3)],'k')
%            plot3(cos(0:0.05:(2*pi)),sin(0:0.05:(2*pi)),zeros(length(0:0.05:(2*pi)),1),':','Color','b')
%            plot3([-1 1],[0 0],[0 0],':','Color','b')
%            quiver3(Lwt(1),Lwt(2),Lwt(3),F_L(1),F_L(2),F_L(3),'r')
%            quiver3(Rwt(1),Rwt(2),Rwt(3),F_R(1),F_R(2),F_R(3),'r')
%        end
%        hold off
%        axis equal
%        xlabel('x')
%        ylabel('y')
%        zlabel('z')
% 
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(1,b),t_man(b),M_a_L(1,b),t_man(b),M_a_R(1,b),t_man(b),M_a_L(1,b)+M_a_R(1,b)+M_I(1,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(1,b)),t_man(b),ones(length(b),1)*mean(M_a_L(1,b)+M_a_R(1,b)),t_man(b),ones(length(b),1)*mean(M_a_L(1,b)+M_a_R(1,b)+M_I(1,b)))
%        end
%        hold off
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(2,b),t_man(b),M_a_L(2,b),t_man(b),M_a_R(2,b),t_man(b),M_a_L(2,b)+M_a_R(2,b)+M_I(2,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(2,b)),t_man(b),ones(length(b),1)*mean(M_a_L(2,b)+M_a_R(2,b)),t_man(b),ones(length(b),1)*mean(M_a_L(2,b)+M_a_R(2,b)+M_I(2,b)))
%        end
%        hold off
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(3,b),t_man(b),M_a_L(3,b),t_man(b),M_a_R(3,b),t_man(b),M_a_L(3,b)+M_a_R(3,b)+M_I(3,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(3,b)),t_man(b),ones(length(b),1)*mean(M_a_L(3,b)+M_a_R(3,b)),t_man(b),ones(length(b),1)*mean(M_a_L(3,b)+M_a_R(3,b)+M_I(3,b)))
%        end
%        hold off
%        
%        pause
%        
%        
%     end       
    
    


end

