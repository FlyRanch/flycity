function PLOT_generator2( settings, pathDB )

    % File that creates and saves plots. The plots are divided into three
    % groups, the first group consists of plots which evaluate the
    % smoothing process. The second group generates plots which display all
    % the wing kinematic variables. The third group generates all
    % trajectory related plots.
    
    addpath(char(settings.path_names(1)))
    alldirs=dir;
    
    addpath(char(settings.path_names(5)))
    
    for i=1:size(pathDB.x,2)
        
       close all;
       clearvars -except alldirs settings pathDB i

       start = find(isnan(pathDB.x(:,i))==0, 1 );
       stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
       
       % Plot path and add body reference frame after each 100 frames
   
       %Calculate the postions of the left and right wing hinges w.r.t. the bodyframe
       load_name = strcat(alldirs(i+2).name, '/flytracks/ManualFit_flytracks');
       load(load_name);

%         BL = ManualFit.params.bodyscale*(ManualFit.params.bodylen+ManualFit.params.headlen);
%         
%         wing_l = (abs(max(ManualFit.params.wing(:,1)))+abs(min(ManualFit.params.wing(:,1))))*ManualFit.params.wingscale;

        BL = pathDB.body_l(i);
        wing_l = pathDB.wing_l(i);
    
        % Read pathDB
    
        t = pathDB.t(start:stop,1);

        x_r = pathDB.x(start:stop,i);
        y_r = pathDB.y(start:stop,i);
        z_r = pathDB.z(start:stop,i);

        x_f = pathDB.x_filt(start:stop,i);
        y_f = pathDB.y_filt(start:stop,i);
        z_f = pathDB.z_filt(start:stop,i);

        u_f = pathDB.u_filt(start:stop,i);
        v_f = pathDB.v_filt(start:stop,i);
        w_f = pathDB.w_filt(start:stop,i);

        ax_f = pathDB.ax_filt(start:stop,i);
        ay_f = pathDB.ay_filt(start:stop,i);
        az_f = pathDB.az_filt(start:stop,i);

        qb1_r = pathDB.qb1(start:stop,i);
        qb2_r = pathDB.qb2(start:stop,i);
        qb3_r = pathDB.qb3(start:stop,i);
        qb4_r = pathDB.qb4(start:stop,i);

        qL1_r = pathDB.qL1(start:stop,i);
        qL2_r = pathDB.qL2(start:stop,i);
        qL3_r = pathDB.qL3(start:stop,i);
        qL4_r = pathDB.qL4(start:stop,i);

        qR1_r = pathDB.qR1(start:stop,i);
        qR2_r = pathDB.qR2(start:stop,i);
        qR3_r = pathDB.qR3(start:stop,i);
        qR4_r = pathDB.qR4(start:stop,i);

        qb1_f = pathDB.qb1_filt(start:stop,i);
        qb2_f = pathDB.qb2_filt(start:stop,i);
        qb3_f = pathDB.qb3_filt(start:stop,i);
        qb4_f = pathDB.qb4_filt(start:stop,i);

        wbx_f = pathDB.b_omega1(start:stop,i);
        wby_f = pathDB.b_omega2(start:stop,i);
        wbz_f = pathDB.b_omega3(start:stop,i);

        qL1_f1 = pathDB.qL1_filt1(start:stop,i);
        qL2_f1 = pathDB.qL2_filt1(start:stop,i);
        qL3_f1 = pathDB.qL3_filt1(start:stop,i);
        qL4_f1 = pathDB.qL4_filt1(start:stop,i);

        qR1_f1 = pathDB.qR1_filt1(start:stop,i);
        qR2_f1 = pathDB.qR2_filt1(start:stop,i);
        qR3_f1 = pathDB.qR3_filt1(start:stop,i);
        qR4_f1 = pathDB.qR4_filt1(start:stop,i);

        wLx_f2 = pathDB.omega1_L(start:stop,i);
        wLy_f2 = pathDB.omega2_L(start:stop,i);
        wLz_f2 = pathDB.omega3_L(start:stop,i);

        qL1_f2 = pathDB.qL1_filt2(start:stop,i);
        qL2_f2 = pathDB.qL2_filt2(start:stop,i);
        qL3_f2 = pathDB.qL3_filt2(start:stop,i);
        qL4_f2 = pathDB.qL4_filt2(start:stop,i);

        wRx_f2 = pathDB.omega1_R(start:stop,i);
        wRy_f2 = pathDB.omega2_R(start:stop,i);
        wRz_f2 = pathDB.omega3_R(start:stop,i);

        qR1_f2 = pathDB.qR1_filt2(start:stop,i);
        qR2_f2 = pathDB.qR2_filt2(start:stop,i);
        qR3_f2 = pathDB.qR3_filt2(start:stop,i);
        qR4_f2 = pathDB.qR4_filt2(start:stop,i);
        
        alfa_b = pathDB.b_alfa(start:stop,i);
        beta_b = pathDB.b_beta(start:stop,i);
        
        u_wing_L = pathDB.u_wing_L(start:stop,:,i);
        v_wing_L = pathDB.v_wing_L(start:stop,:,i);
        w_wing_L = pathDB.w_wing_L(start:stop,:,i);
        
        u_wing_R = pathDB.u_wing_R(start:stop,:,i);
        v_wing_R = pathDB.v_wing_R(start:stop,:,i);
        w_wing_R = pathDB.w_wing_R(start:stop,:,i);
        
        u_wing_rot_L(:,:,1) = pathDB.u_wing_rot_L(start:stop,:,i);
        u_wing_rot_L(:,:,2) = pathDB.v_wing_rot_L(start:stop,:,i);
        u_wing_rot_L(:,:,3) = pathDB.w_wing_rot_L(start:stop,:,i);
    
        u_wing_rot_R(:,:,1) = pathDB.u_wing_rot_R(start:stop,:,i);
        u_wing_rot_R(:,:,2) = pathDB.v_wing_rot_R(start:stop,:,i);
        u_wing_rot_R(:,:,3) = pathDB.w_wing_rot_R(start:stop,:,i);
    
        u_trans_L(:,:,1) = pathDB.u_wing_trans_L(start:stop,:,i);
        u_trans_L(:,:,2) = pathDB.v_wing_trans_L(start:stop,:,i);
        u_trans_L(:,:,3) = pathDB.w_wing_trans_L(start:stop,:,i);
    
        u_trans_R(:,:,1) = pathDB.u_wing_trans_R(start:stop,:,i);
        u_trans_R(:,:,2) = pathDB.v_wing_trans_R(start:stop,:,i);
        u_trans_R(:,:,3) = pathDB.w_wing_trans_R(start:stop,:,i);
    
        u_body_rot_L(:,:,1) = pathDB.u_wing_body_rot_L(start:stop,:,i);
        u_body_rot_L(:,:,2) = pathDB.v_wing_body_rot_L(start:stop,:,i);
        u_body_rot_L(:,:,3) = pathDB.w_wing_body_rot_L(start:stop,:,i);
    
        u_body_rot_R(:,:,1) = pathDB.u_wing_body_rot_R(start:stop,:,i);
        u_body_rot_R(:,:,2) = pathDB.v_wing_body_rot_R(start:stop,:,i);
        u_body_rot_R(:,:,3) = pathDB.w_wing_body_rot_R(start:stop,:,i);
        
        u_body = pathDB.u_body(start:stop,i);
        v_body = pathDB.v_body(start:stop,i);
        w_body = pathDB.w_body(start:stop,i);
        
        ax_body = pathDB.ax_body(start:stop,i);
        ay_body = pathDB.ay_body(start:stop,i);
        az_body = pathDB.az_body(start:stop,i);
        
        alfa_L = pathDB.alfa_L(start:stop,:,i);
        beta_L = pathDB.beta_L(start:stop,:,i);
        
        alfa_R = pathDB.alfa_R(start:stop,:,i);
        beta_R = pathDB.beta_R(start:stop,:,i);
        
        Lwingtip = pathDB.Lwingtip(start:stop,:,i);
        Rwingtip = pathDB.Rwingtip(start:stop,:,i);
        
        phi_L = pathDB.phi_L(start:stop,i);
        theta_L = pathDB.theta_L(start:stop,i);
        eta_L = pathDB.eta_L(start:stop,i);
        
        phi_R = pathDB.phi_R(start:stop,i);
        theta_R = pathDB.theta_R(start:stop,i);
        eta_R = pathDB.eta_R(start:stop,i);
        
        k = 1;
        
        while isnan(pathDB.L_wingbeat_loc(k,1,i)) == 0
          
            k = k+1;
            
        end
        
        Lwingbeat_loc = pathDB.L_wingbeat_loc(1:k-1,:,i);
        Rwingbeat_loc = pathDB.R_wingbeat_loc(1:k-1,:,i);
        
        % Determine whether plots will be saved or not [0 = not, 1 = save]
        %save_on_off = settings.save_on_off;
        
       % Movie plots:
        
       %Wing_kinematics_3D_stationary(settings,pathDB,i,save_on_off)

       %Wing_kinematics_3D_moving(settings,pathDB,i,save_on_off)
       
       %Trajectory_plus_orientation_movie(settings,pathDB,i,save_on_off)
              
       %strokeplane_force_movie(settings,pathDB,i,save_on_off)
        
       %arm_force_movie(settings,pathDB,i,save_on_off)

       %Force_3D_stationary(settings,pathDB,i,save_on_off)

       %Force_3D_moving(settings,pathDB,i,save_on_off)


        
        clear start stop
        
    end
    
    save_on_off = settings.save_on_off;
    
    Tracker_movie(settings, pathDB, 13)
    
%     Wing_kinematics_3D_moving(settings,pathDB,13,save_on_off)
%        
%     Trajectory_plus_orientation_movie(settings,pathDB,13,save_on_off)

%    strokeplane_force_movie(settings,pathDB,13,save_on_off)

end