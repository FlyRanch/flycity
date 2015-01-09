function [temp] = Wing_ref_frame_saveWBkin(settings,pathDB)

    savefile = 'WBkin.mat';
    
    % nr_of_wing_sections:
    
    nr_wing_sec = 10;
        
    u_wing_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_trans_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_trans_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_trans_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_trans_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_trans_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_trans_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_body_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_body_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_body_rot_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    u_wing_body_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    v_wing_body_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    w_wing_body_rot_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    alfa_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    beta_L = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    alfa_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    beta_R = nan(size(pathDB.x,1), nr_wing_sec+2, size(pathDB.x,2));
    
    Lwingtip = nan(size(pathDB.x,1),3,size(pathDB.x,2));
    Rwingtip = nan(size(pathDB.x,1),3,size(pathDB.x,2));
    
    Lwingbeat_loc = nan(150,2,size(pathDB.x,2));
    Rwingbeat_loc = nan(150,2,size(pathDB.x,2));
    
    wingbeat_f = nan(size(pathDB.x));
    
    phi_L = nan(size(pathDB.x));
    theta_L = nan(size(pathDB.x));
    eta_L = nan(size(pathDB.x));
    
    phi_R = nan(size(pathDB.x));
    theta_R = nan(size(pathDB.x));
    eta_R = nan(size(pathDB.x));
    
    
    for i=1:size(pathDB.qb1_filt,2)
    
        % start and stop point for the measurements
        start = find(isnan(pathDB.qb1(:,i))==0, 1 );
        stop = find(isnan(pathDB.qb1(:,i))==0, 1, 'last' );
        
        wing_l = pathDB.wing_l(i);
        
        joint_loc_L = pathDB.joint_pos_L(:,i);

        joint_loc_R = pathDB.joint_pos_R(:,i);
      
        wing_loc_L = wing_l.*[zeros(1,nr_wing_sec+2); -1 -1+(0.5/nr_wing_sec):(1/nr_wing_sec):-(0.5/nr_wing_sec) 0; zeros(1,nr_wing_sec+2)];

        wing_loc_R = wing_l.*[zeros(1,nr_wing_sec+2); 1 1-(0.5/nr_wing_sec):-(1/nr_wing_sec):(0.5/nr_wing_sec) 0; zeros(1,nr_wing_sec+2)];

        for j = start:stop
            
            omega_b = [pathDB.b_omega1(j,i); pathDB.b_omega2(j,i); pathDB.b_omega3(j,i)];
            
            omega_L = [pathDB.omega1_L(j,i); pathDB.omega2_L(j,i); pathDB.omega3_L(j,i)];
            
            omega_R = [pathDB.omega1_R(j,i); pathDB.omega2_R(j,i); pathDB.omega3_R(j,i)];
            
            qL = [pathDB.qL1_filt2(j,i); pathDB.qL2_filt2(j,i); pathDB.qL3_filt2(j,i); pathDB.qL4_filt2(j,i)];
            
            qR = [pathDB.qR1_filt2(j,i); pathDB.qR2_filt2(j,i); pathDB.qR3_filt2(j,i); pathDB.qR4_filt2(j,i)];
            
            U_body = [pathDB.u_body(j,i); pathDB.v_body(j,i); pathDB.w_body(j,i)];
            
            % Calculate the wing velocity at different locations on the wing:
            
            [ U_wing_L, U_wing_R , U_wing_rot_L, U_wing_rot_R, U_trans_L, U_trans_R, U_body_rot_L, U_body_rot_R] = Wing_vel(wing_loc_L, wing_loc_R, omega_b, omega_L, omega_R, joint_loc_L, joint_loc_R, qL, qR, U_body);
            
            u_wing_L(j,:,i) = U_wing_L(1,:);
            v_wing_L(j,:,i) = U_wing_L(2,:);
            w_wing_L(j,:,i) = U_wing_L(3,:);
            
            u_wing_R(j,:,i) = U_wing_R(1,:);
            v_wing_R(j,:,i) = U_wing_R(2,:);
            w_wing_R(j,:,i) = U_wing_R(3,:);  
            
            u_wing_rot_L(j,:,i) = U_wing_rot_L(1,:);
            v_wing_rot_L(j,:,i) = U_wing_rot_L(2,:);
            w_wing_rot_L(j,:,i) = U_wing_rot_L(3,:);

            u_wing_rot_R(j,:,i) = U_wing_rot_R(1,:);
            v_wing_rot_R(j,:,i) = U_wing_rot_R(2,:);
            w_wing_rot_R(j,:,i) = U_wing_rot_R(3,:);

            u_wing_trans_L(j,:,i) = U_trans_L(1,:);
            v_wing_trans_L(j,:,i) = U_trans_L(2,:);
            w_wing_trans_L(j,:,i) = U_trans_L(3,:);

            u_wing_trans_R(j,:,i) = U_trans_R(1,:);
            v_wing_trans_R(j,:,i) = U_trans_R(2,:);
            w_wing_trans_R(j,:,i) = U_trans_R(3,:);

            u_wing_body_rot_L(j,:,i) = U_body_rot_L(1,:);
            v_wing_body_rot_L(j,:,i) = U_body_rot_L(2,:);
            w_wing_body_rot_L(j,:,i) = U_body_rot_L(3,:);

            u_wing_body_rot_R(j,:,i) = U_body_rot_R(1,:);
            v_wing_body_rot_R(j,:,i) = U_body_rot_R(2,:);
            w_wing_body_rot_R(j,:,i) = U_body_rot_R(3,:);
            
            clear omega_b omega_L omega_R qL qR U_body
            
            
            
            % Calculate the angle of attack and angle of side-slip at wing.
            
            [ Alfa_L, Beta_L, Alfa_R, Beta_R ] = wing_alfa_beta( u_wing_L(j,:,i), v_wing_L(j,:,i), w_wing_L(j,:,i), u_wing_R(j,:,i), v_wing_R(j,:,i), w_wing_R(j,:,i) );
            
            alfa_L(j,:,i) = Alfa_L;
            beta_L(j,:,i) = Beta_L;

            alfa_R(j,:,i) = Alfa_R;
            beta_R(j,:,i) = Beta_R;    

        
        end
        
        
        % Calculate the duration of the upstroke and the downstroke:
            
            settings.sequence_names(i)
        
            [ L_wingtip_path, R_wingtip_path, L_wingbeat_loc, R_wingbeat_loc] = Wingbeat(pathDB.qL1_filt2(start:stop,i),pathDB.qL2_filt2(start:stop,i),pathDB.qL3_filt2(start:stop,i),pathDB.qL4_filt2(start:stop,i), ...
                                                                                pathDB.qR1_filt2(start:stop,i),pathDB.qR2_filt2(start:stop,i),pathDB.qR3_filt2(start:stop,i),pathDB.qR4_filt2(start:stop,i), wing_l);
                                                                            
            
        
            Lwingtip(start:stop,:,i) = L_wingtip_path;
            Rwingtip(start:stop,:,i) = R_wingtip_path;
            Lwingbeat_loc(1:size(L_wingbeat_loc,1),1,i) = L_wingbeat_loc(:,1);
            Lwingbeat_loc(1:size(L_wingbeat_loc,1),2,i) = L_wingbeat_loc(:,2);
            Rwingbeat_loc(1:size(R_wingbeat_loc,1),1,i) = R_wingbeat_loc(:,1);
            Rwingbeat_loc(1:size(R_wingbeat_loc,1),2,i) = R_wingbeat_loc(:,2);
            
            
         % Calculate phi, theta and eta of the strokeplane and calculculate
         % the wingbeatfrequency per set of down- and up-strokes. The
         % wingbeatfrequency will be analyzed starting from the first
         % downstroke:
         
         % Calculate the stroke plane angles:
            
         [ L_phi, L_theta, L_eta, R_phi, R_theta, R_eta ] = wing_stroke_plane( pathDB.qL1_filt2(start:stop,i), pathDB.qL2_filt2(start:stop,i), pathDB.qL3_filt2(start:stop,i), pathDB.qL4_filt2(start:stop,i), ...
                                                                                pathDB.qR1_filt2(start:stop,i), pathDB.qR2_filt2(start:stop,i), pathDB.qR3_filt2(start:stop,i), pathDB.qR4_filt2(start:stop,i),...
                                                                                settings);
         
        phi_L(start:stop,i) = L_phi;
        theta_L(start:stop,i) = L_theta;
        eta_L(start:stop,i) = L_eta;
    
        phi_R(start:stop,i) = R_phi;
        theta_R(start:stop,i) = R_theta;
        eta_R(start:stop,i) = R_eta;
         
        
    end
    
    temp.u_wing_L = u_wing_L;
    temp.v_wing_L = v_wing_L;
    temp.w_wing_L = w_wing_L;
    
    temp.u_wing_R = u_wing_R;
    temp.v_wing_R = v_wing_R;
    temp.w_wing_R = w_wing_R;
    
    temp.u_wing_rot_L = u_wing_rot_L;
    temp.v_wing_rot_L = v_wing_rot_L;
    temp.w_wing_rot_L = w_wing_rot_L;
    
    temp.u_wing_rot_R = u_wing_rot_R;
    temp.v_wing_rot_R = v_wing_rot_R;
    temp.w_wing_rot_R = w_wing_rot_R;
    
    temp.u_wing_trans_L = u_wing_trans_L;
    temp.v_wing_trans_L = v_wing_trans_L;
    temp.w_wing_trans_L = w_wing_trans_L;
    
    temp.u_wing_trans_R = u_wing_trans_R;
    temp.v_wing_trans_R = v_wing_trans_R;
    temp.w_wing_trans_R = w_wing_trans_R;
    
    temp.u_wing_body_rot_L = u_wing_body_rot_L;
    temp.v_wing_body_rot_L = v_wing_body_rot_L;
    temp.w_wing_body_rot_L = w_wing_body_rot_L;
    
    temp.u_wing_body_rot_R = u_wing_body_rot_R;
    temp.v_wing_body_rot_R = v_wing_body_rot_R;
    temp.w_wing_body_rot_R = w_wing_body_rot_R;
    
    temp.alfa_L = alfa_L;
    temp.beta_L = beta_L;
    
    temp.alfa_R = alfa_R;
    temp.beta_R = beta_R;
    
    temp.Lwingtip = Lwingtip;
    temp.Rwingtip = Rwingtip;

    temp.Lwingbeat_loc = Lwingbeat_loc;
    temp.Rwingbeat_loc = Rwingbeat_loc;
    
    temp.phi_L = phi_L;
    temp.theta_L = theta_L;
    temp.eta_L = eta_L;

    temp.phi_R = phi_R;
    temp.theta_R = theta_R;
    temp.eta_R = eta_R;

    save(savefile,'u_wing_L', 'v_wing_L', 'w_wing_L', 'u_wing_R', 'v_wing_R', 'w_wing_R', 'u_wing_rot_L', 'v_wing_rot_L', 'w_wing_rot_L', ...
        'u_wing_rot_R', 'v_wing_rot_R', 'w_wing_rot_R', 'u_wing_trans_L', 'v_wing_trans_L', 'w_wing_trans_L', 'u_wing_trans_R', ...
        'v_wing_trans_R', 'w_wing_trans_R', 'u_wing_body_rot_L', 'v_wing_body_rot_L', 'w_wing_body_rot_L', 'u_wing_body_rot_R', 'v_wing_body_rot_R', ...
        'w_wing_body_rot_R', 'alfa_L', 'beta_L', 'alfa_R', 'beta_R', 'Lwingtip', 'Rwingtip', 'Lwingbeat_loc', 'Rwingbeat_loc','phi_L','theta_L', ...
        'eta_L','phi_R','theta_R','eta_R')
    