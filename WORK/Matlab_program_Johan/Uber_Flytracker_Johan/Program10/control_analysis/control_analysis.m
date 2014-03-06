function [a_dev_man,a_dev_glob] = control_analysis( pathDB )


    n_pol_theta = length(pathDB.a_glob.theta1)-1;
    n_pol_eta = length(pathDB.a_glob.eta1)-1;
    n_pol_phi = length(pathDB.a_glob.phi1)-1;
    
    roll_maneuvers = pathDB.all_maneuvers.roll_maneuvers;
    pitch_maneuvers = pathDB.all_maneuvers.pitch_maneuvers;
    yaw_maneuvers = pathDB.all_maneuvers.yaw_maneuvers;
    ax_maneuvers = pathDB.all_maneuvers.ax_maneuvers;
    ay_maneuvers = pathDB.all_maneuvers.ay_maneuvers;
    az_maneuvers = pathDB.all_maneuvers.az_maneuvers;
    pure_roll = pathDB.all_maneuvers.pure_roll;
    pure_pitch = pathDB.all_maneuvers.pure_pitch;
    pure_yaw = pathDB.all_maneuvers.pure_yaw;
    pure_ax = pathDB.all_maneuvers.pure_ax;
    pure_ay = pathDB.all_maneuvers.pure_ay;
    pure_az = pathDB.all_maneuvers.pure_az;
    
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
    
    [ F_ax_theta, F_ax_eta, F_ax_phi] = FM_matrix2( ax_maneuvers, ax_man_name, nr_ax_man, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );
    [ F_ay_theta, F_ay_eta, F_ay_phi] = FM_matrix2( ay_maneuvers, ay_man_name, nr_ay_man, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );
    [ F_az_theta, F_az_eta, F_az_phi] = FM_matrix2( az_maneuvers, az_man_name, nr_az_man, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );

    [ F_pure_ax_theta, F_pure_ax_eta, F_pure_ax_phi] = FM_matrix2( pure_ax, pure_ax_name, nr_pure_ax, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );
    [ F_pure_ay_theta, F_pure_ay_eta, F_pure_ay_phi] = FM_matrix2( pure_ay, pure_ay_name, nr_pure_ay, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );
    [ F_pure_az_theta, F_pure_az_eta, F_pure_az_phi] = FM_matrix2( pure_az, pure_az_name, nr_pure_az, n_pol_theta, n_pol_eta, n_pol_phi, 1, 0 );
    
    [ M_roll_theta, M_roll_eta, M_roll_phi ]    = FM_matrix2( roll_maneuvers, roll_man_name, nr_roll_man, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );
    [ M_pitch_theta, M_pitch_eta, M_pitch_phi ] = FM_matrix2( pitch_maneuvers, pitch_man_name, nr_pitch_man, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );
    [ M_yaw_theta, M_yaw_eta, M_yaw_phi ]       = FM_matrix2( yaw_maneuvers, yaw_man_name, nr_yaw_man, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );
    
    [ M_pure_roll_theta, M_pure_roll_eta, M_pure_roll_phi ]     = FM_matrix2( pure_roll, pure_roll_name, nr_pure_roll, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );
    [ M_pure_pitch_theta, M_pure_pitch_eta, M_pure_pitch_phi ]  = FM_matrix2( pure_pitch, pure_pitch_name, nr_pure_pitch, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );
    [ M_pure_yaw_theta, M_pure_yaw_eta, M_pure_yaw_phi ]        = FM_matrix2( pure_yaw, pure_yaw_name, nr_pure_yaw, n_pol_theta, n_pol_eta, n_pol_phi, 0, 1 );

    
    %----------------------------------------------------------------------
    
    
    size(F_ax_theta)
    size(F_ay_theta)
    size(F_az_theta)
    size(F_pure_ax_theta)
    size(F_pure_ay_theta)
    size(F_pure_az_theta)
    
    size(M_roll_theta)
    size(M_pitch_theta)
    size(M_yaw_theta)
    size(M_pure_roll_theta)
    size(M_pure_pitch_theta)
    size(M_pure_yaw_theta)
    

    
    %----------------------------------------------------------------------
    

    % Non-dimensionalize the moments and forces:

    non_dim_ax_F_L = (F_ax_theta(:,1)-F_ax_theta(:,10))./(F_ax_theta(:,8).^2.*F_ax_theta(:,7).*F_ax_theta(:,18));
    non_dim_ay_F_L = (F_ay_theta(:,2)-F_ay_theta(:,11))./(F_ay_theta(:,8).^2.*F_ay_theta(:,7).*F_ay_theta(:,18));
    non_dim_az_F_L = (F_az_theta(:,3)-F_az_theta(:,12))./(F_az_theta(:,8).^2.*F_az_theta(:,7).*F_az_theta(:,18));
    non_dim_ax_F_R = (F_ax_theta(:,4)-F_ax_theta(:,13))./(F_ax_theta(:,8).^2.*F_ax_theta(:,7).*F_ax_theta(:,18));
    non_dim_ay_F_R = (F_ay_theta(:,5)-F_ay_theta(:,14))./(F_ay_theta(:,8).^2.*F_ay_theta(:,7).*F_ay_theta(:,18));
    non_dim_az_F_R = (F_az_theta(:,6)-F_az_theta(:,15))./(F_az_theta(:,8).^2.*F_az_theta(:,7).*F_az_theta(:,18));

    non_dim_roll_M_L  = (M_roll_theta(:,1)-M_roll_theta(:,10))./(M_roll_theta(:,8).^2.*M_roll_theta(:,7).^2.*M_roll_theta(:,18));
    non_dim_pitch_M_L = (M_pitch_theta(:,2)-M_pitch_theta(:,11))./(M_pitch_theta(:,8).^2.*M_pitch_theta(:,7).^2.*M_pitch_theta(:,18));
    non_dim_yaw_M_L   = (M_yaw_theta(:,3)-M_yaw_theta(:,12))./(M_yaw_theta(:,8).^2.*M_yaw_theta(:,7).^2.*M_yaw_theta(:,18));
    non_dim_roll_M_R  = (M_roll_theta(:,4)-M_roll_theta(:,13))./(M_roll_theta(:,8).^2.*M_roll_theta(:,7).^2.*M_roll_theta(:,18));
    non_dim_pitch_M_R = (M_pitch_theta(:,5)-M_pitch_theta(:,14))./(M_pitch_theta(:,8).^2.*M_pitch_theta(:,7).^2.*M_pitch_theta(:,18));
    non_dim_yaw_M_R   = (M_yaw_theta(:,6)-M_yaw_theta(:,15))./(M_yaw_theta(:,8).^2.*M_yaw_theta(:,7).^2.*M_yaw_theta(:,18));    
    
    non_dim_ax_pure_F_L = (F_pure_ax_theta(:,1)-F_pure_ax_theta(:,10))./(F_pure_ax_theta(:,8).^2.*F_pure_ax_theta(:,7).*F_pure_ax_theta(:,18));
    non_dim_ay_pure_F_L = (F_pure_ay_theta(:,2)-F_pure_ay_theta(:,11))./(F_pure_ay_theta(:,8).^2.*F_pure_ay_theta(:,7).*F_pure_ay_theta(:,18));
    non_dim_az_pure_F_L = (F_pure_az_theta(:,3)-F_pure_az_theta(:,12))./(F_pure_az_theta(:,8).^2.*F_pure_az_theta(:,7).*F_pure_az_theta(:,18));
    non_dim_ax_pure_F_R = (F_pure_ax_theta(:,4)-F_pure_ax_theta(:,13))./(F_pure_ax_theta(:,8).^2.*F_pure_ax_theta(:,7).*F_pure_ax_theta(:,18));
    non_dim_ay_pure_F_R = (F_pure_ay_theta(:,5)-F_pure_ay_theta(:,14))./(F_pure_ay_theta(:,8).^2.*F_pure_ay_theta(:,7).*F_pure_ay_theta(:,18));
    non_dim_az_pure_F_R = (F_pure_az_theta(:,6)-F_pure_az_theta(:,15))./(F_pure_az_theta(:,8).^2.*F_pure_az_theta(:,7).*F_pure_az_theta(:,18));    

    non_dim_roll_pure_M_L  = (M_pure_roll_theta(:,1)-M_pure_roll_theta(:,10))./(M_pure_roll_theta(:,8).^2.*M_pure_roll_theta(:,7).^2.*M_pure_roll_theta(:,18));
    non_dim_pitch_pure_M_L = (M_pure_pitch_theta(:,2)-M_pure_pitch_theta(:,11))./(M_pure_pitch_theta(:,8).^2.*M_pure_pitch_theta(:,7).^2.*M_pure_pitch_theta(:,18));
    non_dim_yaw_pure_M_L   = (M_pure_yaw_theta(:,3)-M_pure_yaw_theta(:,12))./(M_pure_yaw_theta(:,8).^2.*M_pure_yaw_theta(:,7).^2.*M_pure_yaw_theta(:,18));
    non_dim_roll_pure_M_R  = (M_pure_roll_theta(:,4)-M_pure_roll_theta(:,13))./(M_pure_roll_theta(:,8).^2.*M_pure_roll_theta(:,7).^2.*M_pure_roll_theta(:,18));
    non_dim_pitch_pure_M_R = (M_pure_pitch_theta(:,5)-M_pure_pitch_theta(:,14))./(M_pure_pitch_theta(:,8).^2.*M_pure_pitch_theta(:,7).^2.*M_pure_pitch_theta(:,18));
    non_dim_yaw_pure_M_R   = (M_pure_yaw_theta(:,6)-M_pure_yaw_theta(:,15))./(M_pure_yaw_theta(:,8).^2.*M_pure_yaw_theta(:,7).^2.*M_pure_yaw_theta(:,18));
    
    %----------------------------------------------------------------------

        
    n_a_theta = (2*(n_pol_theta+1));
    n_a_eta = (2*(n_pol_eta+1));
    n_a_phi = (2*(n_pol_phi+1));
    
    % average the symmetric motion and mirror the assymetric motion:

    non_dim_ax_F_theta_sym      = zeros(length(non_dim_ax_F_L),2+2*n_a_theta);
    non_dim_ay_F_theta_asym     = zeros(length(non_dim_ay_F_L),2+2*n_a_theta);
    non_dim_az_F_theta_sym      = zeros(length(non_dim_az_F_L),2+2*n_a_theta);
    
    non_dim_ax_F_eta_sym      = zeros(length(non_dim_ax_F_L),2+2*n_a_eta);
    non_dim_ay_F_eta_asym     = zeros(length(non_dim_ay_F_L),2+2*n_a_eta);
    non_dim_az_F_eta_sym      = zeros(length(non_dim_az_F_L),2+2*n_a_eta);
    
    non_dim_ax_F_phi_sym      = zeros(length(non_dim_ax_F_L),2+2*n_a_phi);
    non_dim_ay_F_phi_asym     = zeros(length(non_dim_ay_F_L),2+2*n_a_phi);
    non_dim_az_F_phi_sym      = zeros(length(non_dim_az_F_L),2+2*n_a_phi);    

    non_dim_ax_pure_F_theta_sym      = zeros(length(non_dim_ax_pure_F_L),2+2*n_a_theta);
    non_dim_ay_pure_F_theta_asym     = zeros(length(non_dim_ay_pure_F_L),2+2*n_a_theta);
    non_dim_az_pure_F_theta_sym      = zeros(length(non_dim_az_pure_F_L),2+2*n_a_theta);
    
    non_dim_ax_pure_F_eta_sym      = zeros(length(non_dim_ax_pure_F_L),2+2*n_a_eta);
    non_dim_ay_pure_F_eta_asym     = zeros(length(non_dim_ay_pure_F_L),2+2*n_a_eta);
    non_dim_az_pure_F_eta_sym      = zeros(length(non_dim_az_pure_F_L),2+2*n_a_eta);
    
    non_dim_ax_pure_F_phi_sym      = zeros(length(non_dim_ax_pure_F_L),2+2*n_a_phi);
    non_dim_ay_pure_F_phi_asym     = zeros(length(non_dim_ay_pure_F_L),2+2*n_a_phi);
    non_dim_az_pure_F_phi_sym      = zeros(length(non_dim_az_pure_F_L),2+2*n_a_phi);
    
    
    non_dim_roll_M_theta_asym  = zeros(length(non_dim_roll_M_L),2+2*n_a_theta);
    non_dim_pitch_M_theta_sym  = zeros(length(non_dim_pitch_M_L),2+2*n_a_theta);
    non_dim_yaw_M_theta_asym   = zeros(length(non_dim_yaw_M_L),2+2*n_a_theta);
    
    non_dim_roll_M_eta_asym    = zeros(length(non_dim_roll_M_L),2+2*n_a_eta);
    non_dim_pitch_M_eta_sym    = zeros(length(non_dim_pitch_M_L),2+2*n_a_eta);
    non_dim_yaw_M_eta_asym     = zeros(length(non_dim_yaw_M_L),2+2*n_a_eta);
    
    non_dim_roll_M_phi_asym    = zeros(length(non_dim_roll_M_L),2+2*n_a_phi);
    non_dim_pitch_M_phi_sym    = zeros(length(non_dim_pitch_M_L),2+2*n_a_phi);
    non_dim_yaw_M_phi_asym     = zeros(length(non_dim_yaw_M_L),2+2*n_a_phi);
    
    non_dim_roll_pure_M_theta_asym  = zeros(length(non_dim_roll_pure_M_L),2+2*n_a_theta);
    non_dim_pitch_pure_M_theta_sym  = zeros(length(non_dim_pitch_pure_M_L),2+2*n_a_theta);
    non_dim_yaw_pure_M_theta_asym   = zeros(length(non_dim_yaw_pure_M_L),2+2*n_a_theta);
    
    non_dim_roll_pure_M_eta_asym    = zeros(length(non_dim_roll_pure_M_L),2+2*n_a_eta);
    non_dim_pitch_pure_M_eta_sym    = zeros(length(non_dim_pitch_pure_M_L),2+2*n_a_eta);
    non_dim_yaw_pure_M_eta_asym     = zeros(length(non_dim_yaw_pure_M_L),2+2*n_a_eta);
    
    non_dim_roll_pure_M_phi_asym    = zeros(length(non_dim_roll_pure_M_L),2+2*n_a_phi);
    non_dim_pitch_pure_M_phi_sym    = zeros(length(non_dim_pitch_pure_M_L),2+2*n_a_phi);
    non_dim_yaw_pure_M_phi_asym     = zeros(length(non_dim_yaw_pure_M_L),2+2*n_a_phi);
 
    
    non_dim_freq = zeros(length(non_dim_az_F_L),8);
    
    for i = 1:length(non_dim_az_F_L)
        
        non_dim_az_F_L_t(1) = (F_az_theta(i,1)-F_az_theta(i,10))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));
        non_dim_az_F_L_t(2) = (F_az_theta(i,2)-F_az_theta(i,11))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));
        non_dim_az_F_L_t(3) = (F_az_theta(i,3)-F_az_theta(i,12))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));
        non_dim_az_F_R_t(1) = (F_az_theta(i,4)-F_az_theta(i,13))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));
        non_dim_az_F_R_t(2) = (F_az_theta(i,5)-F_az_theta(i,14))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));
        non_dim_az_F_R_t(3) = (F_az_theta(i,6)-F_az_theta(i,15))./(F_az_theta(i,8).^2.*F_az_theta(i,7).*F_az_theta(i,18));       
        
        non_dim_freq(i,:) = [ F_az_theta(i,8) F_az_theta(i,9) non_dim_az_F_L_t non_dim_az_F_R_t];
        
    end

    
    for i = 1:length(non_dim_ax_F_L)
        
        non_dim_ax_F_theta_sym(i,:)    = [non_dim_ax_F_L(i)  non_dim_ax_F_R(i)  F_ax_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_ax_F_eta_sym(i,:)      = [non_dim_ax_F_L(i)  non_dim_ax_F_R(i)  F_ax_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_ax_F_phi_sym(i,:)      = [non_dim_ax_F_L(i)  non_dim_ax_F_R(i)  F_ax_phi(i,18+(1:(2*n_a_phi)))];
   
    end    
    
    for i = 1:length(non_dim_ay_F_L)
        
        if (non_dim_ay_F_L(i)+non_dim_ay_F_R(i)) >= 0
    
            non_dim_ay_F_theta_asym(i,:)   = [non_dim_ay_F_L(i)    non_dim_ay_F_R(i)   F_ay_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_ay_F_eta_asym(i,:)     = [non_dim_ay_F_L(i)    non_dim_ay_F_R(i)   F_ay_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_ay_F_phi_asym(i,:)     = [non_dim_ay_F_L(i)    non_dim_ay_F_R(i)   F_ay_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_ay_F_L(i)+non_dim_ay_F_R(i)) < 0
            
            non_dim_ay_F_theta_asym(i,:)   = [-non_dim_ay_F_R(i)    -non_dim_ay_F_L(i)   F_ay_theta(i,18+n_a_theta+(1:(n_a_theta))) F_ay_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_ay_F_eta_asym(i,:)     = [-non_dim_ay_F_R(i)    -non_dim_ay_F_L(i)   F_ay_eta(i,18+n_a_eta+(1:(n_a_eta)))       F_ay_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_ay_F_phi_asym(i,:)     = [-non_dim_ay_F_R(i)    -non_dim_ay_F_L(i)   F_ay_phi(i,18+n_a_phi+(1:(n_a_phi)))       F_ay_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end
    
    for i = 1:length(non_dim_az_F_L)
        
        non_dim_az_F_theta_sym(i,:)    = [non_dim_az_F_L(i)  non_dim_az_F_R(i)  F_az_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_az_F_eta_sym(i,:)      = [non_dim_az_F_L(i)  non_dim_az_F_R(i)  F_az_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_az_F_phi_sym(i,:)      = [non_dim_az_F_L(i)  non_dim_az_F_R(i)  F_az_phi(i,18+(1:(2*n_a_phi)))];
   
    end  
    

    for i = 1:length(non_dim_ax_pure_F_L)
        
        non_dim_ax_pure_F_theta_sym(i,:)    = [non_dim_ax_pure_F_L(i)  non_dim_ax_pure_F_R(i)  F_pure_ax_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_ax_pure_F_eta_sym(i,:)      = [non_dim_ax_pure_F_L(i)  non_dim_ax_pure_F_R(i)  F_pure_ax_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_ax_pure_F_phi_sym(i,:)      = [non_dim_ax_pure_F_L(i)  non_dim_ax_pure_F_R(i)  F_pure_ax_phi(i,18+(1:(2*n_a_phi)))];
   
    end    
    
    for i = 1:length(non_dim_ay_pure_F_L)
        
        if (non_dim_ay_pure_F_L(i)+non_dim_ay_pure_F_R(i)) >= 0
    
            non_dim_ay_pure_F_theta_asym(i,:)   = [non_dim_ay_pure_F_L(i) non_dim_ay_pure_F_R(i) F_pure_ay_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_ay_pure_F_eta_asym(i,:)     = [non_dim_ay_pure_F_L(i) non_dim_ay_pure_F_R(i) F_pure_ay_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_ay_pure_F_phi_asym(i,:)     = [non_dim_ay_pure_F_L(i) non_dim_ay_pure_F_R(i) F_pure_ay_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_ay_pure_F_L(i)+non_dim_ay_pure_F_R(i)) < 0
            
            non_dim_ay_pure_F_theta_asym(i,:)   = [-non_dim_ay_pure_F_R(i)   -non_dim_ay_pure_F_L(i) F_pure_ay_theta(i,18+n_a_theta+(1:(n_a_theta))) F_pure_ay_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_ay_pure_F_eta_asym(i,:)     = [-non_dim_ay_pure_F_R(i)   -non_dim_ay_pure_F_L(i) F_pure_ay_eta(i,18+n_a_eta+(1:(n_a_eta)))       F_pure_ay_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_ay_pure_F_phi_asym(i,:)     = [-non_dim_ay_pure_F_R(i)   -non_dim_ay_pure_F_L(i) F_pure_ay_phi(i,18+n_a_phi+(1:(n_a_phi)))       F_pure_ay_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end
    
    for i = 1:length(non_dim_az_pure_F_L)
        
        non_dim_az_pure_F_theta_sym(i,:)    = [non_dim_az_pure_F_L(i)  non_dim_az_pure_F_R(i)    F_pure_az_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_az_pure_F_eta_sym(i,:)      = [non_dim_az_pure_F_L(i)  non_dim_az_pure_F_R(i)    F_pure_az_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_az_pure_F_phi_sym(i,:)      = [non_dim_az_pure_F_L(i)  non_dim_az_pure_F_R(i)    F_pure_az_phi(i,18+(1:(2*n_a_phi)))];
   
    end  
    
    
    for i = 1:length(non_dim_roll_M_L)
        
        if (non_dim_roll_M_L(i)+non_dim_roll_M_R(i)) >= 0
    
            non_dim_roll_M_theta_asym(i,:)   = [non_dim_roll_M_L(i) non_dim_roll_M_R(i) M_roll_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_roll_M_eta_asym(i,:)     = [non_dim_roll_M_L(i) non_dim_roll_M_R(i) M_roll_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_roll_M_phi_asym(i,:)     = [non_dim_roll_M_L(i) non_dim_roll_M_R(i) M_roll_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_roll_M_L(i)+non_dim_roll_M_R(i)) < 0
            
            non_dim_roll_M_theta_asym(i,:)   = [-non_dim_roll_M_R(i)   -non_dim_roll_M_L(i) M_roll_theta(i,18+n_a_theta+(1:(n_a_theta))) M_roll_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_roll_M_eta_asym(i,:)     = [-non_dim_roll_M_R(i)   -non_dim_roll_M_L(i) M_roll_eta(i,18+n_a_eta+(1:(n_a_eta)))       M_roll_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_roll_M_phi_asym(i,:)     = [-non_dim_roll_M_R(i)   -non_dim_roll_M_L(i) M_roll_phi(i,18+n_a_phi+(1:(n_a_phi)))       M_roll_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end    
    
    
    for i = 1:length(non_dim_pitch_M_L)

        
        non_dim_pitch_M_theta_sym(i,:)    = [non_dim_pitch_M_L(i)  non_dim_pitch_M_R(i)  M_pitch_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_pitch_M_eta_sym(i,:)      = [non_dim_pitch_M_L(i)  non_dim_pitch_M_R(i)  M_pitch_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_pitch_M_phi_sym(i,:)      = [non_dim_pitch_M_L(i)  non_dim_pitch_M_R(i)  M_pitch_phi(i,18+(1:(2*n_a_phi)))];

   
    end
    
    for i = 1:length(non_dim_yaw_M_L)
        
        if (non_dim_yaw_M_L(i)+non_dim_yaw_M_R(i)) >= 0
    
            non_dim_yaw_M_theta_asym(i,:)   = [non_dim_yaw_M_L(i) non_dim_yaw_M_R(i) M_yaw_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_yaw_M_eta_asym(i,:)     = [non_dim_yaw_M_L(i) non_dim_yaw_M_R(i) M_yaw_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_yaw_M_phi_asym(i,:)     = [non_dim_yaw_M_L(i) non_dim_yaw_M_R(i) M_yaw_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_yaw_M_L(i)+non_dim_yaw_M_R(i)) < 0
            
            non_dim_yaw_M_theta_asym(i,:)   = [-non_dim_yaw_M_R(i)   -non_dim_yaw_M_L(i) M_yaw_theta(i,18+n_a_theta+(1:(n_a_theta))) M_yaw_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_yaw_M_eta_asym(i,:)     = [-non_dim_yaw_M_R(i)   -non_dim_yaw_M_L(i) M_yaw_eta(i,18+n_a_eta+(1:(n_a_eta)))       M_yaw_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_yaw_M_phi_asym(i,:)     = [-non_dim_yaw_M_R(i)   -non_dim_yaw_M_L(i) M_yaw_phi(i,18+n_a_phi+(1:(n_a_phi)))       M_yaw_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end     

    

    
    
    for i = 1:length(non_dim_roll_pure_M_L)
        
        if (non_dim_roll_pure_M_L(i)+non_dim_roll_pure_M_R(i)) >= 0
    
            non_dim_roll_pure_M_theta_asym(i,:)   = [non_dim_roll_pure_M_L(i) non_dim_roll_pure_M_R(i) M_pure_roll_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_roll_pure_M_eta_asym(i,:)     = [non_dim_roll_pure_M_L(i) non_dim_roll_pure_M_R(i) M_pure_roll_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_roll_pure_M_phi_asym(i,:)     = [non_dim_roll_pure_M_L(i) non_dim_roll_pure_M_R(i) M_pure_roll_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_roll_pure_M_L(i)+non_dim_roll_pure_M_R(i)) < 0
            
            non_dim_roll_pure_M_theta_asym(i,:)   = [-non_dim_roll_pure_M_R(i) -non_dim_roll_pure_M_L(i) M_pure_roll_theta(i,18+n_a_theta+(1:(n_a_theta))) M_pure_roll_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_roll_pure_M_eta_asym(i,:)     = [-non_dim_roll_pure_M_R(i) -non_dim_roll_pure_M_L(i) M_pure_roll_eta(i,18+n_a_eta+(1:(n_a_eta)))       M_pure_roll_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_roll_pure_M_phi_asym(i,:)     = [-non_dim_roll_pure_M_R(i) -non_dim_roll_pure_M_L(i) M_pure_roll_phi(i,18+n_a_phi+(1:(n_a_phi)))       M_pure_roll_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end    
    
    
    for i = 1:length(non_dim_pitch_pure_M_L)
        
        non_dim_pitch_pure_M_theta_sym(i,:)    = [non_dim_pitch_pure_M_L(i)  non_dim_pitch_pure_M_R(i)  M_pure_pitch_theta(i,18+(1:(2*n_a_theta)))];
        non_dim_pitch_pure_M_eta_sym(i,:)      = [non_dim_pitch_pure_M_L(i)  non_dim_pitch_pure_M_R(i)  M_pure_pitch_eta(i,18+(1:(2*n_a_eta)))];
        non_dim_pitch_pure_M_phi_sym(i,:)      = [non_dim_pitch_pure_M_L(i)  non_dim_pitch_pure_M_R(i)  M_pure_pitch_phi(i,18+(1:(2*n_a_phi)))];
   
    end
    
    for i = 1:length(non_dim_yaw_pure_M_L)
        
        if (non_dim_yaw_pure_M_L(i)+non_dim_yaw_pure_M_R(i)) >= 0
    
            non_dim_yaw_pure_M_theta_asym(i,:)   = [non_dim_yaw_pure_M_L(i) non_dim_yaw_pure_M_R(i) M_pure_yaw_theta(i,18+(1:(2*n_a_theta)))];
            non_dim_yaw_pure_M_eta_asym(i,:)     = [non_dim_yaw_pure_M_L(i) non_dim_yaw_pure_M_R(i) M_pure_yaw_eta(i,18+(1:(2*n_a_eta)))];
            non_dim_yaw_pure_M_phi_asym(i,:)     = [non_dim_yaw_pure_M_L(i) non_dim_yaw_pure_M_R(i) M_pure_yaw_phi(i,18+(1:(2*n_a_phi)))];
        
        elseif (non_dim_yaw_pure_M_L(i)+non_dim_yaw_pure_M_R(i)) < 0
            
            non_dim_yaw_pure_M_theta_asym(i,:)   = [-non_dim_yaw_pure_M_R(i)   -non_dim_yaw_pure_M_L(i) M_pure_yaw_theta(i,18+n_a_theta+(1:(n_a_theta))) M_pure_yaw_theta(i,18+(1:(n_a_theta)))  ];
            non_dim_yaw_pure_M_eta_asym(i,:)     = [-non_dim_yaw_pure_M_R(i)   -non_dim_yaw_pure_M_L(i) M_pure_yaw_eta(i,18+n_a_eta+(1:(n_a_eta)))       M_pure_yaw_eta(i,18+(1:(n_a_eta)))      ];
            non_dim_yaw_pure_M_phi_asym(i,:)     = [-non_dim_yaw_pure_M_R(i)   -non_dim_yaw_pure_M_L(i) M_pure_yaw_phi(i,18+n_a_phi+(1:(n_a_phi)))       M_pure_yaw_phi(i,18+(1:(n_a_phi)))      ];            
            
        end
    
    end     
    
    
    max(non_dim_ax_F_L+non_dim_ax_F_R)
    max(non_dim_ay_F_L+non_dim_ay_F_R)
    max(non_dim_az_F_L+non_dim_az_F_R)
    
    max(non_dim_roll_pure_M_L+non_dim_roll_pure_M_R)
    max(non_dim_pitch_pure_M_L+non_dim_pitch_pure_M_R)
    max(non_dim_yaw_pure_M_L+non_dim_yaw_pure_M_R)
    
    min(non_dim_ax_F_L+non_dim_ax_F_R)
    min(non_dim_ay_F_L+non_dim_ay_F_R)
    min(non_dim_az_F_L+non_dim_az_F_R)
    
    min(non_dim_roll_pure_M_L+non_dim_roll_pure_M_R)
    min(non_dim_pitch_pure_M_L+non_dim_pitch_pure_M_R)
    min(non_dim_yaw_pure_M_L+non_dim_yaw_pure_M_R)
    
    %----------------------------------------------------------------------
    
    
    dev_fit = {};
    
    [ dev_fit.ax_theta_p, dev_fit.ax_theta_n ] = sym_motion_regression( n_a_theta, non_dim_ax_F_theta_sym);
    [ dev_fit.ax_eta_p, dev_fit.ax_eta_n ] = sym_motion_regression( n_a_eta, non_dim_ax_F_eta_sym);
    [ dev_fit.ax_phi_p, dev_fit.ax_phi_n ] = sym_motion_regression( n_a_phi, non_dim_ax_F_phi_sym);
    
    [ dev_fit.ay_theta_p, dev_fit.ay_theta_n, dev_fit.ay_theta_pn ] = asym_motion_regression( n_a_theta, non_dim_ay_F_theta_asym );
    [ dev_fit.ay_eta_p, dev_fit.ay_eta_n, dev_fit.ay_eta_pn ] = asym_motion_regression( n_a_eta, non_dim_ay_F_eta_asym );
    [ dev_fit.ay_phi_p, dev_fit.ay_phi_n, dev_fit.ay_phi_pn ] = asym_motion_regression( n_a_phi, non_dim_ay_F_phi_asym );

    [ dev_fit.az_theta_p, dev_fit.az_theta_n ] = sym_motion_regression( n_a_theta, non_dim_az_F_theta_sym);
    [ dev_fit.az_eta_p, dev_fit.az_eta_n ] = sym_motion_regression( n_a_eta, non_dim_az_F_eta_sym);
    [ dev_fit.az_phi_p, dev_fit.az_phi_n ] = sym_motion_regression( n_a_phi, non_dim_az_F_phi_sym);
    
    [ dev_fit.roll_theta_p, dev_fit.roll_theta_n, dev_fit.roll_theta_pn ] = asym_motion_regression( n_a_theta, non_dim_roll_M_theta_asym );
    [ dev_fit.roll_eta_p, dev_fit.roll_eta_n, dev_fit.roll_eta_pn ] = asym_motion_regression( n_a_eta, non_dim_roll_M_eta_asym );
    [ dev_fit.roll_phi_p, dev_fit.roll_phi_n, dev_fit.roll_phi_pn ] = asym_motion_regression( n_a_phi, non_dim_roll_M_phi_asym );
    
    [ dev_fit.pitch_theta_p, dev_fit.pitch_theta_n ] = sym_motion_regression( n_a_theta, non_dim_pitch_M_theta_sym);
    [ dev_fit.pitch_eta_p, dev_fit.pitch_eta_n ] = sym_motion_regression( n_a_eta, non_dim_pitch_M_eta_sym);
    [ dev_fit.pitch_phi_p, dev_fit.pitch_phi_n ] = sym_motion_regression( n_a_phi, non_dim_pitch_M_phi_sym);
    
    [ dev_fit.yaw_theta_p, dev_fit.yaw_theta_n, dev_fit.yaw_theta_pn ] = asym_motion_regression( n_a_theta, non_dim_yaw_M_theta_asym );
    [ dev_fit.yaw_eta_p, dev_fit.yaw_eta_n, dev_fit.yaw_eta_pn ] = asym_motion_regression( n_a_eta, non_dim_yaw_M_eta_asym );
    [ dev_fit.yaw_phi_p, dev_fit.yaw_phi_n, dev_fit.yaw_phi_pn ] = asym_motion_regression( n_a_phi, non_dim_yaw_M_phi_asym );
    
    
    freq_regression( non_dim_freq )


    % plot the fits:
    
%     plot_sym_motion( n_a_theta, non_dim_ax_F_theta_sym, dev_fit.ax_theta_p, dev_fit.ax_theta_n )
%     plot_sym_motion( n_a_eta, non_dim_ax_F_eta_sym, dev_fit.ax_eta_p, dev_fit.ax_eta_n )
%     plot_sym_motion( n_a_phi, non_dim_ax_F_phi_sym, dev_fit.ax_phi_p, dev_fit.ax_phi_n )
% 
%     pause
%     
%     close all
%     
%     plot_asym_motion( n_a_theta, non_dim_ay_F_theta_asym, dev_fit.ay_theta_p, dev_fit.ay_theta_n, dev_fit.ay_theta_pn )
%     plot_asym_motion( n_a_eta, non_dim_ay_F_eta_asym, dev_fit.ay_eta_p, dev_fit.ay_eta_n, dev_fit.ay_eta_pn )
%     plot_asym_motion( n_a_phi, non_dim_ay_F_phi_asym, dev_fit.ay_phi_p, dev_fit.ay_phi_n, dev_fit.ay_phi_pn )
%     
%     pause
%     
%     close all
%     
%     plot_sym_motion( n_a_theta, non_dim_az_F_theta_sym, dev_fit.az_theta_p, dev_fit.az_theta_n )
%     plot_sym_motion( n_a_eta, non_dim_az_F_eta_sym, dev_fit.az_eta_p, dev_fit.az_eta_n )
%     plot_sym_motion( n_a_phi, non_dim_az_F_phi_sym, dev_fit.az_phi_p, dev_fit.az_phi_n )
%     
%     pause
%     
%     close all
%     
%     plot_asym_motion( n_a_theta, non_dim_roll_M_theta_asym, dev_fit.roll_theta_p, dev_fit.roll_theta_n, dev_fit.roll_theta_pn )
%     plot_asym_motion( n_a_eta, non_dim_roll_M_eta_asym, dev_fit.roll_eta_p, dev_fit.roll_eta_n, dev_fit.roll_eta_pn )
%     plot_asym_motion( n_a_phi, non_dim_roll_M_phi_asym, dev_fit.roll_phi_p, dev_fit.roll_phi_n, dev_fit.roll_phi_pn )
%     
%     pause
%     
%     close all
%     
%     plot_sym_motion( n_a_theta, non_dim_pitch_M_theta_sym, dev_fit.pitch_theta_p, dev_fit.pitch_theta_n )
%     plot_sym_motion( n_a_eta, non_dim_pitch_M_eta_sym, dev_fit.pitch_eta_p, dev_fit.pitch_eta_n )
%     plot_sym_motion( n_a_phi, non_dim_pitch_M_phi_sym, dev_fit.pitch_phi_p, dev_fit.pitch_phi_n )
%     
%     pause
%     
%     close all
%     
%     plot_asym_motion( n_a_theta, non_dim_yaw_M_theta_asym, dev_fit.yaw_theta_p, dev_fit.yaw_theta_n, dev_fit.yaw_theta_pn )
%     plot_asym_motion( n_a_eta, non_dim_yaw_M_eta_asym, dev_fit.yaw_eta_p, dev_fit.yaw_eta_n, dev_fit.yaw_eta_pn )
%     plot_asym_motion( n_a_phi, non_dim_yaw_M_phi_asym, dev_fit.yaw_phi_p, dev_fit.yaw_phi_n, dev_fit.yaw_phi_pn )
%     
%     pause
%     
%     close all



    % Create d_a_fit as a function of force or moment:
    
    cor_coeff = 0;
    
    order_dev_theta = n_pol_theta+1;
    order_dev_eta = n_pol_eta+1;
    order_dev_phi = n_pol_phi+1;
    
    a_dev_man = {};

    a_dev_man.ax_theta_forward  = a_dev_fit( n_pol_theta, dev_fit.ax_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.ax_theta_back     = a_dev_fit( n_pol_theta, dev_fit.ax_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.ax_eta_forward    = a_dev_fit( n_pol_eta, dev_fit.ax_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.ax_eta_back       = a_dev_fit( n_pol_eta, dev_fit.ax_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.ax_phi_forward    = a_dev_fit( n_pol_phi, dev_fit.ax_phi_p, order_dev_phi, cor_coeff );
    a_dev_man.ax_phi_back       = a_dev_fit( n_pol_phi, dev_fit.ax_phi_n, order_dev_phi, cor_coeff );
    
    a_dev_man.ay_theta_left     = a_dev_fit( n_pol_theta, dev_fit.ay_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.ay_theta_right    = a_dev_fit( n_pol_theta, dev_fit.ay_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.ay_eta_left       = a_dev_fit( n_pol_eta, dev_fit.ay_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.ay_eta_right      = a_dev_fit( n_pol_eta, dev_fit.ay_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.ay_phi_left       = a_dev_fit( n_pol_phi, dev_fit.ay_phi_p, order_dev_phi, cor_coeff );
    a_dev_man.ay_phi_right      = a_dev_fit( n_pol_phi, dev_fit.ay_phi_n, order_dev_phi, cor_coeff );
    
    a_dev_man.az_theta_down     = a_dev_fit( n_pol_theta, dev_fit.az_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.az_theta_up       = a_dev_fit( n_pol_theta, dev_fit.az_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.az_eta_down       = a_dev_fit( n_pol_eta, dev_fit.az_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.az_eta_up         = a_dev_fit( n_pol_eta, dev_fit.az_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.az_phi_down    	= a_dev_fit( n_pol_phi, dev_fit.az_phi_p, order_dev_phi, cor_coeff );
    a_dev_man.az_phi_up         = a_dev_fit( n_pol_phi, dev_fit.az_phi_n, order_dev_phi, cor_coeff );
    
    a_dev_man.roll_theta_left   = a_dev_fit( n_pol_theta, dev_fit.roll_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.roll_theta_right  = a_dev_fit( n_pol_theta, dev_fit.roll_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.roll_eta_left     = a_dev_fit( n_pol_eta, dev_fit.roll_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.roll_eta_right    = a_dev_fit( n_pol_eta, dev_fit.roll_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.roll_phi_left     = a_dev_fit( n_pol_phi, dev_fit.roll_phi_p, order_dev_phi, cor_coeff );
    a_dev_man.roll_phi_right    = a_dev_fit( n_pol_phi, dev_fit.roll_phi_n, order_dev_phi, cor_coeff );
    
    a_dev_man.pitch_theta_forward   = a_dev_fit( n_pol_theta, dev_fit.pitch_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.pitch_theta_back      = a_dev_fit( n_pol_theta, dev_fit.pitch_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.pitch_eta_forward     = a_dev_fit( n_pol_eta, dev_fit.pitch_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.pitch_eta_back        = a_dev_fit( n_pol_eta, dev_fit.pitch_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.pitch_phi_forward     = a_dev_fit( n_pol_phi, dev_fit.pitch_phi_n, order_dev_phi, cor_coeff );
    a_dev_man.pitch_phi_back        = a_dev_fit( n_pol_phi, dev_fit.pitch_phi_p, order_dev_phi, cor_coeff );
        
    a_dev_man.yaw_theta_left    = a_dev_fit( n_pol_theta, dev_fit.yaw_theta_p, order_dev_theta, cor_coeff );
    a_dev_man.yaw_theta_right   = a_dev_fit( n_pol_theta, dev_fit.yaw_theta_n, order_dev_theta, cor_coeff );
    a_dev_man.yaw_eta_left      = a_dev_fit( n_pol_eta, dev_fit.yaw_eta_p, order_dev_eta, cor_coeff );
    a_dev_man.yaw_eta_right     = a_dev_fit( n_pol_eta, dev_fit.yaw_eta_n, order_dev_eta, cor_coeff );
    a_dev_man.yaw_phi_left      = a_dev_fit( n_pol_phi, dev_fit.yaw_phi_p, order_dev_phi, cor_coeff );
    a_dev_man.yaw_phi_right     = a_dev_fit( n_pol_phi, dev_fit.yaw_phi_n, order_dev_phi, cor_coeff );
    
    
    
    [ t, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, 0.55, 200, 0, 1, 0 );
    [ ~, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, 0.55, 200, 0, 1, 0 );
    [ ~, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, 0.55, 200, 0, 1, 0 );
    
    
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.ax_theta_forward,t,X_eta*a_dev_man.ax_eta_forward,t,X_phi*a_dev_man.ax_phi_forward)
    title('Forward acceleration in x-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.ax_theta_back,t,X_eta*a_dev_man.ax_eta_back,t,X_phi*a_dev_man.ax_phi_back)
    title('Backward acceleration in x-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off
    
    
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.ay_theta_left,t,X_eta*a_dev_man.ay_eta_left,t,X_phi*a_dev_man.ay_phi_left)
    title('Left wingbeat, positive acceleration in y-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.ay_theta_right,t,X_eta*a_dev_man.ay_eta_right,t,X_phi*a_dev_man.ay_phi_right)
    title('Right wingbeat, positive acceleration in y-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off    
    
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.az_theta_down,t,X_eta*a_dev_man.az_eta_down,t,X_phi*a_dev_man.az_phi_down)
    title(' Downward acceleration in z-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.az_theta_up,t,X_eta*a_dev_man.az_eta_up,t,X_phi*a_dev_man.az_phi_up)
    title('Upward acceleration in z-direction')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off 
        
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.roll_theta_left,t,X_eta*a_dev_man.roll_eta_left,t,X_phi*a_dev_man.roll_phi_left)
    title('Left wingbeat, positive roll motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.roll_theta_right,t,X_eta*a_dev_man.roll_eta_right,t,X_phi*a_dev_man.roll_phi_right)
    title('Right wingbeat, positive roll motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off
    
    
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.pitch_theta_forward,t,X_eta*a_dev_man.pitch_eta_forward,t,X_phi*a_dev_man.pitch_phi_forward)
    title('Pitch up motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.pitch_theta_back,t,X_eta*a_dev_man.pitch_eta_back,t,X_phi*a_dev_man.pitch_phi_back)
    title('Pitch down motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off
    
    
    figure()
    hold on
    subplot(2,1,1); plot(t,X_theta*a_dev_man.yaw_theta_left,t,X_eta*a_dev_man.yaw_eta_left,t,X_phi*a_dev_man.yaw_phi_left)
    title('Left wingbeat, positive yaw motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    subplot(2,1,2); plot(t,X_theta*a_dev_man.yaw_theta_right,t,X_eta*a_dev_man.yaw_eta_right,t,X_phi*a_dev_man.yaw_phi_right)
    title('Right wingbeat, positive yaw motion')
    legend('\theta','\eta','\phi')
    xlabel('non dimensional wing beat time')
    ylabel('delta angle [rad]')
    hold off
    
    
    
    a_dev_glob = {};
    
    f_glob = pathDB.f_glob;
    
    down_up_glob = pathDB.down_up_glob;
    
    w_l_glob = mean(pathDB.wing_l);
    
    a_glob_theta = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]; 
    a_glob_eta = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]; 
    a_glob_phi = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]; 
    
    a_dev_glob.ax_theta_forward  = a_man_wb_FM( n_pol_theta, a_dev_man.ax_theta_forward, down_up_glob );
    a_dev_glob.ax_theta_back     = a_man_wb_FM( n_pol_theta, a_dev_man.ax_theta_back, down_up_glob );
    a_dev_glob.ax_eta_forward    = a_man_wb_FM( n_pol_eta, a_dev_man.ax_eta_forward, down_up_glob );
    a_dev_glob.ax_eta_back       = a_man_wb_FM( n_pol_eta, a_dev_man.ax_eta_back, down_up_glob );
    a_dev_glob.ax_phi_forward    = a_man_wb_FM( n_pol_phi, a_dev_man.ax_phi_forward, down_up_glob );
    a_dev_glob.ax_phi_back       = a_man_wb_FM( n_pol_phi, a_dev_man.ax_phi_back, down_up_glob );
    
    a_dev_glob.ay_theta_left     = a_man_wb_FM( n_pol_theta, a_dev_man.ay_theta_left, down_up_glob );
    a_dev_glob.ay_theta_right    = a_man_wb_FM( n_pol_theta, a_dev_man.ay_theta_right, down_up_glob );
    a_dev_glob.ay_eta_left       = a_man_wb_FM( n_pol_eta, a_dev_man.ay_eta_left, down_up_glob );
    a_dev_glob.ay_eta_right      = a_man_wb_FM( n_pol_eta, a_dev_man.ay_eta_right, down_up_glob );
    a_dev_glob.ay_phi_left       = a_man_wb_FM( n_pol_phi, a_dev_man.ay_phi_left, down_up_glob );
    a_dev_glob.ay_phi_right      = a_man_wb_FM( n_pol_phi, a_dev_man.ay_phi_right, down_up_glob );
    
    a_dev_glob.az_theta_down     = a_man_wb_FM( n_pol_theta, a_dev_man.az_theta_down, down_up_glob );
    a_dev_glob.az_theta_up       = a_man_wb_FM( n_pol_theta, a_dev_man.az_theta_up, down_up_glob );
    a_dev_glob.az_eta_down       = a_man_wb_FM( n_pol_eta, a_dev_man.az_eta_down, down_up_glob );
    a_dev_glob.az_eta_up         = a_man_wb_FM( n_pol_eta, a_dev_man.az_eta_up, down_up_glob );
    a_dev_glob.az_phi_down    	 = a_man_wb_FM( n_pol_phi, a_dev_man.az_phi_down, down_up_glob );
    a_dev_glob.az_phi_up         = a_man_wb_FM( n_pol_phi, a_dev_man.az_phi_up, down_up_glob );
    
    a_dev_glob.roll_theta_left   = a_man_wb_FM( n_pol_theta, a_dev_man.roll_theta_left, down_up_glob );
    a_dev_glob.roll_theta_right  = a_man_wb_FM( n_pol_theta, a_dev_man.roll_theta_right, down_up_glob );
    a_dev_glob.roll_eta_left     = a_man_wb_FM( n_pol_eta, a_dev_man.roll_eta_left, down_up_glob );
    a_dev_glob.roll_eta_right    = a_man_wb_FM( n_pol_eta, a_dev_man.roll_eta_right, down_up_glob );
    a_dev_glob.roll_phi_left     = a_man_wb_FM( n_pol_phi, a_dev_man.roll_phi_left, down_up_glob );
    a_dev_glob.roll_phi_right    = a_man_wb_FM( n_pol_phi, a_dev_man.roll_phi_right, down_up_glob );
    
    a_dev_glob.pitch_theta_forward   = a_man_wb_FM( n_pol_theta, a_dev_man.pitch_theta_forward, down_up_glob );
    a_dev_glob.pitch_theta_back      = a_man_wb_FM( n_pol_theta, a_dev_man.pitch_theta_back, down_up_glob );
    a_dev_glob.pitch_eta_forward     = a_man_wb_FM( n_pol_eta, a_dev_man.pitch_eta_forward, down_up_glob );
    a_dev_glob.pitch_eta_back        = a_man_wb_FM( n_pol_eta, a_dev_man.pitch_eta_back, down_up_glob );
    a_dev_glob.pitch_phi_forward     = a_man_wb_FM( n_pol_phi, a_dev_man.pitch_phi_forward, down_up_glob );
    a_dev_glob.pitch_phi_back        = a_man_wb_FM( n_pol_phi, a_dev_man.pitch_phi_back, down_up_glob );
        
    a_dev_glob.yaw_theta_left    = a_man_wb_FM( n_pol_theta, a_dev_man.yaw_theta_left, down_up_glob );
    a_dev_glob.yaw_theta_right   = a_man_wb_FM( n_pol_theta, a_dev_man.yaw_theta_right, down_up_glob );
    a_dev_glob.yaw_eta_left      = a_man_wb_FM( n_pol_eta, a_dev_man.yaw_eta_left, down_up_glob );
    a_dev_glob.yaw_eta_right     = a_man_wb_FM( n_pol_eta, a_dev_man.yaw_eta_right, down_up_glob );
    a_dev_glob.yaw_phi_left      = a_man_wb_FM( n_pol_phi, a_dev_man.yaw_phi_left, down_up_glob );
    a_dev_glob.yaw_phi_right     = a_man_wb_FM( n_pol_phi, a_dev_man.yaw_phi_right, down_up_glob );
    
    
    F_ax_forward_range = 0.0198;
    F_ax_back_range = -0.0362;
    
    F_ay_range = 0.0281;
    
    F_az_down_range = 0.0646;
    F_az_up_range = -0.0662;
    
    M_roll_range = 0.0193;
    
    M_pitch_forward_range = -1.752e-4;
    M_pitch_back_range = 0.0138;
    
    M_yaw_range = 0.0108;
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.ax_theta_forward, a_dev_glob.ax_eta_forward, a_dev_glob.ax_phi_forward, F_ax_forward_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.ax_theta_back, a_dev_glob.ax_eta_back, a_dev_glob.ax_phi_back, F_ax_back_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    legend('forward','backward')
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.ay_theta_left, a_dev_glob.ay_eta_left, a_dev_glob.ay_phi_left, F_ay_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.ay_theta_right, a_dev_glob.ay_eta_right, a_dev_glob.ay_phi_right, F_ay_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    legend('left wing','right wing')
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.az_theta_down, a_dev_glob.az_eta_down, a_dev_glob.az_phi_down, F_az_down_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.az_theta_up, a_dev_glob.az_eta_up, a_dev_glob.az_phi_up, F_az_up_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.roll_theta_left, a_dev_glob.roll_eta_left, a_dev_glob.roll_phi_left, M_roll_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.roll_theta_right, a_dev_glob.roll_eta_right, a_dev_glob.roll_phi_right, M_roll_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.pitch_theta_forward, a_dev_glob.pitch_eta_forward, a_dev_glob.pitch_phi_forward, M_pitch_forward_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.pitch_theta_back, a_dev_glob.pitch_eta_back, a_dev_glob.pitch_phi_back, M_pitch_back_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    
    figure()
    hold on
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.yaw_theta_left, a_dev_glob.yaw_eta_left, a_dev_glob.yaw_phi_left, M_yaw_range, f_glob, down_up_glob, 8 ,'r')
    plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_glob_theta, a_glob_eta, a_glob_phi, a_dev_glob.yaw_theta_right, a_dev_glob.yaw_eta_right, a_dev_glob.yaw_phi_right, M_yaw_range, f_glob, down_up_glob, 8 ,'b')
    hold off
    



    
end

