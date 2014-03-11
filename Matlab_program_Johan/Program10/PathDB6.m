function PathDB6( settings, pathDB)


    % Program that calculates and saves the mean value and the standard
    % deviation of a number of variables characterizing the body and wing
    % kinematics per 4 wingbeats. Every sequence of 4 wingbeats overlaps
    % with 2 wingbeats of the previous and 2 wingbeats of the following
    % sequences of 4 wingbeats. This data is used to generate histograms
    % which can show a possible relation between the flight parameters and
    % wing and body kinematics.
    
    savefile = 'pathDB6.mat';
    
    % Create nan matrices:
    
    n = 150;
    
    % Body:
    
    phi_body_mean = nan(n,2,size(pathDB.x,2));
    theta_body_mean = nan(n,2,size(pathDB.x,2));
    xsi_body_mean = nan(n,2,size(pathDB.x,2));
    
    phi_body_sd = nan(n,2,size(pathDB.x,2));
    theta_body_sd = nan(n,2,size(pathDB.x,2));
    xsi_body_sd = nan(n,2,size(pathDB.x,2));
    
    omegax_body_mean = nan(n,2,size(pathDB.x,2));
    omegay_body_mean = nan(n,2,size(pathDB.x,2));
    omegaz_body_mean = nan(n,2,size(pathDB.x,2));
    Omega_body_mean = nan(n,2,size(pathDB.x,2));
    
    omegax_body_sd = nan(n,2,size(pathDB.x,2));
    omegay_body_sd = nan(n,2,size(pathDB.x,2));
    omegaz_body_sd = nan(n,2,size(pathDB.x,2));
    Omega_body_sd = nan(n,2,size(pathDB.x,2));
    
    alfa_body_mean = nan(n,2,size(pathDB.x,2)); 
    beta_body_mean = nan(n,2,size(pathDB.x,2));
    
    alfa_body_sd = nan(n,2,size(pathDB.x,2));
    beta_body_sd = nan(n,2,size(pathDB.x,2));
    
    u_body_mean = nan(n,2,size(pathDB.x,2));
    v_body_mean = nan(n,2,size(pathDB.x,2));
    w_body_mean = nan(n,2,size(pathDB.x,2));
    U_body_mean = nan(n,2,size(pathDB.x,2));
    
    u_body_sd = nan(n,2,size(pathDB.x,2)); 
    v_body_sd = nan(n,2,size(pathDB.x,2));
    w_body_sd = nan(n,2,size(pathDB.x,2));
    U_body_sd = nan(n,2,size(pathDB.x,2));
    
    ax_body_mean = nan(n,2,size(pathDB.x,2));
    ay_body_mean = nan(n,2,size(pathDB.x,2));
    az_body_mean = nan(n,2,size(pathDB.x,2));
    a_body_mean = nan(n,2,size(pathDB.x,2));
    
    ax_body_sd = nan(n,2,size(pathDB.x,2));
    ay_body_sd = nan(n,2,size(pathDB.x,2));
    az_body_sd = nan(n,2,size(pathDB.x,2));
    a_body_sd = nan(n,2,size(pathDB.x,2));
    
    % Left wing:
    
    phi_L_up_mean = nan(n,2,size(pathDB.x,2)); 
    theta_L_up_mean = nan(n,2,size(pathDB.x,2));
    eta_L_up_mean = nan(n,2,size(pathDB.x,2));
    
    phi_L_up_sd = nan(n,2,size(pathDB.x,2));
    theta_L_up_sd = nan(n,2,size(pathDB.x,2));
    eta_L_up_sd = nan(n,2,size(pathDB.x,2));
    
    phi_L_down_mean = nan(n,2,size(pathDB.x,2));
    theta_L_down_mean = nan(n,2,size(pathDB.x,2));
    eta_L_down_mean = nan(n,2,size(pathDB.x,2));
    
    phi_L_down_sd = nan(n,2,size(pathDB.x,2));
    theta_L_down_sd = nan(n,2,size(pathDB.x,2));
    eta_L_down_sd = nan(n,2,size(pathDB.x,2));
    
    alfa_L_up_mean = nan(n,2,size(pathDB.x,2)); 
    beta_L_up_mean = nan(n,2,size(pathDB.x,2));
    
    alfa_L_up_sd = nan(n,2,size(pathDB.x,2));
    beta_L_up_sd = nan(n,2,size(pathDB.x,2));
    
    alfa_L_down_mean = nan(n,2,size(pathDB.x,2));
    beta_L_down_mean = nan(n,2,size(pathDB.x,2));
    
    alfa_L_down_sd = nan(n,2,size(pathDB.x,2));
    beta_L_down_sd = nan(n,2,size(pathDB.x,2));
    
    u_L_up_mean = nan(n,2,size(pathDB.x,2));
    v_L_up_mean = nan(n,2,size(pathDB.x,2));
    w_L_up_mean = nan(n,2,size(pathDB.x,2));
    U_L_up_mean = nan(n,2,size(pathDB.x,2));
    
    u_L_up_sd = nan(n,2,size(pathDB.x,2));
    v_L_up_sd = nan(n,2,size(pathDB.x,2));
    w_L_up_sd = nan(n,2,size(pathDB.x,2));
    U_L_up_sd = nan(n,2,size(pathDB.x,2));
    
    u_L_down_mean = nan(n,2,size(pathDB.x,2));
    v_L_down_mean = nan(n,2,size(pathDB.x,2));
    w_L_down_mean = nan(n,2,size(pathDB.x,2));
    U_L_down_mean = nan(n,2,size(pathDB.x,2));
    
    u_L_down_sd = nan(n,2,size(pathDB.x,2)); 
    v_L_down_sd = nan(n,2,size(pathDB.x,2));
    w_L_down_sd = nan(n,2,size(pathDB.x,2));
    U_L_down_sd = nan(n,2,size(pathDB.x,2));
    
    omegax_L_up_mean = nan(n,2,size(pathDB.x,2));
    omegay_L_up_mean = nan(n,2,size(pathDB.x,2));
    omegaz_L_up_mean = nan(n,2,size(pathDB.x,2));
    Omega_L_up_mean = nan(n,2,size(pathDB.x,2));
    
    omegax_L_up_sd = nan(n,2,size(pathDB.x,2));
    omegay_L_up_sd = nan(n,2,size(pathDB.x,2));
    omegaz_L_up_sd = nan(n,2,size(pathDB.x,2));
    Omega_L_up_sd = nan(n,2,size(pathDB.x,2));
    
    omegax_L_down_mean = nan(n,2,size(pathDB.x,2));
    omegay_L_down_mean = nan(n,2,size(pathDB.x,2));
    omegaz_L_down_mean = nan(n,2,size(pathDB.x,2));
    Omega_L_down_mean = nan(n,2,size(pathDB.x,2));
    
    omegax_L_down_sd = nan(n,2,size(pathDB.x,2));
    omegay_L_down_sd = nan(n,2,size(pathDB.x,2));
    omegaz_L_down_sd = nan(n,2,size(pathDB.x,2));
    Omega_L_down_sd = nan(n,2,size(pathDB.x,2));
    
    
    % Right wing:
    
    phi_R_up_mean = nan(n,2,size(pathDB.x,2));
    theta_R_up_mean = nan(n,2,size(pathDB.x,2));
    eta_R_up_mean = nan(n,2,size(pathDB.x,2));
    
    phi_R_up_sd = nan(n,2,size(pathDB.x,2));
    theta_R_up_sd = nan(n,2,size(pathDB.x,2));
    eta_R_up_sd = nan(n,2,size(pathDB.x,2));
    
    phi_R_down_mean = nan(n,2,size(pathDB.x,2));
    theta_R_down_mean = nan(n,2,size(pathDB.x,2));
    eta_R_down_mean = nan(n,2,size(pathDB.x,2));
    
    phi_R_down_sd = nan(n,2,size(pathDB.x,2));
    theta_R_down_sd = nan(n,2,size(pathDB.x,2));
    eta_R_down_sd = nan(n,2,size(pathDB.x,2));
    
    alfa_R_up_mean = nan(n,2,size(pathDB.x,2));
    beta_R_up_mean = nan(n,2,size(pathDB.x,2));
    
    alfa_R_up_sd = nan(n,2,size(pathDB.x,2));
    beta_R_up_sd = nan(n,2,size(pathDB.x,2));
    
    alfa_R_down_mean = nan(n,2,size(pathDB.x,2));
    beta_R_down_mean = nan(n,2,size(pathDB.x,2));
    
    alfa_R_down_sd = nan(n,2,size(pathDB.x,2));
    beta_R_down_sd = nan(n,2,size(pathDB.x,2));
    
    u_R_up_mean = nan(n,2,size(pathDB.x,2));
    v_R_up_mean = nan(n,2,size(pathDB.x,2));
    w_R_up_mean = nan(n,2,size(pathDB.x,2));
    U_R_up_mean = nan(n,2,size(pathDB.x,2));
    
    u_R_up_sd = nan(n,2,size(pathDB.x,2));
    v_R_up_sd = nan(n,2,size(pathDB.x,2));
    w_R_up_sd = nan(n,2,size(pathDB.x,2));
    U_R_up_sd = nan(n,2,size(pathDB.x,2));
    
    u_R_down_mean = nan(n,2,size(pathDB.x,2));
    v_R_down_mean = nan(n,2,size(pathDB.x,2));
    w_R_down_mean = nan(n,2,size(pathDB.x,2));
    U_R_down_mean = nan(n,2,size(pathDB.x,2));
    
    u_R_down_sd = nan(n,2,size(pathDB.x,2));
    v_R_down_sd = nan(n,2,size(pathDB.x,2));
    w_R_down_sd = nan(n,2,size(pathDB.x,2));
    U_R_down_sd = nan(n,2,size(pathDB.x,2));
    
    omegax_R_up_mean = nan(n,2,size(pathDB.x,2));
    omegay_R_up_mean = nan(n,2,size(pathDB.x,2));
    omegaz_R_up_mean = nan(n,2,size(pathDB.x,2));
    Omega_R_up_mean = nan(n,2,size(pathDB.x,2));
    
    omegax_R_up_sd = nan(n,2,size(pathDB.x,2));
    omegay_R_up_sd = nan(n,2,size(pathDB.x,2));
    omegaz_R_up_sd = nan(n,2,size(pathDB.x,2));
    Omega_R_up_sd = nan(n,2,size(pathDB.x,2));
    
    omegax_R_down_mean = nan(n,2,size(pathDB.x,2));
    omegay_R_down_mean = nan(n,2,size(pathDB.x,2));
    omegaz_R_down_mean = nan(n,2,size(pathDB.x,2));
    Omega_R_down_mean = nan(n,2,size(pathDB.x,2));
    
    omegax_R_down_sd = nan(n,2,size(pathDB.x,2));
    omegay_R_down_sd = nan(n,2,size(pathDB.x,2));
    omegaz_R_down_sd = nan(n,2,size(pathDB.x,2));
    Omega_R_down_sd = nan(n,2,size(pathDB.x,2));
    
    
    for i=1:size(pathDB.x,2)
        
        
    start_meas = find(isnan(pathDB.x(:,i))==0, 1 );
    stop_meas = find(isnan(pathDB.x(:,i))==0, 1 ,'last');
    
    % Every movie sequence wil start with a downstroke:
          
    % start and stop point for downstrokes, left wing
    start_down_L = pathDB.L_wingbeat_loc(1,2,i);
    stop_down_L = find(isnan(pathDB.L_wingbeat_loc(:,2,i))==0, 1, 'last' );
    
    
    % start and stop point for upstrokes, left wing
    if start_down_L < pathDB.L_wingbeat_loc(1,1,i)
        
        start_up_L = 0;
        stop_up_L = stop_down_L;
        
    elseif start_down_L > pathDB.L_wingbeat_loc(1,1,i)
        
        start_up_L = 1;
        stop_up_L = stop_down_L;
        
    end
    
    % start and stop point for downstrokes, right wing
    start_down_R = pathDB.R_wingbeat_loc(1,2,i);
    stop_down_R = find(isnan(pathDB.R_wingbeat_loc(:,2,i))==0, 1, 'last' );
    
    % start and stop point for upstrokes, right wing
    if start_down_R < pathDB.R_wingbeat_loc(1,1,i)
        
        start_up_R = 0;
        stop_up_R = stop_down_R;
        
    elseif start_down_R > pathDB.R_wingbeat_loc(1,1,i)
        
        start_up_R = 1;
        stop_up_R = stop_down_R;
        
    end

    
    N_L = stop_up_L-start_up_L;
    
    N_R = stop_up_R-start_up_R;

    
    
    if N_L ~= N_R
        
        'number of wingbeats left is not equal to the number of wingbeats right'
        
    end
    
    % Analyze the body kinematic parameters for sets of 4 wingbeats (in
    % order to easy implement the code it is assumed that the body follows
    % the left wingbeats, this introduces a temporal error, however in
    % general the deviation between the left and right wingbeats is 1/7500
    % of a second so the error involved will remain small):
    
    %N_seq = ((N_L-mod((N_L-2),2)-2)/2)-1;
    
    N_seq = N_L;
    

    
    for j = 1:1:N_seq-1

        
        a = start_meas+pathDB.L_wingbeat_loc(j,2,i)-1;
        b = start_meas+pathDB.L_wingbeat_loc(j+1,2,i)-2;

        % Body:
        
%         phi_body_mean(j,1,i) = Euler_mean(pathDB.b_roll(a:b,i));
%         theta_body_mean(j,1,i) = Euler_mean(pathDB.b_pitch(a:b,i));
%         xsi_body_mean(j,1,i) = Euler_mean(pathDB.b_yaw(a:b,i));
%         
%         phi_body_sd(j,1,i) = Euler_sd(pathDB.b_roll(a:b,i));
%         theta_body_sd(j,1,i) = Euler_sd(pathDB.b_pitch(a:b,i));
%         xsi_body_sd(j,1,i) = Euler_sd(pathDB.b_yaw(a:b,i));

        [phi_body_mean(j,1,i), theta_body_mean(j,1,i), xsi_body_mean(j,1,i)] = Euler_mean(pathDB.qb1_filt(a:b,i), pathDB.qb2_filt(a:b,i), pathDB.qb3_filt(a:b,i), pathDB.qb4_filt(a:b,i));
        
        phi_body_sd(j,1,i) = std(pathDB.b_roll(a:b,i));
        theta_body_sd(j,1,i) = std(pathDB.b_pitch(a:b,i));
        xsi_body_sd(j,1,i) = std(pathDB.b_yaw(a:b,i));
        
        phi_body_mean(j,2,i) = i*1000+j;
        theta_body_mean(j,2,i) = i*1000+j;
        xsi_body_mean(j,2,i) = i*1000+j;
        
        phi_body_sd(j,2,i) = i*1000+j;
        theta_body_sd(j,2,i) = i*1000+j;
        xsi_body_sd(j,2,i) = i*1000+j;
        
        omegax_body_mean(j,1,i) = mean(pathDB.b_omega1(a:b,i));
        omegay_body_mean(j,1,i) = mean(pathDB.b_omega2(a:b,i));
        omegaz_body_mean(j,1,i) = mean(pathDB.b_omega3(a:b,i));
        Omega_body_mean(j,1,i) = mean(sqrt(pathDB.b_omega1(a:b,i).^2+pathDB.b_omega2(a:b,i).^2+pathDB.b_omega3(a:b,i).^2));
    
        omegax_body_sd(j,1,i) = std(pathDB.b_omega1(a:b,i));
        omegay_body_sd(j,1,i) = std(pathDB.b_omega2(a:b,i));
        omegaz_body_sd(j,1,i) = std(pathDB.b_omega3(a:b,i));
        Omega_body_sd(j,1,i) = std(sqrt(pathDB.b_omega1(a:b,i).^2+pathDB.b_omega2(a:b,i).^2+pathDB.b_omega3(a:b,i).^2));
        
        omegax_body_mean(j,2,i) = i*1000+j;
        omegay_body_mean(j,2,i) = i*1000+j;
        omegaz_body_mean(j,2,i) = i*1000+j;
        Omega_body_mean(j,2,i) = i*1000+j;
    
        omegax_body_sd(j,2,i) = i*1000+j;
        omegay_body_sd(j,2,i) = i*1000+j;
        omegaz_body_sd(j,2,i) = i*1000+j;
        Omega_body_sd(j,2,i) = i*1000+j;
        
        alfa_body_mean(j,1,i) = mean(pathDB.b_alfa(a:b,i));
        beta_body_mean(j,1,i) = mean(pathDB.b_beta(a:b,i));
    
        alfa_body_sd(j,1,i) = std(pathDB.b_alfa(a:b,i));
        beta_body_sd(j,1,i) = std(pathDB.b_beta(a:b,i));
        
        alfa_body_mean(j,2,i) = i*1000+j;
        beta_body_mean(j,2,i) = i*1000+j;
    
        alfa_body_sd(j,2,i) = i*1000+j;
        beta_body_sd(j,2,i) = i*1000+j;
        
        u_body_mean(j,1,i) = mean(pathDB.u_body(a:b,i));
        v_body_mean(j,1,i) = mean(pathDB.v_body(a:b,i));
        w_body_mean(j,1,i) = mean(pathDB.w_body(a:b,i));
        U_body_mean(j,1,i) = mean(sqrt(pathDB.u_body(a:b,i).^2+pathDB.v_body(a:b,i).^2+pathDB.w_body(a:b,i).^2));
    
        u_body_sd(j,1,i) = std(pathDB.u_body(a:b,i));
        v_body_sd(j,1,i) = std(pathDB.v_body(a:b,i));
        w_body_sd(j,1,i) = std(pathDB.w_body(a:b,i));
        U_body_sd(j,1,i) = std(sqrt(pathDB.u_body(a:b,i).^2+pathDB.v_body(a:b,i).^2+pathDB.w_body(a:b,i).^2));
        
        u_body_mean(j,2,i) = i*1000+j;
        v_body_mean(j,2,i) = i*1000+j;
        w_body_mean(j,2,i) = i*1000+j;
        U_body_mean(j,2,i) = i*1000+j;
    
        u_body_sd(j,2,i) = i*1000+j;
        v_body_sd(j,2,i) = i*1000+j;
        w_body_sd(j,2,i) = i*1000+j;
        U_body_sd(j,2,i) = i*1000+j;
        
        ax_body_mean(j,1,i) = mean(pathDB.ax_body(a:b,i));
        ay_body_mean(j,1,i) = mean(pathDB.ay_body(a:b,i));
        az_body_mean(j,1,i) = mean(pathDB.az_body(a:b,i));
        a_body_mean(j,1,i) = mean(sqrt(pathDB.ax_body(a:b,i).^2+pathDB.ay_body(a:b,i).^2+pathDB.az_body(a:b,i).^2));
    
        ax_body_sd(j,1,i) = std(pathDB.ax_body(a:b,i));
        ay_body_sd(j,1,i) = std(pathDB.ay_body(a:b,i));
        az_body_sd(j,1,i) = std(pathDB.az_body(a:b,i));
        a_body_sd(j,1,i) = std(sqrt(pathDB.ax_body(a:b,i).^2+pathDB.ay_body(a:b,i).^2+pathDB.az_body(a:b,i).^2));
        
        ax_body_mean(j,2,i) = i*1000+j;
        ay_body_mean(j,2,i) = i*1000+j;
        az_body_mean(j,2,i) = i*1000+j;
        a_body_mean(j,2,i) = i*1000+j;
    
        ax_body_sd(j,2,i) = i*1000+j;
        ay_body_sd(j,2,i) = i*1000+j;
        az_body_sd(j,2,i) = i*1000+j;
        a_body_sd(j,2,i) = i*1000+j;
       
        clear a b
        
    end

    
    N_seq_L = N_L;
   
    
    
    for j = 1:1:N_seq_L-1
    
        a_up = start_meas+pathDB.L_wingbeat_loc(start_up_L+j,1,i)-1;
        b_up = start_meas+pathDB.L_wingbeat_loc(j+1,2,i)-2;

        
        
        a_down = start_meas+pathDB.L_wingbeat_loc(j,2,i)-1;
        b_down = start_meas+pathDB.L_wingbeat_loc(start_up_L+j,1,i)-2;
        
        
        % Cut off the pronation and suppination phase of the wingstrokes
        
        a_up2 = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)+1;
        b_up2 = start_meas+pathDB.R_wingbeat_loc(j+1,2,i)-4;
        
        a_down2 = start_meas+pathDB.R_wingbeat_loc(j,2,i)+1;
        b_down2 = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)-4;

    
        % Left wing:
        [phi_L_up_mean(j,1,i), theta_L_up_mean(j,1,i), eta_L_up_mean(j,1,i)] = Euler_mean(pathDB.qL1_filt2(a_up:b_up,i), pathDB.qL2_filt2(a_up:b_up,i), pathDB.qL3_filt2(a_up:b_up,i), pathDB.qL4_filt2(a_up:b_up,i));
        
        phi_L_up_sd(j,1,i) = std(pathDB.phi_L(a_up:b_up,i));
        theta_L_up_sd(j,1,i) = std(pathDB.theta_L(a_up:b_up,i));
        eta_L_up_sd(j,1,i) = std(pathDB.eta_L(a_up:b_up,i));
        
        phi_L_up_mean(j,2,i) = i*1000+j;
        theta_L_up_mean(j,2,i) = i*1000+j;
        eta_L_up_mean(j,2,i) = i*1000+j;

        phi_L_up_sd(j,2,i) = i*1000+j;
        theta_L_up_sd(j,2,i) = i*1000+j;
        eta_L_up_sd(j,2,i) = i*1000+j;
        
        [phi_L_down_mean(j,1,i), theta_L_down_mean(j,1,i), eta_L_down_mean(j,1,i)] = Euler_mean(pathDB.qL1_filt2(a_down:b_down,i), pathDB.qL2_filt2(a_down:b_down,i), pathDB.qL3_filt2(a_down:b_down,i), pathDB.qL4_filt2(a_down:b_down,i));

        phi_L_down_sd(j,1,i) = std(pathDB.phi_L(a_down:b_down,i));
        theta_L_down_sd(j,1,i) = std(pathDB.theta_L(a_down:b_down,i));
        eta_L_down_sd(j,1,i) = std(pathDB.eta_L(a_down:b_down,i));
        
        phi_L_down_mean(j,2,i) = i*1000+j;
        theta_L_down_mean(j,2,i) = i*1000+j;
        eta_L_down_mean(j,2,i) = i*1000+j;

        phi_L_down_sd(j,2,i) = i*1000+j;
        theta_L_down_sd(j,2,i) = i*1000+j;
        eta_L_down_sd(j,2,i) = i*1000+j;

%         alfa_L_up_mean(j,1,i) = mean(pathDB.alfa_L(a_up:b_up,2,i));
%         beta_L_up_mean(j,1,i) = mean(pathDB.beta_L(a_up:b_up,2,i));
%         
%         alfa_L_up_sd(j,1,i) = std(pathDB.alfa_L(a_up:b_up,2,i));
%         beta_L_up_sd(j,1,i) = std(pathDB.beta_L(a_up:b_up,2,i));
%         
%         alfa_L_up_mean(j,2,i) = i*1000+j;
%         beta_L_up_mean(j,2,i) = i*1000+j;
% 
%         alfa_L_up_sd(j,2,i) = i*1000+j;
%         beta_L_up_sd(j,2,i) = i*1000+j;
% 
%         alfa_L_down_mean(j,1,i) = mean(pathDB.alfa_L(a_down:b_down,2,i));
%         beta_L_down_mean(j,1,i) = mean(pathDB.beta_L(a_down:b_down,2,i));
% 
%         alfa_L_down_sd(j,1,i) = std(pathDB.alfa_L(a_down:b_down,2,i));
%         beta_L_down_sd(j,1,i) = std(pathDB.beta_L(a_down:b_down,2,i));
%         
%         alfa_L_down_mean(j,2,i) = i*1000+j;
%         beta_L_down_mean(j,2,i) = i*1000+j;
% 
%         alfa_L_down_sd(j,2,i) = i*1000+j;
%         beta_L_down_sd(j,2,i) = i*1000+j;
% 
%         u_L_up_mean(j,1,i) = mean(pathDB.u_wing_L(a_up:b_up,2,i));
%         v_L_up_mean(j,1,i) = mean(pathDB.v_wing_L(a_up:b_up,2,i));
%         w_L_up_mean(j,1,i) = mean(pathDB.w_wing_L(a_up:b_up,2,i));
%         U_L_up_mean(j,1,i) = mean(sqrt(pathDB.u_wing_L(a_up:b_up,2,i).^2 + ...
%                                        pathDB.v_wing_L(a_up:b_up,2,i).^2 + ...
%                                        pathDB.w_wing_L(a_up:b_up,2,i).^2));
% 
%         u_L_up_sd(j,1,i) = std(pathDB.u_wing_L(a_up:b_up,2,i));
%         v_L_up_sd(j,1,i) = std(pathDB.v_wing_L(a_up:b_up,2,i));
%         w_L_up_sd(j,1,i) = std(pathDB.w_wing_L(a_up:b_up,2,i));
%         U_L_up_sd(j,1,i) = std(sqrt(pathDB.u_wing_L(a_up:b_up,2,i).^2 + ...
%                                     pathDB.v_wing_L(a_up:b_up,2,i).^2 + ...
%                                     pathDB.w_wing_L(a_up:b_up,2,i).^2));
%         
%         u_L_up_mean(j,2,i) = i*1000+j;
%         v_L_up_mean(j,2,i) = i*1000+j;
%         w_L_up_mean(j,2,i) = i*1000+j;
%         U_L_up_mean(j,2,i) = i*1000+j;
% 
%         u_L_up_sd(j,2,i) = i*1000+j;
%         v_L_up_sd(j,2,i) = i*1000+j;
%         w_L_up_sd(j,2,i) = i*1000+j;
%         U_L_up_sd(j,2,i) = i*1000+j;
%         
%         u_L_down_mean(j,1,i) = mean(pathDB.u_wing_L(a_down:b_down,2,i));
%         v_L_down_mean(j,1,i) = mean(pathDB.v_wing_L(a_down:b_down,2,i));
%         w_L_down_mean(j,1,i) = mean(pathDB.w_wing_L(a_down:b_down,2,i));
%         U_L_down_mean(j,1,i) = mean(sqrt(pathDB.u_wing_L(a_down:b_down,2,i).^2 + ...
%                                          pathDB.v_wing_L(a_down:b_down,2,i).^2 + ...
%                                          pathDB.w_wing_L(a_down:b_down,2,i).^2));
% 
%         u_L_down_sd(j,1,i) = std(pathDB.u_wing_L(a_down:b_down,2,i));
%         v_L_down_sd(j,1,i) = std(pathDB.v_wing_L(a_down:b_down,2,i));
%         w_L_down_sd(j,1,i) = std(pathDB.w_wing_L(a_down:b_down,2,i));
%         U_L_down_sd(j,1,i) = std(sqrt( pathDB.u_wing_L(a_down:b_down,2,i).^2 + ...
%                                        pathDB.v_wing_L(a_down:b_down,2,i).^2 + ...
%                                        pathDB.w_wing_L(a_down:b_down,2,i).^2));
%         
%         u_L_down_mean(j,2,i) = i*1000+j;
%         v_L_down_mean(j,2,i) = i*1000+j;
%         w_L_down_mean(j,2,i) = i*1000+j;
%         U_L_down_mean(j,2,i) = i*1000+j;
% 
%         u_L_down_sd(j,2,i) = i*1000+j;
%         v_L_down_sd(j,2,i) = i*1000+j;
%         w_L_down_sd(j,2,i) = i*1000+j;
%         U_L_down_sd(j,2,i) = i*1000+j;
% 
%         omegax_L_up_mean(j,1,i) = mean(pathDB.omega1_L(a_up:b_up,i));
%         omegay_L_up_mean(j,1,i) = mean(pathDB.omega2_L(a_up:b_up,i));
%         omegaz_L_up_mean(j,1,i) = mean(pathDB.omega3_L(a_up:b_up,i));
%         Omega_L_up_mean(j,1,i) = mean(sqrt(pathDB.omega1_L(a_up:b_up,i).^2 + ...
%                                            pathDB.omega2_L(a_up:b_up,i).^2 + ...
%                                            pathDB.omega3_L(a_up:b_up,i).^2));
%                                         
%         omegax_L_up_sd(j,1,i) = std(pathDB.omega1_L(a_up:b_up,i));
%         omegay_L_up_sd(j,1,i) = std(pathDB.omega2_L(a_up:b_up,i));
%         omegaz_L_up_sd(j,1,i) = std(pathDB.omega3_L(a_up:b_up,i));
%         Omega_L_up_sd(j,1,i) = std(sqrt(pathDB.omega1_L(a_up:b_up,i).^2 + ...
%                                         pathDB.omega2_L(a_up:b_up,i).^2 + ...
%                                         pathDB.omega3_L(a_up:b_up,i).^2));

        alfa_L_up_mean(j,1,i) = mean(pathDB.alfa_L(a_up2:b_up2,2,i));
        beta_L_up_mean(j,1,i) = mean(pathDB.beta_L(a_up2:b_up2,2,i));
        
        alfa_L_up_sd(j,1,i) = std(pathDB.alfa_L(a_up2:b_up2,2,i));
        beta_L_up_sd(j,1,i) = std(pathDB.beta_L(a_up2:b_up2,2,i));
        
        alfa_L_up_mean(j,2,i) = i*1000+j;
        beta_L_up_mean(j,2,i) = i*1000+j;

        alfa_L_up_sd(j,2,i) = i*1000+j;
        beta_L_up_sd(j,2,i) = i*1000+j;

        alfa_L_down_mean(j,1,i) = mean(pathDB.alfa_L(a_down2:b_down2,2,i));
        beta_L_down_mean(j,1,i) = mean(pathDB.beta_L(a_down2:b_down2,2,i));

        alfa_L_down_sd(j,1,i) = std(pathDB.alfa_L(a_down2:b_down2,2,i));
        beta_L_down_sd(j,1,i) = std(pathDB.beta_L(a_down2:b_down2,2,i));
        
        alfa_L_down_mean(j,2,i) = i*1000+j;
        beta_L_down_mean(j,2,i) = i*1000+j;

        alfa_L_down_sd(j,2,i) = i*1000+j;
        beta_L_down_sd(j,2,i) = i*1000+j;

        u_L_up_mean(j,1,i) = mean(pathDB.u_wing_L(a_up2:b_up2,2,i));
        v_L_up_mean(j,1,i) = mean(pathDB.v_wing_L(a_up2:b_up2,2,i));
        w_L_up_mean(j,1,i) = mean(pathDB.w_wing_L(a_up2:b_up2,2,i));
        U_L_up_mean(j,1,i) = mean(sqrt(pathDB.u_wing_L(a_up2:b_up2,2,i).^2 + ...
                                       pathDB.v_wing_L(a_up2:b_up2,2,i).^2 + ...
                                       pathDB.w_wing_L(a_up2:b_up2,2,i).^2));

        u_L_up_sd(j,1,i) = std(pathDB.u_wing_L(a_up2:b_up2,2,i));
        v_L_up_sd(j,1,i) = std(pathDB.v_wing_L(a_up2:b_up2,2,i));
        w_L_up_sd(j,1,i) = std(pathDB.w_wing_L(a_up2:b_up2,2,i));
        U_L_up_sd(j,1,i) = std(sqrt(pathDB.u_wing_L(a_up2:b_up2,2,i).^2 + ...
                                    pathDB.v_wing_L(a_up2:b_up2,2,i).^2 + ...
                                    pathDB.w_wing_L(a_up2:b_up2,2,i).^2));
        
        u_L_up_mean(j,2,i) = i*1000+j;
        v_L_up_mean(j,2,i) = i*1000+j;
        w_L_up_mean(j,2,i) = i*1000+j;
        U_L_up_mean(j,2,i) = i*1000+j;

        u_L_up_sd(j,2,i) = i*1000+j;
        v_L_up_sd(j,2,i) = i*1000+j;
        w_L_up_sd(j,2,i) = i*1000+j;
        U_L_up_sd(j,2,i) = i*1000+j;
        
        u_L_down_mean(j,1,i) = mean(pathDB.u_wing_L(a_down2:b_down2,2,i));
        v_L_down_mean(j,1,i) = mean(pathDB.v_wing_L(a_down2:b_down2,2,i));
        w_L_down_mean(j,1,i) = mean(pathDB.w_wing_L(a_down2:b_down2,2,i));
        U_L_down_mean(j,1,i) = mean(sqrt(pathDB.u_wing_L(a_down2:b_down2,2,i).^2 + ...
                                         pathDB.v_wing_L(a_down2:b_down2,2,i).^2 + ...
                                         pathDB.w_wing_L(a_down2:b_down2,2,i).^2));

        u_L_down_sd(j,1,i) = std(pathDB.u_wing_L(a_down2:b_down2,2,i));
        v_L_down_sd(j,1,i) = std(pathDB.v_wing_L(a_down2:b_down2,2,i));
        w_L_down_sd(j,1,i) = std(pathDB.w_wing_L(a_down2:b_down2,2,i));
        U_L_down_sd(j,1,i) = std(sqrt( pathDB.u_wing_L(a_down2:b_down2,2,i).^2 + ...
                                       pathDB.v_wing_L(a_down2:b_down2,2,i).^2 + ...
                                       pathDB.w_wing_L(a_down2:b_down2,2,i).^2));
        
        u_L_down_mean(j,2,i) = i*1000+j;
        v_L_down_mean(j,2,i) = i*1000+j;
        w_L_down_mean(j,2,i) = i*1000+j;
        U_L_down_mean(j,2,i) = i*1000+j;

        u_L_down_sd(j,2,i) = i*1000+j;
        v_L_down_sd(j,2,i) = i*1000+j;
        w_L_down_sd(j,2,i) = i*1000+j;
        U_L_down_sd(j,2,i) = i*1000+j;

        omegax_L_up_mean(j,1,i) = mean(pathDB.omega1_L(a_up:b_up,i));
        omegay_L_up_mean(j,1,i) = mean(pathDB.omega2_L(a_up:b_up,i));
        omegaz_L_up_mean(j,1,i) = mean(pathDB.omega3_L(a_up:b_up,i));
        Omega_L_up_mean(j,1,i) = mean(sqrt(pathDB.omega1_L(a_up:b_up,i).^2 + ...
                                           pathDB.omega2_L(a_up:b_up,i).^2 + ...
                                           pathDB.omega3_L(a_up:b_up,i).^2));
                                        
        omegax_L_up_sd(j,1,i) = std(pathDB.omega1_L(a_up:b_up,i));
        omegay_L_up_sd(j,1,i) = std(pathDB.omega2_L(a_up:b_up,i));
        omegaz_L_up_sd(j,1,i) = std(pathDB.omega3_L(a_up:b_up,i));
        Omega_L_up_sd(j,1,i) = std(sqrt(pathDB.omega1_L(a_up:b_up,i).^2 + ...
                                        pathDB.omega2_L(a_up:b_up,i).^2 + ...
                                        pathDB.omega3_L(a_up:b_up,i).^2));
        
        omegax_L_up_mean(j,2,i) = i*1000+j;
        omegay_L_up_mean(j,2,i) = i*1000+j;
        omegaz_L_up_mean(j,2,i) = i*1000+j;
        Omega_L_up_mean(j,2,i) = i*1000+j;

        omegax_L_up_sd(j,2,i) = i*1000+j;
        omegay_L_up_sd(j,2,i) = i*1000+j;
        omegaz_L_up_sd(j,2,i) = i*1000+j;
        Omega_L_up_sd(j,2,i) = i*1000+j;
        
        omegax_L_down_mean(j,1,i) = mean(pathDB.omega1_L(a_down:b_down,i));
        omegay_L_down_mean(j,1,i) = mean(pathDB.omega2_L(a_down:b_down,i));
        omegaz_L_down_mean(j,1,i) = mean(pathDB.omega3_L(a_down:b_down,i));
        Omega_L_down_mean(j,1,i) = mean(sqrt(pathDB.omega1_L(a_down:b_down,i).^2 + ...
                                             pathDB.omega2_L(a_down:b_down,i).^2 + ...
                                             pathDB.omega3_L(a_down:b_down,i).^2));

        omegax_L_down_sd(j,1,i) = std(pathDB.omega1_L(a_down:b_down,i));
        omegay_L_down_sd(j,1,i) = std(pathDB.omega2_L(a_down:b_down,i));
        omegaz_L_down_sd(j,1,i) = std(pathDB.omega3_L(a_down:b_down,i));
        Omega_L_down_sd(j,1,i) = std(sqrt(pathDB.omega1_L(a_down:b_down,i).^2 + ...
                                          pathDB.omega2_L(a_down:b_down,i).^2 + ...
                                          pathDB.omega3_L(a_down:b_down,i).^2));
        
        omegax_L_down_mean(j,2,i) = i*1000+j;
        omegay_L_down_mean(j,2,i) = i*1000+j;
        omegaz_L_down_mean(j,2,i) = i*1000+j;
        Omega_L_down_mean(j,2,i) = i*1000+j;

        omegax_L_down_sd(j,2,i) = i*1000+j;
        omegay_L_down_sd(j,2,i) = i*1000+j;
        omegaz_L_down_sd(j,2,i) = i*1000+j;
        Omega_L_down_sd(j,2,i) = i*1000+j;
        
        
        
        clear a_up b_up a_down b_down a_up b_up2  a_down b_down2
    end

    N_seq_R = N_R;
    
    for j = 1:1:N_seq_R-1
        
        a_up = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)-1;
        b_up = start_meas+pathDB.R_wingbeat_loc(j+1,2,i)-2;

        
        
        a_down = start_meas+pathDB.R_wingbeat_loc(j,2,i)-1;
        b_down = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)-2;
        
        
        
        % Cut off the pronation and suppination phase of the wingstrokes
        
        a_up2 = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)+1;
        b_up2 = start_meas+pathDB.R_wingbeat_loc(j+1,2,i)-4;
        
        a_down2 = start_meas+pathDB.R_wingbeat_loc(j,2,i)+1;
        b_down2 = start_meas+pathDB.R_wingbeat_loc(start_up_R+j,1,i)-4;

    
        % Right wing:
        
        phi_R_up_mean(j,1,i) = mean(pathDB.phi_R(a_up:b_up,i));
        theta_R_up_mean(j,1,i) = mean(pathDB.theta_R(a_up:b_up,i));
        eta_R_up_mean(j,1,i) = mean(pathDB.eta_R(a_up:b_up,i));
        
        [phi_R_up_mean(j,1,i), theta_R_up_mean(j,1,i), eta_R_up_mean(j,1,i)] = Euler_mean(pathDB.qR1_filt2(a_up:b_up,i), pathDB.qR2_filt2(a_up:b_up,i), pathDB.qR3_filt2(a_up:b_up,i), pathDB.qR4_filt2(a_up:b_up,i));

        phi_R_up_sd(j,1,i) = std(pathDB.phi_R(a_up:b_up,i));
        theta_R_up_sd(j,1,i) = std(pathDB.theta_R(a_up:b_up,i));
        eta_R_up_sd(j,1,i) = std(pathDB.eta_R(a_up:b_up,i));
        
        phi_R_up_mean(j,2,i) = i*1000+j;
        theta_R_up_mean(j,2,i) = i*1000+j;
        eta_R_up_mean(j,2,i) = i*1000+j;

        phi_R_up_sd(j,2,i) = i*1000+j;
        theta_R_up_sd(j,2,i) = i*1000+j;
        eta_R_up_sd(j,2,i) = i*1000+j;
        
        [phi_R_down_mean(j,1,i), theta_R_down_mean(j,1,i), eta_R_down_mean(j,1,i)] = Euler_mean(pathDB.qR1_filt2(a_down:b_down,i), pathDB.qR2_filt2(a_down:b_down,i), pathDB.qR3_filt2(a_down:b_down,i), pathDB.qR4_filt2(a_down:b_down,i));

        phi_R_down_sd(j,1,i) = std(pathDB.phi_R(a_down:b_down,i));
        theta_R_down_sd(j,1,i) = std(pathDB.theta_R(a_down:b_down,i));
        eta_R_down_sd(j,1,i) = std(pathDB.eta_R(a_down:b_down,i));
        
        phi_R_down_mean(j,2,i) = i*1000+j;
        theta_R_down_mean(j,2,i) = i*1000+j;
        eta_R_down_mean(j,2,i) = i*1000+j;

        phi_R_down_sd(j,2,i) = i*1000+j;
        theta_R_down_sd(j,2,i) = i*1000+j;
        eta_R_down_sd(j,2,i) = i*1000+j;

%         alfa_R_up_mean(j,1,i) = mean(pathDB.alfa_R(a_up:b_up,2,i));
%         beta_R_up_mean(j,1,i) = mean(pathDB.beta_R(a_up:b_up,2,i));
% 
%         alfa_R_up_sd(j,1,i) = std(pathDB.alfa_R(a_up:b_up,2,i));
%         beta_R_up_sd(j,1,i) = std(pathDB.beta_R(a_up:b_up,2,i));
%         
%         alfa_R_up_mean(j,2,i) = i*1000+j;
%         beta_R_up_mean(j,2,i) = i*1000+j;
% 
%         alfa_R_up_sd(j,2,i) = i*1000+j;
%         beta_R_up_sd(j,2,i) = i*1000+j;
% 
%         alfa_R_down_mean(j,1,i) = mean(pathDB.alfa_R(a_down:b_down,2,i));
%         beta_R_down_mean(j,1,i) = mean(pathDB.beta_R(a_down:b_down,2,i));
% 
%         alfa_R_down_sd(j,1,i) = std(pathDB.alfa_R(a_down:b_down,2,i));
%         beta_R_down_sd(j,1,i) = std(pathDB.beta_R(a_down:b_down,2,i));
%         
%         alfa_R_down_mean(j,2,i) = i*1000+j;
%         beta_R_down_mean(j,2,i) = i*1000+j;
% 
%         alfa_R_down_sd(j,2,i) = i*1000+j;
%         beta_R_down_sd(j,2,i) = i*1000+j;
% 
%         u_R_up_mean(j,1,i) = mean(pathDB.u_wing_R(a_up:b_up,2,i));
%         v_R_up_mean(j,1,i) = mean(pathDB.v_wing_R(a_up:b_up,2,i));
%         w_R_up_mean(j,1,i) = mean(pathDB.w_wing_R(a_up:b_up,2,i));
%         U_R_up_mean(j,1,i) = mean(sqrt(pathDB.u_wing_R(a_up:b_up,2,i).^2 + ...
%                                        pathDB.v_wing_R(a_up:b_up,2,i).^2 + ...
%                                        pathDB.w_wing_R(a_up:b_up,2,i).^2));
% 
%         u_R_up_sd(j,1,i) = std(pathDB.u_wing_R(a_up:b_up,2,i));
%         v_R_up_sd(j,1,i) = std(pathDB.v_wing_R(a_up:b_up,2,i));
%         w_R_up_sd(j,1,i) = std(pathDB.w_wing_R(a_up:b_up,2,i));
%         U_R_up_sd(j,1,i) = std(sqrt(pathDB.u_wing_R(a_up:b_up,2,i).^2 + ...
%                                     pathDB.v_wing_R(a_up:b_up,2,i).^2 + ...
%                                     pathDB.w_wing_R(a_up:b_up,2,i).^2));
%         
%         u_R_up_mean(j,2,i) = i*1000+j;
%         v_R_up_mean(j,2,i) = i*1000+j;
%         w_R_up_mean(j,2,i) = i*1000+j;
%         U_R_up_mean(j,2,i) = i*1000+j;
% 
%         u_R_up_sd(j,2,i) = i*1000+j;
%         v_R_up_sd(j,2,i) = i*1000+j;
%         w_R_up_sd(j,2,i) = i*1000+j;
%         U_R_up_sd(j,2,i) = i*1000+j;
%         
%         u_R_down_mean(j,1,i) = mean(pathDB.u_wing_R(a_down:b_down,2,i));
%         v_R_down_mean(j,1,i) = mean(pathDB.v_wing_R(a_down:b_down,2,i));
%         w_R_down_mean(j,1,i) = mean(pathDB.w_wing_R(a_down:b_down,2,i));
%         U_R_down_mean(j,1,i) = mean(sqrt(pathDB.u_wing_R(a_down:b_down,2,i).^2 + ...
%                                          pathDB.v_wing_R(a_down:b_down,2,i).^2 + ...
%                                          pathDB.w_wing_R(a_down:b_down,2,i).^2));
% 
%         u_R_down_sd(j,1,i) = std(pathDB.u_wing_R(a_down:b_down,2,i));
%         v_R_down_sd(j,1,i) = std(pathDB.v_wing_R(a_down:b_down,2,i));
%         w_R_down_sd(j,1,i) = std(pathDB.w_wing_R(a_down:b_down,2,i));
%         U_R_down_sd(j,1,i) = std(sqrt(pathDB.u_wing_R(a_down:b_down,2,i).^2 + ...
%                                       pathDB.v_wing_R(a_down:b_down,2,i).^2 + ...
%                                       pathDB.w_wing_R(a_down:b_down,2,i).^2));
%         
%         u_R_down_mean(j,2,i) = i*1000+j;
%         v_R_down_mean(j,2,i) = i*1000+j;
%         w_R_down_mean(j,2,i) = i*1000+j;
%         U_R_down_mean(j,2,i) = i*1000+j;
% 
%         u_R_down_sd(j,2,i) = i*1000+j;
%         v_R_down_sd(j,2,i) = i*1000+j;
%         w_R_down_sd(j,2,i) = i*1000+j;
%         U_R_down_sd(j,2,i) = i*1000+j;

        alfa_R_up_mean(j,1,i) = mean(pathDB.alfa_R(a_up2:b_up2,2,i));
        beta_R_up_mean(j,1,i) = mean(pathDB.beta_R(a_up2:b_up2,2,i));

        alfa_R_up_sd(j,1,i) = std(pathDB.alfa_R(a_up2:b_up2,2,i));
        beta_R_up_sd(j,1,i) = std(pathDB.beta_R(a_up2:b_up2,2,i));
        
        alfa_R_up_mean(j,2,i) = i*1000+j;
        beta_R_up_mean(j,2,i) = i*1000+j;

        alfa_R_up_sd(j,2,i) = i*1000+j;
        beta_R_up_sd(j,2,i) = i*1000+j;

        alfa_R_down_mean(j,1,i) = mean(pathDB.alfa_R(a_down2:b_down2,2,i));
        beta_R_down_mean(j,1,i) = mean(pathDB.beta_R(a_down2:b_down2,2,i));

        alfa_R_down_sd(j,1,i) = std(pathDB.alfa_R(a_down2:b_down2,2,i));
        beta_R_down_sd(j,1,i) = std(pathDB.beta_R(a_down2:b_down2,2,i));
        
        alfa_R_down_mean(j,2,i) = i*1000+j;
        beta_R_down_mean(j,2,i) = i*1000+j;

        alfa_R_down_sd(j,2,i) = i*1000+j;
        beta_R_down_sd(j,2,i) = i*1000+j;

        u_R_up_mean(j,1,i) = mean(pathDB.u_wing_R(a_up2:b_up2,2,i));
        v_R_up_mean(j,1,i) = mean(pathDB.v_wing_R(a_up2:b_up2,2,i));
        w_R_up_mean(j,1,i) = mean(pathDB.w_wing_R(a_up2:b_up2,2,i));
        U_R_up_mean(j,1,i) = mean(sqrt(pathDB.u_wing_R(a_up2:b_up2,2,i).^2 + ...
                                       pathDB.v_wing_R(a_up2:b_up2,2,i).^2 + ...
                                       pathDB.w_wing_R(a_up2:b_up2,2,i).^2));

        u_R_up_sd(j,1,i) = std(pathDB.u_wing_R(a_up2:b_up2,2,i));
        v_R_up_sd(j,1,i) = std(pathDB.v_wing_R(a_up2:b_up2,2,i));
        w_R_up_sd(j,1,i) = std(pathDB.w_wing_R(a_up2:b_up2,2,i));
        U_R_up_sd(j,1,i) = std(sqrt(pathDB.u_wing_R(a_up2:b_up2,2,i).^2 + ...
                                    pathDB.v_wing_R(a_up2:b_up2,2,i).^2 + ...
                                    pathDB.w_wing_R(a_up2:b_up2,2,i).^2));
        
        u_R_up_mean(j,2,i) = i*1000+j;
        v_R_up_mean(j,2,i) = i*1000+j;
        w_R_up_mean(j,2,i) = i*1000+j;
        U_R_up_mean(j,2,i) = i*1000+j;

        u_R_up_sd(j,2,i) = i*1000+j;
        v_R_up_sd(j,2,i) = i*1000+j;
        w_R_up_sd(j,2,i) = i*1000+j;
        U_R_up_sd(j,2,i) = i*1000+j;
        
        u_R_down_mean(j,1,i) = mean(pathDB.u_wing_R(a_down2:b_down2,2,i));
        v_R_down_mean(j,1,i) = mean(pathDB.v_wing_R(a_down2:b_down2,2,i));
        w_R_down_mean(j,1,i) = mean(pathDB.w_wing_R(a_down2:b_down2,2,i));
        U_R_down_mean(j,1,i) = mean(sqrt(pathDB.u_wing_R(a_down2:b_down2,2,i).^2 + ...
                                         pathDB.v_wing_R(a_down2:b_down2,2,i).^2 + ...
                                         pathDB.w_wing_R(a_down2:b_down2,2,i).^2));

        u_R_down_sd(j,1,i) = std(pathDB.u_wing_R(a_down2:b_down2,2,i));
        v_R_down_sd(j,1,i) = std(pathDB.v_wing_R(a_down2:b_down2,2,i));
        w_R_down_sd(j,1,i) = std(pathDB.w_wing_R(a_down2:b_down2,2,i));
        U_R_down_sd(j,1,i) = std(sqrt(pathDB.u_wing_R(a_down2:b_down2,2,i).^2 + ...
                                      pathDB.v_wing_R(a_down2:b_down2,2,i).^2 + ...
                                      pathDB.w_wing_R(a_down2:b_down2,2,i).^2));
        
        u_R_down_mean(j,2,i) = i*1000+j;
        v_R_down_mean(j,2,i) = i*1000+j;
        w_R_down_mean(j,2,i) = i*1000+j;
        U_R_down_mean(j,2,i) = i*1000+j;

        u_R_down_sd(j,2,i) = i*1000+j;
        v_R_down_sd(j,2,i) = i*1000+j;
        w_R_down_sd(j,2,i) = i*1000+j;
        U_R_down_sd(j,2,i) = i*1000+j;

        omegax_R_up_mean(j,1,i) = mean(pathDB.omega1_R(a_up:b_up,i));
        omegay_R_up_mean(j,1,i) = mean(pathDB.omega2_R(a_up:b_up,i));
        omegaz_R_up_mean(j,1,i) = mean(pathDB.omega3_R(a_up:b_up,i));
        Omega_R_up_mean(j,1,i) = mean(sqrt(pathDB.omega1_R(a_up:b_up,i).^2 + ...
                                           pathDB.omega2_R(a_up:b_up,i).^2 + ...
                                           pathDB.omega3_R(a_up:b_up,i).^2));
                                        
        omegax_R_up_sd(j,1,i) = std(pathDB.omega1_R(a_up:b_up,i));
        omegay_R_up_sd(j,1,i) = std(pathDB.omega2_R(a_up:b_up,i));
        omegaz_R_up_sd(j,1,i) = std(pathDB.omega3_R(a_up:b_up,i));
        Omega_R_up_sd(j,1,i) = std(sqrt(pathDB.omega1_R(a_up:b_up,i).^2 + ...
                                        pathDB.omega2_R(a_up:b_up,i).^2 + ...
                                        pathDB.omega3_R(a_up:b_up,i).^2));
        
        omegax_R_up_mean(j,2,i) = i*1000+j;
        omegay_R_up_mean(j,2,i) = i*1000+j;
        omegaz_R_up_mean(j,2,i) = i*1000+j;
        Omega_R_up_mean(j,2,i) = i*1000+j;

        omegax_R_up_sd(j,2,i) = i*1000+j;
        omegay_R_up_sd(j,2,i) = i*1000+j;
        omegaz_R_up_sd(j,2,i) = i*1000+j;
        Omega_R_up_sd(j,2,i) = i*1000+j;
        
        omegax_R_down_mean(j,1,i) = mean(pathDB.omega1_R(a_down:b_down,i));
        omegay_R_down_mean(j,1,i) = mean(pathDB.omega2_R(a_down:b_down,i));
        omegaz_R_down_mean(j,1,i) = mean(pathDB.omega3_R(a_down:b_down,i));
        Omega_R_down_mean(j,1,i) = mean(sqrt(pathDB.omega1_R(a_down:b_down,i).^2 + ...
                                             pathDB.omega2_R(a_down:b_down,i).^2 + ...
                                             pathDB.omega3_R(a_down:b_down,i).^2));

        omegax_R_down_sd(j,1,i) = std(pathDB.omega1_R(a_down:b_down,i));
        omegay_R_down_sd(j,1,i) = std(pathDB.omega2_R(a_down:b_down,i));
        omegaz_R_down_sd(j,1,i) = std(pathDB.omega3_R(a_down:b_down,i));
        Omega_R_down_sd(j,1,i) = std(sqrt(pathDB.omega1_R(a_down:b_down,i).^2 + ...
                                          pathDB.omega2_R(a_down:b_down,i).^2 + ...
                                          pathDB.omega3_R(a_down:b_down,i).^2));
        
        omegax_R_down_mean(j,2,i) = i*1000+j;
        omegay_R_down_mean(j,2,i) = i*1000+j;
        omegaz_R_down_mean(j,2,i) = i*1000+j;
        Omega_R_down_mean(j,2,i) = i*1000+j;

        omegax_R_down_sd(j,2,i) = i*1000+j;
        omegay_R_down_sd(j,2,i) = i*1000+j;
        omegaz_R_down_sd(j,2,i) = i*1000+j;
        Omega_R_down_sd(j,2,i) = i*1000+j;
        
        
        
        clear a_up b_up  a_down b_down a_up b_up2  a_down b_down2
        
    end

    


    end



save(savefile,'phi_body_mean','theta_body_mean','xsi_body_mean','phi_body_sd','theta_body_sd','xsi_body_sd','omegax_body_mean','omegay_body_mean','omegaz_body_mean','Omega_body_mean','omegax_body_sd', ...
            'omegay_body_sd','omegaz_body_sd','Omega_body_sd','alfa_body_mean','beta_body_mean','alfa_body_sd','beta_body_sd','u_body_mean','v_body_mean','w_body_mean','U_body_mean','u_body_sd', ...
            'v_body_sd','w_body_sd','U_body_sd','ax_body_mean','ay_body_mean','az_body_mean','a_body_mean','ax_body_sd','ay_body_sd','az_body_sd','a_body_sd','phi_L_up_mean','theta_L_up_mean', ...
            'eta_L_up_mean','phi_L_up_sd','theta_L_up_sd','eta_L_up_sd','phi_L_down_mean','theta_L_down_mean','eta_L_down_mean','phi_L_down_sd','theta_L_down_sd','eta_L_down_sd','alfa_L_up_mean', ...
            'beta_L_up_mean','alfa_L_up_sd','beta_L_up_sd','alfa_L_down_mean','beta_L_down_mean','alfa_L_down_sd','beta_L_down_sd','u_L_up_mean','v_L_up_mean','w_L_up_mean','U_L_up_mean','u_L_up_sd', ...
            'v_L_up_sd','w_L_up_sd','U_L_up_sd','u_L_down_mean','v_L_down_mean','w_L_down_mean','U_L_down_mean','u_L_down_sd','v_L_down_sd','w_L_down_sd','U_L_down_sd','omegax_L_up_mean', ...
            'omegay_L_up_mean','omegaz_L_up_mean','Omega_L_up_mean','omegax_L_up_sd','omegay_L_up_sd','omegaz_L_up_sd','Omega_L_up_sd','omegax_L_down_mean','omegay_L_down_mean','omegaz_L_down_mean', ...
            'Omega_L_down_mean','omegax_L_down_sd','omegay_L_down_sd','omegaz_L_down_sd','Omega_L_down_sd','phi_R_up_mean','theta_R_up_mean','eta_R_up_mean','phi_R_up_sd','theta_R_up_sd','eta_R_up_sd', ...
            'phi_R_down_mean','theta_R_down_mean','eta_R_down_mean','phi_R_down_sd','theta_R_down_sd','eta_R_down_sd','alfa_R_up_mean','beta_R_up_mean','alfa_R_up_sd','beta_R_up_sd','alfa_R_down_mean', ...
            'beta_R_down_mean','alfa_R_down_sd','beta_R_down_sd','u_R_up_mean','v_R_up_mean','w_R_up_mean','U_R_up_mean','u_R_up_sd','v_R_up_sd','w_R_up_sd','U_R_up_sd','u_R_down_mean','v_R_down_mean', ...
            'w_R_down_mean','U_R_down_mean','u_R_down_sd','v_R_down_sd','w_R_down_sd','U_R_down_sd','omegax_R_up_mean','omegay_R_up_mean','omegaz_R_up_mean','Omega_R_up_mean','omegax_R_up_sd', ...
            'omegay_R_up_sd','omegaz_R_up_sd','Omega_R_up_sd','omegax_R_down_mean','omegay_R_down_mean','omegaz_R_down_mean','Omega_R_down_mean','omegax_R_down_sd','omegay_R_down_sd', ...
            'omegaz_R_down_sd','Omega_R_down_sd')

end



%     N_L = stop_up_L-start_up_L-1;
%     
%     N_R = stop_up_R-start_up_R-1;
% 
%     
%     
%     if N_L ~= N_R
%         
%         'number of wingbeats left is not equal to the number of wingbeats right'
%         
%     end
%     
%     % Analyze the body kinematic parameters for sets of 4 wingbeats (in
%     % order to easy implement the code it is assumed that the body follows
%     % the left wingbeats, this introduces a temporal error, however in
%     % general the deviation between the left and right wingbeats is 1/7500
%     % of a second so the error involved will remain small):
%     
%     %N_seq = ((N_L-mod((N_L-2),2)-2)/2)-1;
%     
%     N_seq = ((N_L-mod(N_L,4))/4)+((N_L-2-start_up_L-mod((N_L-2-start_up_L),4))/4);
%     
% 
%     
%     for j = 1:1:N_seq
% 
%         
%         a = start_meas+pathDB.L_wingbeat_loc(j*2-1,2,i)-1;
%         b = start_meas+pathDB.L_wingbeat_loc(j*2+3,2,i)-2;
% 
%         % Body:
%         
%         phi_body_mean(j,1,i) = mean(pathDB.b_roll(a:b,i));
%         theta_body_mean(j,1,i) = mean(pathDB.b_pitch(a:b,i));
%         xsi_body_mean(j,1,i) = mean(pathDB.b_yaw(a:b,i));
%         
%         phi_body_sd(j,1,i) = std(pathDB.b_roll(a:b,i));
%         theta_body_sd(j,1,i) = std(pathDB.b_pitch(a:b,i));
%         xsi_body_sd(j,1,i) = std(pathDB.b_yaw(a:b,i));
%         
%         phi_body_mean(j,2,i) = i*1000+j;
%         theta_body_mean(j,2,i) = i*1000+j;
%         xsi_body_mean(j,2,i) = i*1000+j;
%         
%         phi_body_sd(j,2,i) = i*1000+j;
%         theta_body_sd(j,2,i) = i*1000+j;
%         xsi_body_sd(j,2,i) = i*1000+j;
%         
%         omegax_body_mean(j,1,i) = mean(pathDB.b_omega1(a:b,i));
%         omegay_body_mean(j,1,i) = mean(pathDB.b_omega2(a:b,i));
%         omegaz_body_mean(j,1,i) = mean(pathDB.b_omega3(a:b,i));
%         Omega_body_mean(j,1,i) = mean(sqrt(pathDB.b_omega1(a:b,i).^2+pathDB.b_omega2(a:b,i).^2+pathDB.b_omega3(a:b,i).^2));
%     
%         omegax_body_sd(j,1,i) = std(pathDB.b_omega1(a:b,i));
%         omegay_body_sd(j,1,i) = std(pathDB.b_omega2(a:b,i));
%         omegaz_body_sd(j,1,i) = std(pathDB.b_omega3(a:b,i));
%         Omega_body_sd(j,1,i) = std(sqrt(pathDB.b_omega1(a:b,i).^2+pathDB.b_omega2(a:b,i).^2+pathDB.b_omega3(a:b,i).^2));
%         
%         omegax_body_mean(j,2,i) = i*1000+j;
%         omegay_body_mean(j,2,i) = i*1000+j;
%         omegaz_body_mean(j,2,i) = i*1000+j;
%         Omega_body_mean(j,2,i) = i*1000+j;
%     
%         omegax_body_sd(j,2,i) = i*1000+j;
%         omegay_body_sd(j,2,i) = i*1000+j;
%         omegaz_body_sd(j,2,i) = i*1000+j;
%         Omega_body_sd(j,2,i) = i*1000+j;
%         
%         alfa_body_mean(j,1,i) = mean(pathDB.b_alfa(a:b,i));
%         beta_body_mean(j,1,i) = mean(pathDB.b_beta(a:b,i));
%     
%         alfa_body_sd(j,1,i) = std(pathDB.b_alfa(a:b,i));
%         beta_body_sd(j,1,i) = std(pathDB.b_beta(a:b,i));
%         
%         alfa_body_mean(j,2,i) = i*1000+j;
%         beta_body_mean(j,2,i) = i*1000+j;
%     
%         alfa_body_sd(j,2,i) = i*1000+j;
%         beta_body_sd(j,2,i) = i*1000+j;
%         
%         u_body_mean(j,1,i) = mean(pathDB.u_body(a:b,i));
%         v_body_mean(j,1,i) = mean(pathDB.v_body(a:b,i));
%         w_body_mean(j,1,i) = mean(pathDB.w_body(a:b,i));
%         U_body_mean(j,1,i) = mean(sqrt(pathDB.u_body(a:b,i).^2+pathDB.v_body(a:b,i).^2+pathDB.w_body(a:b,i).^2));
%     
%         u_body_sd(j,1,i) = std(pathDB.u_body(a:b,i));
%         v_body_sd(j,1,i) = std(pathDB.v_body(a:b,i));
%         w_body_sd(j,1,i) = std(pathDB.w_body(a:b,i));
%         U_body_sd(j,1,i) = std(sqrt(pathDB.u_body(a:b,i).^2+pathDB.v_body(a:b,i).^2+pathDB.w_body(a:b,i).^2));
%         
%         u_body_mean(j,2,i) = i*1000+j;
%         v_body_mean(j,2,i) = i*1000+j;
%         w_body_mean(j,2,i) = i*1000+j;
%         U_body_mean(j,2,i) = i*1000+j;
%     
%         u_body_sd(j,2,i) = i*1000+j;
%         v_body_sd(j,2,i) = i*1000+j;
%         w_body_sd(j,2,i) = i*1000+j;
%         U_body_sd(j,2,i) = i*1000+j;
%         
%         ax_body_mean(j,1,i) = mean(pathDB.ax_body(a:b,i));
%         ay_body_mean(j,1,i) = mean(pathDB.ay_body(a:b,i));
%         az_body_mean(j,1,i) = mean(pathDB.az_body(a:b,i));
%         a_body_mean(j,1,i) = mean(sqrt(pathDB.ax_body(a:b,i).^2+pathDB.ay_body(a:b,i).^2+pathDB.az_body(a:b,i).^2));
%     
%         ax_body_sd(j,1,i) = std(pathDB.ax_body(a:b,i));
%         ay_body_sd(j,1,i) = std(pathDB.ay_body(a:b,i));
%         az_body_sd(j,1,i) = std(pathDB.az_body(a:b,i));
%         a_body_sd(j,1,i) = std(sqrt(pathDB.ax_body(a:b,i).^2+pathDB.ay_body(a:b,i).^2+pathDB.az_body(a:b,i).^2));
%         
%         ax_body_mean(j,2,i) = i*1000+j;
%         ay_body_mean(j,2,i) = i*1000+j;
%         az_body_mean(j,2,i) = i*1000+j;
%         a_body_mean(j,2,i) = i*1000+j;
%     
%         ax_body_sd(j,2,i) = i*1000+j;
%         ay_body_sd(j,2,i) = i*1000+j;
%         az_body_sd(j,2,i) = i*1000+j;
%         a_body_sd(j,2,i) = i*1000+j;
%        
%         clear a b
%         
%     end
% 
%     
%     N_seq_L = ((N_L-mod(N_L,4))/4)+((N_L-2-start_up_L-mod((N_L-2-start_up_L),4))/4);
%    
%     
%     
%     for j = 1:1:N_seq_L
%     
%         a_up = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2-1,1,i)-1;
%         b_up = start_meas+pathDB.L_wingbeat_loc(j*2,2,i)-2;
%         
%         c_up = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2,1,i)-1;
%         d_up = start_meas+pathDB.L_wingbeat_loc(j*2+1,2,i)-2;
%         
%         e_up = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2+1,1,i)-1;
%         f_up = start_meas+pathDB.L_wingbeat_loc(j*2+2,2,i)-2;
%         
%         g_up = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2+2,1,i)-1;
%         h_up = start_meas+pathDB.L_wingbeat_loc(j*2+3,2,i)-2;
%         
%         
%         a_down = start_meas+pathDB.L_wingbeat_loc(j*2-1,2,i)-1;
%         b_down = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2-1,1,i)-2;
%        
%         c_down = start_meas+pathDB.L_wingbeat_loc(j*2,2,i)-1;
%         d_down = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2,1,i)-2;
%         
%         e_down = start_meas+pathDB.L_wingbeat_loc(j*2+1,2,i)-1;
%         f_down = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2+1,1,i)-2;
%         
%         g_down = start_meas+pathDB.L_wingbeat_loc(j*2+2,2,i)-1;
%         h_down = start_meas+pathDB.L_wingbeat_loc(start_up_L+j*2+2,1,i)-2;
%     
%         % Left wing:
%         
%         phi_L_up_mean(j,1,i) = mean([pathDB.phi_L(a_up:b_up,i); pathDB.phi_L(c_up:d_up,i); pathDB.phi_L(e_up:f_up,i); pathDB.phi_L(g_up:h_up,i)]);
%         theta_L_up_mean(j,1,i) = mean([pathDB.theta_L(a_up:b_up,i); pathDB.theta_L(c_up:d_up,i); pathDB.theta_L(e_up:f_up,i); pathDB.theta_L(g_up:h_up,i)]);
%         eta_L_up_mean(j,1,i) = mean([pathDB.eta_L(a_up:b_up,i); pathDB.eta_L(c_up:d_up,i); pathDB.eta_L(e_up:f_up,i); pathDB.eta_L(g_up:h_up,i)]);
%         
%         phi_L_up_sd(j,1,i) = std([pathDB.phi_L(a_up:b_up,i); pathDB.phi_L(c_up:d_up,i); pathDB.phi_L(e_up:f_up,i); pathDB.phi_L(g_up:h_up,i)]);
%         theta_L_up_sd(j,1,i) = std([pathDB.theta_L(a_up:b_up,i); pathDB.theta_L(c_up:d_up,i); pathDB.theta_L(e_up:f_up,i); pathDB.theta_L(g_up:h_up,i)]);
%         eta_L_up_sd(j,1,i) = std([pathDB.eta_L(a_up:b_up,i); pathDB.eta_L(c_up:d_up,i); pathDB.eta_L(e_up:f_up,i); pathDB.eta_L(g_up:h_up,i)]);
%         
%         phi_L_up_mean(j,2,i) = i*1000+j;
%         theta_L_up_mean(j,2,i) = i*1000+j;
%         eta_L_up_mean(j,2,i) = i*1000+j;
% 
%         phi_L_up_sd(j,2,i) = i*1000+j;
%         theta_L_up_sd(j,2,i) = i*1000+j;
%         eta_L_up_sd(j,2,i) = i*1000+j;
% 
%         phi_L_down_mean(j,1,i) = mean([pathDB.phi_L(a_down:b_down,i); pathDB.phi_L(c_down:d_down,i); pathDB.phi_L(e_down:f_down,i); pathDB.phi_L(g_down:h_down,i)]);
%         theta_L_down_mean(j,1,i) = mean([pathDB.theta_L(a_down:b_down,i); pathDB.theta_L(c_down:d_down,i); pathDB.theta_L(e_down:f_down,i); pathDB.theta_L(g_down:h_down,i)]);
%         eta_L_down_mean(j,1,i) = mean([pathDB.eta_L(a_down:b_down,i); pathDB.eta_L(c_down:d_down,i); pathDB.eta_L(e_down:f_down,i); pathDB.eta_L(g_down:h_down,i)]);
% 
%         phi_L_down_sd(j,1,i) = std([pathDB.phi_L(a_down:b_down,i); pathDB.phi_L(c_down:d_down,i); pathDB.phi_L(e_down:f_down,i); pathDB.phi_L(g_down:h_down,i)]);
%         theta_L_down_sd(j,1,i) = std([pathDB.theta_L(a_down:b_down,i); pathDB.theta_L(c_down:d_down,i); pathDB.theta_L(e_down:f_down,i); pathDB.theta_L(g_down:h_down,i)]);
%         eta_L_down_sd(j,1,i) = std([pathDB.eta_L(a_down:b_down,i); pathDB.eta_L(c_down:d_down,i); pathDB.eta_L(e_down:f_down,i); pathDB.eta_L(g_down:h_down,i)]);
%         
%         phi_L_down_mean(j,2,i) = i*1000+j;
%         theta_L_down_mean(j,2,i) = i*1000+j;
%         eta_L_down_mean(j,2,i) = i*1000+j;
% 
%         phi_L_down_sd(j,2,i) = i*1000+j;
%         theta_L_down_sd(j,2,i) = i*1000+j;
%         eta_L_down_sd(j,2,i) = i*1000+j;
% 
%         alfa_L_up_mean(j,1,i) = mean([pathDB.alfa_L(a_up:b_up,2,i); pathDB.alfa_L(c_up:d_up,2,i); pathDB.alfa_L(e_up:f_up,2,i); pathDB.alfa_L(g_up:h_up,2,i)]);
%         beta_L_up_mean(j,1,i) = mean([pathDB.beta_L(a_up:b_up,2,i); pathDB.beta_L(c_up:d_up,2,i); pathDB.beta_L(e_up:f_up,2,i); pathDB.beta_L(g_up:h_up,2,i)]);
%         
%         alfa_L_up_sd(j,1,i) = std([pathDB.alfa_L(a_up:b_up,2,i); pathDB.alfa_L(c_up:d_up,2,i); pathDB.alfa_L(e_up:f_up,2,i); pathDB.alfa_L(g_up:h_up,2,i)]);
%         beta_L_up_sd(j,1,i) = std([pathDB.beta_L(a_up:b_up,2,i); pathDB.beta_L(c_up:d_up,2,i); pathDB.beta_L(e_up:f_up,2,i); pathDB.beta_L(g_up:h_up,2,i)]);
%         
%         alfa_L_up_mean(j,2,i) = i*1000+j;
%         beta_L_up_mean(j,2,i) = i*1000+j;
% 
%         alfa_L_up_sd(j,2,i) = i*1000+j;
%         beta_L_up_sd(j,2,i) = i*1000+j;
% 
%         alfa_L_down_mean(j,1,i) = mean([pathDB.alfa_L(a_down:b_down,2,i); pathDB.alfa_L(c_down:d_down,2,i); pathDB.alfa_L(e_down:f_down,2,i); pathDB.alfa_L(g_down:h_down,2,i)]);
%         beta_L_down_mean(j,1,i) = mean([pathDB.beta_L(a_down:b_down,2,i); pathDB.beta_L(c_down:d_down,2,i); pathDB.beta_L(e_down:f_down,2,i); pathDB.beta_L(g_down:h_down,2,i)]);
% 
%         alfa_L_down_sd(j,1,i) = std([pathDB.alfa_L(a_down:b_down,2,i); pathDB.alfa_L(c_down:d_down,2,i); pathDB.alfa_L(e_down:f_down,2,i); pathDB.alfa_L(g_down:h_down,2,i)]);
%         beta_L_down_sd(j,1,i) = std([pathDB.beta_L(a_down:b_down,2,i); pathDB.beta_L(c_down:d_down,2,i); pathDB.beta_L(e_down:f_down,2,i); pathDB.beta_L(g_down:h_down,2,i)]);
%         
%         alfa_L_down_mean(j,2,i) = i*1000+j;
%         beta_L_down_mean(j,2,i) = i*1000+j;
% 
%         alfa_L_down_sd(j,2,i) = i*1000+j;
%         beta_L_down_sd(j,2,i) = i*1000+j;
% 
%         u_L_up_mean(j,1,i) = mean([pathDB.u_wing_L(a_up:b_up,2,i); pathDB.u_wing_L(c_up:d_up,2,i); pathDB.u_wing_L(e_up:f_up,2,i); pathDB.u_wing_L(g_up:h_up,2,i)]);
%         v_L_up_mean(j,1,i) = mean([pathDB.v_wing_L(a_up:b_up,2,i); pathDB.v_wing_L(c_up:d_up,2,i); pathDB.v_wing_L(e_up:f_up,2,i); pathDB.v_wing_L(g_up:h_up,2,i)]);
%         w_L_up_mean(j,1,i) = mean([pathDB.w_wing_L(a_up:b_up,2,i); pathDB.w_wing_L(c_up:d_up,2,i); pathDB.w_wing_L(e_up:f_up,2,i); pathDB.w_wing_L(g_up:h_up,2,i)]);
%         U_L_up_mean(j,1,i) = mean(sqrt( [pathDB.u_wing_L(a_up:b_up,2,i); pathDB.u_wing_L(c_up:d_up,2,i); pathDB.u_wing_L(e_up:f_up,2,i); pathDB.u_wing_L(g_up:h_up,2,i)].^2 + ...
%                                         [pathDB.v_wing_L(a_up:b_up,2,i); pathDB.v_wing_L(c_up:d_up,2,i); pathDB.v_wing_L(e_up:f_up,2,i); pathDB.v_wing_L(g_up:h_up,2,i)].^2 + ...
%                                         [pathDB.w_wing_L(a_up:b_up,2,i); pathDB.w_wing_L(c_up:d_up,2,i); pathDB.w_wing_L(e_up:f_up,2,i); pathDB.w_wing_L(g_up:h_up,2,i)].^2));
% 
%         u_L_up_sd(j,1,i) = std([pathDB.u_wing_L(a_up:b_up,2,i); pathDB.u_wing_L(c_up:d_up,2,i); pathDB.u_wing_L(e_up:f_up,2,i); pathDB.u_wing_L(g_up:h_up,2,i)]);
%         v_L_up_sd(j,1,i) = std([pathDB.v_wing_L(a_up:b_up,2,i); pathDB.v_wing_L(c_up:d_up,2,i); pathDB.v_wing_L(e_up:f_up,2,i); pathDB.v_wing_L(g_up:h_up,2,i)]);
%         w_L_up_sd(j,1,i) = std([pathDB.w_wing_L(a_up:b_up,2,i); pathDB.w_wing_L(c_up:d_up,2,i); pathDB.w_wing_L(e_up:f_up,2,i); pathDB.w_wing_L(g_up:h_up,2,i)]);
%         U_L_up_sd(j,1,i) = std(sqrt( [pathDB.u_wing_L(a_up:b_up,2,i); pathDB.u_wing_L(c_up:d_up,2,i); pathDB.u_wing_L(e_up:f_up,2,i); pathDB.u_wing_L(g_up:h_up,2,i)].^2 + ...
%                                      [pathDB.v_wing_L(a_up:b_up,2,i); pathDB.v_wing_L(c_up:d_up,2,i); pathDB.v_wing_L(e_up:f_up,2,i); pathDB.v_wing_L(g_up:h_up,2,i)].^2 + ...
%                                      [pathDB.w_wing_L(a_up:b_up,2,i); pathDB.w_wing_L(c_up:d_up,2,i); pathDB.w_wing_L(e_up:f_up,2,i); pathDB.w_wing_L(g_up:h_up,2,i)].^2));
%         
%         u_L_up_mean(j,2,i) = i*1000+j;
%         v_L_up_mean(j,2,i) = i*1000+j;
%         w_L_up_mean(j,2,i) = i*1000+j;
%         U_L_up_mean(j,2,i) = i*1000+j;
% 
%         u_L_up_sd(j,2,i) = i*1000+j;
%         v_L_up_sd(j,2,i) = i*1000+j;
%         w_L_up_sd(j,2,i) = i*1000+j;
%         U_L_up_sd(j,2,i) = i*1000+j;
%         
%         u_L_down_mean(j,1,i) = mean([pathDB.u_wing_L(a_down:b_down,2,i); pathDB.u_wing_L(c_down:d_down,2,i); pathDB.u_wing_L(e_down:f_down,2,i); pathDB.u_wing_L(g_down:h_down,2,i)]);
%         v_L_down_mean(j,1,i) = mean([pathDB.v_wing_L(a_down:b_down,2,i); pathDB.v_wing_L(c_down:d_down,2,i); pathDB.v_wing_L(e_down:f_down,2,i); pathDB.v_wing_L(g_down:h_down,2,i)]);
%         w_L_down_mean(j,1,i) = mean([pathDB.w_wing_L(a_down:b_down,2,i); pathDB.w_wing_L(c_down:d_down,2,i); pathDB.w_wing_L(e_down:f_down,2,i); pathDB.w_wing_L(g_down:h_down,2,i)]);
%         U_L_down_mean(j,1,i) = mean(sqrt( [pathDB.u_wing_L(a_down:b_down,2,i); pathDB.u_wing_L(c_down:d_down,2,i); pathDB.u_wing_L(e_down:f_down,2,i); pathDB.u_wing_L(g_down:h_down,2,i)].^2 + ...
%                                           [pathDB.v_wing_L(a_down:b_down,2,i); pathDB.v_wing_L(c_down:d_down,2,i); pathDB.v_wing_L(e_down:f_down,2,i); pathDB.v_wing_L(g_down:h_down,2,i)].^2 + ...
%                                           [pathDB.w_wing_L(a_down:b_down,2,i); pathDB.w_wing_L(c_down:d_down,2,i); pathDB.w_wing_L(e_down:f_down,2,i); pathDB.w_wing_L(g_down:h_down,2,i)].^2));
% 
%         u_L_down_sd(j,1,i) = std([pathDB.u_wing_L(a_down:b_down,2,i); pathDB.u_wing_L(c_down:d_down,2,i); pathDB.u_wing_L(e_down:f_down,2,i); pathDB.u_wing_L(g_down:h_down,2,i)]);
%         v_L_down_sd(j,1,i) = std([pathDB.v_wing_L(a_down:b_down,2,i); pathDB.v_wing_L(c_down:d_down,2,i); pathDB.v_wing_L(e_down:f_down,2,i); pathDB.v_wing_L(g_down:h_down,2,i)]);
%         w_L_down_sd(j,1,i) = std([pathDB.w_wing_L(a_down:b_down,2,i); pathDB.w_wing_L(c_down:d_down,2,i); pathDB.w_wing_L(e_down:f_down,2,i); pathDB.w_wing_L(g_down:h_down,2,i)]);
%         U_L_down_sd(j,1,i) = std(sqrt( [pathDB.u_wing_L(a_down:b_down,2,i); pathDB.u_wing_L(c_down:d_down,2,i); pathDB.u_wing_L(e_down:f_down,2,i); pathDB.u_wing_L(g_down:h_down,2,i)].^2 + ...
%                                        [pathDB.v_wing_L(a_down:b_down,2,i); pathDB.v_wing_L(c_down:d_down,2,i); pathDB.v_wing_L(e_down:f_down,2,i); pathDB.v_wing_L(g_down:h_down,2,i)].^2 + ...
%                                        [pathDB.w_wing_L(a_down:b_down,2,i); pathDB.w_wing_L(c_down:d_down,2,i); pathDB.w_wing_L(e_down:f_down,2,i); pathDB.w_wing_L(g_down:h_down,2,i)].^2));
%         
%         u_L_down_mean(j,2,i) = i*1000+j;
%         v_L_down_mean(j,2,i) = i*1000+j;
%         w_L_down_mean(j,2,i) = i*1000+j;
%         U_L_down_mean(j,2,i) = i*1000+j;
% 
%         u_L_down_sd(j,2,i) = i*1000+j;
%         v_L_down_sd(j,2,i) = i*1000+j;
%         w_L_down_sd(j,2,i) = i*1000+j;
%         U_L_down_sd(j,2,i) = i*1000+j;
% 
%         omegax_L_up_mean(j,1,i) = mean([pathDB.omega1_L(a_up:b_up,i); pathDB.omega1_L(c_up:d_up,i); pathDB.omega1_L(e_up:f_up,i); pathDB.omega1_L(g_up:h_up,i)]);
%         omegay_L_up_mean(j,1,i) = mean([pathDB.omega2_L(a_up:b_up,i); pathDB.omega2_L(c_up:d_up,i); pathDB.omega2_L(e_up:f_up,i); pathDB.omega2_L(g_up:h_up,i)]);
%         omegaz_L_up_mean(j,1,i) = mean([pathDB.omega3_L(a_up:b_up,i); pathDB.omega3_L(c_up:d_up,i); pathDB.omega3_L(e_up:f_up,i); pathDB.omega3_L(g_up:h_up,i)]);
%         Omega_L_up_mean(j,1,i) = mean(sqrt( [pathDB.omega1_L(a_up:b_up,i); pathDB.omega1_L(c_up:d_up,i); pathDB.omega1_L(e_up:f_up,i); pathDB.omega1_L(g_up:h_up,i)].^2 + ...
%                                             [pathDB.omega2_L(a_up:b_up,i); pathDB.omega2_L(c_up:d_up,i); pathDB.omega2_L(e_up:f_up,i); pathDB.omega2_L(g_up:h_up,i)].^2 + ...
%                                             [pathDB.omega3_L(a_up:b_up,i); pathDB.omega3_L(c_up:d_up,i); pathDB.omega3_L(e_up:f_up,i); pathDB.omega3_L(g_up:h_up,i)].^2));
%                                         
%         omegax_L_up_sd(j,1,i) = std([pathDB.omega1_L(a_up:b_up,i); pathDB.omega1_L(c_up:d_up,i); pathDB.omega1_L(e_up:f_up,i); pathDB.omega1_L(g_up:h_up,i)]);
%         omegay_L_up_sd(j,1,i) = std([pathDB.omega2_L(a_up:b_up,i); pathDB.omega2_L(c_up:d_up,i); pathDB.omega2_L(e_up:f_up,i); pathDB.omega2_L(g_up:h_up,i)]);
%         omegaz_L_up_sd(j,1,i) = std([pathDB.omega3_L(a_up:b_up,i); pathDB.omega3_L(c_up:d_up,i); pathDB.omega3_L(e_up:f_up,i); pathDB.omega3_L(g_up:h_up,i)]);
%         Omega_L_up_sd(j,1,i) = std(sqrt([pathDB.omega1_L(a_up:b_up,i); pathDB.omega1_L(c_up:d_up,i); pathDB.omega1_L(e_up:f_up,i); pathDB.omega1_L(g_up:h_up,i)].^2 + ...
%                                         [pathDB.omega2_L(a_up:b_up,i); pathDB.omega2_L(c_up:d_up,i); pathDB.omega2_L(e_up:f_up,i); pathDB.omega2_L(g_up:h_up,i)].^2 + ...
%                                         [pathDB.omega3_L(a_up:b_up,i); pathDB.omega3_L(c_up:d_up,i); pathDB.omega3_L(e_up:f_up,i); pathDB.omega3_L(g_up:h_up,i)].^2));
%         
%         omegax_L_up_mean(j,2,i) = i*1000+j;
%         omegay_L_up_mean(j,2,i) = i*1000+j;
%         omegaz_L_up_mean(j,2,i) = i*1000+j;
%         Omega_L_up_mean(j,2,i) = i*1000+j;
% 
%         omegax_L_up_sd(j,2,i) = i*1000+j;
%         omegay_L_up_sd(j,2,i) = i*1000+j;
%         omegaz_L_up_sd(j,2,i) = i*1000+j;
%         Omega_L_up_sd(j,2,i) = i*1000+j;
%         
%         omegax_L_down_mean(j,1,i) = mean([pathDB.omega1_L(a_down:b_down,i); pathDB.omega1_L(c_down:d_down,i); pathDB.omega1_L(e_down:f_down,i); pathDB.omega1_L(g_down:h_down,i)]);
%         omegay_L_down_mean(j,1,i) = mean([pathDB.omega2_L(a_down:b_down,i); pathDB.omega2_L(c_down:d_down,i); pathDB.omega2_L(e_down:f_down,i); pathDB.omega2_L(g_down:h_down,i)]);
%         omegaz_L_down_mean(j,1,i) = mean([pathDB.omega3_L(a_down:b_down,i); pathDB.omega3_L(c_down:d_down,i); pathDB.omega3_L(e_down:f_down,i); pathDB.omega3_L(g_down:h_down,i)]);
%         Omega_L_down_mean(j,1,i) = mean(sqrt( [pathDB.omega1_L(a_down:b_down,i); pathDB.omega1_L(c_down:d_down,i); pathDB.omega1_L(e_down:f_down,i); pathDB.omega1_L(g_down:h_down,i)].^2 + ...
%                                               [pathDB.omega2_L(a_down:b_down,i); pathDB.omega2_L(c_down:d_down,i); pathDB.omega2_L(e_down:f_down,i); pathDB.omega2_L(g_down:h_down,i)].^2 + ...
%                                               [pathDB.omega3_L(a_down:b_down,i); pathDB.omega3_L(c_down:d_down,i); pathDB.omega3_L(e_down:f_down,i); pathDB.omega3_L(g_down:h_down,i)].^2));
% 
%         omegax_L_down_sd(j,1,i) = std([pathDB.omega1_L(a_down:b_down,i); pathDB.omega1_L(c_down:d_down,i); pathDB.omega1_L(e_down:f_down,i); pathDB.omega1_L(g_down:h_down,i)]);
%         omegay_L_down_sd(j,1,i) = std([pathDB.omega2_L(a_down:b_down,i); pathDB.omega2_L(c_down:d_down,i); pathDB.omega2_L(e_down:f_down,i); pathDB.omega2_L(g_down:h_down,i)]);
%         omegaz_L_down_sd(j,1,i) = std([pathDB.omega3_L(a_down:b_down,i); pathDB.omega3_L(c_down:d_down,i); pathDB.omega3_L(e_down:f_down,i); pathDB.omega3_L(g_down:h_down,i)]);
%         Omega_L_down_sd(j,1,i) = std(sqrt( [pathDB.omega1_L(a_down:b_down,i); pathDB.omega1_L(c_down:d_down,i); pathDB.omega1_L(e_down:f_down,i); pathDB.omega1_L(g_down:h_down,i)].^2 + ...
%                                            [pathDB.omega2_L(a_down:b_down,i); pathDB.omega2_L(c_down:d_down,i); pathDB.omega2_L(e_down:f_down,i); pathDB.omega2_L(g_down:h_down,i)].^2 + ...
%                                            [pathDB.omega3_L(a_down:b_down,i); pathDB.omega3_L(c_down:d_down,i); pathDB.omega3_L(e_down:f_down,i); pathDB.omega3_L(g_down:h_down,i)].^2));
%         
%         omegax_L_down_mean(j,2,i) = i*1000+j;
%         omegay_L_down_mean(j,2,i) = i*1000+j;
%         omegaz_L_down_mean(j,2,i) = i*1000+j;
%         Omega_L_down_mean(j,2,i) = i*1000+j;
% 
%         omegax_L_down_sd(j,2,i) = i*1000+j;
%         omegay_L_down_sd(j,2,i) = i*1000+j;
%         omegaz_L_down_sd(j,2,i) = i*1000+j;
%         Omega_L_down_sd(j,2,i) = i*1000+j;
%         
%         
%         
%         clear a_up b_up c_up d_up e_up f_up g_up h_up a_down b_down c_down d_down e_down f_down g_down h_down
%         
%     end
% 
%     N_seq_R = ((N_R-mod(N_R,4))/4)+((N_R-2-start_up_R-mod((N_R-2-start_up_R),4))/4);
%     
%     for j = 1:1:N_seq_R
%         
%         a_up = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2-1,1,i)-1;
%         b_up = start_meas+pathDB.R_wingbeat_loc(j*2,2,i)-2;
%         
%         c_up = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2,1,i)-1;
%         d_up = start_meas+pathDB.R_wingbeat_loc(j*2+1,2,i)-2;
%         
%         e_up = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2+1,1,i)-1;
%         f_up = start_meas+pathDB.R_wingbeat_loc(j*2+2,2,i)-2;
%         
%         g_up = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2+2,1,i)-1;
%         h_up = start_meas+pathDB.R_wingbeat_loc(j*2+3,2,i)-2;
%         
%         
%         a_down = start_meas+pathDB.R_wingbeat_loc(j*2-1,2,i)-1;
%         b_down = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2-1,1,i)-2;
%        
%         c_down = start_meas+pathDB.R_wingbeat_loc(j*2,2,i)-1;
%         d_down = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2,1,i)-2;
%         
%         e_down = start_meas+pathDB.R_wingbeat_loc(j*2+1,2,i)-1;
%         f_down = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2+1,1,i)-2;
%         
%         g_down = start_meas+pathDB.R_wingbeat_loc(j*2+2,2,i)-1;
%         h_down = start_meas+pathDB.R_wingbeat_loc(start_up_R+j*2+2,1,i)-2;
%     
%         % Left wing:
%         
%         phi_R_up_mean(j,1,i) = mean([pathDB.phi_R(a_up:b_up,i); pathDB.phi_R(c_up:d_up,i); pathDB.phi_R(e_up:f_up,i); pathDB.phi_R(g_up:h_up,i)]);
%         theta_R_up_mean(j,1,i) = mean([pathDB.theta_R(a_up:b_up,i); pathDB.theta_R(c_up:d_up,i); pathDB.theta_R(e_up:f_up,i); pathDB.theta_R(g_up:h_up,i)]);
%         eta_R_up_mean(j,1,i) = mean([pathDB.eta_R(a_up:b_up,i); pathDB.eta_R(c_up:d_up,i); pathDB.eta_R(e_up:f_up,i); pathDB.eta_R(g_up:h_up,i)]);
% 
%         phi_R_up_sd(j,1,i) = std([pathDB.phi_R(a_up:b_up,i); pathDB.phi_R(c_up:d_up,i); pathDB.phi_R(e_up:f_up,i); pathDB.phi_R(g_up:h_up,i)]);
%         theta_R_up_sd(j,1,i) = std([pathDB.theta_R(a_up:b_up,i); pathDB.theta_R(c_up:d_up,i); pathDB.theta_R(e_up:f_up,i); pathDB.theta_R(g_up:h_up,i)]);
%         eta_R_up_sd(j,1,i) = std([pathDB.eta_R(a_up:b_up,i); pathDB.eta_R(c_up:d_up,i); pathDB.eta_R(e_up:f_up,i); pathDB.eta_R(g_up:h_up,i)]);
%         
%         phi_R_up_mean(j,2,i) = i*1000+j;
%         theta_R_up_mean(j,2,i) = i*1000+j;
%         eta_R_up_mean(j,2,i) = i*1000+j;
% 
%         phi_R_up_sd(j,2,i) = i*1000+j;
%         theta_R_up_sd(j,2,i) = i*1000+j;
%         eta_R_up_sd(j,2,i) = i*1000+j;
% 
%         phi_R_down_mean(j,1,i) = mean([pathDB.phi_R(a_down:b_down,i); pathDB.phi_R(c_down:d_down,i); pathDB.phi_R(e_down:f_down,i); pathDB.phi_R(g_down:h_down,i)]);
%         theta_R_down_mean(j,1,i) = mean([pathDB.theta_R(a_down:b_down,i); pathDB.theta_R(c_down:d_down,i); pathDB.theta_R(e_down:f_down,i); pathDB.theta_R(g_down:h_down,i)]);
%         eta_R_down_mean(j,1,i) = mean([pathDB.eta_R(a_down:b_down,i); pathDB.eta_R(c_down:d_down,i); pathDB.eta_R(e_down:f_down,i); pathDB.eta_R(g_down:h_down,i)]);
% 
%         phi_R_down_sd(j,1,i) = std([pathDB.phi_R(a_down:b_down,i); pathDB.phi_R(c_down:d_down,i); pathDB.phi_R(e_down:f_down,i); pathDB.phi_R(g_down:h_down,i)]);
%         theta_R_down_sd(j,1,i) = std([pathDB.theta_R(a_down:b_down,i); pathDB.theta_R(c_down:d_down,i); pathDB.theta_R(e_down:f_down,i); pathDB.theta_R(g_down:h_down,i)]);
%         eta_R_down_sd(j,1,i) = std([pathDB.eta_R(a_down:b_down,i); pathDB.eta_R(c_down:d_down,i); pathDB.eta_R(e_down:f_down,i); pathDB.eta_R(g_down:h_down,i)]);
%         
%         phi_R_down_mean(j,2,i) = i*1000+j;
%         theta_R_down_mean(j,2,i) = i*1000+j;
%         eta_R_down_mean(j,2,i) = i*1000+j;
% 
%         phi_R_down_sd(j,2,i) = i*1000+j;
%         theta_R_down_sd(j,2,i) = i*1000+j;
%         eta_R_down_sd(j,2,i) = i*1000+j;
% 
%         alfa_R_up_mean(j,1,i) = mean([pathDB.alfa_R(a_up:b_up,2,i); pathDB.alfa_R(c_up:d_up,2,i); pathDB.alfa_R(e_up:f_up,2,i); pathDB.alfa_R(g_up:h_up,2,i)]);
%         beta_R_up_mean(j,1,i) = mean([pathDB.beta_R(a_up:b_up,2,i); pathDB.beta_R(c_up:d_up,2,i); pathDB.beta_R(e_up:f_up,2,i); pathDB.beta_R(g_up:h_up,2,i)]);
% 
%         alfa_R_up_sd(j,1,i) = std([pathDB.alfa_R(a_up:b_up,2,i); pathDB.alfa_R(c_up:d_up,2,i); pathDB.alfa_R(e_up:f_up,2,i); pathDB.alfa_R(g_up:h_up,2,i)]);
%         beta_R_up_sd(j,1,i) = std([pathDB.beta_R(a_up:b_up,2,i); pathDB.beta_R(c_up:d_up,2,i); pathDB.beta_R(e_up:f_up,2,i); pathDB.beta_R(g_up:h_up,2,i)]);
%         
%         alfa_R_up_mean(j,2,i) = i*1000+j;
%         beta_R_up_mean(j,2,i) = i*1000+j;
% 
%         alfa_R_up_sd(j,2,i) = i*1000+j;
%         beta_R_up_sd(j,2,i) = i*1000+j;
% 
%         alfa_R_down_mean(j,1,i) = mean([pathDB.alfa_R(a_down:b_down,2,i); pathDB.alfa_R(c_down:d_down,2,i); pathDB.alfa_R(e_down:f_down,2,i); pathDB.alfa_R(g_down:h_down,2,i)]);
%         beta_R_down_mean(j,1,i) = mean([pathDB.beta_R(a_down:b_down,2,i); pathDB.beta_R(c_down:d_down,2,i); pathDB.beta_R(e_down:f_down,2,i); pathDB.beta_R(g_down:h_down,2,i)]);
% 
%         alfa_R_down_sd(j,1,i) = std([pathDB.alfa_R(a_down:b_down,2,i); pathDB.alfa_R(c_down:d_down,2,i); pathDB.alfa_R(e_down:f_down,2,i); pathDB.alfa_R(g_down:h_down,2,i)]);
%         beta_R_down_sd(j,1,i) = std([pathDB.beta_R(a_down:b_down,2,i); pathDB.beta_R(c_down:d_down,2,i); pathDB.beta_R(e_down:f_down,2,i); pathDB.beta_R(g_down:h_down,2,i)]);
%         
%         alfa_R_down_mean(j,2,i) = i*1000+j;
%         beta_R_down_mean(j,2,i) = i*1000+j;
% 
%         alfa_R_down_sd(j,2,i) = i*1000+j;
%         beta_R_down_sd(j,2,i) = i*1000+j;
% 
%         u_R_up_mean(j,1,i) = mean([pathDB.u_wing_R(a_up:b_up,2,i); pathDB.u_wing_R(c_up:d_up,2,i); pathDB.u_wing_R(e_up:f_up,2,i); pathDB.u_wing_R(g_up:h_up,2,i)]);
%         v_R_up_mean(j,1,i) = mean([pathDB.v_wing_R(a_up:b_up,2,i); pathDB.v_wing_R(c_up:d_up,2,i); pathDB.v_wing_R(e_up:f_up,2,i); pathDB.v_wing_R(g_up:h_up,2,i)]);
%         w_R_up_mean(j,1,i) = mean([pathDB.w_wing_R(a_up:b_up,2,i); pathDB.w_wing_R(c_up:d_up,2,i); pathDB.w_wing_R(e_up:f_up,2,i); pathDB.w_wing_R(g_up:h_up,2,i)]);
%         U_R_up_mean(j,1,i) = mean(sqrt( [pathDB.u_wing_R(a_up:b_up,2,i); pathDB.u_wing_R(c_up:d_up,2,i); pathDB.u_wing_R(e_up:f_up,2,i); pathDB.u_wing_R(g_up:h_up,2,i)].^2 + ...
%                                         [pathDB.v_wing_R(a_up:b_up,2,i); pathDB.v_wing_R(c_up:d_up,2,i); pathDB.v_wing_R(e_up:f_up,2,i); pathDB.v_wing_R(g_up:h_up,2,i)].^2 + ...
%                                         [pathDB.w_wing_R(a_up:b_up,2,i); pathDB.w_wing_R(c_up:d_up,2,i); pathDB.w_wing_R(e_up:f_up,2,i); pathDB.w_wing_R(g_up:h_up,2,i)].^2));
% 
%         u_R_up_sd(j,1,i) = std([pathDB.u_wing_R(a_up:b_up,2,i); pathDB.u_wing_R(c_up:d_up,2,i); pathDB.u_wing_R(e_up:f_up,2,i); pathDB.u_wing_R(g_up:h_up,2,i)]);
%         v_R_up_sd(j,1,i) = std([pathDB.v_wing_R(a_up:b_up,2,i); pathDB.v_wing_R(c_up:d_up,2,i); pathDB.v_wing_R(e_up:f_up,2,i); pathDB.v_wing_R(g_up:h_up,2,i)]);
%         w_R_up_sd(j,1,i) = std([pathDB.w_wing_R(a_up:b_up,2,i); pathDB.w_wing_R(c_up:d_up,2,i); pathDB.w_wing_R(e_up:f_up,2,i); pathDB.w_wing_R(g_up:h_up,2,i)]);
%         U_R_up_sd(j,1,i) = std(sqrt( [pathDB.u_wing_R(a_up:b_up,2,i); pathDB.u_wing_R(c_up:d_up,2,i); pathDB.u_wing_R(e_up:f_up,2,i); pathDB.u_wing_R(g_up:h_up,2,i)].^2 + ...
%                                      [pathDB.v_wing_R(a_up:b_up,2,i); pathDB.v_wing_R(c_up:d_up,2,i); pathDB.v_wing_R(e_up:f_up,2,i); pathDB.v_wing_R(g_up:h_up,2,i)].^2 + ...
%                                      [pathDB.w_wing_R(a_up:b_up,2,i); pathDB.w_wing_R(c_up:d_up,2,i); pathDB.w_wing_R(e_up:f_up,2,i); pathDB.w_wing_R(g_up:h_up,2,i)].^2));
%         
%         u_R_up_mean(j,2,i) = i*1000+j;
%         v_R_up_mean(j,2,i) = i*1000+j;
%         w_R_up_mean(j,2,i) = i*1000+j;
%         U_R_up_mean(j,2,i) = i*1000+j;
% 
%         u_R_up_sd(j,2,i) = i*1000+j;
%         v_R_up_sd(j,2,i) = i*1000+j;
%         w_R_up_sd(j,2,i) = i*1000+j;
%         U_R_up_sd(j,2,i) = i*1000+j;
%         
%         u_R_down_mean(j,1,i) = mean([pathDB.u_wing_R(a_down:b_down,2,i); pathDB.u_wing_R(c_down:d_down,2,i); pathDB.u_wing_R(e_down:f_down,2,i); pathDB.u_wing_R(g_down:h_down,2,i)]);
%         v_R_down_mean(j,1,i) = mean([pathDB.v_wing_R(a_down:b_down,2,i); pathDB.v_wing_R(c_down:d_down,2,i); pathDB.v_wing_R(e_down:f_down,2,i); pathDB.v_wing_R(g_down:h_down,2,i)]);
%         w_R_down_mean(j,1,i) = mean([pathDB.w_wing_R(a_down:b_down,2,i); pathDB.w_wing_R(c_down:d_down,2,i); pathDB.w_wing_R(e_down:f_down,2,i); pathDB.w_wing_R(g_down:h_down,2,i)]);
%         U_R_down_mean(j,1,i) = mean(sqrt( [pathDB.u_wing_R(a_down:b_down,2,i); pathDB.u_wing_R(c_down:d_down,2,i); pathDB.u_wing_R(e_down:f_down,2,i); pathDB.u_wing_R(g_down:h_down,2,i)].^2 + ...
%                                           [pathDB.v_wing_R(a_down:b_down,2,i); pathDB.v_wing_R(c_down:d_down,2,i); pathDB.v_wing_R(e_down:f_down,2,i); pathDB.v_wing_R(g_down:h_down,2,i)].^2 + ...
%                                           [pathDB.w_wing_R(a_down:b_down,2,i); pathDB.w_wing_R(c_down:d_down,2,i); pathDB.w_wing_R(e_down:f_down,2,i); pathDB.w_wing_R(g_down:h_down,2,i)].^2));
% 
%         u_R_down_sd(j,1,i) = std([pathDB.u_wing_R(a_down:b_down,2,i); pathDB.u_wing_R(c_down:d_down,2,i); pathDB.u_wing_R(e_down:f_down,2,i); pathDB.u_wing_R(g_down:h_down,2,i)]);
%         v_R_down_sd(j,1,i) = std([pathDB.v_wing_R(a_down:b_down,2,i); pathDB.v_wing_R(c_down:d_down,2,i); pathDB.v_wing_R(e_down:f_down,2,i); pathDB.v_wing_R(g_down:h_down,2,i)]);
%         w_R_down_sd(j,1,i) = std([pathDB.w_wing_R(a_down:b_down,2,i); pathDB.w_wing_R(c_down:d_down,2,i); pathDB.w_wing_R(e_down:f_down,2,i); pathDB.w_wing_R(g_down:h_down,2,i)]);
%         U_R_down_sd(j,1,i) = std(sqrt( [pathDB.u_wing_R(a_down:b_down,2,i); pathDB.u_wing_R(c_down:d_down,2,i); pathDB.u_wing_R(e_down:f_down,2,i); pathDB.u_wing_R(g_down:h_down,2,i)].^2 + ...
%                                        [pathDB.v_wing_R(a_down:b_down,2,i); pathDB.v_wing_R(c_down:d_down,2,i); pathDB.v_wing_R(e_down:f_down,2,i); pathDB.v_wing_R(g_down:h_down,2,i)].^2 + ...
%                                        [pathDB.w_wing_R(a_down:b_down,2,i); pathDB.w_wing_R(c_down:d_down,2,i); pathDB.w_wing_R(e_down:f_down,2,i); pathDB.w_wing_R(g_down:h_down,2,i)].^2));
%         
%         u_R_down_mean(j,2,i) = i*1000+j;
%         v_R_down_mean(j,2,i) = i*1000+j;
%         w_R_down_mean(j,2,i) = i*1000+j;
%         U_R_down_mean(j,2,i) = i*1000+j;
% 
%         u_R_down_sd(j,2,i) = i*1000+j;
%         v_R_down_sd(j,2,i) = i*1000+j;
%         w_R_down_sd(j,2,i) = i*1000+j;
%         U_R_down_sd(j,2,i) = i*1000+j;
% 
%         omegax_R_up_mean(j,1,i) = mean([pathDB.omega1_R(a_up:b_up,i); pathDB.omega1_R(c_up:d_up,i); pathDB.omega1_R(e_up:f_up,i); pathDB.omega1_R(g_up:h_up,i)]);
%         omegay_R_up_mean(j,1,i) = mean([pathDB.omega2_R(a_up:b_up,i); pathDB.omega2_R(c_up:d_up,i); pathDB.omega2_R(e_up:f_up,i); pathDB.omega2_R(g_up:h_up,i)]);
%         omegaz_R_up_mean(j,1,i) = mean([pathDB.omega3_R(a_up:b_up,i); pathDB.omega3_R(c_up:d_up,i); pathDB.omega3_R(e_up:f_up,i); pathDB.omega3_R(g_up:h_up,i)]);
%         Omega_R_up_mean(j,1,i) = mean(sqrt( [pathDB.omega1_R(a_up:b_up,i); pathDB.omega1_R(c_up:d_up,i); pathDB.omega1_R(e_up:f_up,i); pathDB.omega1_R(g_up:h_up,i)].^2 + ...
%                                             [pathDB.omega2_R(a_up:b_up,i); pathDB.omega2_R(c_up:d_up,i); pathDB.omega2_R(e_up:f_up,i); pathDB.omega2_R(g_up:h_up,i)].^2 + ...
%                                             [pathDB.omega3_R(a_up:b_up,i); pathDB.omega3_R(c_up:d_up,i); pathDB.omega3_R(e_up:f_up,i); pathDB.omega3_R(g_up:h_up,i)].^2));
%                                         
%         omegax_R_up_sd(j,1,i) = std([pathDB.omega1_R(a_up:b_up,i); pathDB.omega1_R(c_up:d_up,i); pathDB.omega1_R(e_up:f_up,i); pathDB.omega1_R(g_up:h_up,i)]);
%         omegay_R_up_sd(j,1,i) = std([pathDB.omega2_R(a_up:b_up,i); pathDB.omega2_R(c_up:d_up,i); pathDB.omega2_R(e_up:f_up,i); pathDB.omega2_R(g_up:h_up,i)]);
%         omegaz_R_up_sd(j,1,i) = std([pathDB.omega3_R(a_up:b_up,i); pathDB.omega3_R(c_up:d_up,i); pathDB.omega3_R(e_up:f_up,i); pathDB.omega3_R(g_up:h_up,i)]);
%         Omega_R_up_sd(j,1,i) = std(sqrt([pathDB.omega1_R(a_up:b_up,i); pathDB.omega1_R(c_up:d_up,i); pathDB.omega1_R(e_up:f_up,i); pathDB.omega1_R(g_up:h_up,i)].^2 + ...
%                                         [pathDB.omega2_R(a_up:b_up,i); pathDB.omega2_R(c_up:d_up,i); pathDB.omega2_R(e_up:f_up,i); pathDB.omega2_R(g_up:h_up,i)].^2 + ...
%                                         [pathDB.omega3_R(a_up:b_up,i); pathDB.omega3_R(c_up:d_up,i); pathDB.omega3_R(e_up:f_up,i); pathDB.omega3_R(g_up:h_up,i)].^2));
%         
%         omegax_R_up_mean(j,2,i) = i*1000+j;
%         omegay_R_up_mean(j,2,i) = i*1000+j;
%         omegaz_R_up_mean(j,2,i) = i*1000+j;
%         Omega_R_up_mean(j,2,i) = i*1000+j;
% 
%         omegax_R_up_sd(j,2,i) = i*1000+j;
%         omegay_R_up_sd(j,2,i) = i*1000+j;
%         omegaz_R_up_sd(j,2,i) = i*1000+j;
%         Omega_R_up_sd(j,2,i) = i*1000+j;
%         
%         omegax_R_down_mean(j,1,i) = mean([pathDB.omega1_R(a_down:b_down,i); pathDB.omega1_R(c_down:d_down,i); pathDB.omega1_R(e_down:f_down,i); pathDB.omega1_R(g_down:h_down,i)]);
%         omegay_R_down_mean(j,1,i) = mean([pathDB.omega2_R(a_down:b_down,i); pathDB.omega2_R(c_down:d_down,i); pathDB.omega2_R(e_down:f_down,i); pathDB.omega2_R(g_down:h_down,i)]);
%         omegaz_R_down_mean(j,1,i) = mean([pathDB.omega3_R(a_down:b_down,i); pathDB.omega3_R(c_down:d_down,i); pathDB.omega3_R(e_down:f_down,i); pathDB.omega3_R(g_down:h_down,i)]);
%         Omega_R_down_mean(j,1,i) = mean(sqrt( [pathDB.omega1_R(a_down:b_down,i); pathDB.omega1_R(c_down:d_down,i); pathDB.omega1_R(e_down:f_down,i); pathDB.omega1_R(g_down:h_down,i)].^2 + ...
%                                               [pathDB.omega2_R(a_down:b_down,i); pathDB.omega2_R(c_down:d_down,i); pathDB.omega2_R(e_down:f_down,i); pathDB.omega2_R(g_down:h_down,i)].^2 + ...
%                                               [pathDB.omega3_R(a_down:b_down,i); pathDB.omega3_R(c_down:d_down,i); pathDB.omega3_R(e_down:f_down,i); pathDB.omega3_R(g_down:h_down,i)].^2));
% 
%         omegax_R_down_sd(j,1,i) = std([pathDB.omega1_R(a_down:b_down,i); pathDB.omega1_R(c_down:d_down,i); pathDB.omega1_R(e_down:f_down,i); pathDB.omega1_R(g_down:h_down,i)]);
%         omegay_R_down_sd(j,1,i) = std([pathDB.omega2_R(a_down:b_down,i); pathDB.omega2_R(c_down:d_down,i); pathDB.omega2_R(e_down:f_down,i); pathDB.omega2_R(g_down:h_down,i)]);
%         omegaz_R_down_sd(j,1,i) = std([pathDB.omega3_R(a_down:b_down,i); pathDB.omega3_R(c_down:d_down,i); pathDB.omega3_R(e_down:f_down,i); pathDB.omega3_R(g_down:h_down,i)]);
%         Omega_R_down_sd(j,1,i) = std(sqrt( [pathDB.omega1_R(a_down:b_down,i); pathDB.omega1_R(c_down:d_down,i); pathDB.omega1_R(e_down:f_down,i); pathDB.omega1_R(g_down:h_down,i)].^2 + ...
%                                            [pathDB.omega2_R(a_down:b_down,i); pathDB.omega2_R(c_down:d_down,i); pathDB.omega2_R(e_down:f_down,i); pathDB.omega2_R(g_down:h_down,i)].^2 + ...
%                                            [pathDB.omega3_R(a_down:b_down,i); pathDB.omega3_R(c_down:d_down,i); pathDB.omega3_R(e_down:f_down,i); pathDB.omega3_R(g_down:h_down,i)].^2));
%         
%         omegax_R_down_mean(j,2,i) = i*1000+j;
%         omegay_R_down_mean(j,2,i) = i*1000+j;
%         omegaz_R_down_mean(j,2,i) = i*1000+j;
%         Omega_R_down_mean(j,2,i) = i*1000+j;
% 
%         omegax_R_down_sd(j,2,i) = i*1000+j;
%         omegay_R_down_sd(j,2,i) = i*1000+j;
%         omegaz_R_down_sd(j,2,i) = i*1000+j;
%         Omega_R_down_sd(j,2,i) = i*1000+j;
%         
%         
%         
%         clear a_up b_up c_up d_up e_up f_up g_up h_up a_down b_down c_down d_down e_down f_down g_down h_down
%         
%     end
