function Robofly_man_tests( settings, pathDB )

   savefile = 'Robofly_maneuvering_tests.mat';

   nr_of_wb = 7; % number of wingbeats per sequence
    
   n_pol_theta = (length(pathDB.a_glob.theta1)-1);
   n_pol_eta = (length(pathDB.a_glob.eta1)-1);
   n_pol_phi = (length(pathDB.a_glob.phi1)-1);
   
   
   % Global average wingbeats---------------------------------------------
    
   
   n_wb_glob = nr_of_wb;
       
    wing_l_glob = mean(pathDB.wing_l);
    down_up_glob = pathDB.down_up_glob;
    f_glob = pathDB.f_glob;
    
   % FX forward wingbeats---------------------------------------------
    
%    FX_forward_range = [ 0:(0.02/19):0.02 0:(0.02/19):0.02 0:(0.02/19):0.02 0:(0.02/19):0.02  ];

    FX_forward_range = [ 0:(0.0145/19):0.0145 0:(0.0145/19):0.0145 0:(0.0145/19):0.0145 0:(0.0145/19):0.0145  ];
   
   n_test_FX_forward = length(FX_forward_range);
    
   n_wb_FX_forward = nr_of_wb;
    
   a_theta_FX_forward_L = nan((n_pol_theta+1)*2,n_wb_FX_forward,n_test_FX_forward);
   a_eta_FX_forward_L = nan((n_pol_eta+1)*2,n_wb_FX_forward,n_test_FX_forward);
   a_phi_FX_forward_L = nan((n_pol_phi+1)*2,n_wb_FX_forward,n_test_FX_forward);
    
    
   a_theta_FX_forward_R = nan((n_pol_theta+1)*2,n_wb_FX_forward,n_test_FX_forward);
   a_eta_FX_forward_R = nan((n_pol_eta+1)*2,n_wb_FX_forward,n_test_FX_forward);
   a_phi_FX_forward_R = nan((n_pol_phi+1)*2,n_wb_FX_forward,n_test_FX_forward);
    
    
    for i = 1:n_test_FX_forward
        
        a_dev_theta_L   = pathDB.control_glob.ax_theta_forward.*FX_forward_range(i);
        a_dev_eta_L     = pathDB.control_glob.ax_eta_forward.*FX_forward_range(i);
        a_dev_phi_L     = pathDB.control_glob.ax_phi_forward.*FX_forward_range(i);
        a_dev_theta_R   = pathDB.control_glob.ax_theta_forward.*FX_forward_range(i);
        a_dev_eta_R     = pathDB.control_glob.ax_eta_forward.*FX_forward_range(i);
        a_dev_phi_R     = pathDB.control_glob.ax_phi_forward.*FX_forward_range(i);
        
        if i > 0 && i <= (n_test_FX_forward/4)
        
            for j = 1:n_wb_FX_forward

                a_theta_FX_forward_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_FX_forward_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_FX_forward_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_FX_forward_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_FX_forward_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_FX_forward_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end
        
        elseif i > (n_test_FX_forward/4) && i <= (2*n_test_FX_forward/4)
            
            for j = 1:n_wb_FX_forward

                a_theta_FX_forward_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_FX_forward_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_FX_forward_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_FX_forward_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_FX_forward_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_FX_forward_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end
            
        elseif i > (2*n_test_FX_forward/4) && i <= (3*n_test_FX_forward/4)
            
            for j = 1:n_wb_FX_forward

                a_theta_FX_forward_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_FX_forward_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_FX_forward_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_FX_forward_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_FX_forward_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_FX_forward_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end
                
        elseif i > (3*n_test_FX_forward/4) && i <= (4*n_test_FX_forward/4)
            
            for j = 1:n_wb_FX_forward

                a_theta_FX_forward_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_FX_forward_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_FX_forward_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_FX_forward_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_FX_forward_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_FX_forward_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end
                
        end
        
    end
    
    wing_l_FX_forward = wing_l_glob;
    down_up_FX_forward = down_up_glob;
    f_FX_forward = f_glob;
    t_FX_forward = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    theta_FX_forward_L = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    eta_FX_forward_L = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    phi_FX_forward_L = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    theta_FX_forward_R = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    eta_FX_forward_R = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    phi_FX_forward_R = zeros(n_wb_FX_forward*100+1,n_test_FX_forward);
    
        
    for i = 1:n_test_FX_forward
 
        for j = 1:n_wb_FX_forward
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_FX_forward, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_FX_forward, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_FX_forward, 101, (j-1), j, 0 );

            t_FX_forward(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_FX_forward_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FX_forward_L(:,j,i);
            eta_FX_forward_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FX_forward_L(:,j,i);
            phi_FX_forward_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FX_forward_L(:,j,i);
            theta_FX_forward_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FX_forward_R(:,j,i);
            eta_FX_forward_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FX_forward_R(:,j,i);
            phi_FX_forward_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FX_forward_R(:,j,i);

        end
        
           
    end
    
    
    
   % FX back wingbeats---------------------------------------------
    
%    FX_back_range = [ 0:(-0.036/19):-0.036 0:(-0.036/19):-0.036 0:(-0.036/19):-0.036 0:(-0.036/19):-0.036 ];
 
   FX_back_range = [ 0:(-0.0265/19):-0.0265 0:(-0.0265/19):-0.0265 0:(-0.0265/19):-0.0265 0:(-0.0265/19):-0.0265 ];

   n_test_FX_back = length(FX_back_range);
    
   n_wb_FX_back = nr_of_wb;
    
   a_theta_FX_back_L = nan((n_pol_theta+1)*2,n_wb_FX_back,n_test_FX_back);
   a_eta_FX_back_L = nan((n_pol_eta+1)*2,n_wb_FX_back,n_test_FX_back);
   a_phi_FX_back_L = nan((n_pol_phi+1)*2,n_wb_FX_back,n_test_FX_back);
    
    
   a_theta_FX_back_R = nan((n_pol_theta+1)*2,n_wb_FX_back,n_test_FX_back);
   a_eta_FX_back_R = nan((n_pol_eta+1)*2,n_wb_FX_back,n_test_FX_back);
   a_phi_FX_back_R = nan((n_pol_phi+1)*2,n_wb_FX_back,n_test_FX_back);
    
    
    for i = 1:n_test_FX_back
        
        a_dev_theta_L   = pathDB.control_glob.ax_theta_back.*FX_back_range(i);
        a_dev_eta_L     = pathDB.control_glob.ax_eta_back.*FX_back_range(i);
        a_dev_phi_L     = pathDB.control_glob.ax_phi_back.*FX_back_range(i);
        a_dev_theta_R   = pathDB.control_glob.ax_theta_back.*FX_back_range(i);
        a_dev_eta_R     = pathDB.control_glob.ax_eta_back.*FX_back_range(i);
        a_dev_phi_R     = pathDB.control_glob.ax_phi_back.*FX_back_range(i);
        
        if i > 0 && i <= (n_test_FX_back/4)
        
            for j = 1:n_wb_FX_back
        
            a_theta_FX_back_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_FX_back_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_FX_back_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_FX_back_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_FX_back_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_FX_back_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
            end
        
        elseif i > (n_test_FX_back/4) && i <= (2*n_test_FX_back/4)
            
            for j = 1:n_wb_FX_back
        
            a_theta_FX_back_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_FX_back_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FX_back_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_FX_back_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_FX_back_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FX_back_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
            end
            
        elseif i > (2*n_test_FX_back/4) && i <= (3*n_test_FX_back/4)
            
            for j = 1:n_wb_FX_back
        
            a_theta_FX_back_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FX_back_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_FX_back_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_FX_back_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FX_back_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_FX_back_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
            end
                
        elseif i > (3*n_test_FX_back/4) && i <= (4*n_test_FX_back/4)
            
            for j = 1:n_wb_FX_back
        
            a_theta_FX_back_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FX_back_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FX_back_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_FX_back_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FX_back_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FX_back_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
            end
                
        end
        
    end
    
    wing_l_FX_back = wing_l_glob;
    down_up_FX_back = down_up_glob;
    f_FX_back = f_glob;
    t_FX_back = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    theta_FX_back_L = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    eta_FX_back_L = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    phi_FX_back_L = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    theta_FX_back_R = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    eta_FX_back_R = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    phi_FX_back_R = zeros(n_wb_FX_back*100+1,n_test_FX_back);
    
        
    for i = 1:n_test_FX_back
 
        for j = 1:n_wb_FX_back
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_FX_back, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_FX_back, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_FX_back, 101, (j-1), j, 0 );

            t_FX_back(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_FX_back_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FX_back_L(:,j,i);
            eta_FX_back_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FX_back_L(:,j,i);
            phi_FX_back_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FX_back_L(:,j,i);
            theta_FX_back_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FX_back_R(:,j,i);
            eta_FX_back_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FX_back_R(:,j,i);
            phi_FX_back_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FX_back_R(:,j,i);

        end
    
    end
    

   % FY back wingbeats---------------------------------------------
    
%    FY_range = [ 0:(0.028/19):0.028 0:(0.028/19):0.028 0:(0.028/19):0.028 0:(0.028/19):0.028 ];

    FY_range = [ 0:(0.0206/19):0.0206 0:(0.0206/19):0.0206 0:(0.0206/19):0.0206 0:(0.0206/19):0.0206 ];
   
   n_test_FY = length(FY_range);
    
   n_wb_FY = nr_of_wb;
    
   a_theta_FY_L = nan((n_pol_theta+1)*2,n_wb_FY,n_test_FY);
   a_eta_FY_L = nan((n_pol_eta+1)*2,n_wb_FY,n_test_FY);
   a_phi_FY_L = nan((n_pol_phi+1)*2,n_wb_FY,n_test_FY);
    
    
   a_theta_FY_R = nan((n_pol_theta+1)*2,n_wb_FY,n_test_FY);
   a_eta_FY_R = nan((n_pol_eta+1)*2,n_wb_FY,n_test_FY);
   a_phi_FY_R = nan((n_pol_phi+1)*2,n_wb_FY,n_test_FY);
    
    
    for i = 1:n_test_FY
        
       
        if i > 0 && i <= (n_test_FY/4)

                    a_dev_theta_L   = pathDB.control_glob.ay_theta_left.*FY_range(i);
                    a_dev_eta_L     = pathDB.control_glob.ay_eta_left.*FY_range(i);
                    a_dev_phi_L     = pathDB.control_glob.ay_phi_left.*FY_range(i);
                    a_dev_theta_R   = pathDB.control_glob.ay_theta_right.*FY_range(i);
                    a_dev_eta_R     = pathDB.control_glob.ay_eta_right.*FY_range(i);
                    a_dev_phi_R     = pathDB.control_glob.ay_phi_right.*FY_range(i);

                    for j = 1:n_wb_FY

                        a_theta_FY_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                        a_eta_FY_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                        a_phi_FY_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                        a_theta_FY_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                        a_eta_FY_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                        a_phi_FY_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

                    end
        
        elseif i > (n_test_FY/4) && i <= (2*n_test_FY/4)

                    a_dev_theta_L   = pathDB.control_glob.ay_theta_left.*FY_range(i);
                    a_dev_eta_L     = pathDB.control_glob.ay_eta_left.*FY_range(i);
                    a_dev_phi_L     = pathDB.control_glob.ay_phi_left.*FY_range(i);
                    a_dev_theta_R   = pathDB.control_glob.ay_theta_right.*FY_range(i);
                    a_dev_eta_R     = pathDB.control_glob.ay_eta_right.*FY_range(i);
                    a_dev_phi_R     = pathDB.control_glob.ay_phi_right.*FY_range(i);

                    for j = 1:n_wb_FY

                        a_theta_FY_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                        a_eta_FY_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                        a_phi_FY_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                        a_theta_FY_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                        a_eta_FY_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                        a_phi_FY_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

                    end         
            
        elseif i > (2*n_test_FY/4) && i <= (3*n_test_FY/4)
            
                    a_dev_theta_L   = pathDB.control_glob.ay_theta_left.*FY_range(i);
                    a_dev_eta_L     = pathDB.control_glob.ay_eta_left.*FY_range(i);
                    a_dev_phi_L     = pathDB.control_glob.ay_phi_left.*FY_range(i);
                    a_dev_theta_R   = pathDB.control_glob.ay_theta_right.*FY_range(i);
                    a_dev_eta_R     = pathDB.control_glob.ay_eta_right.*FY_range(i);
                    a_dev_phi_R     = pathDB.control_glob.ay_phi_right.*FY_range(i);

                    for j = 1:n_wb_FY

                        a_theta_FY_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                        a_eta_FY_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                        a_phi_FY_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                        a_theta_FY_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                        a_eta_FY_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                        a_phi_FY_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

                    end         
                
        elseif i > (3*n_test_FY/4) && i <= (4*n_test_FY/4)

                    a_dev_theta_L   = pathDB.control_glob.ay_theta_left.*FY_range(i);
                    a_dev_eta_L     = pathDB.control_glob.ay_eta_left.*FY_range(i);
                    a_dev_phi_L     = pathDB.control_glob.ay_phi_left.*FY_range(i);
                    a_dev_theta_R   = pathDB.control_glob.ay_theta_right.*FY_range(i);
                    a_dev_eta_R     = pathDB.control_glob.ay_eta_right.*FY_range(i);
                    a_dev_phi_R     = pathDB.control_glob.ay_phi_right.*FY_range(i);

                    for j = 1:n_wb_FY

                        a_theta_FY_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                        a_eta_FY_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                        a_phi_FY_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                        a_theta_FY_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                        a_eta_FY_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                        a_phi_FY_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

                    end           
                
        end
        
    end
    
    wing_l_FY = wing_l_glob;
    down_up_FY = down_up_glob;
    f_FY = f_glob;
    t_FY = zeros(n_wb_FY*100+1,n_test_FY);
    theta_FY_L = zeros(n_wb_FY*100+1,n_test_FY);
    eta_FY_L = zeros(n_wb_FY*100+1,n_test_FY);
    phi_FY_L = zeros(n_wb_FY*100+1,n_test_FY);
    theta_FY_R = zeros(n_wb_FY*100+1,n_test_FY);
    eta_FY_R = zeros(n_wb_FY*100+1,n_test_FY);
    phi_FY_R = zeros(n_wb_FY*100+1,n_test_FY);
    
        
    for i = 1:n_test_FY
 
        for j = 1:n_wb_FY
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_FY, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_FY, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_FY, 101, (j-1), j, 0 );

            t_FY(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_FY_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FY_L(:,j,i);
            eta_FY_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FY_L(:,j,i);
            phi_FY_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FY_L(:,j,i);
            theta_FY_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FY_R(:,j,i);
            eta_FY_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FY_R(:,j,i);
            phi_FY_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FY_R(:,j,i);

        end
    
    end
    
    
   % FZ up wingbeats---------------------------------------------
    
%    FZ_up_range = [ 0:(-0.066/19):-0.066 0:(-0.066/19):-0.066 0:(-0.066/19):-0.066 0:(-0.066/19):-0.066 ];

% this time 71% of the maximum force production

   FZ_up_range = [ 0:(-0.0313/19):-0.0313 0:(-0.0313/19):-0.0313 0:(-0.0313/19):-0.0313 0:(-0.0313/19):-0.0313 ];
   
   n_test_FZ_up = length(FZ_up_range);
    
   n_wb_FZ_up = nr_of_wb;
    
   a_theta_FZ_up_L = nan((n_pol_theta+1)*2,n_wb_FZ_up,n_test_FZ_up);
   a_eta_FZ_up_L = nan((n_pol_eta+1)*2,n_wb_FZ_up,n_test_FZ_up);
   a_phi_FZ_up_L = nan((n_pol_phi+1)*2,n_wb_FZ_up,n_test_FZ_up);
    
    
   a_theta_FZ_up_R = nan((n_pol_theta+1)*2,n_wb_FZ_up,n_test_FZ_up);
   a_eta_FZ_up_R = nan((n_pol_eta+1)*2,n_wb_FZ_up,n_test_FZ_up);
   a_phi_FZ_up_R = nan((n_pol_phi+1)*2,n_wb_FZ_up,n_test_FZ_up);
    
    
    for i = 1:n_test_FZ_up
        
        a_dev_theta_L   = pathDB.control_glob.az_theta_up.*FZ_up_range(i);
        a_dev_eta_L     = pathDB.control_glob.az_eta_up.*FZ_up_range(i);
        a_dev_phi_L     = pathDB.control_glob.az_phi_up.*FZ_up_range(i);
        a_dev_theta_R   = pathDB.control_glob.az_theta_up.*FZ_up_range(i);
        a_dev_eta_R     = pathDB.control_glob.az_eta_up.*FZ_up_range(i);
        a_dev_phi_R     = pathDB.control_glob.az_phi_up.*FZ_up_range(i);
        
        
        if i > 0 && i <= (n_test_FZ_up/4)
        
        for j = 1:n_wb_FZ_up
        
            a_theta_FZ_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_FZ_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_FZ_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_FZ_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_FZ_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_FZ_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
        end            
        
        elseif i > (n_test_FZ_up/4) && i <= (2*n_test_FZ_up/4)
            
        for j = 1:n_wb_FZ_up
        
            a_theta_FZ_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_FZ_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FZ_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_FZ_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_FZ_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FZ_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
        end            
            
        elseif i > (2*n_test_FZ_up/4) && i <= (3*n_test_FZ_up/4)
            
         for j = 1:n_wb_FZ_up
        
            a_theta_FZ_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FZ_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_FZ_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_FZ_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FZ_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_FZ_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
        end           
                
        elseif i > (3*n_test_FZ_up/4) && i <= (4*n_test_FZ_up/4)
            
         for j = 1:n_wb_FZ_up
        
            a_theta_FZ_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FZ_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FZ_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_FZ_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_FZ_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_FZ_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
        end           
                
        end
        
    end
    
    wing_l_FZ_up = wing_l_glob;
    down_up_FZ_up = down_up_glob;
    f_FZ_up = f_glob;
    t_FZ_up = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    theta_FZ_up_L = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    eta_FZ_up_L = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    phi_FZ_up_L = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    theta_FZ_up_R = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    eta_FZ_up_R = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    phi_FZ_up_R = zeros(n_wb_FZ_up*100+1,n_test_FZ_up);
    
        
    for i = 1:n_test_FZ_up
 
        for j = 1:n_wb_FZ_up
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_FZ_up, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_FZ_up, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_FZ_up, 101, (j-1), j, 0 );

            t_FZ_up(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_FZ_up_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FZ_up_L(:,j,i);
            eta_FZ_up_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FZ_up_L(:,j,i);
            phi_FZ_up_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FZ_up_L(:,j,i);
            theta_FZ_up_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_FZ_up_R(:,j,i);
            eta_FZ_up_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_FZ_up_R(:,j,i);
            phi_FZ_up_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_FZ_up_R(:,j,i);

        end
    
    end
    
    
   % MX wingbeats---------------------------------------------
    
   % Robofly can only handle deviation angles up to 30 degrees, to make
   % sure the robofly still runs we set the maximum roll moment to 1.2
   % times the maximum measured treshold instead of 1.5
   
%    MX_range = [ 0:(0.015/19):0.015 0:(0.015/19):0.015 0:(0.015/19):0.015 0:(0.015/19):0.015 ];

    MX_range = [ 0:(0.0142/19):0.0142 0:(0.0142/19):0.0142 0:(0.0142/19):0.0142 0:(0.0142/19):0.0142 ];
   
   n_test_MX = length(MX_range);
    
   n_wb_MX = nr_of_wb;
    
   a_theta_MX_L = nan((n_pol_theta+1)*2,n_wb_MX,n_test_MX);
   a_eta_MX_L = nan((n_pol_eta+1)*2,n_wb_MX,n_test_MX);
   a_phi_MX_L = nan((n_pol_phi+1)*2,n_wb_MX,n_test_MX);
    
    
   a_theta_MX_R = nan((n_pol_theta+1)*2,n_wb_MX,n_test_MX);
   a_eta_MX_R = nan((n_pol_eta+1)*2,n_wb_MX,n_test_MX);
   a_phi_MX_R = nan((n_pol_phi+1)*2,n_wb_MX,n_test_MX);
    
    
    for i = 1:n_test_MX

        if i <= (n_test_MX/4)
        
            a_dev_theta_L   = pathDB.control_glob.roll_theta_left.*MX_range(i);
            a_dev_eta_L     = pathDB.control_glob.roll_eta_left.*MX_range(i);
            a_dev_phi_L     = pathDB.control_glob.roll_phi_left.*MX_range(i);
            a_dev_theta_R   = pathDB.control_glob.roll_theta_right.*MX_range(i);
            a_dev_eta_R     = pathDB.control_glob.roll_eta_right.*MX_range(i);
            a_dev_phi_R     = pathDB.control_glob.roll_phi_right.*MX_range(i);

            for j = 1:n_wb_MX

                a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end       
            
        
        elseif i > (n_test_MX/4) && i <= (2*n_test_MX/4)
        
            a_dev_theta_L   = pathDB.control_glob.roll_theta_left.*MX_range(i);
            a_dev_eta_L     = pathDB.control_glob.roll_eta_left.*MX_range(i);
            a_dev_phi_L     = pathDB.control_glob.roll_phi_left.*MX_range(i);
            a_dev_theta_R   = pathDB.control_glob.roll_theta_right.*MX_range(i);
            a_dev_eta_R     = pathDB.control_glob.roll_eta_right.*MX_range(i);
            a_dev_phi_R     = pathDB.control_glob.roll_phi_right.*MX_range(i);

            for j = 1:n_wb_MX

                a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end        
            
        elseif i > (2*n_test_MX/4) && i <= (3*n_test_MX/4)
        
            a_dev_theta_L   = pathDB.control_glob.roll_theta_left.*MX_range(i);
            a_dev_eta_L     = pathDB.control_glob.roll_eta_left.*MX_range(i);
            a_dev_phi_L     = pathDB.control_glob.roll_phi_left.*MX_range(i);
            a_dev_theta_R   = pathDB.control_glob.roll_theta_right.*MX_range(i);
            a_dev_eta_R     = pathDB.control_glob.roll_eta_right.*MX_range(i);
            a_dev_phi_R     = pathDB.control_glob.roll_phi_right.*MX_range(i);

            for j = 1:n_wb_MX

                a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end         
                
        elseif i > (3*n_test_MX/4) && i <= (4*n_test_MX/4)
        
            a_dev_theta_L   = pathDB.control_glob.roll_theta_left.*MX_range(i);
            a_dev_eta_L     = pathDB.control_glob.roll_eta_left.*MX_range(i);
            a_dev_phi_L     = pathDB.control_glob.roll_phi_left.*MX_range(i);
            a_dev_theta_R   = pathDB.control_glob.roll_theta_right.*MX_range(i);
            a_dev_eta_R     = pathDB.control_glob.roll_eta_right.*MX_range(i);
            a_dev_phi_R     = pathDB.control_glob.roll_phi_right.*MX_range(i);

            for j = 1:n_wb_MX

                a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end         
                
        end
        
    end
    
    wing_l_MX = wing_l_glob;
    down_up_MX = down_up_glob;
    f_MX = f_glob;
    t_MX = zeros(n_wb_MX*100+1,n_test_MX);
    theta_MX_L = zeros(n_wb_MX*100+1,n_test_MX);
    eta_MX_L = zeros(n_wb_MX*100+1,n_test_MX);
    phi_MX_L = zeros(n_wb_MX*100+1,n_test_MX);
    theta_MX_R = zeros(n_wb_MX*100+1,n_test_MX);
    eta_MX_R = zeros(n_wb_MX*100+1,n_test_MX);
    phi_MX_R = zeros(n_wb_MX*100+1,n_test_MX);
    
        
    for i = 1:n_test_MX
 
        for j = 1:n_wb_MX
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_MX, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_MX, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_MX, 101, (j-1), j, 0 );

            t_MX(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_MX_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MX_L(:,j,i);
            eta_MX_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MX_L(:,j,i);
            phi_MX_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MX_L(:,j,i);
            theta_MX_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MX_R(:,j,i);
            eta_MX_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MX_R(:,j,i);
            phi_MX_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MX_R(:,j,i);

        end
    
    end
    
    
    
   % MY back wingbeats---------------------------------------------
    
%    MY_up_range = [ 0:(0.014/19):0.014 0:(0.014/19):0.014 0:(0.014/19):0.014 0:(0.014/19):0.014 ];

    MY_up_range = [ 0:(0.0101/19):0.0101 0:(0.0101/19):0.0101 0:(0.0101/19):0.0101 0:(0.0101/19):0.0101 ];
   
   n_test_MY_up = length(MY_up_range);
    
   n_wb_MY_up = nr_of_wb;
    
   a_theta_MY_up_L = nan((n_pol_theta+1)*2,n_wb_MY_up,n_test_MY_up);
   a_eta_MY_up_L = nan((n_pol_eta+1)*2,n_wb_MY_up,n_test_MY_up);
   a_phi_MY_up_L = nan((n_pol_phi+1)*2,n_wb_MY_up,n_test_MY_up);
    
    
   a_theta_MY_up_R = nan((n_pol_theta+1)*2,n_wb_MY_up,n_test_MY_up);
   a_eta_MY_up_R = nan((n_pol_eta+1)*2,n_wb_MY_up,n_test_MY_up);
   a_phi_MY_up_R = nan((n_pol_phi+1)*2,n_wb_MY_up,n_test_MY_up);
    
    
    for i = 1:n_test_MY_up
        
        a_dev_theta_L   = pathDB.control_glob.pitch_theta_back.*MY_up_range(i);
        a_dev_eta_L     = pathDB.control_glob.pitch_eta_back.*MY_up_range(i);
        a_dev_phi_L     = pathDB.control_glob.pitch_phi_back.*MY_up_range(i);
        a_dev_theta_R   = pathDB.control_glob.pitch_theta_back.*MY_up_range(i);
        a_dev_eta_R     = pathDB.control_glob.pitch_eta_back.*MY_up_range(i);
        a_dev_phi_R     = pathDB.control_glob.pitch_phi_back.*MY_up_range(i);
        
                
        if i > 0 && i <= (n_test_MY_up/4)
        
            for j = 1:n_wb_MY_up
        
            a_theta_MY_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_MY_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_MY_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_MY_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_MY_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_MY_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
            end
        
        elseif i > (n_test_MY_up/4) && i <= (2*n_test_MY_up/4)
            
            for j = 1:n_wb_MY_up
        
            a_theta_MY_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
            a_eta_MY_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_MY_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_MY_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
            a_eta_MY_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_MY_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
            end
            
        elseif i > (2*n_test_MY_up/4) && i <= (3*n_test_MY_up/4)
            
            for j = 1:n_wb_MY_up
        
            a_theta_MY_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_MY_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
            a_phi_MY_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
            a_theta_MY_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_MY_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
            a_phi_MY_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
            end
                
        elseif i > (3*n_test_MY_up/4) && i <= (4*n_test_MY_up/4)
            
            for j = 1:n_wb_MY_up
        
            a_theta_MY_up_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_MY_up_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_MY_up_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
            a_theta_MY_up_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
            a_eta_MY_up_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
            a_phi_MY_up_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;
        
            end
                
        end
        
    end
    
    wing_l_MY_up = wing_l_glob;
    down_up_MY_up = down_up_glob;
    f_MY_up = f_glob;
    t_MY_up = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    theta_MY_up_L = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    eta_MY_up_L = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    phi_MY_up_L = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    theta_MY_up_R = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    eta_MY_up_R = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    phi_MY_up_R = zeros(n_wb_MY_up*100+1,n_test_MY_up);
    
        
    for i = 1:n_test_MY_up
 
        for j = 1:n_wb_MY_up
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_MY_up, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_MY_up, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_MY_up, 101, (j-1), j, 0 );

            t_MY_up(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_MY_up_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MY_up_L(:,j,i);
            eta_MY_up_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MY_up_L(:,j,i);
            phi_MY_up_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MY_up_L(:,j,i);
            theta_MY_up_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MY_up_R(:,j,i);
            eta_MY_up_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MY_up_R(:,j,i);
            phi_MY_up_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MY_up_R(:,j,i);

        end
        
    
    end 
    
    
    
   % MZ wingbeats---------------------------------------------
    
%    MZ_range = [ 0:(0.011/19):0.011 0:(0.011/19):0.011 0:(0.011/19):0.011 0:(0.011/19):0.011 ];

    MZ_range = [ 0:(0.0079/19):0.0079 0:(0.0079/19):0.0079 0:(0.0079/19):0.0079 0:(0.0079/19):0.0079 ];
   
   n_test_MZ = length(MZ_range);
    
   n_wb_MZ = nr_of_wb;
    
   a_theta_MZ_L = nan((n_pol_theta+1)*2,n_wb_MZ,n_test_MZ);
   a_eta_MZ_L = nan((n_pol_eta+1)*2,n_wb_MZ,n_test_MZ);
   a_phi_MZ_L = nan((n_pol_phi+1)*2,n_wb_MZ,n_test_MZ);
    
    
   a_theta_MZ_R = nan((n_pol_theta+1)*2,n_wb_MZ,n_test_MZ);
   a_eta_MZ_R = nan((n_pol_eta+1)*2,n_wb_MZ,n_test_MZ);
   a_phi_MZ_R = nan((n_pol_phi+1)*2,n_wb_MZ,n_test_MZ);
    
    
    for i = 1:n_test_MZ
        
            if i <= (n_test_MZ/4)
        
            a_dev_theta_L   = pathDB.control_glob.yaw_theta_left.*MZ_range(i);
            a_dev_eta_L     = pathDB.control_glob.yaw_eta_left.*MZ_range(i);
            a_dev_phi_L     = pathDB.control_glob.yaw_phi_left.*MZ_range(i);
            a_dev_theta_R   = pathDB.control_glob.yaw_theta_right.*MZ_range(i);
            a_dev_eta_R     = pathDB.control_glob.yaw_eta_right.*MZ_range(i);
            a_dev_phi_R     = pathDB.control_glob.yaw_phi_right.*MZ_range(i);

            for j = 1:n_wb_MZ

                a_theta_MZ_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_MZ_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_MZ_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_MZ_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_MZ_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_MZ_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end
        
        elseif i > (n_test_MZ/4) && i <= (2*n_test_MZ/4)
        
            a_dev_theta_L   = pathDB.control_glob.yaw_theta_left.*MZ_range(i);
            a_dev_eta_L     = pathDB.control_glob.yaw_eta_left.*MZ_range(i);
            a_dev_phi_L     = pathDB.control_glob.yaw_phi_left.*MZ_range(i);
            a_dev_theta_R   = pathDB.control_glob.yaw_theta_right.*MZ_range(i);
            a_dev_eta_R     = pathDB.control_glob.yaw_eta_right.*MZ_range(i);
            a_dev_phi_R     = pathDB.control_glob.yaw_phi_right.*MZ_range(i);

            for j = 1:n_wb_MZ

                a_theta_MZ_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                a_eta_MZ_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MZ_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_MZ_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                a_eta_MZ_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MZ_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end
            
        elseif i > (2*n_test_MZ/4) && i <= (3*n_test_MZ/4)
        
            a_dev_theta_L   = pathDB.control_glob.yaw_theta_left.*MZ_range(i);
            a_dev_eta_L     = pathDB.control_glob.yaw_eta_left.*MZ_range(i);
            a_dev_phi_L     = pathDB.control_glob.yaw_phi_left.*MZ_range(i);
            a_dev_theta_R   = pathDB.control_glob.yaw_theta_right.*MZ_range(i);
            a_dev_eta_R     = pathDB.control_glob.yaw_eta_right.*MZ_range(i);
            a_dev_phi_R     = pathDB.control_glob.yaw_phi_right.*MZ_range(i);

            for j = 1:n_wb_MZ

                a_theta_MZ_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MZ_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                a_phi_MZ_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                a_theta_MZ_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MZ_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
                a_phi_MZ_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            end
                
        elseif i > (3*n_test_MZ/4) && i <= (4*n_test_MZ/4)
        
            a_dev_theta_L   = pathDB.control_glob.yaw_theta_left.*MZ_range(i);
            a_dev_eta_L     = pathDB.control_glob.yaw_eta_left.*MZ_range(i);
            a_dev_phi_L     = pathDB.control_glob.yaw_phi_left.*MZ_range(i);
            a_dev_theta_R   = pathDB.control_glob.yaw_theta_right.*MZ_range(i);
            a_dev_eta_R     = pathDB.control_glob.yaw_eta_right.*MZ_range(i);
            a_dev_phi_R     = pathDB.control_glob.yaw_phi_right.*MZ_range(i);

            for j = 1:n_wb_MZ

                a_theta_MZ_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MZ_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MZ_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                a_theta_MZ_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                a_eta_MZ_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                a_phi_MZ_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_R;

            end
                
        end
        
    end
    
    wing_l_MZ = wing_l_glob;
    down_up_MZ = down_up_glob;
    f_MZ = f_glob;
    t_MZ = zeros(n_wb_MZ*100+1,n_test_MZ);
    theta_MZ_L = zeros(n_wb_MZ*100+1,n_test_MZ);
    eta_MZ_L = zeros(n_wb_MZ*100+1,n_test_MZ);
    phi_MZ_L = zeros(n_wb_MZ*100+1,n_test_MZ);
    theta_MZ_R = zeros(n_wb_MZ*100+1,n_test_MZ);
    eta_MZ_R = zeros(n_wb_MZ*100+1,n_test_MZ);
    phi_MZ_R = zeros(n_wb_MZ*100+1,n_test_MZ);
    
        
    for i = 1:n_test_MZ
 
        for j = 1:n_wb_MZ
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_MZ, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_MZ, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_MZ, 101, (j-1), j, 0 );

            t_MZ(((j-1)*100+1):(j*100+1),i) = t_test;
            theta_MZ_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MZ_L(:,j,i);
            eta_MZ_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MZ_L(:,j,i);
            phi_MZ_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MZ_L(:,j,i);
            theta_MZ_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_MZ_R(:,j,i);
            eta_MZ_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_MZ_R(:,j,i);
            phi_MZ_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_MZ_R(:,j,i);

        end
    
    end
    
    
    
    %----------------------------------------------------------------------
       
   
   
   FX_forward_test = {};
   
   FX_forward_test.nr_of_wb = n_wb_FX_forward;
   FX_forward_test.nr_of_tests = n_test_FX_forward;
   FX_forward_test.a_theta_L = a_theta_FX_forward_L;
   FX_forward_test.a_theta_R = a_theta_FX_forward_R;
   FX_forward_test.a_eta_L = a_eta_FX_forward_L;
   FX_forward_test.a_eta_R = a_eta_FX_forward_R;
   FX_forward_test.a_phi_L = a_phi_FX_forward_L;
   FX_forward_test.a_phi_R = a_phi_FX_forward_R;
   FX_forward_test.freq = f_FX_forward;
   FX_forward_test.down_up = down_up_FX_forward;
   FX_forward_test.wing_length = wing_l_FX_forward;
   FX_forward_test.time_ref = t_FX_forward;
   FX_forward_test.theta_L_ref = theta_FX_forward_L;
   FX_forward_test.theta_R_ref = theta_FX_forward_R;
   FX_forward_test.eta_L_ref = eta_FX_forward_L;
   FX_forward_test.eta_R_ref = eta_FX_forward_R;
   FX_forward_test.phi_L_ref = phi_FX_forward_L;
   FX_forward_test.phi_R_ref = phi_FX_forward_R;
   FX_forward_test.n_pol_theta = n_pol_theta;
   FX_forward_test.n_pol_eta = n_pol_eta;
   FX_forward_test.n_pol_phi = n_pol_phi;

   figure()
   plot(FX_forward_test.time_ref,FX_forward_test.phi_L_ref)
   
   figure()
   plot(FX_forward_test.time_ref,FX_forward_test.phi_R_ref)
   
   FX_back_test = {};
   
   FX_back_test.nr_of_wb = n_wb_FX_back;
   FX_back_test.nr_of_tests = n_test_FX_back;
   FX_back_test.a_theta_L = a_theta_FX_back_L;
   FX_back_test.a_theta_R = a_theta_FX_back_R;
   FX_back_test.a_eta_L = a_eta_FX_back_L;
   FX_back_test.a_eta_R = a_eta_FX_back_R;
   FX_back_test.a_phi_L = a_phi_FX_back_L;
   FX_back_test.a_phi_R = a_phi_FX_back_R;
   FX_back_test.freq = f_FX_back;
   FX_back_test.down_up = down_up_FX_back;
   FX_back_test.wing_length = wing_l_FX_back;
   FX_back_test.time_ref = t_FX_back;
   FX_back_test.theta_L_ref = theta_FX_back_L;
   FX_back_test.theta_R_ref = theta_FX_back_R;
   FX_back_test.eta_L_ref = eta_FX_back_L;
   FX_back_test.eta_R_ref = eta_FX_back_R;
   FX_back_test.phi_L_ref = phi_FX_back_L;
   FX_back_test.phi_R_ref = phi_FX_back_R;
   FX_back_test.n_pol_theta = n_pol_theta;
   FX_back_test.n_pol_eta = n_pol_eta;
   FX_back_test.n_pol_phi = n_pol_phi;

   figure()
   plot(FX_back_test.time_ref,FX_back_test.phi_L_ref)
   
   figure()
   plot(FX_back_test.time_ref,FX_back_test.phi_R_ref)
   
   FY_test = {};
   
   FY_test.nr_of_wb = n_wb_FY;
   FY_test.nr_of_tests = n_test_FY;
   FY_test.a_theta_L = a_theta_FY_L;
   FY_test.a_theta_R = a_theta_FY_R;
   FY_test.a_eta_L = a_eta_FY_L;
   FY_test.a_eta_R = a_eta_FY_R;
   FY_test.a_phi_L = a_phi_FY_L;
   FY_test.a_phi_R = a_phi_FY_R;
   FY_test.freq = f_FY;
   FY_test.down_up = down_up_FY;
   FY_test.wing_length = wing_l_FY;
   FY_test.time_ref = t_FY;
   FY_test.theta_L_ref = theta_FY_L;
   FY_test.theta_R_ref = theta_FY_R;
   FY_test.eta_L_ref = eta_FY_L;
   FY_test.eta_R_ref = eta_FY_R;
   FY_test.phi_L_ref = phi_FY_L;
   FY_test.phi_R_ref = phi_FY_R;
   FY_test.n_pol_theta = n_pol_theta;
   FY_test.n_pol_eta = n_pol_eta;
   FY_test.n_pol_phi = n_pol_phi;
   
   figure()
   plot(FY_test.time_ref,FY_test.phi_L_ref)
   
   figure()
   plot(FY_test.time_ref,FY_test.phi_R_ref)

   
   FZ_up_test = {};
   
   FZ_up_test.nr_of_wb = n_wb_FZ_up;
   FZ_up_test.nr_of_tests = n_test_FZ_up;
   FZ_up_test.a_theta_L = a_theta_FZ_up_L;
   FZ_up_test.a_theta_R = a_theta_FZ_up_R;
   FZ_up_test.a_eta_L = a_eta_FZ_up_L;
   FZ_up_test.a_eta_R = a_eta_FZ_up_R;
   FZ_up_test.a_phi_L = a_phi_FZ_up_L;
   FZ_up_test.a_phi_R = a_phi_FZ_up_R;
   FZ_up_test.freq = f_FZ_up;
   FZ_up_test.down_up = down_up_FZ_up;
   FZ_up_test.wing_length = wing_l_FZ_up;
   FZ_up_test.time_ref = t_FZ_up;
   FZ_up_test.theta_L_ref = theta_FZ_up_L;
   FZ_up_test.theta_R_ref = theta_FZ_up_R;
   FZ_up_test.eta_L_ref = eta_FZ_up_L;
   FZ_up_test.eta_R_ref = eta_FZ_up_R;
   FZ_up_test.phi_L_ref = phi_FZ_up_L;
   FZ_up_test.phi_R_ref = phi_FZ_up_R;
   FZ_up_test.n_pol_theta = n_pol_theta;
   FZ_up_test.n_pol_eta = n_pol_eta;
   FZ_up_test.n_pol_phi = n_pol_phi;
   
   figure()
   plot(FZ_up_test.time_ref,FZ_up_test.phi_L_ref)
   
   figure()
   plot(FZ_up_test.time_ref,FZ_up_test.phi_R_ref)
   
   MX_test = {};
   
   MX_test.nr_of_wb = n_wb_MX;
   MX_test.nr_of_tests = n_test_MX;
   MX_test.a_theta_L = a_theta_MX_L;
   MX_test.a_theta_R = a_theta_MX_R;
   MX_test.a_eta_L = a_eta_MX_L;
   MX_test.a_eta_R = a_eta_MX_R;
   MX_test.a_phi_L = a_phi_MX_L;
   MX_test.a_phi_R = a_phi_MX_R;
   MX_test.freq = f_MX;
   MX_test.down_up = down_up_MX;
   MX_test.wing_length = wing_l_MX;
   MX_test.time_ref = t_MX;
   MX_test.theta_L_ref = theta_MX_L;
   MX_test.theta_R_ref = theta_MX_R;
   MX_test.eta_L_ref = eta_MX_L;
   MX_test.eta_R_ref = eta_MX_R;
   MX_test.phi_L_ref = phi_MX_L;
   MX_test.phi_R_ref = phi_MX_R;
   MX_test.n_pol_theta = n_pol_theta;
   MX_test.n_pol_eta = n_pol_eta;
   MX_test.n_pol_phi = n_pol_phi;
   
   figure()
   plot(MX_test.time_ref,MX_test.phi_L_ref)
   
   figure()
   plot(MX_test.time_ref,MX_test.phi_R_ref)

   MY_up_test = {};
   
   MY_up_test.nr_of_wb = n_wb_MY_up;
   MY_up_test.nr_of_tests = n_test_MY_up;
   MY_up_test.a_theta_L = a_theta_MY_up_L;
   MY_up_test.a_theta_R = a_theta_MY_up_R;
   MY_up_test.a_eta_L = a_eta_MY_up_L;
   MY_up_test.a_eta_R = a_eta_MY_up_R;
   MY_up_test.a_phi_L = a_phi_MY_up_L;
   MY_up_test.a_phi_R = a_phi_MY_up_R;
   MY_up_test.freq = f_MY_up;
   MY_up_test.down_up = down_up_MY_up;
   MY_up_test.wing_length = wing_l_MY_up;
   MY_up_test.time_ref = t_MY_up;
   MY_up_test.theta_L_ref = theta_MY_up_L;
   MY_up_test.theta_R_ref = theta_MY_up_R;
   MY_up_test.eta_L_ref = eta_MY_up_L;
   MY_up_test.eta_R_ref = eta_MY_up_R;
   MY_up_test.phi_L_ref = phi_MY_up_L;
   MY_up_test.phi_R_ref = phi_MY_up_R;
   MY_up_test.n_pol_theta = n_pol_theta;
   MY_up_test.n_pol_eta = n_pol_eta;
   MY_up_test.n_pol_phi = n_pol_phi;
    
   figure()
   plot(MY_up_test.time_ref,MY_up_test.phi_L_ref)
   
   figure()
   plot(MY_up_test.time_ref,MY_up_test.phi_R_ref)
   
   
   MZ_test = {};
   
   MZ_test.nr_of_wb = n_wb_MZ;
   MZ_test.nr_of_tests = n_test_MZ;
   MZ_test.a_theta_L = a_theta_MZ_L;
   MZ_test.a_theta_R = a_theta_MZ_R;
   MZ_test.a_eta_L = a_eta_MZ_L;
   MZ_test.a_eta_R = a_eta_MZ_R;
   MZ_test.a_phi_L = a_phi_MZ_L;
   MZ_test.a_phi_R = a_phi_MZ_R;
   MZ_test.freq = f_MZ;
   MZ_test.down_up = down_up_MZ;
   MZ_test.wing_length = wing_l_MZ;
   MZ_test.time_ref = t_MZ;
   MZ_test.theta_L_ref = theta_MZ_L;
   MZ_test.theta_R_ref = theta_MZ_R;
   MZ_test.eta_L_ref = eta_MZ_L;
   MZ_test.eta_R_ref = eta_MZ_R;
   MZ_test.phi_L_ref = phi_MZ_L;
   MZ_test.phi_R_ref = phi_MZ_R;
   MZ_test.n_pol_theta = n_pol_theta;
   MZ_test.n_pol_eta = n_pol_eta;
   MZ_test.n_pol_phi = n_pol_phi;  
   
   figure()
   plot(MZ_test.time_ref,MZ_test.phi_L_ref)
   
   figure()
   plot(MZ_test.time_ref,MZ_test.phi_R_ref)
    
   
   save(savefile,'FX_forward_test','FX_back_test','FY_test','FZ_up_test','MX_test','MY_up_test','MZ_test')

    
end