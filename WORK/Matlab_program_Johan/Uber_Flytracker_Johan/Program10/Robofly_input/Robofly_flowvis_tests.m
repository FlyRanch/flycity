function Robofly_flowvis_tests( settings, pathDB )

   savefile = 'Robofly_flow_vis_tests.mat';

   nr_of_wb = 13; % number of wingbeats per sequence
    
   n_pol_theta = (length(pathDB.a_glob.theta1)-1);
   n_pol_eta = (length(pathDB.a_glob.eta1)-1);
   n_pol_phi = (length(pathDB.a_glob.phi1)-1);
   
   
   % Global average wingbeats---------------------------------------------
    
   
   n_wb_glob = nr_of_wb;
       
    wing_l_glob = mean(pathDB.wing_l);
    down_up_glob = pathDB.down_up_glob;
    f_glob = pathDB.f_glob;






   % MX wingbeats---------------------------------------------

    MX_range = [ 0:(0.0142/4):0.0142 ];
   
   n_test_MX = length(MX_range);
    
   n_wb_MX = nr_of_wb;
    
   a_theta_MX_L = nan((n_pol_theta+1)*2,n_wb_MX,n_test_MX);
   a_eta_MX_L = nan((n_pol_eta+1)*2,n_wb_MX,n_test_MX);
   a_phi_MX_L = nan((n_pol_phi+1)*2,n_wb_MX,n_test_MX);
    
    
   a_theta_MX_R = nan((n_pol_theta+1)*2,n_wb_MX,n_test_MX);
   a_eta_MX_R = nan((n_pol_eta+1)*2,n_wb_MX,n_test_MX);
   a_phi_MX_R = nan((n_pol_phi+1)*2,n_wb_MX,n_test_MX);
    
    
    for i = 1:n_test_MX
        
            a_dev_theta_L   = pathDB.control_glob.roll_theta_left.*MX_range(i);
            a_dev_eta_L     = pathDB.control_glob.roll_eta_left.*MX_range(i);
            a_dev_phi_L     = pathDB.control_glob.roll_phi_left.*MX_range(i);
            a_dev_theta_R   = pathDB.control_glob.roll_theta_right.*MX_range(i);
            a_dev_eta_R     = pathDB.control_glob.roll_eta_right.*MX_range(i);
            a_dev_phi_R     = pathDB.control_glob.roll_phi_right.*MX_range(i);
            
            a_trans_theta_L = transition_fit( a_dev_theta_L, n_pol_theta, down_up_glob );
            a_trans_eta_L   = transition_fit( a_dev_eta_L, n_pol_eta, down_up_glob );
            a_trans_phi_L   = transition_fit( a_dev_phi_L, n_pol_phi, down_up_glob );
            
            a_trans_theta_R = transition_fit( a_dev_theta_R, n_pol_theta, down_up_glob );
            a_trans_eta_R   = transition_fit( a_dev_eta_R, n_pol_eta, down_up_glob );
            a_trans_phi_R   = transition_fit( a_dev_phi_R, n_pol_phi, down_up_glob );


            for j = 1:n_wb_MX
                
                if j <= 6
                    
                    a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                    a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                    a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                    a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];
                    a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];
                    a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
                    
                elseif j == 7
                    
                    a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_trans_theta_L;
                    a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_trans_eta_L;
                    a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_trans_phi_L;
                    a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_trans_theta_R;
                    a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_trans_eta_R;
                    a_phi_MX_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_trans_phi_R;
                    
                elseif j > 7

                    a_theta_MX_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_L;
                    a_eta_MX_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_L;
                    a_phi_MX_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2]+a_dev_phi_L;
                    a_theta_MX_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2]+a_dev_theta_R;
                    a_eta_MX_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2]+a_dev_eta_R;
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
   plot(MX_test.time_ref,MX_test.theta_L_ref,MX_test.time_ref,MX_test.theta_R_ref)
   
   figure()
   plot(MX_test.time_ref,MX_test.eta_L_ref,MX_test.time_ref,MX_test.eta_R_ref)
   
   figure()
   plot(MX_test.time_ref,MX_test.phi_L_ref,MX_test.time_ref,MX_test.phi_R_ref)
    
   save(savefile,'MX_test')
   
end
