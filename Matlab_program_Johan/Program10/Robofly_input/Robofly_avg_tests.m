function Robofly_avg_tests( settings, pathDB )

    savefile = 'Robofly_avg_tests.mat';

    %----------------------------------------------------------------------
    %
    % Generate Robofly_avg_tests.mat containing the following data:
    %
    %   test_id
    %   
    %   number of simulated wingbeats per tests
    %
    %   vector of Legendre fitting coefficients for theta_L, eta_L, phi_L
    %   and theta_R, eta_R, phi_R for each simulated wingbeat
    %
    %   Robofly flap frequency per simulated wingbeat
    %
    %----------------------------------------------------------------------
    
    
    nr_of_wb = 7; % number of wingbeats per sequence
    
    avg_seqs = [1 4 6 8 15]; % Selection of sequences which are representatives of sequence averaged wingbeats.
    
    n_avg_seqs = length(avg_seqs);
    
    n_pol_theta = (length(pathDB.a_glob.theta1)-1);
    
    n_pol_eta = (length(pathDB.a_glob.eta1)-1);
    
    n_pol_phi = (length(pathDB.a_glob.phi1)-1);
    
        
    % Global average wingbeats---------------------------------------------
    
    n_test_glob = 3;
    
    n_wb_glob = nr_of_wb;
    
    a_theta_glob_L = nan((n_pol_theta+1)*2,n_wb_glob,n_test_glob);
    
    a_eta_glob_L = nan((n_pol_eta+1)*2,n_wb_glob,n_test_glob);
    
    a_phi_glob_L = nan((n_pol_phi+1)*2,n_wb_glob,n_test_glob);
    
    
    a_theta_glob_R = nan((n_pol_theta+1)*2,n_wb_glob,n_test_glob);
    
    a_eta_glob_R = nan((n_pol_eta+1)*2,n_wb_glob,n_test_glob);
    
    a_phi_glob_R = nan((n_pol_phi+1)*2,n_wb_glob,n_test_glob);
    
    
    for i = 1:n_test_glob
        
        for j = 1:n_wb_glob
        
            a_theta_glob_L(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];

            a_eta_glob_L(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];

            a_phi_glob_L(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];

            a_theta_glob_R(:,j,i) = [pathDB.a_glob.theta1; pathDB.a_glob.theta2];

            a_eta_glob_R(:,j,i) = [pathDB.a_glob.eta1; pathDB.a_glob.eta2];

            a_phi_glob_R(:,j,i) = [pathDB.a_glob.phi1; pathDB.a_glob.phi2];
        
        end
        
    end
    
    wing_l_glob = mean(pathDB.wing_l);
    
    down_up_glob = pathDB.down_up_glob;
    
    f_glob = pathDB.f_glob;
    
    t_glob = zeros(nr_of_wb*100+1,n_test_glob);
    
    theta_glob_L = zeros(nr_of_wb*100+1,n_test_glob);
    
    eta_glob_L = zeros(nr_of_wb*100+1,n_test_glob);
    
    phi_glob_L = zeros(nr_of_wb*100+1,n_test_glob);
    
    theta_glob_R = zeros(nr_of_wb*100+1,n_test_glob);
    
    eta_glob_R = zeros(nr_of_wb*100+1,n_test_glob);
    
    phi_glob_R = zeros(nr_of_wb*100+1,n_test_glob);
    
        
    for i = 1:n_test_glob
 
        for j = 1:nr_of_wb
            
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_glob, 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_glob, 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_glob, 101, (j-1), j, 0 );

            t_glob(((j-1)*100+1):(j*100+1),i) = t_test;

            theta_glob_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_glob_L(:,j,i);

            eta_glob_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_glob_L(:,j,i);

            phi_glob_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_glob_L(:,j,i);

            theta_glob_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_glob_R(:,j,i);

            eta_glob_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_glob_R(:,j,i);

            phi_glob_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_glob_R(:,j,i);

        end
    
    end
    

    %----------------------------------------------------------------------
    
        
    a_avg_temp = pathDB.a_avg;
    
    a_avg_names = fieldnames(a_avg_temp);
    
    % down-up ratio:
    
    down_up_avg_names = fieldnames(pathDB.down_up_avg);
        
    down_up_avg = nan(1,n_avg_seqs);
    
    for i = 1:n_avg_seqs
    
        down_up_avg(i) = pathDB.down_up_avg.(char(down_up_avg_names(avg_seqs(i))));
    
    end
    
    % Average frequency:
    
    f_avg_names = fieldnames(pathDB.f_avg);
    
    f_avg = nan(1,n_avg_seqs);
    
    for i = 1:n_avg_seqs
    
        f_avg(i) = pathDB.f_avg.(char(f_avg_names(avg_seqs(i))));
    
    end
    
    % Wing lenght
    
    wing_lengths = pathDB.wing_l(avg_seqs);
    
    %----------------------------------------------------------------------
    
    
    % Sequence averaged symmetric wingbeats--------------------------------
    
    n_test_sym = n_avg_seqs;
    
    n_wb_sym = nr_of_wb;
    
    a_theta_sym_L = nan((n_pol_theta+1)*2,n_wb_sym,n_test_sym);
    
    a_eta_sym_L = nan((n_pol_eta+1)*2,n_wb_sym,n_test_sym);
   
    a_phi_sym_L = nan((n_pol_phi+1)*2,n_wb_sym,n_test_sym);
    
    
    a_theta_sym_R = nan((n_pol_theta+1)*2,n_wb_sym,n_test_sym);
    
    a_eta_sym_R = nan((n_pol_eta+1)*2,n_wb_sym,n_test_sym);
    
    a_phi_sym_R = nan((n_pol_phi+1)*2,n_wb_sym,n_test_sym);
    
    
    for i = 1:n_test_sym
                
        a_t = a_avg_temp.(char(a_avg_names(avg_seqs(i))));
        
        for j = 1:n_wb_sym
            
            a_theta_sym_L(:,j,i) = [a_t.theta_LR1; a_t.theta_LR2];
            
            a_eta_sym_L(:,j,i) = [a_t.eta_LR1; a_t.eta_LR2];
            
            a_phi_sym_L(:,j,i) = [a_t.phi_LR1; a_t.phi_LR2];
            
            a_theta_sym_R(:,j,i) = [a_t.theta_LR1; a_t.theta_LR2];
            
            a_eta_sym_R(:,j,i) = [a_t.eta_LR1; a_t.eta_LR2];
            
            a_phi_sym_R(:,j,i) = [a_t.phi_LR1; a_t.phi_LR2];
            
        end
        
    end
    
    wing_l_sym = wing_lengths;
    
    down_up_sym = zeros(n_wb_sym,n_test_sym);
    
    f_sym = zeros(n_wb_sym,n_test_sym);
    
    for i = 1:n_test_sym
        
        for j = 1:n_wb_sym
           
            down_up_sym(j,i) = down_up_avg(i);
            
            f_sym(j,i) = f_avg(i);
            
        end
        
    end
    
    t_sym = zeros(n_wb_sym*100+1,n_test_sym);
    
    theta_sym_L = zeros(n_wb_sym*100+1,n_test_sym);
    
    eta_sym_L = zeros(n_wb_sym*100+1,n_test_sym);
    
    phi_sym_L = zeros(n_wb_sym*100+1,n_test_sym);
    
    theta_sym_R = zeros(n_wb_sym*100+1,n_test_sym);
    
    eta_sym_R = zeros(n_wb_sym*100+1,n_test_sym);
    
    phi_sym_R = zeros(n_wb_sym*100+1,n_test_sym);
    
        
    for i = 1:n_test_sym
 
        for j = 1:n_wb_sym
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg(i), 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg(i), 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg(i), 101, (j-1), j, 0 );

        t_sym(((j-1)*100+1):(j*100+1),i) = t_test;

        theta_sym_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_sym_L(:,j,i);

        eta_sym_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_sym_L(:,j,i);

        phi_sym_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_sym_L(:,j,i);

        theta_sym_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_sym_R(:,j,i);

        eta_sym_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_sym_R(:,j,i);

        phi_sym_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_sym_R(:,j,i);

        end
    
    end

%     clmp = colormap(jet(n_avg_seqs));
%         
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),theta_sym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),eta_sym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),phi_sym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),theta_sym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),eta_sym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_sym(:,i),phi_sym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
    
    
    

    
    
    %----------------------------------------------------------------------
    
    
    
    % Sequence averaged asymmetric wingbeats--------------------------------
    
    n_test_asym = n_avg_seqs;
    
    n_wb_asym = nr_of_wb;    
    
    a_theta_asym_L = nan((n_pol_theta+1)*2,n_wb_asym,n_test_asym);
    
    a_eta_asym_L = nan((n_pol_eta+1)*2,n_wb_asym,n_test_asym);
   
    a_phi_asym_L = nan((n_pol_phi+1)*2,n_wb_asym,n_test_asym);
    
    
    a_theta_asym_R = nan((n_pol_theta+1)*2,n_wb_asym,n_test_asym);
    
    a_eta_asym_R = nan((n_pol_eta+1)*2,n_wb_asym,n_test_asym);
    
    a_phi_asym_R = nan((n_pol_phi+1)*2,n_wb_asym,n_test_asym);
    
    
    for i = 1:n_test_asym
                
        a_t = a_avg_temp.(char(a_avg_names(avg_seqs(i))));
        
        for j = 1:n_wb_asym
            
            a_theta_asym_L(:,j,i) = [a_t.theta_L1; a_t.theta_L2];
            
            a_eta_asym_L(:,j,i) = [a_t.eta_L1; a_t.eta_L2];
            
            a_phi_asym_L(:,j,i) = [a_t.phi_L1; a_t.phi_L2];
            
            a_theta_asym_R(:,j,i) = [a_t.theta_R1; a_t.theta_R2];
            
            a_eta_asym_R(:,j,i) = [a_t.eta_R1; a_t.eta_R2];
            
            a_phi_asym_R(:,j,i) = [a_t.phi_R1; a_t.phi_R2];
            
        end
        
    end
    
    wing_l_asym = wing_lengths;
    
    down_up_asym = zeros(n_wb_asym,n_test_asym);
    
    f_asym = zeros(n_wb_asym,n_test_asym);
    
    for i = 1:n_test_asym
        
        for j = 1:n_wb_asym
           
            down_up_asym(j,i) = down_up_avg(i);
            
            f_asym(j,i) = f_avg(i);
            
        end
        
    end
    
    
    t_asym = zeros(n_wb_asym*100+1,n_test_asym);
    
    theta_asym_L = zeros(n_wb_asym*100+1,n_test_asym);
    
    eta_asym_L = zeros(n_wb_asym*100+1,n_test_asym);
    
    phi_asym_L = zeros(n_wb_asym*100+1,n_test_asym);
    
    theta_asym_R = zeros(n_wb_asym*100+1,n_test_asym);
    
    eta_asym_R = zeros(n_wb_asym*100+1,n_test_asym);
    
    phi_asym_R = zeros(n_wb_asym*100+1,n_test_asym);
    
        
    for i = 1:n_test_asym
 
        for j = 1:n_wb_asym
            [ ~, X_test_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg(i), 101, (j-1), j, 0 );
            [ ~, X_test_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg(i), 101, (j-1), j, 0 );
            [ t_test, X_test_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg(i), 101, (j-1), j, 0 );

        t_asym(((j-1)*100+1):(j*100+1),i) = t_test;

        theta_asym_L(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_asym_L(:,j,i);

        eta_asym_L(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_asym_L(:,j,i);

        phi_asym_L(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_asym_L(:,j,i);

        theta_asym_R(((j-1)*100+1):(j*100+1),i) = X_test_theta*a_theta_asym_R(:,j,i);

        eta_asym_R(((j-1)*100+1):(j*100+1),i) = X_test_eta*a_eta_asym_R(:,j,i);

        phi_asym_R(((j-1)*100+1):(j*100+1),i) = X_test_phi*a_phi_asym_R(:,j,i);

        end
    
    end
        
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),theta_asym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),eta_asym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),phi_asym_L(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),theta_asym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),eta_asym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
%     
%     figure()
%     hold on
%     for i = 1:n_avg_seqs
%         plot(t_asym(:,i),phi_asym_R(:,i),'Color',clmp(i,:))
%     end
%     hold off
    
    
    %----------------------------------------------------------------------
    
    
    
    % Combine the tests in structures:
    
    glob_test = {};
    
    sym_test = {};
    
    asym_test1 = {};
    
    asym_test2 = {};
    
    
    
    glob_test.nr_of_wb = n_wb_glob;
    glob_test.nr_of_tests = n_test_glob;
    glob_test.a_theta_L = a_theta_glob_L;
    glob_test.a_theta_R = a_theta_glob_R;
    glob_test.a_eta_L = a_eta_glob_L;
    glob_test.a_eta_R = a_eta_glob_R;
    glob_test.a_phi_L = a_phi_glob_L;
    glob_test.a_phi_R = a_phi_glob_R;
    glob_test.freq = f_glob;
    glob_test.down_up = down_up_glob;
    glob_test.wing_length = wing_l_glob;
    glob_test.time_ref = t_glob;
    glob_test.theta_L_ref = theta_glob_L;
    glob_test.theta_R_ref = theta_glob_R;
    glob_test.eta_L_ref = eta_glob_L;
    glob_test.eta_R_ref = eta_glob_R;
    glob_test.phi_L_ref = phi_glob_L;
    glob_test.phi_R_ref = phi_glob_R;
    glob_test.n_pol_theta = n_pol_theta;
    glob_test.n_pol_eta = n_pol_eta;
    glob_test.n_pol_phi = n_pol_phi;
        
    sym_test.nr_of_wb = n_wb_sym;
    sym_test.nr_of_tests = n_test_sym;
    sym_test.a_theta_L = a_theta_sym_L;
    sym_test.a_theta_R = a_theta_sym_R;
    sym_test.a_eta_L = a_eta_sym_L;
    sym_test.a_eta_R = a_eta_sym_R;
    sym_test.a_phi_L = a_phi_sym_L;
    sym_test.a_phi_R = a_phi_sym_R;
    sym_test.freq = f_sym;
    sym_test.down_up = down_up_sym;
    sym_test.wing_length = wing_l_sym;
    sym_test.time_ref = t_sym;
    sym_test.theta_L_ref = theta_sym_L;
    sym_test.theta_R_ref = theta_sym_R;
    sym_test.eta_L_ref = eta_sym_L;
    sym_test.eta_R_ref = eta_sym_R;
    sym_test.phi_L_ref = phi_sym_L;
    sym_test.phi_R_ref = phi_sym_R;
    sym_test.n_pol_theta = n_pol_theta;
    sym_test.n_pol_eta = n_pol_eta;
    sym_test.n_pol_phi = n_pol_phi;

    
    asym_test1.nr_of_wb = n_wb_asym;
    asym_test1.nr_of_tests = n_test_asym;
    asym_test1.a_theta_L = a_theta_asym_L;
    asym_test1.a_theta_R = a_theta_asym_R;
    asym_test1.a_eta_L = a_eta_asym_L;
    asym_test1.a_eta_R = a_eta_asym_R;
    asym_test1.a_phi_L = a_phi_asym_L;
    asym_test1.a_phi_R = a_phi_asym_R;
    asym_test1.freq = f_asym;
    asym_test1.down_up = down_up_asym;
    asym_test1.wing_length = wing_l_asym;
    asym_test1.time_ref = t_asym;
    asym_test1.theta_L_ref = theta_asym_L;
    asym_test1.theta_R_ref = theta_asym_R;
    asym_test1.eta_L_ref = eta_asym_L;
    asym_test1.eta_R_ref = eta_asym_R;
    asym_test1.phi_L_ref = phi_asym_L;
    asym_test1.phi_R_ref = phi_asym_R;
    asym_test1.n_pol_theta = n_pol_theta;
    asym_test1.n_pol_eta = n_pol_eta;
    asym_test1.n_pol_phi = n_pol_phi;
    
    asym_test2.nr_of_wb = n_wb_asym;
    asym_test2.nr_of_tests = n_test_asym;
    asym_test2.a_theta_L = a_theta_asym_R;
    asym_test2.a_theta_R = a_theta_asym_L;
    asym_test2.a_eta_L = a_eta_asym_R;
    asym_test2.a_eta_R = a_eta_asym_L;
    asym_test2.a_phi_L = a_phi_asym_R;
    asym_test2.a_phi_R = a_phi_asym_L;
    asym_test2.freq = f_asym;
    asym_test2.down_up = down_up_asym;
    asym_test2.wing_length = wing_l_asym;
    asym_test2.time_ref = t_asym;
    asym_test2.theta_L_ref = theta_asym_R;
    asym_test2.theta_R_ref = theta_asym_L;
    asym_test2.eta_L_ref = eta_asym_R;
    asym_test2.eta_R_ref = eta_asym_L;
    asym_test2.phi_L_ref = phi_asym_R;
    asym_test2.phi_R_ref = phi_asym_L;
    asym_test2.n_pol_theta = n_pol_theta;
    asym_test2.n_pol_eta = n_pol_eta;
    asym_test2.n_pol_phi = n_pol_phi;

    save(savefile,'glob_test','sym_test','asym_test1','asym_test2')


end

