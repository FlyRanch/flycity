function [ a_glob ] = Global_wingbeat( a_avg )


    % Determine the global wingbeat:
    
    a_theta     = [];
    a_eta       = [];
    a_phi       = [];
    f           = [];
    down_up     = [];
    m_fly       = [];
    wing_l      = [];
    
    for i = 1:size(a_avg.theta_LR,2)
        
        if isnan(a_avg.theta_LR(1,i)) == 0
        
        a_theta     = [a_theta a_avg.theta_LR(:,i)];
        a_eta       = [a_eta a_avg.eta_LR(:,i)];
        a_phi       = [a_phi a_avg.phi_LR(:,i)];
        f           = [f a_avg.f(i)];
        down_up     = [down_up a_avg.down_up(i)];
        m_fly       = [m_fly a_avg.m_fly(i)];
        wing_l      = [wing_l a_avg.wing_l(i)];
        
        end
        
    end
    
    nr_seq          = size(a_theta,2);
    
    f_glob          = mean(f);
    down_up_glob    = mean(down_up);
    m_fly_glob      = mean(m_fly);
    wing_l_glob     = mean(wing_l);
    
    n_pol_theta     = (size(a_theta,1)-2)/2;
    n_pol_eta       = (size(a_eta,1)-2)/2;
    n_pol_phi       = (size(a_phi,1)-2)/2;
    
    % Weights:
    
    W_theta_1       = ones(1,nr_seq)*(1/nr_seq);
    W_theta_2       = ones(1,nr_seq)*(1/nr_seq);
    W_eta_1         = ones(1,nr_seq)*(1/nr_seq);
    W_eta_2         = ones(1,nr_seq)*(1/nr_seq);
    W_phi_1         = ones(1,nr_seq)*(1/nr_seq);
    W_phi_2         = ones(1,nr_seq)*(1/nr_seq);
    

    % Use function average_fit to obtain the average for the left and
    % right wing combined:

    [a_avg_theta_1, a_avg_theta_2] = average_fit(   a_theta((1:(n_pol_theta+1)),:) , ...
                                                        a_theta((n_pol_theta+2):(2*(n_pol_theta+1)),:), ...
                                                         n_pol_theta,W_theta_1,W_theta_2,down_up_glob);
                                                        
    [a_avg_eta_1, a_avg_eta_2] = average_fit(   a_eta((1:(n_pol_eta+1)),:) , ...
                                                    a_eta((n_pol_eta+2):(2*(n_pol_eta+1)),:), ...
                                                    n_pol_eta,W_eta_1,W_eta_2,down_up_glob);
                                                         
    [a_avg_phi_1, a_avg_phi_2] = average_fit(   a_phi((1:(n_pol_phi+1)),:) , ...
                                                    a_phi((n_pol_phi+2):(2*(n_pol_phi+1)),:), ...
                                                    n_pol_phi,W_phi_1,W_phi_2,down_up_glob);

        
    % Store average coefficients in a_avg     
        
    a_glob.theta    = [a_avg_theta_1; a_avg_theta_2];
    a_glob.eta      = [a_avg_eta_1; a_avg_eta_2];
    a_glob.phi      = [a_avg_phi_1; a_avg_phi_2];
    a_glob.f        = f_glob;
    a_glob.down_up  = down_up_glob;
    a_glob.m_fly    = m_fly_glob;
    a_glob.wing_l   = wing_l_glob;
    
end

