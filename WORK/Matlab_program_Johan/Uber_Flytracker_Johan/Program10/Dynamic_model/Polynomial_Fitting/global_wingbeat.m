function [ a_glob, f_glob, down_up_glob ] = global_wingbeat( a_avg_tot, f_avg_tot, down_up_tot, n_pol_theta, n_pol_eta, n_pol_phi )


    % Determine global average wingbeat from sequence averaged wingbeats.
    
    nr_of_seq = length(fieldnames(a_avg_tot));
        
    a_avg_theta_LR1 = [];
    a_avg_theta_LR2 = [];
    
    a_avg_eta_LR1 = [];
    a_avg_eta_LR2 = [];
    
    a_avg_phi_LR1 = [];
    a_avg_phi_LR2 = [];
    
    down_up = [];
    
    f_avg = [];
    
    a_avg_tot_names = fieldnames(a_avg_tot);
    
    down_up_tot_names = fieldnames(down_up_tot);
    
    f_avg_names = fieldnames(f_avg_tot);
    
    for i = 1:nr_of_seq
        
        a_temp = a_avg_tot.(char(a_avg_tot_names(i)));
        
        down_up_temp = down_up_tot.(char(down_up_tot_names(i)));
        
        f_avg_temp = f_avg_tot.(char(f_avg_names(i)));
        
        if isstruct(a_temp)==1
            
            a_avg_theta_LR1 = [a_avg_theta_LR1 a_temp.theta_LR1];
            a_avg_theta_LR2 = [a_avg_theta_LR2 a_temp.theta_LR2];
            
            a_avg_eta_LR1 = [a_avg_eta_LR1 a_temp.eta_LR1];
            a_avg_eta_LR2 = [a_avg_eta_LR2 a_temp.eta_LR2];
            
            a_avg_phi_LR1 = [a_avg_phi_LR1 a_temp.phi_LR1];
            a_avg_phi_LR2 = [a_avg_phi_LR2 a_temp.phi_LR2];
            
            down_up = [ down_up; down_up_temp ];
            
            f_avg = [f_avg; f_avg_temp];
            
        end
        
    end
    
    a_glob = {};
    
    nr_of_avg = size(a_avg_theta_LR1,2);

    
    [ a_theta1, a_theta2] = average_fit(a_avg_theta_LR1,a_avg_theta_LR2,n_pol_theta,ones(1,nr_of_avg)*(1/nr_of_avg),ones(1,nr_of_avg)*(1/nr_of_avg),mean(down_up));
    [ a_eta1, a_eta2] = average_fit(a_avg_eta_LR1,a_avg_eta_LR2,n_pol_eta,ones(1,nr_of_avg)*(1/nr_of_avg),ones(1,nr_of_avg)*(1/nr_of_avg),mean(down_up));
    [ a_phi1, a_phi2] = average_fit(a_avg_phi_LR1,a_avg_phi_LR2,n_pol_phi,ones(1,nr_of_avg)*(1/nr_of_avg),ones(1,nr_of_avg)*(1/nr_of_avg),mean(down_up));
    
    a_glob.theta1 = a_theta1;
    a_glob.theta2 = a_theta2;
    a_glob.eta1 = a_eta1;
    a_glob.eta2 = a_eta2;
    a_glob.phi1 = a_phi1;
    a_glob.phi2 = a_phi2;
    
    f_glob = mean(f_avg);
    
    down_up_glob = mean(down_up);

end

