function [a_sym, a_dev] = Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,down_up,down_up_avg)


    % Determine the deviation between the average left wingbeat and the
    % actual left wingbeats and the deviation between the average right
    % wingbeat and the actual right wingbeats:
    
    nr_wb = length(a_fit.theta_L1(1,:));
    
    a_dev = {};
        
    a_dev_theta_L1 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_theta_L2 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_eta_L1 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_eta_L2 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_phi_L1 = zeros(n_pol_phi+1,nr_wb);
    
    a_dev_phi_L2 = zeros(n_pol_phi+1,nr_wb);
    
    
    a_dev_theta_R1 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_theta_R2 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_eta_R1 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_eta_R2 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_phi_R1 = zeros(n_pol_phi+1,nr_wb);
    
    a_dev_phi_R2 = zeros(n_pol_phi+1,nr_wb);
    
    
    
    for i = 1:nr_wb
        
        a_dev_theta_L1(:,i) = a_fit.theta_L1(:,i)-a_avg.theta_L1;

        a_dev_theta_L2(:,i) = a_fit.theta_L2(:,i)-a_avg.theta_L2;

        a_dev_eta_L1(:,i) = a_fit.eta_L1(:,i)-a_avg.eta_L1;

        a_dev_eta_L2(:,i) = a_fit.eta_L2(:,i)-a_avg.eta_L2;

        a_dev_phi_L1(:,i) = a_fit.phi_L1(:,i)-a_avg.phi_L1;
        
        a_dev_phi_L2(:,i) = a_fit.phi_L2(:,i)-a_avg.phi_L2;


        a_dev_theta_R1(:,i) = a_fit.theta_R1(:,i)-a_avg.theta_R1;

        a_dev_theta_R2(:,i) = a_fit.theta_R2(:,i)-a_avg.theta_R2;

        a_dev_eta_R1(:,i) = a_fit.eta_R1(:,i)-a_avg.eta_R1;

        a_dev_eta_R2(:,i) = a_fit.eta_R2(:,i)-a_avg.eta_R2;

        a_dev_phi_R1(:,i) = a_fit.phi_R1(:,i)-a_avg.phi_R1;

        a_dev_phi_R2(:,i) = a_fit.phi_R2(:,i)-a_avg.phi_R2;   
        
    end
    
    a_dev.theta_L1 = a_dev_theta_L1;
    
    a_dev.theta_L2 = a_dev_theta_L2;
    
    a_dev.eta_L1 = a_dev_eta_L1;
    
    a_dev.eta_L2 = a_dev_eta_L2;
    
    a_dev.phi_L1 = a_dev_phi_L1;
    
    a_dev.phi_L2 = a_dev_phi_L2;
    
     
    a_dev.theta_R1 = a_dev_theta_R1;
    
    a_dev.theta_R2 = a_dev_theta_R2;
    
    a_dev.eta_R1 = a_dev_eta_R1;
    
    a_dev.eta_R2 = a_dev_eta_R2;
    
    a_dev.phi_R1 = a_dev_phi_R1;
    
    a_dev.phi_R2 = a_dev_phi_R2;
    
    
    
    % Add up the deviation to the average symmetric wingbeat to create a
    % 'symmetric' wingbeat:
    
    a_sym = {};
    
    a_sym_theta_L1 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_theta_L2 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_eta_L1 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_eta_L2 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_phi_L1 = zeros(n_pol_phi+1,nr_wb);
    
    a_sym_phi_L2 = zeros(n_pol_phi+1,nr_wb);
    
    
    a_sym_theta_R1 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_theta_R2 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_eta_R1 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_eta_R2 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_phi_R1 = zeros(n_pol_phi+1,nr_wb);
    
    a_sym_phi_R2 = zeros(n_pol_phi+1,nr_wb);    
    
    
    for i = 1:nr_wb
        
        a_sym_theta_L1(:,i) = a_avg.theta_LR1+a_dev_theta_L1(:,i);

        a_sym_theta_L2(:,i) = a_avg.theta_LR2+a_dev_theta_L2(:,i);

        a_sym_eta_L1(:,i) = a_avg.eta_LR1+a_dev_eta_L1(:,i);

        a_sym_eta_L2(:,i) = a_avg.eta_LR2+a_dev_eta_L2(:,i);

        a_sym_phi_L1(:,i) = a_avg.phi_LR1+a_dev_phi_L1(:,i);
        
        a_sym_phi_L2(:,i) = a_avg.phi_LR2+a_dev_phi_L2(:,i);


        a_sym_theta_R1(:,i) = a_avg.theta_LR1+a_dev_theta_L1(:,i);

        a_sym_theta_R2(:,i) = a_avg.theta_LR2+a_dev_theta_L2(:,i);

        a_sym_eta_R1(:,i) = a_avg.eta_LR1+a_dev_eta_L1(:,i);

        a_sym_eta_R2(:,i) = a_avg.eta_LR2+a_dev_eta_L2(:,i);

        a_sym_phi_R1(:,i) = a_avg.phi_LR1+a_dev_phi_L1(:,i);

        a_sym_phi_R2(:,i) = a_avg.phi_LR2+a_dev_phi_L2(:,i);   
        
    end    
    
    a_sym.theta_L1 = a_sym_theta_L1;
    
    a_sym.theta_L2 = a_sym_theta_L2;
    
    a_sym.eta_L1 = a_sym_eta_L1;
    
    a_sym.eta_L2 = a_sym_eta_L2;
    
    a_sym.phi_L1 = a_sym_phi_L1;
    
    a_sym.phi_L2 = a_sym_phi_L2;
    
    
    a_sym.theta_R1 = a_sym_theta_R1;
    
    a_sym.theta_R2 = a_sym_theta_R2;
    
    a_sym.eta_R1 = a_sym_eta_R1;
    
    a_sym.eta_R2 = a_sym_eta_R2;
    
    a_sym.phi_R1 = a_sym_phi_R1;
    
    a_sym.phi_R2 = a_sym_phi_R2;


end

