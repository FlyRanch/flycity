function maneuver_wb = select_maneuver( dev_data, trigger_wb, nr_wb, c_var )


    % Select maneuvering wingbeats based on their deviation from the
    % average wingbeat. 
    
    % The selection procedure is performed in the following way:
    %
    %   1: Determine the variance of the symmetric and asymmetric deviation
    %   of the left and right wings with respect to the average wingbeat
    %   before the stimulus starts.
    %
    %   2: Determine for each wingbeat if their absolute mean exceeds c_var
    %   times the variance of step 1.
    %
    %
    %----------------------------------------------------------------------

    
    % Step 1:
    
    
    % Symmetric:
    
    var_sym_theta = var((dev_data.theta_L(201:(100*trigger_wb))+dev_data.theta_R(201:(100*trigger_wb)))/2);
    
    var_sym_eta = var((dev_data.eta_L(201:(100*trigger_wb))+dev_data.eta_R(201:(100*trigger_wb)))/2);
    
    var_sym_phi = var((dev_data.phi_L(201:(100*trigger_wb))+dev_data.phi_R(201:(100*trigger_wb)))/2);

    
    % Asymmetric:
    
    var_asym_theta = var(dev_data.theta_L(201:(100*trigger_wb))-dev_data.theta_R(201:(100*trigger_wb)));
    
    var_asym_eta = var(dev_data.eta_L(201:(100*trigger_wb))-dev_data.eta_R(201:(100*trigger_wb)));
    
    var_asym_phi = var(dev_data.phi_L(201:(100*trigger_wb))-dev_data.phi_R(201:(100*trigger_wb)));
    
    % Step 2:
    
    sym_theta_mean = zeros(nr_wb,1);
    
    sym_eta_mean = zeros(nr_wb,1);
    
    sym_phi_mean = zeros(nr_wb,1);
    
    asym_theta_mean = zeros(nr_wb,1);
    
    asym_eta_mean = zeros(nr_wb,1);
    
    asym_phi_mean = zeros(nr_wb,1);
    
    for k = 1:nr_wb
        
        sym_theta_mean(k) = mean((dev_data.theta_L(((k-1)*100+1):(k*100))+dev_data.theta_R(((k-1)*100+1):(k*100)))/2);
        
        sym_eta_mean(k) = mean((dev_data.eta_L(((k-1)*100+1):(k*100))+dev_data.eta_R(((k-1)*100+1):(k*100)))/2);
        
        sym_phi_mean(k) = mean((dev_data.phi_L(((k-1)*100+1):(k*100))+dev_data.phi_R(((k-1)*100+1):(k*100)))/2);
        
        asym_theta_mean(k) = mean(dev_data.theta_L(((k-1)*100+1):(k*100))-dev_data.theta_R(((k-1)*100+1):(k*100)));
        
        asym_eta_mean(k) = mean(dev_data.eta_L(((k-1)*100+1):(k*100))-dev_data.eta_R(((k-1)*100+1):(k*100)));
        
        asym_phi_mean(k) = mean(dev_data.phi_L(((k-1)*100+1):(k*100))-dev_data.phi_R(((k-1)*100+1):(k*100)));
        
    end
    
    Man_wb = zeros(nr_wb,1);
    
    Man_wb_sym_theta = zeros(nr_wb,1);
    
    Man_wb_sym_eta = zeros(nr_wb,1);
    
    Man_wb_sym_phi = zeros(nr_wb,1);
    
    Man_wb_asym_theta = zeros(nr_wb,1);
    
    Man_wb_asym_eta = zeros(nr_wb,1);
    
    Man_wb_asym_phi = zeros(nr_wb,1);
    
    for k = 1:nr_wb

        if abs(sym_theta_mean(k)) > (c_var*sqrt(var_sym_theta))
            
            Man_wb(k) = 1;
            
            Man_wb_sym_theta(k) = 1;
            
        elseif abs(sym_eta_mean(k)) > (c_var*sqrt(var_sym_eta))
            
            Man_wb(k) = 1;
            
            Man_wb_sym_eta(k) = 1;
            
        elseif abs(sym_phi_mean(k)) > (c_var*sqrt(var_sym_phi))
            
            Man_wb(k) = 1;
            
            Man_wb_sym_phi(k) = 1;
            
        elseif abs(asym_theta_mean(k)) > (c_var*sqrt(var_asym_theta))
            
            Man_wb(k) = 1;
            
            Man_wb_asym_theta(k) = 1;
            
        elseif abs(asym_eta_mean(k)) > (c_var*sqrt(var_asym_eta))
            
            Man_wb(k) = 1;
            
            Man_wb_asym_eta(k) = 1;
            
        elseif abs(asym_phi_mean(k)) > (c_var*sqrt(var_asym_phi))
            
            Man_wb(k) = 1;
            
            Man_wb_asym_phi(k) = 1;
            
        end
        
    end
    
    maneuver_wb = {};
    
    
    maneuver_wb.id = find(Man_wb);

    maneuver_wb.sym_theta_id = find(Man_wb_sym_theta);
    
    maneuver_wb.sym_eta_id = find(Man_wb_sym_eta);
    
    maneuver_wb.sym_phi_id = find(Man_wb_sym_phi);
    
    maneuver_wb.asym_theta_id = find(Man_wb_asym_theta);
    
    maneuver_wb.asym_eta_id = find(Man_wb_asym_eta);
    
    maneuver_wb.asym_phi_id = find(Man_wb_asym_phi);
    
    
    
    
   
    
%     figure()
%     plot(radtodeg(sym_theta_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_sym_theta*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_sym_theta*ones(nr_wb,1)).^0.5),'r')
%     hold off
%     
%     figure()
%     plot(radtodeg(sym_eta_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_sym_eta*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_sym_eta*ones(nr_wb,1)).^0.5),'r')
%     hold off
%     
%     figure()
%     plot(radtodeg(sym_phi_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_sym_phi*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_sym_phi*ones(nr_wb,1)).^0.5),'r')
%     hold off
%     
%     figure()
%     plot(radtodeg(asym_theta_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_asym_theta*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_asym_theta*ones(nr_wb,1)).^0.5),'r')
%     
%     figure()
%     plot(radtodeg(asym_eta_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_asym_eta*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_asym_eta*ones(nr_wb,1)).^0.5),'r')
%     hold off
%     
%     figure()
%     plot(radtodeg(asym_phi_mean),'o')
%     hold on
%     plot(radtodeg(c_var*(var_asym_phi*ones(nr_wb,1)).^0.5),'r')
%     plot(radtodeg(-c_var*(var_asym_phi*ones(nr_wb,1)).^0.5),'r')
%     hold off
%     
%     
%     
%     
%     
%     
%     pause
%     
%     close all


    


end

