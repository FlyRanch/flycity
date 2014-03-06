function [a_fit,a_avg,f_avg,down_up,trigger_wb,down_up_ratio] = Standard_wingbeat( settings, pathDB, seq_nr, order_pol_theta, order_pol_eta, order_pol_phi )


    % Create a standard wingbeat:
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    t = pathDB.t(start:stop);
    
    dt = pathDB.t(2)-pathDB.t(1);

    eta_L = pathDB.eta_L(start:stop,seq_nr);

    theta_L = pathDB.theta_L(start:stop,seq_nr);

    phi_L = pathDB.phi_L(start:stop,seq_nr);

    eta_R = pathDB.eta_R(start:stop,seq_nr);

    theta_R = pathDB.theta_R(start:stop,seq_nr);

    phi_R = pathDB.phi_R(start:stop,seq_nr);
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    
    
    % Wingbeat locations

    wb_loc = zeros(nr_wb,2);
    
    for k = 1:nr_wb

        if k == nr_wb

            wb_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );

            wb_loc(k,:) = [ pathDB.wingbeat_time(nr_wb,1,seq_nr) pathDB.wingbeat_time(nr_wb,wb_end,seq_nr) ];

        else

            wb_loc(k,:) = [ pathDB.wingbeat_time(k,1,seq_nr) pathDB.wingbeat_time(k+1,1,seq_nr) ];

        end

    end


    wb_loc12 = zeros(2*nr_wb,2);

    for k = 1:nr_wb
       
        down_start = pathDB.down_time_L(k,1,seq_nr);
        
        up_start = pathDB.up_time_L(k,1,seq_nr);
        
        up_end = pathDB.up_time_L(k,1,seq_nr)+find(isnan(pathDB.up_time_L(k,:,seq_nr))==0,1,'last');
        
        
        wb_loc12(2*k-1,:) = [down_start up_start];
        
        wb_loc12(2*k,:) = [up_start up_end];

    end
            

    
    n_pol_theta = order_pol_theta;
    n_pol_eta = order_pol_eta;
    n_pol_phi = order_pol_phi;
    
    
    a_fit_theta_L = Piecewise_polynomial_fit(theta_L,n_pol_theta,wb_loc12);
    
    a_fit_eta_L = Piecewise_polynomial_fit(eta_L,n_pol_eta,wb_loc12);

    a_fit_phi_L  = Piecewise_polynomial_fit(phi_L,n_pol_phi,wb_loc12);

    
    a_fit_theta_L1 = a_fit_theta_L(:,1:2:(end-1));

    a_fit_eta_L1 = a_fit_eta_L(:,1:2:(end-1));

    a_fit_phi_L1 = a_fit_phi_L(:,1:2:(end-1));
    
    
    a_fit_theta_L2 = a_fit_theta_L(:,2:2:end);

    a_fit_eta_L2 = a_fit_eta_L(:,2:2:end);

    a_fit_phi_L2 = a_fit_phi_L(:,2:2:end);
    
    
    a_fit_theta_R = Piecewise_polynomial_fit(theta_R,n_pol_theta,wb_loc12);

    a_fit_eta_R = Piecewise_polynomial_fit(eta_R,n_pol_eta,wb_loc12);

    a_fit_phi_R = Piecewise_polynomial_fit(phi_R,n_pol_phi,wb_loc12);
    
    
    a_fit_theta_R1 = a_fit_theta_R(:,1:2:(end-1));

    a_fit_eta_R1 = a_fit_eta_R(:,1:2:(end-1));

    a_fit_phi_R1 = a_fit_phi_R(:,1:2:(end-1));
    
    
    a_fit_theta_R2 = a_fit_theta_R(:,2:2:end);

    a_fit_eta_R2 = a_fit_eta_R(:,2:2:end);

    a_fit_phi_R2 = a_fit_phi_R(:,2:2:end);    
    
    
    a_fit_f = 1./((wb_loc(:,2)-wb_loc(:,1)-1)*dt);
    
       
    
    
    % Store data in structure a_fit:
    
    a_fit = {};
    
    a_avg = {};
    
    a_fit.theta_L1 = a_fit_theta_L1;
    a_fit.theta_L2 = a_fit_theta_L2;
    a_fit.eta_L1 = a_fit_eta_L1;
    a_fit.eta_L2 = a_fit_eta_L2;
    a_fit.phi_L1 = a_fit_phi_L1;
    a_fit.phi_L2 = a_fit_phi_L2;

    a_fit.theta_R1 = a_fit_theta_R1;
    a_fit.theta_R2 = a_fit_theta_R2;
    a_fit.eta_R1 = a_fit_eta_R1;
    a_fit.eta_R2 = a_fit_eta_R2;
    a_fit.phi_R1 = a_fit_phi_R1;
    a_fit.phi_R2 = a_fit_phi_R2;
    
    a_fit.wb_loc_12 = wb_loc12;
    
    a_fit.f = a_fit_f;

    
    % find the wingbeat which includes the trigger:
    
    trigger_wb = 0;
    
    seq_trigger_shift = { '20130205_S0004'; ...
                          '20130205_S0008'; ...
                          '20130206_S0003'; ...
                          '20130206_S0006'; ...
                          '20130206_S0007'; ...
                          '20130206_S0008'; ...
                          '20130207_S0002'; ...
                          '20130208_S0005' };
    
    trig_shift = sum(strcmp(char(settings.sequence_names(seq_nr)),seq_trigger_shift)==1);
                          
    if trig_shift ~= 0
        
        trigger_frame = 2683-start+1;
    
    else
        
        trigger_frame = 2795-start+1;
        
    end
    
    for j = 1:nr_wb
        
        if wb_loc(j,1)<trigger_frame && wb_loc(j,2)>trigger_frame
            
            trigger_wb = j;
            
        elseif wb_loc(j,2) == trigger_frame
        
            trigger_wb = j;
            
        end
    end
    
    f_avg =0;
    
    down_up = 0;
    
    % Determine average ratio between down and upstroke:

    down_up_ratio = zeros(nr_wb,1);

        for m = 1:nr_wb
            
                down_length_L = length(find(isnan(pathDB.down_time_L(m,:,seq_nr))==0));

                up_length_L = length(find(isnan(pathDB.up_time_L(m,:,seq_nr))==0));

                down_length_R = length(find(isnan(pathDB.down_time_R(m,:,seq_nr))==0));

                up_length_R = length(find(isnan(pathDB.up_time_R(m,:,seq_nr))==0));

                down_up_ratio(m) = 0.5*((down_length_L/(down_length_L+up_length_L))+(down_length_R/(down_length_R+up_length_R)));

        end
    
    
    % Find the average wingbeat:
    
    if trigger_wb > 3
        
        [ W_theta_L1, W_theta_L2, W_theta_R1, W_theta_R2, W_theta_LR1, W_theta_LR2 ] = weight(a_fit.theta_L1, a_fit.theta_L2, a_fit.theta_R1, a_fit.theta_R2, n_pol_theta, trigger_wb);
        
        [ W_eta_L1, W_eta_L2, W_eta_R1, W_eta_R2, W_eta_LR1, W_eta_LR2 ] = weight(a_fit.eta_L1, a_fit.eta_L2, a_fit.eta_R1, a_fit.eta_R2, n_pol_eta, trigger_wb);
        
        [ W_phi_L1, W_phi_L2, W_phi_R1, W_phi_R2, W_phi_LR1, W_phi_LR2 ] = weight(a_fit.phi_L1, a_fit.phi_L2, a_fit.phi_R1, a_fit.phi_R2, n_pol_phi, trigger_wb);
        

        down_up = mean(down_up_ratio(3:trigger_wb));


        [a_avg_theta_L1, a_avg_theta_L2] = average_fit(a_fit.theta_L1(:,3:trigger_wb),a_fit.theta_L2(:,3:trigger_wb),n_pol_theta,W_theta_L1,W_theta_L2,down_up);
        
        [a_avg_theta_R1, a_avg_theta_R2] = average_fit(a_fit.theta_R1(:,3:trigger_wb),a_fit.theta_R2(:,3:trigger_wb),n_pol_theta,W_theta_R1,W_theta_R2,down_up);
        
        [a_avg_eta_L1, a_avg_eta_L2] = average_fit(a_fit.eta_L1(:,3:trigger_wb),a_fit.eta_L2(:,3:trigger_wb),n_pol_eta,W_eta_L1,W_eta_L2,down_up);
        
        [a_avg_eta_R1, a_avg_eta_R2] = average_fit(a_fit.eta_R1(:,3:trigger_wb),a_fit.eta_R2(:,3:trigger_wb),n_pol_eta,W_eta_R1,W_eta_R2,down_up);
        
        [a_avg_phi_L1, a_avg_phi_L2] = average_fit(a_fit.phi_L1(:,3:trigger_wb),a_fit.phi_L2(:,3:trigger_wb),n_pol_phi,W_phi_L1,W_phi_L2,down_up);
        
        [a_avg_phi_R1, a_avg_phi_R2] = average_fit(a_fit.phi_R1(:,3:trigger_wb),a_fit.phi_R2(:,3:trigger_wb),n_pol_phi,W_phi_R1,W_phi_R2,down_up);     


        [a_avg_theta_LR1, a_avg_theta_LR2] = average_fit([a_fit.theta_L1(:,3:trigger_wb) a_fit.theta_R1(:,3:trigger_wb)],[a_fit.theta_L2(:,3:trigger_wb) a_fit.theta_R2(:,3:trigger_wb)],n_pol_theta,W_theta_LR1,W_theta_LR2,down_up);

        [a_avg_eta_LR1, a_avg_eta_LR2] = average_fit([a_fit.eta_L1(:,3:trigger_wb) a_fit.eta_R1(:,3:trigger_wb)],[a_fit.eta_L2(:,3:trigger_wb) a_fit.eta_R2(:,3:trigger_wb)],n_pol_eta,W_eta_LR1,W_eta_LR2,down_up);

        [a_avg_phi_LR1, a_avg_phi_LR2] = average_fit([a_fit.phi_L1(:,3:trigger_wb) a_fit.phi_R1(:,3:trigger_wb)],[a_fit.phi_L2(:,3:trigger_wb) a_fit.phi_R2(:,3:trigger_wb)],n_pol_phi,W_phi_LR1,W_phi_LR2,down_up);   

        
        % Store average coefficients in a_avg
        
        a_avg.theta_L1 = a_avg_theta_L1;
        a_avg.theta_L2 = a_avg_theta_L2;
        a_avg.eta_L1 = a_avg_eta_L1;
        a_avg.eta_L2 = a_avg_eta_L2;
        a_avg.phi_L1 = a_avg_phi_L1;
        a_avg.phi_L2 = a_avg_phi_L2;
        

        a_avg.theta_R1 = a_avg_theta_R1;
        a_avg.theta_R2 = a_avg_theta_R2;
        a_avg.eta_R1 = a_avg_eta_R1;
        a_avg.eta_R2 = a_avg_eta_R2;
        a_avg.phi_R1 = a_avg_phi_R1;
        a_avg.phi_R2 = a_avg_phi_R2;        
        
        
        a_avg.theta_LR1 = a_avg_theta_LR1;
        a_avg.theta_LR2 = a_avg_theta_LR2;
        a_avg.eta_LR1 = a_avg_eta_LR1;
        a_avg.eta_LR2 = a_avg_eta_LR2;
        a_avg.phi_LR1 = a_avg_phi_LR1;
        a_avg.phi_LR2 = a_avg_phi_LR2;  
        
        
        
        % Determine average wingbeat frequency:

        f_avg = 1/(((wb_loc(trigger_wb,2)-wb_loc(3,1))*dt)/(trigger_wb-2));
  
                
    end
    
end