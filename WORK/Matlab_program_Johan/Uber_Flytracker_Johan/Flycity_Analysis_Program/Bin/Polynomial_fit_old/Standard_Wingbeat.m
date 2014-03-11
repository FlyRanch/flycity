function [ a_fit, a_avg, a_dev ] = Standard_Wingbeat( settings, pathDB, seq_nr )

    % Obtain local polynomial fit and sequence averaged fit:
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    start   = settings.start_stop(seq_nr,1);
    stop    = settings.start_stop(seq_nr,2);
    
    t = pathDB.t;
    dt = pathDB.dt;

    theta_L     = pathDB.wing_kin.theta_L(seq_nr,:);
    eta_L       = pathDB.wing_kin.eta_L(seq_nr,:);
    phi_L       = pathDB.wing_kin.phi_L(seq_nr,:);
    theta_R     = pathDB.wing_kin.theta_R(seq_nr,:);
    eta_R       = pathDB.wing_kin.eta_R(seq_nr,:);
    phi_R       = pathDB.wing_kin.phi_R(seq_nr,:);
    
    nr_wb = pathDB.wingbeats.nr_of_wb(seq_nr);
    
    wingbeat_loc    = pathDB.wingbeats.wingbeat_loc(1:nr_wb,:,seq_nr);
    downstroke_loc  = pathDB.wingbeats.downstroke_loc(1:nr_wb,:,seq_nr);
    upstroke_loc    = pathDB.wingbeats.upstroke_loc(1:nr_wb,:,seq_nr);
    
    a_fit.theta_L   = nan(2*(n_pol_theta+1),200);
    a_fit.eta_L     = nan(2*(n_pol_eta+1),200);
    a_fit.phi_L     = nan(2*(n_pol_phi+1),200);
    a_fit.theta_R   = nan(2*(n_pol_theta+1),200);
    a_fit.eta_R     = nan(2*(n_pol_eta+1),200);
    a_fit.phi_R     = nan(2*(n_pol_phi+1),200);
    
    wb_loc = nan(nr_wb*2,2);
    
    for k = 1:nr_wb
        
        wb_loc(2*k-1,:) = [ downstroke_loc(k,1) downstroke_loc(k,2)+1 ];
        wb_loc(2*k,:)   = [ upstroke_loc(k,1) upstroke_loc(k,2)+1 ];
        
    end
    
    a_theta_L = Piecewise_polynomial_fit(theta_L',n_pol_theta,wb_loc);
    a_eta_L = Piecewise_polynomial_fit(eta_L',n_pol_eta,wb_loc);
    a_phi_L = Piecewise_polynomial_fit(phi_L',n_pol_phi,wb_loc);
    a_theta_R = Piecewise_polynomial_fit(theta_R',n_pol_theta,wb_loc);
    a_eta_R = Piecewise_polynomial_fit(eta_R',n_pol_eta,wb_loc);
    a_phi_R = Piecewise_polynomial_fit(phi_R',n_pol_phi,wb_loc);
    
    a_fit.theta_L(:,1:nr_wb)   = [a_theta_L(:,1:2:end); a_theta_L(:,2:2:end)];
    a_fit.eta_L(:,1:nr_wb)     = [a_eta_L(:,1:2:end); a_eta_L(:,2:2:end)];
    a_fit.phi_L(:,1:nr_wb)     = [a_phi_L(:,1:2:end); a_phi_L(:,2:2:end)];
    a_fit.theta_R(:,1:nr_wb)   = [a_theta_R(:,1:2:end); a_theta_R(:,2:2:end)];
    a_fit.eta_R(:,1:nr_wb)     = [a_eta_R(:,1:2:end); a_eta_R(:,2:2:end)];
    a_fit.phi_R(:,1:nr_wb)     = [a_phi_R(:,1:2:end); a_phi_R(:,2:2:end)];
    
    % Determine downstroke upstroke ratio:
    
    a_fit.down_up = (downstroke_loc(:,2)-downstroke_loc(:,1)+1)./(wingbeat_loc(:,2)-wingbeat_loc(:,1)+1);
    
    % Determine the wingbeat frequency:
    
    a_fit.f = 1./((wingbeat_loc(:,2)-wingbeat_loc(:,1)+1)*dt);


    % Find the average wingbeat:
    
    % find the wingbeat which includes the trigger:
    
    trigger_wb = 0;
    
    seq_trigger_shift = settings.trigger_shift;
    
    trig_shift = sum(strcmp(char(settings.sequence_names(seq_nr)),seq_trigger_shift)==1);
                          
    if trig_shift ~= 0
        trigger_frame = 2683;
    else
        trigger_frame = 2795;
    end
    
    for j = 1:nr_wb
        
        if wingbeat_loc(j,1)<=trigger_frame && wingbeat_loc(j,2)>trigger_frame
            
            trigger_wb = j;
            
        elseif wingbeat_loc(j,2) == trigger_frame
        
            trigger_wb = j;
            
        end
    end
    
    % Find the average wingbeat:
    
    a_avg.theta_L    = nan((n_pol_theta+1)*2,1);
    a_avg.eta_L      = nan((n_pol_eta+1)*2,1);
    a_avg.phi_L      = nan((n_pol_phi+1)*2,1);
    
    a_avg.theta_R    = nan((n_pol_theta+1)*2,1);
    a_avg.eta_R      = nan((n_pol_eta+1)*2,1);
    a_avg.phi_R      = nan((n_pol_phi+1)*2,1);
    
    a_avg.theta_LR   = nan((n_pol_theta+1)*2,1);
    a_avg.eta_LR     = nan((n_pol_eta+1)*2,1);
    a_avg.phi_LR     = nan((n_pol_phi+1)*2,1);
    
    a_avg.down_up    = nan;
    a_avg.f          = nan;
    a_avg.trigger_wb = trigger_wb;
    
    a_dev.theta_L = nan(2*(n_pol_theta+1),200);
    a_dev.eta_L = nan(2*(n_pol_eta+1),200);
    a_dev.phi_L = nan(2*(n_pol_phi+1),200);
    
    a_dev.theta_R = nan(2*(n_pol_theta+1),200);
    a_dev.eta_R = nan(2*(n_pol_eta+1),200);
    a_dev.phi_R = nan(2*(n_pol_phi+1),200);
    
    if trigger_wb >= 5  % minimum of 3 wingbeats in order to determine the average
        
        a_avg.down_up   = mean(a_fit.down_up(3:trigger_wb));  % exclude the first 2 wingbeats from analysis
    
        a_avg.f         = mean(a_fit.f(3:trigger_wb));  % exclude the first 2 wingbeats from analysis
    
        
        % Determine weights for the weighted average:
        
        [ W_theta_L1, W_theta_L2, W_theta_R1, W_theta_R2, W_theta_LR1, W_theta_LR2 ] = weight(  a_fit.theta_L(1:(n_pol_theta+1),1:nr_wb), ...
                                                                                                a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),1:nr_wb), ...
                                                                                                a_fit.theta_R(1:(n_pol_theta+1),1:nr_wb), ... 
                                                                                                a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),1:nr_wb), ... 
                                                                                                n_pol_theta, trigger_wb);
        
        [ W_eta_L1, W_eta_L2, W_eta_R1, W_eta_R2, W_eta_LR1, W_eta_LR2 ] = weight(  a_fit.eta_L(1:(n_pol_eta+1),1:nr_wb), ... 
                                                                                    a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),1:nr_wb), ... 
                                                                                    a_fit.eta_R(1:(n_pol_eta+1),1:nr_wb), ... 
                                                                                    a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),1:nr_wb), ...
                                                                                    n_pol_eta, trigger_wb);
                                                                                
        [ W_phi_L1, W_phi_L2, W_phi_R1, W_phi_R2, W_phi_LR1, W_phi_LR2 ] = weight(  a_fit.phi_L(1:(n_pol_phi+1),1:nr_wb), ... 
                                                                                    a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),1:nr_wb), ...
                                                                                    a_fit.phi_R(1:(n_pol_phi+1),1:nr_wb), ...
                                                                                    a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),1:nr_wb), ...
                                                                                    n_pol_phi, trigger_wb);
        
        % Use function average_fit to obtain the average for the left and
        % right wing seperately:

        [a_avg_theta_L1, a_avg_theta_L2] = average_fit( a_fit.theta_L((1:(n_pol_theta+1)),3:trigger_wb), ...
                                                        a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb), ...
                                                        n_pol_theta,W_theta_L1,W_theta_L2,a_avg.down_up);
                                                    
        [a_avg_theta_R1, a_avg_theta_R2] = average_fit( a_fit.theta_R((1:(n_pol_theta+1)),3:trigger_wb), ...
                                                        a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb), ...
                                                        n_pol_theta,W_theta_R1,W_theta_R2,a_avg.down_up);
                                                    
        [a_avg_eta_L1, a_avg_eta_L2] = average_fit( a_fit.eta_L(1:(n_pol_eta+1),3:trigger_wb), ...
                                                    a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb), ...
                                                    n_pol_eta,W_eta_L1,W_eta_L2,a_avg.down_up);
                                                
        [a_avg_eta_R1, a_avg_eta_R2] = average_fit( a_fit.eta_R(1:(n_pol_eta+1),3:trigger_wb), ...
                                                    a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb), ...
                                                    n_pol_eta,W_eta_R1,W_eta_R2,a_avg.down_up);
                                                
        [a_avg_phi_L1, a_avg_phi_L2] = average_fit( a_fit.phi_L(1:(n_pol_phi+1),3:trigger_wb), ...
                                                    a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb), ...
                                                    n_pol_phi,W_phi_L1,W_phi_L2,a_avg.down_up);
                                                
        [a_avg_phi_R1, a_avg_phi_R2] = average_fit( a_fit.phi_R(1:(n_pol_phi+1),3:trigger_wb), ...
                                                    a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb), ...
                                                    n_pol_phi,W_phi_R1,W_phi_R2,a_avg.down_up);     

        % Use function average_fit to obtain the average for the left and
        % right wing combined:

        [a_avg_theta_LR1, a_avg_theta_LR2] = average_fit(   [a_fit.theta_L((1:(n_pol_theta+1)),3:trigger_wb) ...
                                                             a_fit.theta_R((1:(n_pol_theta+1)),3:trigger_wb)], ...
                                                            [a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb) ... 
                                                             a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb)], ...
                                                             n_pol_theta,W_theta_LR1,W_theta_LR2,a_avg.down_up);
                                                        
        [a_avg_eta_LR1, a_avg_eta_LR2] = average_fit(   [a_fit.eta_L((1:(n_pol_eta+1)),3:trigger_wb) ...
                                                         a_fit.eta_R((1:(n_pol_eta+1)),3:trigger_wb)], ...
                                                        [a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb) ... 
                                                         a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb)], ...
                                                         n_pol_eta,W_eta_LR1,W_eta_LR2,a_avg.down_up);
                                                         
        [a_avg_phi_LR1, a_avg_phi_LR2] = average_fit(   [a_fit.phi_L((1:(n_pol_phi+1)),3:trigger_wb) ...
                                                         a_fit.phi_R((1:(n_pol_phi+1)),3:trigger_wb)], ...
                                                        [a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb) ... 
                                                         a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb)], ...
                                                         n_pol_phi,W_phi_LR1,W_phi_LR2,a_avg.down_up);

        
        % Store average coefficients in a_avg
        
        a_avg.theta_L   = [a_avg_theta_L1; a_avg_theta_L2];
        a_avg.eta_L     = [a_avg_eta_L1; a_avg_eta_L2];
        a_avg.phi_L     = [a_avg_phi_L1; a_avg_phi_L2];

        a_avg.theta_R   = [a_avg_theta_R1; a_avg_theta_R2];
        a_avg.eta_R     = [a_avg_eta_R1; a_avg_eta_R2];
        a_avg.phi_R     = [a_avg_phi_R1; a_avg_phi_R2];     
        
        a_avg.theta_LR = [a_avg_theta_LR1; a_avg_theta_LR2];
        a_avg.eta_LR   = [a_avg_eta_LR1; a_avg_eta_LR2];
        a_avg.phi_LR   = [a_avg_phi_LR1; a_avg_phi_LR2];
        
        % Determine the deviation from the average wingbeat for the left and
        % right wings seperately:
        
        for k = 1:nr_wb

            a_dev.theta_L(:,k)   = a_fit.theta_L(:,k)-a_avg.theta_L;
            a_dev.eta_L(:,k)     = a_fit.eta_L(:,k)-a_avg.eta_L;
            a_dev.phi_L(:,k)     = a_fit.phi_L(:,k)-a_avg.phi_L;

            a_dev.theta_R(:,k)   = a_fit.theta_R(:,k)-a_avg.theta_R;
            a_dev.eta_R(:,k)     = a_fit.eta_R(:,k)-a_avg.eta_R;
            a_dev.phi_R(:,k)     = a_fit.phi_R(:,k)-a_avg.phi_R;
            
        end
                
    end    
    
        
    
%     figure()
%     hold on
%     plot(t,theta_L)
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_fit.theta_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,eta_L)
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_fit.eta_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,phi_L)
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_fit.phi_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,theta_R)
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_fit.theta_R(:,k),'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,eta_R)
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_fit.eta_R(:,k),'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,phi_R)
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_fit.phi_R(:,k),'g')
%     end
%     hold off
%     
%     
%     if isnan(a_avg.down_up)==0
%     
%     
%     figure()
%     hold on
%     plot(t,theta_L)
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_avg.theta_L,'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,eta_L)
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_avg.eta_L,'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,phi_L)
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_avg.phi_L,'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,theta_R)
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_avg.theta_R,'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,eta_R)
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_avg.eta_R,'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     plot(t,phi_R)
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_avg.phi_R,'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_dev.theta_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_dev.eta_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_dev.phi_L(:,k),'r')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_theta*a_dev.theta_R(:,k),'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_eta*a_dev.eta_R(:,k),'g')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_wb, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k),100,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
%     plot(t_wb,X_phi*a_dev.phi_R(:,k),'g')
%     end
%     hold off
%     
%     
%     end

end

