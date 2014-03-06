function [a_fit_L,a_fit_R,a_avg_L,a_avg_R,a_avg_LR,f_avg,down_up,trigger_wb] = Standard_wingbeat( settings, pathDB, seq_nr, order_pol_theta, order_pol_eta, order_pol_phi )



    a_fit_L = {};
    
    a_fit_R =  {};
    
    a_avg_L =  {};

    a_avg_R =  {};
    
    a_avg_LR =  {};
    
    f_avg = 0;
    
    

    % Program to determine the average steady-state wingbeat:

    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
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

    n_pol_theta = order_pol_theta;
    n_pol_eta = order_pol_eta;
    n_pol_phi = order_pol_phi;

    a_fit_theta_L = Piecewise_polynomial_fit(theta_L,n_pol_theta,wb_loc);

    a_fit_eta_L = Piecewise_polynomial_fit(eta_L,n_pol_eta,wb_loc);

    a_fit_phi_L = Piecewise_polynomial_fit(phi_L,n_pol_phi,wb_loc);


    a_fit_theta_R = Piecewise_polynomial_fit(theta_R,n_pol_theta,wb_loc);

    a_fit_eta_R = Piecewise_polynomial_fit(eta_R,n_pol_eta,wb_loc);

    a_fit_phi_R = Piecewise_polynomial_fit(phi_R,n_pol_phi,wb_loc);
    
    
    
    a_fit_L.theta = a_fit_theta_L;
    
    a_fit_L.eta = a_fit_eta_L;
    
    a_fit_L.phi = a_fit_phi_L;
    
    a_fit_R.theta = a_fit_theta_R;
    
    a_fit_R.eta = a_fit_eta_R;
    
    a_fit_R.phi = a_fit_phi_R;
%    
%     
%     
%     a_fit = Piecewise_Bernstein_polynomial_fit(theta_L,100,n_pol_theta,wb_loc,a_fit_theta_L)
    

    % Determine the wingbeat average
    
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
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
    if trigger_wb > 3
    
    % Find for the variance between every wingbeat in the interval:
    % wb_nr==3 until wb_nr==trigger_frame:
    
    var_theta_L = zeros(trigger_wb-2,trigger_wb-2,n_pol_theta+1);
    
    var_theta_R = zeros(trigger_wb-2,trigger_wb-2,n_pol_theta+1);
    
    var_eta_L = zeros(trigger_wb-2,trigger_wb-2,n_pol_eta+1);
    
    var_eta_R = zeros(trigger_wb-2,trigger_wb-2,n_pol_eta+1);
    
    var_phi_L = zeros(trigger_wb-2,trigger_wb-2,n_pol_phi+1);
    
    var_phi_R = zeros(trigger_wb-2,trigger_wb-2,n_pol_phi+1);
    
    var_theta_LR = zeros(trigger_wb-2,trigger_wb-2,n_pol_theta+1);
      
    var_eta_LR = zeros(trigger_wb-2,trigger_wb-2,n_pol_eta+1);
    
    var_phi_LR = zeros(trigger_wb-2,trigger_wb-2,n_pol_phi+1);

    
    for m = 1:(trigger_wb-2)
            
        for n = 1:(trigger_wb-2)
            
            if m ~= n
            
                var_theta_L(m,n,:) = var([a_fit_theta_L(:,m+2)' ; a_fit_theta_L(:,n+2)']);

                var_theta_R(m,n,:) =  var([a_fit_theta_R(:,m+2)' ; a_fit_theta_R(:,n+2)']);

                var_eta_L(m,n,:) = var([a_fit_eta_L(:,m+2)' ; a_fit_eta_L(:,n+2)']);

                var_eta_R(m,n,:) = var([a_fit_eta_R(:,m+2)' ; a_fit_eta_R(:,n+2)']);

                var_phi_L(m,n,:) = var([a_fit_phi_L(:,m+2)' ; a_fit_phi_L(:,n+2)']);

                var_phi_R(m,n,:) = var([a_fit_phi_R(:,m+2)' ; a_fit_phi_R(:,n+2)']);
                
                var_theta_LR(m,n,:) = var([a_fit_theta_L(:,m+2)' ; a_fit_theta_R(:,n+2)']);
    
                var_eta_LR(m,n,:) = var([a_fit_eta_L(:,m+2)' ; a_fit_eta_R(:,n+2)']);
    
                var_phi_LR(m,n,:) = var([a_fit_phi_L(:,m+2)' ; a_fit_phi_R(:,n+2)']);
            
            end
            
        end
        
    end    
    
    weight_theta_L = zeros(n_pol_theta+1,trigger_wb-2);
    
    weight_theta_R = zeros(n_pol_theta+1,trigger_wb-2);
    
    for i = 1:(n_pol_theta+1)
        
        weight_theta_L(i,:) = sum(var_theta_L(:,:,i));
        
        weight_theta_R(i,:) = sum(var_theta_R(:,:,i));
        
    end
    
    W_theta_L = (1/sum(1./sum(weight_theta_L.^2)))*(1./sum(weight_theta_L.^2));
    
    W_theta_R = (1/sum(1./sum(weight_theta_R.^2)))*(1./sum(weight_theta_R.^2));
    
    
    weight_eta_L = zeros(n_pol_eta+1,trigger_wb-2);
    
    weight_eta_R = zeros(n_pol_eta+1,trigger_wb-2);
    
    for i = 1:(n_pol_eta+1)
        
        weight_eta_L(i,:) = sum(var_eta_L(:,:,i));
        
        weight_eta_R(i,:) = sum(var_eta_R(:,:,i));
        
    end    
    
    W_eta_L = (1/sum(1./sum(weight_eta_L.^2)))*(1./sum(weight_eta_L.^2));
    
    W_eta_R = (1/sum(1./sum(weight_eta_R.^2)))*(1./sum(weight_eta_R.^2));
    
    
    
    weight_phi_L = zeros(n_pol_phi+1,trigger_wb-2);
    
    weight_phi_R = zeros(n_pol_phi+1,trigger_wb-2);
    
    for i = 1:(n_pol_phi+1)
        
        weight_phi_L(i,:) = sum(var_phi_L(:,:,i));
        
        weight_phi_R(i,:) = sum(var_phi_R(:,:,i));
        
    end        

    W_phi_L = (1/sum(1./sum(weight_phi_L.^2)))*(1./sum(weight_phi_L.^2));
    
    W_phi_R = (1/sum(1./sum(weight_phi_R.^2)))*(1./sum(weight_phi_R.^2));
    
  
    
    wb_loc_mean = zeros(trigger_wb-2,2);
    
    for k = 1:(trigger_wb-2)
       
        if k == (trigger_wb-2)
            
            wb_end = find(isnan(pathDB.wingbeat_time(trigger_wb,:,seq_nr))==0, 1, 'last' );
        
            wb_loc_mean(k,:) = [ pathDB.wingbeat_time(trigger_wb,1,seq_nr) pathDB.wingbeat_time(trigger_wb,wb_end,seq_nr) ];
            
        else
            
            wb_loc_mean(k,:) = [ pathDB.wingbeat_time(k+2,1,seq_nr) pathDB.wingbeat_time(k+3,1,seq_nr) ];
            
        end
        
    end
    
    Y_theta_L = theta_L;
    
    a_avg_theta_L = average_fit(Y_theta_L,n_pol_theta,wb_loc_mean,W_theta_L);
    
    Y_theta_R = theta_R;
    
    a_avg_theta_R = average_fit(Y_theta_R,n_pol_theta,wb_loc_mean,W_theta_R);
    
    
    
    Y_eta_L = eta_L;
    
    a_avg_eta_L = average_fit(Y_eta_L,n_pol_eta,wb_loc_mean,W_eta_L);
    
    Y_eta_R = eta_R;
    
    a_avg_eta_R = average_fit(Y_eta_R,n_pol_eta,wb_loc_mean,W_eta_R);
    
    
    
    Y_phi_L = phi_L;
   
    a_avg_phi_L = average_fit(Y_phi_L,n_pol_phi,wb_loc_mean,W_phi_L);
    
    Y_phi_R = phi_R;
    
    a_avg_phi_R = average_fit(Y_phi_R,n_pol_phi,wb_loc_mean,W_phi_R);
    
    
    
    
    
    % Create average wingbeat----------------------------------------------

    
    
    weight_theta_LR = zeros(n_pol_theta+1,2*(trigger_wb-2));
    
    for i = 1:(n_pol_theta+1)
        
        weight_theta_LR(i,:) = sum([var_theta_L(:,:,i) var_theta_LR(:,:,i); var_theta_LR(:,:,i) var_theta_R(:,:,i)]);
        
    end
    
    W_theta_LR = (1/sum(1./sum(weight_theta_LR.^2)))*(1./sum(weight_theta_LR.^2));
    
    
    weight_eta_LR = zeros(n_pol_eta+1,2*(trigger_wb-2));
    
    for i = 1:(n_pol_eta+1)
        
        weight_eta_LR(i,:) = sum([var_eta_L(:,:,i) var_eta_LR(:,:,i); var_eta_LR(:,:,i) var_eta_R(:,:,i)]);
        
    end
    
    W_eta_LR = (1/sum(1./sum(weight_eta_LR.^2)))*(1./sum(weight_eta_LR.^2)); 
    
    
    weight_phi_LR = zeros(n_pol_phi+1,2*(trigger_wb-2));
    
    for i = 1:(n_pol_phi+1)
        
        weight_phi_LR(i,:) = sum([var_phi_L(:,:,i) var_phi_LR(:,:,i); var_phi_LR(:,:,i) var_phi_R(:,:,i)]);
        
    end
    
    W_phi_LR = (1/sum(1./sum(weight_phi_LR.^2)))*(1./sum(weight_phi_LR.^2));
    
    
    Y_theta_LR = [theta_L(wb_loc_mean(1,1):wb_loc_mean(end,2)); theta_R(wb_loc_mean(1,1):wb_loc_mean(end,2))];
       
    Y_eta_LR = [eta_L(wb_loc_mean(1,1):wb_loc_mean(end,2)); eta_R(wb_loc_mean(1,1):wb_loc_mean(end,2))];
    
    Y_phi_LR = [phi_L(wb_loc_mean(1,1):wb_loc_mean(end,2)); phi_R(wb_loc_mean(1,1):wb_loc_mean(end,2))];
    
 
    wb_loc_LR = zeros(2*(trigger_wb-2),2);
    
    wb_loc_LR = [ [wb_loc_mean(:,1)-wb_loc_mean(1,1)+1; wb_loc_mean(:,1)-2*wb_loc_mean(1,1)+2+wb_loc_mean(end,2)] [wb_loc_mean(:,2)-wb_loc_mean(1,1)+1; wb_loc_mean(:,2)-2*wb_loc_mean(1,1)+2+wb_loc_mean(end,2)] ];
    

    a_avg_theta_LR = average_fit(Y_theta_LR,n_pol_theta,wb_loc_LR,W_theta_LR);

    a_avg_eta_LR = average_fit(Y_eta_LR,n_pol_eta,wb_loc_LR,W_eta_LR);
    
    a_avg_phi_LR = average_fit(Y_phi_LR,n_pol_phi,wb_loc_LR,W_phi_LR);
    
    
    a_avg_L.theta = a_avg_theta_L;
    
    a_avg_L.eta = a_avg_eta_L;
    
    a_avg_L.phi = a_avg_phi_L;
    
    
    a_avg_R.theta = a_avg_theta_R;
    
    a_avg_R.eta = a_avg_eta_R;
    
    a_avg_R.phi = a_avg_phi_R;
    
    
    a_avg_LR.theta = a_avg_theta_LR;
    
    a_avg_LR.eta = a_avg_eta_LR;
    
    a_avg_LR.phi = a_avg_phi_LR;
    
    
    % Determine average wingbeat frequency:
    
    f_avg = 1/(((wb_loc_mean(end,1)-wb_loc_mean(1,1))*dt)/length(wb_loc_mean(:,1)));
    
    
    % Determine average ratio between down and upstroke:
    
    down_up_ratio = zeros(trigger_wb-2,1);
    
    for m = 1:(trigger_wb-2)
        
        if pathDB.L_wingbeat_loc(m+2,2,seq_nr) > pathDB.L_wingbeat_loc(m+2,1,seq_nr);

            down_length_L = pathDB.L_wingbeat_loc(m+3,1,seq_nr)-pathDB.L_wingbeat_loc(m+2,2,seq_nr);

            up_length_L = pathDB.L_wingbeat_loc(m+3,2,seq_nr)-pathDB.L_wingbeat_loc(m+3,1,seq_nr);

            down_length_R = pathDB.R_wingbeat_loc(m+3,1,seq_nr)-pathDB.R_wingbeat_loc(m+2,2,seq_nr);

            up_length_R = pathDB.R_wingbeat_loc(m+3,2,seq_nr)-pathDB.R_wingbeat_loc(m+3,1,seq_nr);

            down_up_ratio(m) = 0.5*((down_length_L/up_length_L)+(down_length_R/up_length_R));
        
        elseif pathDB.L_wingbeat_loc(m+2,2,seq_nr) < pathDB.L_wingbeat_loc(m+2,1,seq_nr)
            
            down_length_L = pathDB.L_wingbeat_loc(m+2,1,seq_nr)-pathDB.L_wingbeat_loc(m+2,2,seq_nr);

            up_length_L = pathDB.L_wingbeat_loc(m+3,2,seq_nr)-pathDB.L_wingbeat_loc(m+2,1,seq_nr);

            down_length_R = pathDB.R_wingbeat_loc(m+2,1,seq_nr)-pathDB.R_wingbeat_loc(m+2,2,seq_nr);

            up_length_R = pathDB.R_wingbeat_loc(m+3,2,seq_nr)-pathDB.R_wingbeat_loc(m+2,1,seq_nr);

            down_up_ratio(m) = 0.5*((down_length_L/up_length_L)+(down_length_R/up_length_R));
            
        end
    
    end
    
    down_up = mean(down_up_ratio);
    
    end
    
    t = -1:0.01:1;
    
    PN_theta = Legendre_polynomial(n_pol_theta,1,t);
    PN_eta = Legendre_polynomial(n_pol_eta,1,t);
    PN_phi = Legendre_polynomial(n_pol_phi,1,t);
    

    
    
%     figure()
%     plot(radtodeg(PN_theta(:,:,1)'*a_avg_theta_L),'r')
%     hold on
%     plot(radtodeg(PN_theta(:,:,1)'*a_avg_theta_R),'g')
%     plot(radtodeg(PN_theta(:,:,1)'*a_avg_theta_LR),'b')
%     hold off
% 
%     figure()
%     plot(radtodeg(PN_eta(:,:,1)'*a_avg_eta_L),'r')
%     hold on
%     plot(radtodeg(PN_eta(:,:,1)'*a_avg_eta_R),'g')
%     plot(radtodeg(PN_eta(:,:,1)'*a_avg_eta_LR),'b')
%     hold off
%     
%     figure()
%     plot(radtodeg(PN_phi(:,:,1)'*a_avg_phi_L),'r')
%     hold on
%     plot(radtodeg(PN_phi(:,:,1)'*a_avg_phi_R),'g')
%     plot(radtodeg(PN_phi(:,:,1)'*a_avg_phi_LR),'b')
%     hold off
%     
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_theta(:,:,1)'*a_fit_theta_L(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_theta(:,:,1)'*a_avg_theta_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_theta(:,:,1)'*a_fit_theta_R(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_theta(:,:,1)'*a_avg_theta_R),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_eta(:,:,1)'*a_fit_eta_L(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_eta(:,:,1)'*a_avg_eta_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_eta(:,:,1)'*a_fit_eta_R(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_eta(:,:,1)'*a_avg_eta_R),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_phi(:,:,1)'*a_fit_phi_L(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_phi(:,:,1)'*a_avg_phi_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(PN_phi(:,:,1)'*a_fit_phi_R(:,k)),'Color',[0.5 0.5 0.5])
%     end
%     plot(radtodeg(PN_phi(:,:,1)'*a_avg_phi_R),'r')
%     hold off
%     
end

