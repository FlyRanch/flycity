function [a_fit,a_avg,f_avg,down_up,trigger_wb,ratio_1_2,ratio_1_2_avg] = Standard_wingbeat( settings, pathDB, seq_nr, order_pol_theta, order_pol_eta, order_pol_phi )


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
    
    wb_loc12 = zeros(2*nr_wb,2);

    for k = 1:nr_wb

        if k == nr_wb

            wb_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );

            wb_loc(k,:) = [ pathDB.wingbeat_time(nr_wb,1,seq_nr) pathDB.wingbeat_time(nr_wb,wb_end,seq_nr) ];

        else

            wb_loc(k,:) = [ pathDB.wingbeat_time(k,1,seq_nr) pathDB.wingbeat_time(k+1,1,seq_nr) ];

        end

    end

    for k = 1:nr_wb

        if k == nr_wb

            wb_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );
            
            frame_half = round((pathDB.wingbeat_time(nr_wb,wb_end,seq_nr)+pathDB.wingbeat_time(nr_wb,1,seq_nr))/2);

            wb_loc1 = [ pathDB.wingbeat_time(nr_wb,1,seq_nr) frame_half ];
            
            wb_loc2 = [frame_half pathDB.wingbeat_time(nr_wb,wb_end,seq_nr)];

        else

           
            frame_half = round((pathDB.wingbeat_time(k,1,seq_nr)+pathDB.wingbeat_time(k+1,1,seq_nr))/2);

            wb_loc1 = [ pathDB.wingbeat_time(k,1,seq_nr) frame_half ];
            
            wb_loc2 = [frame_half pathDB.wingbeat_time(k+1,1,seq_nr)];

        end
        
        wb_loc12(2*k-1,:) = wb_loc1;
        
        wb_loc12(2*k,:) = wb_loc2;

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
    
    down_up_ratio = 0;
    
    down_up = 0;
    
    ratio_1_2 = 0;
    
    ratio_1_2_avg = 0;
    
    
    % Find the average wingbeat:
    
    if trigger_wb > 3
        
        [ W_theta_L1, W_theta_L2, W_theta_R1, W_theta_R2, W_theta_LR1, W_theta_LR2 ] = weight(a_fit.theta_L1, a_fit.theta_L2, a_fit.theta_R1, a_fit.theta_R2, n_pol_theta, trigger_wb);
        
        [ W_eta_L1, W_eta_L2, W_eta_R1, W_eta_R2, W_eta_LR1, W_eta_LR2 ] = weight(a_fit.eta_L1, a_fit.eta_L2, a_fit.eta_R1, a_fit.eta_R2, n_pol_eta, trigger_wb);
        
        [ W_phi_L1, W_phi_L2, W_phi_R1, W_phi_R2, W_phi_LR1, W_phi_LR2 ] = weight(a_fit.phi_L1, a_fit.phi_L2, a_fit.phi_R1, a_fit.phi_R2, n_pol_phi, trigger_wb);
        
        
        wb_loc_mean = zeros(trigger_wb-2,2);
    
        for k = 1:(trigger_wb-2)

            if k == (trigger_wb-2)

                wb_end = find(isnan(pathDB.wingbeat_time(trigger_wb,:,seq_nr))==0, 1, 'last' );

                wb_loc_mean(k,:) = [ pathDB.wingbeat_time(trigger_wb,1,seq_nr) pathDB.wingbeat_time(trigger_wb,wb_end,seq_nr) ];

            else

                wb_loc_mean(k,:) = [ pathDB.wingbeat_time(k+2,1,seq_nr) pathDB.wingbeat_time(k+3,1,seq_nr) ];

            end

            
        end
        
        
        wb_loc_mean12 = zeros(trigger_wb-2,2);
    
        for k = 1:(trigger_wb-2)

            if k == (trigger_wb-2)

                wb_end = find(isnan(pathDB.wingbeat_time(trigger_wb,:,seq_nr))==0, 1, 'last' );
                
                frame_half = round((pathDB.wingbeat_time(trigger_wb,wb_end,seq_nr)+pathDB.wingbeat_time(trigger_wb,1,seq_nr))/2);

                wb_loc_mean1 = [ pathDB.wingbeat_time(trigger_wb,1,seq_nr) frame_half ];
                
                wb_loc_mean2 = [ frame_half pathDB.wingbeat_time(trigger_wb,wb_end,seq_nr) ];

            else
                
                wb_end = find(isnan(pathDB.wingbeat_time(k+2,:,seq_nr))==0, 1, 'last' );

                frame_half = round((pathDB.wingbeat_time(k+2,wb_end,seq_nr)+pathDB.wingbeat_time(k+2,1,seq_nr))/2);
                
                wb_loc_mean1 = [ pathDB.wingbeat_time(k+2,1,seq_nr) frame_half ];
                
                wb_loc_mean2 = [ frame_half pathDB.wingbeat_time(k+3,1,seq_nr) ];

            end
            
            wb_loc_mean12(2*k-1,:) = wb_loc_mean1;
        
            wb_loc_mean12(2*k,:) = wb_loc_mean2;           

        end


%         wb_loc_mean12
%         
% 
%         
%         [a_avg_theta_L1, a_avg_theta_L2] = average_fit(theta_L,n_pol_theta,wb_loc_mean12,W_theta_L1,W_theta_L2);
% 
%         PN_theta = Legendre_polynomial(n_pol_theta,2,-1:0.01:1);
%         PN_eta = Legendre_polynomial(n_pol_eta,2,-1:0.01:1);
%         PN_phi = Legendre_polynomial(n_pol_phi,2,-1:0.01:1);
%         
%         
%         figure()
%         hold on
%         for k = 1:nr_wb
%             m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%             m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%             t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%             t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%             plot(t1,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L1(:,k)),'b')
%             plot(t2,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L2(:,k)),'b')
%             plot(t1,radtodeg(PN_theta(:,:,1)'*a_avg_theta_L1),'g')
%             plot(t2,radtodeg(PN_theta(:,:,1)'*a_avg_theta_L2),'g')
%         end
%         plot(t,radtodeg(theta_L),'r')
%         hold off
%         
%         pause


        
%         [a_avg_theta_L1, a_avg_theta_L2] = average_fit(theta_L,n_pol_theta,wb_loc_mean12,W_theta_L1,W_theta_L2);
%         
%         [a_avg_theta_R1, a_avg_theta_R2] = average_fit(theta_R,n_pol_theta,wb_loc_mean12,W_theta_R1,W_theta_R2);
%         
%         [a_avg_eta_L1, a_avg_eta_L2] = average_fit(eta_L,n_pol_eta,wb_loc_mean12,W_eta_L1,W_eta_L2);
%         
%         [a_avg_eta_R1, a_avg_eta_R2] = average_fit(eta_R,n_pol_eta,wb_loc_mean12,W_eta_R1,W_eta_R2);
%         
%         [a_avg_phi_L1, a_avg_phi_L2] = average_fit(phi_L,n_pol_phi,wb_loc_mean12,W_phi_L1,W_phi_L2);
%         
%         [a_avg_phi_R1, a_avg_phi_R2] = average_fit(phi_R,n_pol_phi,wb_loc_mean12,W_phi_R1,W_phi_R2);

        [a_avg_theta_L1, a_avg_theta_L2] = average_fit(a_fit.theta_L1(:,3:trigger_wb),a_fit.theta_L2(:,3:trigger_wb),n_pol_theta,wb_loc_mean12,W_theta_L1,W_theta_L2);
        
        [a_avg_theta_R1, a_avg_theta_R2] = average_fit(a_fit.theta_R1(:,3:trigger_wb),a_fit.theta_R2(:,3:trigger_wb),n_pol_theta,wb_loc_mean12,W_theta_R1,W_theta_R2);
        
        [a_avg_eta_L1, a_avg_eta_L2] = average_fit(a_fit.eta_L1(:,3:trigger_wb),a_fit.eta_L2(:,3:trigger_wb),n_pol_eta,wb_loc_mean12,W_eta_L1,W_eta_L2);
        
        [a_avg_eta_R1, a_avg_eta_R2] = average_fit(a_fit.eta_R1(:,3:trigger_wb),a_fit.eta_R2(:,3:trigger_wb),n_pol_eta,wb_loc_mean12,W_eta_R1,W_eta_R2);
        
        [a_avg_phi_L1, a_avg_phi_L2] = average_fit(a_fit.phi_L1(:,3:trigger_wb),a_fit.phi_L2(:,3:trigger_wb),n_pol_phi,wb_loc_mean12,W_phi_L1,W_phi_L2);
        
        [a_avg_phi_R1, a_avg_phi_R2] = average_fit(a_fit.phi_R1(:,3:trigger_wb),a_fit.phi_R2(:,3:trigger_wb),n_pol_phi,wb_loc_mean12,W_phi_R1,W_phi_R2);        
        
        
%         Y_theta_LR = [theta_L(wb_loc_mean12(1,1):wb_loc_mean12(end,2)); theta_R(wb_loc_mean12(1,1):wb_loc_mean12(end,2))];
%        
%         Y_eta_LR = [eta_L(wb_loc_mean12(1,1):wb_loc_mean12(end,2)); eta_R(wb_loc_mean12(1,1):wb_loc_mean12(end,2))];
% 
%         Y_phi_LR = [phi_L(wb_loc_mean12(1,1):wb_loc_mean12(end,2)); phi_R(wb_loc_mean12(1,1):wb_loc_mean12(end,2))];


        wb_loc_LR = [wb_loc_mean12(:,1)-wb_loc_mean12(1,1)+1 wb_loc_mean12(:,2)-wb_loc_mean12(1,1)+1; ...
                     wb_loc_mean12(:,1)-2*wb_loc_mean12(1,1)+wb_loc_mean12(end,2)+2 wb_loc_mean12(:,2)-2*wb_loc_mean12(1,1)+wb_loc_mean12(end,2)+2];


        [a_avg_theta_LR1, a_avg_theta_LR2] = average_fit([a_fit.theta_L1(:,3:trigger_wb) a_fit.theta_R1(:,3:trigger_wb)],[a_fit.theta_L2(:,3:trigger_wb) a_fit.theta_R2(:,3:trigger_wb)],n_pol_theta,wb_loc_LR,W_theta_LR1,W_theta_LR2);

        [a_avg_eta_LR1, a_avg_eta_LR2] = average_fit([a_fit.eta_L1(:,3:trigger_wb) a_fit.eta_R1(:,3:trigger_wb)],[a_fit.eta_L2(:,3:trigger_wb) a_fit.eta_R2(:,3:trigger_wb)],n_pol_eta,wb_loc_LR,W_eta_LR1,W_eta_LR2);

        [a_avg_phi_LR1, a_avg_phi_LR2] = average_fit([a_fit.phi_L1(:,3:trigger_wb) a_fit.phi_R1(:,3:trigger_wb)],[a_fit.phi_L2(:,3:trigger_wb) a_fit.phi_R2(:,3:trigger_wb)],n_pol_phi,wb_loc_LR,W_phi_LR1,W_phi_LR2);   
        
        
        
        
        
        
        
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

        f_avg = 1/(((wb_loc_mean(end,1)-wb_loc_mean(1,1))*dt)/length(wb_loc_mean(:,1)));


        % Determine average ratio between down and upstroke:

        down_up_ratio = zeros(trigger_wb-2,1);

        for m = 1:(trigger_wb-2)

%         down_up_ratio = zeros(nr_wb,1);
% 
%         for m = 1:nr_wb

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
        
%         figure()
%         plot(down_up_ratio)
%         
%         pause

        down_up = mean(down_up_ratio);

        
        ratio_1_2 = zeros(nr_wb,1);

        for j = 1:nr_wb

            ratio_1_2(j) = (wb_loc12(j*2-1,2)-wb_loc12(j*2-1,1)+1)/(wb_loc12(j*2,2)-wb_loc12(j*2,1)+wb_loc12(j*2-1,2)-wb_loc12(j*2-1,1)+1);

        end

        ratio_1_2_avg = mean(ratio_1_2(3:trigger_wb));
        
        
                
    end
    
    
    
%     % Plot results:
%     
%     PN_theta = Legendre_polynomial(n_pol_theta,2,-1:0.01:1);
%     PN_eta = Legendre_polynomial(n_pol_eta,2,-1:0.01:1);
%     PN_phi = Legendre_polynomial(n_pol_phi,2,-1:0.01:1);
% 
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L1(:,k)),'b')
%         plot(t2,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L2(:,k)),'b')
%     end
%     plot(t,radtodeg(theta_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L1(:,k)),'b')
%         plot(t2,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L2(:,k)),'b')
%     end
%     plot(t,radtodeg(eta_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_phi(:,:,1)'*a_fit.phi_L1(:,k)),'b')
%         plot(t2,radtodeg(PN_phi(:,:,1)'*a_fit.phi_L2(:,k)),'b')
%     end
%     plot(t,radtodeg(phi_L),'r')
%     hold off
% 
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_theta(:,:,2)'*a_fit.theta_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_theta(:,:,2)'*a_fit.theta_L2(:,k)),'b')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_eta(:,:,2)'*a_fit.eta_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_eta(:,:,2)'*a_fit.eta_L2(:,k)),'b')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_phi(:,:,2)'*a_fit.phi_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_phi(:,:,2)'*a_fit.phi_L2(:,k)),'b')
%     end
%     hold off    
% 
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_theta(:,:,3)'*a_fit.theta_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_theta(:,:,3)'*a_fit.theta_L2(:,k)),'b')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_eta(:,:,3)'*a_fit.eta_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_eta(:,:,3)'*a_fit.eta_L2(:,k)),'b')
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_phi(:,:,3)'*a_fit.phi_L1(:,k)),'r')
%         plot(t2,radtodeg(PN_phi(:,:,3)'*a_fit.phi_L2(:,k)),'b')
%     end
%     hold off    
% 
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L1(:,k)),'b')
%         plot(t2,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L2(:,k)),'b')
%         plot(t1,radtodeg(PN_theta(:,:,1)'*a_avg.theta_L1),'g')
%         plot(t2,radtodeg(PN_theta(:,:,1)'*a_avg.theta_L2),'g')
%     end
%     plot(t,radtodeg(theta_L),'r')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 = wb_loc12(k*2-1,2)-wb_loc12(k*2-1,1);
%         m2 = wb_loc12(k*2,2)-wb_loc12(k*2,1);
%         t1 = t(wb_loc12(k*2-1,1)):((m1*dt)/(200)):t(wb_loc12(k*2-1,2));
%         t2 = t(wb_loc12(k*2,1)):((m2*dt)/(200)):t(wb_loc12(k*2,2));
%         plot(t1,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L1(:,k)),'b')
%         plot(t2,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L2(:,k)),'b')
%         plot(t1,radtodeg(PN_eta(:,:,1)'*a_avg.eta_L1),'g')
%         plot(t2,radtodeg(PN_eta(:,:,1)'*a_avg.eta_L2),'g')
%     end
%     plot(t,radtodeg(eta_L),'r')
%     hold off
%     
%     pause
    
end