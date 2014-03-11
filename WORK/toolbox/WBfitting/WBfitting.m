% function [a_fit, a_avg] = WBfitting()

    n_pol_theta = 12; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 8; % Order of used polynomials 
    
    theta_L = nan(1,50);
    
    theta_R = nan(1,50);
    
    eta_L = nan(1,50);
    
    eta_R = nan(1,50);
    
    phi_L = nan(1,50);
    
    phi_R = nan(1,50);
    
    a_fit = {};
    

    for i = 1:length(pathDB.x(1,:))

    seq_nr = i
        
    % Create a standard wingbeat:
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    t = pathDB.t(start:stop);
    
    dt = pathDB.t(2)-pathDB.t(1);

    theta_L_loc = pathDB.theta_L(start:stop,seq_nr);
    
    theta_R_loc = pathDB.theta_R(start:stop,seq_nr);

    eta_L_loc = pathDB.eta_L(start:stop,seq_nr);
    
    eta_R_loc = pathDB.eta_R(start:stop,seq_nr);

    phi_L_loc = pathDB.phi_L(start:stop,seq_nr);
    
    phi_R_loc = pathDB.phi_R(start:stop,seq_nr);
    
    nr_wb_loc = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    
    
    % Wingbeat locations

    wb_loc = zeros(nr_wb_loc,2);
    
    wb_loc12 = zeros(2*nr_wb_loc,2);

    for k = 1:nr_wb_loc

        if k == nr_wb_loc

            wb_end = find(isnan(pathDB.wingbeat_time(nr_wb_loc,:,seq_nr))==0, 1, 'last' );

            wb_loc(k,:) = [ pathDB.wingbeat_time(nr_wb_loc,1,seq_nr) pathDB.wingbeat_time(nr_wb_loc,wb_end,seq_nr) ];

        else

            wb_loc(k,:) = [ pathDB.wingbeat_time(k,1,seq_nr) pathDB.wingbeat_time(k+1,1,seq_nr) ];

        end

    end  
    
    a_fit_theta_L1 = zeros(nr_wb_loc,n_pol_theta+1);
        
    a_fit_eta_L1 = zeros(nr_wb_loc,n_pol_eta+1);
        
    a_fit_phi_L1 =  zeros(nr_wb_loc,n_pol_phi+1);
        
    a_fit_theta_L2 = zeros(nr_wb_loc,n_pol_theta+1);
        
    a_fit_eta_L2 = zeros(nr_wb_loc,n_pol_eta+1);
        
    a_fit_phi_L2 = zeros(nr_wb_loc,n_pol_phi+1);

    a_fit_theta_R1 = zeros(nr_wb_loc,n_pol_theta+1);
        
    a_fit_eta_R1 = zeros(nr_wb_loc,n_pol_eta+1);
        
    a_fit_phi_R1 = zeros(nr_wb_loc,n_pol_phi+1);
        
    a_fit_theta_R2 = zeros(nr_wb_loc,n_pol_theta+1);
        
    a_fit_eta_R2 = zeros(nr_wb_loc,n_pol_eta+1);
        
    a_fit_phi_R2 = zeros(nr_wb_loc,n_pol_phi+1);   
  
    a_fit_ratio_1_2 = zeros(nr_wb_loc,1);
    
    a_fit_wb1_wb2 = zeros(nr_wb_loc,2);
    
%     theta_L_t = nan(nr_wb_loc,50);
%     
%     theta_R_t = nan(nr_wb_loc,50);
%     
%     eta_L_t = nan(nr_wb_loc,50);
%     
%     eta_R_t = nan(nr_wb_loc,50);
%     
%     phi_L_t = nan(nr_wb_loc,50);
%     
%     phi_R_t = nan(nr_wb_loc,50);
    
    for j = 1:nr_wb_loc
        
        [ a_theta_L1, a_theta_L2, ~ ] = local_fit(theta_L_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_theta);
        
        [ a_eta_L1, a_eta_L2, ~ ] = local_fit(eta_L_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_eta);
        
        [ a_phi_L1, a_phi_L2, ~ ] = local_fit(phi_L_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_phi);

        [ a_theta_R1, a_theta_R2, ~ ] = local_fit(theta_R_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_theta);
        
        [ a_eta_R1, a_eta_R2, ~ ] = local_fit(eta_R_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_eta);
        
        [ a_phi_R1, a_phi_R2, ratio_1_2 ] = local_fit(phi_R_loc(wb_loc(j,1):wb_loc(j,2)),n_pol_phi);
        
        a_fit_theta_L1(j,:) = a_theta_L1;
        
        a_fit_eta_L1(j,:) = a_eta_L1;
        
        a_fit_phi_L1(j,:) = a_phi_L1;
        
        a_fit_theta_L2(j,:) = a_theta_L2;
        
        a_fit_eta_L2(j,:) = a_eta_L2;
        
        a_fit_phi_L2(j,:) = a_phi_L2;

        a_fit_theta_R1(j,:) = a_theta_R1;
        
        a_fit_eta_R1(j,:) = a_eta_R1;
        
        a_fit_phi_R1(j,:) = a_phi_R1;
        
        a_fit_theta_R2(j,:) = a_theta_R2;
        
        a_fit_eta_R2(j,:) = a_eta_R2;
        
        a_fit_phi_R2(j,:) = a_phi_R2;
        
        a_fit_ratio_1_2(j) = ratio_1_2;
        
        a_fit_wb1_wb2(j,1) = ceil((ratio_1_2*(wb_loc(j,2)-wb_loc(j,1)+1))/2);
        
        a_fit_wb1_wb2(j,2) = floor(((2-ratio_1_2)*(wb_loc(j,2)-wb_loc(j,1)+1))/2);
%         
%         wb_loc_end = wb_loc(j,2)-wb_loc(j,1)+1;
%         
%         theta_L_t(j,1:wb_loc_end) = theta_L_loc(wb_loc(j,1):wb_loc(j,2));
%         
%         theta_R_t(j,1:wb_loc_end) = theta_R_loc(wb_loc(j,1):wb_loc(j,2));
%         
%         eta_L_t(j,1:wb_loc_end) = eta_L_loc(wb_loc(j,1):wb_loc(j,2));
%         
%         eta_R_t(j,1:wb_loc_end) = eta_R_loc(wb_loc(j,1):wb_loc(j,2));
%         
%         phi_L_t(j,1:wb_loc_end) = phi_L_loc(wb_loc(j,1):wb_loc(j,2));
%         
%         phi_R_t(j,1:wb_loc_end) = phi_R_loc(wb_loc(j,1):wb_loc(j,2));
        
    end    
    
        if i == 1

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
            
            a_fit.ratio_1_2 = a_fit_ratio_1_2;
            
            a_fit.wb1_wb2 = a_fit_wb1_wb2;
            
%             theta_L = theta_L_t;
% 
%             theta_R = theta_R_t;
% 
%             eta_L = eta_L_t;
% 
%             eta_R = eta_R_t;
% 
%             phi_L = phi_L_t;
% 
%             phi_R = phi_R_t;

        else

            a_fit.theta_L1 = [a_fit.theta_L1; a_fit_theta_L1];
            a_fit.theta_L2 = [a_fit.theta_L2; a_fit_theta_L2];
            a_fit.eta_L1 = [a_fit.eta_L1; a_fit_eta_L1];
            a_fit.eta_L2 = [a_fit.eta_L2; a_fit_eta_L2];
            a_fit.phi_L1 = [a_fit.phi_L1; a_fit_phi_L1];
            a_fit.phi_L2 = [a_fit.phi_L2; a_fit_phi_L2];

            a_fit.theta_R1 = [a_fit.theta_R1; a_fit_theta_R1];
            a_fit.theta_R2 = [a_fit.theta_R2; a_fit_theta_R2];
            a_fit.eta_R1 = [a_fit.eta_R1; a_fit_eta_R1];
            a_fit.eta_R2 = [a_fit.eta_R2; a_fit_eta_R2];
            a_fit.phi_R1 = [a_fit.phi_R1; a_fit_phi_R1];
            a_fit.phi_R2 = [a_fit.phi_R2; a_fit_phi_R2];
            
            a_fit.ratio_1_2 = [a_fit.ratio_1_2; a_fit_ratio_1_2];
            
            a_fit.wb1_wb2 = [a_fit.wb1_wb2; a_fit_wb1_wb2];
%           
%             theta_L = [theta_L; theta_L_t];
% 
%             theta_R = [theta_R; theta_R_t];
% 
%             eta_L = [eta_L; eta_L_t];
% 
%             eta_R = [eta_R; eta_R_t];
% 
%             phi_L = [phi_L; phi_L_t];
% 
%             phi_R = [phi_R; phi_R_t];
            
        end
    
    end
    
    
    % Adjust wb_loc_mean and theta_L etc. then it should work.
        
    wb_loc_mean12 = zeros(2*length(a_fit.wb1_wb2(:,1)),2);
    
    for k = 1:(length(a_fit.wb1_wb2(:,1)))
               
        wb_loc_mean12(2*k-1,:) = [1 a_fit.wb1_wb2(k,1)];
        
        wb_loc_mean12(2*k,:) = [a_fit.wb1_wb2(k,1)+1 (a_fit.wb1_wb2(k,1)+a_fit.wb1_wb2(k,2))];

    end
    


        
        [a_avg_theta_L1, a_avg_theta_L2] = average_fit(a_fit.theta_L1,a_fit.theta_L2,n_pol_theta);
        
        [a_avg_theta_R1, a_avg_theta_R2] = average_fit(a_fit.theta_R1,a_fit.theta_R2,n_pol_theta);
        
        [a_avg_eta_L1, a_avg_eta_L2] = average_fit(a_fit.eta_L1,a_fit.eta_L2,n_pol_eta);
        
        [a_avg_eta_R1, a_avg_eta_R2] = average_fit(a_fit.eta_R1,a_fit.eta_R2,n_pol_eta);
        
        [a_avg_phi_L1, a_avg_phi_L2] = average_fit(a_fit.phi_L1,a_fit.phi_L2,n_pol_phi);
        
        [a_avg_phi_R1, a_avg_phi_R2] = average_fit(a_fit.phi_R1,a_fit.phi_R2,n_pol_phi);        
        

        wb_loc_LR = [wb_loc_mean12; wb_loc_mean12];


        [a_avg_theta_LR1, a_avg_theta_LR2] = average_fit([a_fit.theta_L1; a_fit.theta_R1],[a_fit.theta_L2; a_fit.theta_R2],n_pol_theta);

        [a_avg_eta_LR1, a_avg_eta_LR2] = average_fit([a_fit.eta_L1; a_fit.eta_R1],[a_fit.eta_L2; a_fit.eta_R2],n_pol_eta);

        [a_avg_phi_LR1, a_avg_phi_LR2] = average_fit([a_fit.phi_L1; a_fit.phi_R1],[a_fit.phi_L2; a_fit.phi_R2],n_pol_phi);   
        
        
        
        % Store average coefficients in a_avg
        
        a_avg = {};

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
        
        a_avg.wb1_wb2 = [mean(a_fit.wb1_wb2(:,1)) mean(a_fit.wb1_wb2(:,2))];


        
%         % Determine average wingbeat frequency:
% 
%         f_avg = 1/(((wb_loc_mean(end,1)-wb_loc_mean(1,1))*dt)/length(wb_loc_mean(:,1)));
% 
% 
%         % Determine average ratio between down and upstroke:
% 
%         down_up_ratio = zeros(trigger_wb-2,1);
% 
%         for m = 1:(trigger_wb-2)
% 
%             if pathDB.L_wingbeat_loc(m+2,2,seq_nr) > pathDB.L_wingbeat_loc(m+2,1,seq_nr);
% 
%                 down_length_L = pathDB.L_wingbeat_loc(m+3,1,seq_nr)-pathDB.L_wingbeat_loc(m+2,2,seq_nr);
% 
%                 up_length_L = pathDB.L_wingbeat_loc(m+3,2,seq_nr)-pathDB.L_wingbeat_loc(m+3,1,seq_nr);
% 
%                 down_length_R = pathDB.R_wingbeat_loc(m+3,1,seq_nr)-pathDB.R_wingbeat_loc(m+2,2,seq_nr);
% 
%                 up_length_R = pathDB.R_wingbeat_loc(m+3,2,seq_nr)-pathDB.R_wingbeat_loc(m+3,1,seq_nr);
% 
%                 down_up_ratio(m) = 0.5*((down_length_L/up_length_L)+(down_length_R/up_length_R));
% 
%             elseif pathDB.L_wingbeat_loc(m+2,2,seq_nr) < pathDB.L_wingbeat_loc(m+2,1,seq_nr)
% 
%                 down_length_L = pathDB.L_wingbeat_loc(m+2,1,seq_nr)-pathDB.L_wingbeat_loc(m+2,2,seq_nr);
% 
%                 up_length_L = pathDB.L_wingbeat_loc(m+3,2,seq_nr)-pathDB.L_wingbeat_loc(m+2,1,seq_nr);
% 
%                 down_length_R = pathDB.R_wingbeat_loc(m+2,1,seq_nr)-pathDB.R_wingbeat_loc(m+2,2,seq_nr);
% 
%                 up_length_R = pathDB.R_wingbeat_loc(m+3,2,seq_nr)-pathDB.R_wingbeat_loc(m+2,1,seq_nr);
% 
%                 down_up_ratio(m) = 0.5*((down_length_L/up_length_L)+(down_length_R/up_length_R));
% 
%             end
% 
%         end
% 
%         down_up = mean(down_up_ratio);
% 
%         
%         ratio_1_2 = zeros(nr_wb,1);
% 
%         for j = 1:nr_wb
% 
%             ratio_1_2(j) = (wb_loc12(j*2-1,2)-wb_loc12(j*2-1,1)+1)/(wb_loc12(j*2,2)-wb_loc12(j*2,1)+1);
% 
%         end
% 
%         ratio_1_2_avg = mean(ratio_1_2(3:trigger_wb));
        
        
                

    
    
    
    % Plot results:
    
    PN_theta = Legendre_polynomial(n_pol_theta,2,-1:0.01:1);
    PN_eta = Legendre_polynomial(n_pol_eta,2,-1:0.01:1);
    PN_phi = Legendre_polynomial(n_pol_phi,2,-1:0.01:1);

    nr_wb = length(a_fit.wb1_wb2(:,1))
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(t1,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(t2,radtodeg(PN_theta(:,:,1)'*a_fit.theta_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,radtodeg(PN_theta(:,:,1)'*a_avg.theta_L1),'r')
    plot(t2,radtodeg(PN_theta(:,:,1)'*a_avg.theta_L2),'r')    
    hold off
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(t1,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(t2,radtodeg(PN_eta(:,:,1)'*a_fit.eta_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,radtodeg(PN_eta(:,:,1)'*a_avg.eta_L1),'r')
    plot(t2,radtodeg(PN_eta(:,:,1)'*a_avg.eta_L2),'r')   
    hold off
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(-1:0.005:0,radtodeg(PN_phi(:,:,1)'*a_fit.phi_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(0:0.005:1,radtodeg(PN_phi(:,:,1)'*a_fit.phi_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,radtodeg(PN_phi(:,:,1)'*a_avg.phi_L1),'r')
    plot(t2,radtodeg(PN_phi(:,:,1)'*a_avg.phi_L2),'r')   
    hold off
    


end

