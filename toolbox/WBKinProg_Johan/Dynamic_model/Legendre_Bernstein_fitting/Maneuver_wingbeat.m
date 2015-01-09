function Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,trigger_wb)


    % Determine the deviation between the average left wingbeat and the
    % actual left wingbeats and the deviation between the average right
    % wingbeat and the actual right wingbeats:
    
    nr_wb = length(a_fit_L.theta(1,:));
    
    
    n_pol_theta = n_pol.theta;
    
    n_pol_eta = n_pol.eta;
    
    n_pol_phi = n_pol.phi;
    
    
    dev_theta_L = zeros(201,nr_wb);
    
    dev_eta_L = zeros(201,nr_wb);
    
    dev_phi_L = zeros(201,nr_wb);
    
    dev_theta_R = zeros(201,nr_wb);
    
    dev_eta_R = zeros(201,nr_wb);
    
    dev_phi_R = zeros(201,nr_wb);
    
    
    dev_theta_LR = zeros(201,nr_wb);
    
    dev_eta_LR = zeros(201,nr_wb);
    
    dev_phi_LR = zeros(201,nr_wb);
    
    
    
    theta_L_avg = zeros(201,nr_wb);
    
    theta_R_avg = zeros(201,nr_wb);
    
    theta_avg = zeros(201,nr_wb);

    
    t = -1:0.01:1;
    
    for j = 1:nr_wb
        
        PN_theta = Legendre_polynomial(n_pol_theta,2,t);
        PN_eta = Legendre_polynomial(n_pol_eta,2,t);
        PN_phi = Legendre_polynomial(n_pol_phi,2,t);

        
        dev_eta_L(:,j) = PN_eta(:,:,1)'*a_fit_L.eta(:,j)-PN_eta(:,:,1)'*a_avg_L.eta;
        dev_theta_L(:,j) = PN_theta(:,:,1)'*a_fit_L.theta(:,j)-PN_theta(:,:,1)'*a_avg_L.theta;
        dev_phi_L(:,j) = PN_phi(:,:,1)'*a_fit_L.phi(:,j)-PN_phi(:,:,1)'*a_avg_L.phi;
        
        dev_eta_R(:,j) = PN_eta(:,:,1)'*a_fit_R.eta(:,j)-PN_eta(:,:,1)'*a_avg_R.eta;
        dev_theta_R(:,j) = PN_theta(:,:,1)'*a_fit_R.theta(:,j)-PN_theta(:,:,1)'*a_avg_R.theta;
        dev_phi_R(:,j) = PN_phi(:,:,1)'*a_fit_R.phi(:,j)-PN_phi(:,:,1)'*a_avg_R.phi;


        dev_eta_LR(:,j) = dev_eta_L(:,j)-dev_eta_R(:,j);
        dev_theta_LR(:,j) = dev_theta_L(:,j)-dev_theta_R(:,j);
        dev_phi_LR(:,j) = dev_phi_L(:,j)-dev_phi_R(:,j);
        
        theta_L_avg(:,j) = PN_theta(:,:,1)'*a_avg_LR.theta + dev_theta_L(:,j);
        
        theta_R_avg(:,j) = PN_theta(:,:,1)'*a_avg_LR.theta + dev_theta_R(:,j);

        theta_avg(:,j) = PN_theta(:,:,1)'*a_avg_LR.theta;

    end
    
    figure()
    hold on
    surf(theta_L_avg)
    surf(theta_avg)
    hold off
    
    figure()
    hold on
    surf(theta_R_avg)
    surf(theta_avg)
    hold off
    
    
    
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta left')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta left')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi right')
%     hold off
%     
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta left-right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta left-right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi left-right')
%     hold off
    
    
    

    
    
       
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_dot_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_dot_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta dot left')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_dot_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_dot_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta dot left')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_dot_L(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_dot_L(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi dot right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_dot_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_dot_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta dot right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_dot_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_dot_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta dot right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_dot_R(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_dot_R(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi dot right')
%     hold off
%     
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_eta_dot_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_eta_dot_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation eta dot left-right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_theta_dot_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_theta_dot_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation theta dot left-right')
%     hold off
%     
%     figure()
%     hold on
%     for k = 3:trigger_wb
%         plot(radtodeg(dev_phi_dot_LR(:,k)),'r')
%     end
%     for i = trigger_wb:(nr_wb-2)
%         plot(radtodeg(dev_phi_dot_LR(:,i)),'Color',[0.5 0.5 0.5])
%     end
%     title('deviation phi dot left-right')
%     hold off

    
    
    figure()
    surf(radtodeg(dev_eta_L(:,trigger_wb:(nr_wb-2))))

    figure()
    surf(radtodeg(dev_eta_R(:,trigger_wb:(nr_wb-2))))
    
    figure()
    surf(radtodeg(dev_eta_LR(:,trigger_wb:(nr_wb-2))))

    figure()
    surf(radtodeg(dev_theta_L(:,trigger_wb:(nr_wb-2))))

    figure()
    surf(radtodeg(dev_theta_R(:,trigger_wb:(nr_wb-2))))
    
    figure()
    surf(radtodeg(dev_theta_LR(:,trigger_wb:(nr_wb-2))))

    figure()
    surf(radtodeg(dev_phi_L(:,trigger_wb:(nr_wb-2))))

    figure()
    surf(radtodeg(dev_phi_R(:,trigger_wb:(nr_wb-2))))
    
    figure()
    surf(radtodeg(dev_phi_LR(:,trigger_wb:(nr_wb-2))))
    
    
end

