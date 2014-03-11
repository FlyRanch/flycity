function [ FM_out ] = Robofly_tests( FM_in, a_glob )

   
    nr_steps    = FM_in.nr_steps;
    nr_of_tests = nr_steps*4;
    nr_wb_exp   = FM_in.nr_wb_exp;
    nr_points   = FM_in.nr_points;
    f           = FM_in.f;
    down_up     = FM_in.down_up;
    wing_l      = FM_in.wing_l;
    
    n_pol_theta     = (size(FM_in.a_theta_L,1)-2)/2;
    n_pol_eta       = (size(FM_in.a_eta_L,1)-2)/2;
    n_pol_phi       = (size(FM_in.a_phi_L,1)-2)/2;
    
    a_theta_glob    = a_glob.theta;
    a_eta_glob      = a_glob.eta;
    a_phi_glob      = a_glob.phi;
    
    a_theta_L       = nan((n_pol_theta+1)*2,nr_wb_exp,nr_of_tests);
    a_eta_L         = nan((n_pol_eta+1)*2,nr_wb_exp,nr_of_tests);
    a_phi_L         = nan((n_pol_phi+1)*2,nr_wb_exp,nr_of_tests);
    a_theta_R       = nan((n_pol_theta+1)*2,nr_wb_exp,nr_of_tests);
    a_eta_R         = nan((n_pol_eta+1)*2,nr_wb_exp,nr_of_tests);
    a_phi_R         = nan((n_pol_phi+1)*2,nr_wb_exp,nr_of_tests);
    
    time_ref        = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    
    theta_L_ref     = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    eta_L_ref       = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    phi_L_ref       = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    theta_R_ref     = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    eta_R_ref       = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);
    phi_R_ref       = nan((nr_points-1)*nr_wb_exp+1,nr_of_tests);

    for i = 1:4
        
        % theta, eta and phi on:
        
        k_range = ((i-1)*nr_steps+1):(nr_steps*i);
        
        if i == 1
            
            for j = 1:nr_steps
                
                a_theta_L(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*FM_in.a_theta_L(:,j)')';
                a_eta_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_eta_L(:,j)')';
                a_phi_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_phi_L(:,j)')';
                a_theta_R(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*FM_in.a_theta_R(:,j)')';
                a_eta_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_eta_R(:,j)')';
                a_phi_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_phi_R(:,j)')';
                
            end
            
            time_ref(:,k_range)     = (ones(nr_steps,1)*FM_in.t)';
            
            theta_L_ref(:,k_range)  = FM_in.theta_L';
            eta_L_ref(:,k_range)    = FM_in.eta_L';
            phi_L_ref(:,k_range)    = FM_in.phi_L';
            theta_R_ref(:,k_range)  = FM_in.theta_R';
            eta_R_ref(:,k_range)    = FM_in.eta_R';
            phi_R_ref(:,k_range)    = FM_in.phi_R';
            
        end
        
        % theta only:
        
        if i == 2
            
            for j = 1:nr_steps
                
                a_theta_L(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*FM_in.a_theta_L(:,j)')';
                a_eta_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_eta_glob')';
                a_phi_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_phi_glob')';
                a_theta_R(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*FM_in.a_theta_R(:,j)')';
                a_eta_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_eta_glob')';
                a_phi_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_phi_glob')';
                
            end
            
            time_ref(:,k_range)     = (ones(nr_steps,1)*FM_in.t)';
            
            theta_L_ref(:,k_range)  = FM_in.theta_L';
            eta_L_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.eta_L(1,:))';
            phi_L_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.phi_L(1,:))';
            theta_R_ref(:,k_range)  = FM_in.theta_R';
            eta_R_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.eta_R(1,:))';
            phi_R_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.phi_R(1,:))';
            
        end
        
        % eta only:
        
        if i == 3
            
            for j = 1:nr_steps
                
                a_theta_L(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*a_theta_glob')';
                a_eta_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_eta_L(:,j)')';
                a_phi_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_phi_glob')';
                a_theta_R(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*a_theta_glob')';
                a_eta_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_eta_R(:,j)')';
                a_phi_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_phi_glob')';
                
            end
            
            time_ref(:,k_range)     = (ones(nr_steps,1)*FM_in.t)';
            
            theta_L_ref(:,k_range)  = (ones(nr_steps,1)*FM_in.theta_L(1,:))';
            eta_L_ref(:,k_range)    = FM_in.eta_L';
            phi_L_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.phi_L(1,:))';
            theta_R_ref(:,k_range)  = (ones(nr_steps,1)*FM_in.theta_R(1,:))';
            eta_R_ref(:,k_range)    = FM_in.eta_R';
            phi_R_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.phi_R(1,:))';
            
        end
        
        % phi only:
        
        if i == 4
            
            for j = 1:nr_steps
                
                a_theta_L(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*a_theta_glob')';
                a_eta_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_eta_glob')';
                a_phi_L(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_phi_L(:,j)')';
                a_theta_R(:,:,(i-1)*nr_steps+j)     = (ones(nr_wb_exp,1)*a_theta_glob')';
                a_eta_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*a_eta_glob')';
                a_phi_R(:,:,(i-1)*nr_steps+j)       = (ones(nr_wb_exp,1)*FM_in.a_phi_R(:,j)')';
                
            end
            
            time_ref(:,k_range)     = (ones(nr_steps,1)*FM_in.t)';
            
            theta_L_ref(:,k_range)  = (ones(nr_steps,1)*FM_in.theta_L(1,:))';
            eta_L_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.eta_L(1,:))';
            phi_L_ref(:,k_range)    = FM_in.phi_L';
            theta_R_ref(:,k_range)  = (ones(nr_steps,1)*FM_in.theta_R(1,:))';
            eta_R_ref(:,k_range)    = (ones(nr_steps,1)*FM_in.eta_R(1,:))';
            phi_R_ref(:,k_range)    = FM_in.phi_R';
            
        end
        
    end
    
    FM_out.nr_of_wb     = nr_wb_exp;
    FM_out.nr_of_tests  = nr_of_tests;
    FM_out.a_theta_L    = a_theta_L;
    FM_out.a_eta_L      = a_eta_L;
    FM_out.a_phi_L      = a_phi_L;
    FM_out.a_theta_R    = a_theta_R;
    FM_out.a_eta_R      = a_eta_R;
    FM_out.a_phi_R      = a_phi_R;
    FM_out.freq         = f;
    FM_out.down_up      = down_up;
    FM_out.wing_length  = wing_l;
    FM_out.time_ref     = time_ref;
    FM_out.theta_L_ref  = theta_L_ref;
    FM_out.eta_L_ref    = eta_L_ref;
    FM_out.phi_L_ref    = phi_L_ref;
    FM_out.theta_R_ref  = theta_R_ref;
    FM_out.eta_R_ref    = eta_R_ref;
    FM_out.phi_R_ref    = phi_R_ref;
    FM_out.n_pol_theta  = n_pol_theta;
    FM_out.n_pol_eta    = n_pol_eta;
    FM_out.n_pol_phi    = n_pol_phi;
    
%     % Plot results:
%     
%     for i = 1:4
%         
%         % theta, eta and phi on:
%         
%         if i == 1
%             
%                 figure()
%                 hold on
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,1); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,3); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,5); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\phi_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,2); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,4); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,6); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 hold off
%                 hold off
%                 xlabel('t [s]')
%                 ylabel('\phi_R [deg]')
%             
%         end
%         
%         % theta only:
%         
%         if i == 2
%             
%                 figure()
%                 hold on
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,1); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,3); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,5); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\phi_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,2); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,4); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,6); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 hold off
%                 hold off
%                 xlabel('t [s]')
%                 ylabel('\phi_R [deg]')
%             
%         end
%         
%         % eta only:
%         
%         if i == 3
%             
%                 figure()
%                 hold on
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,1); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,3); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,5); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\phi_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,2); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,4); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,6); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 hold off
%                 hold off
%                 xlabel('t [s]')
%                 ylabel('\phi_R [deg]')
%             
%         end
%         
%         % phi only:
%         
%         if i == 4
%             
%                 figure()
%                 hold on
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,1); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,3); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,5); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_L_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\phi_L [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,2); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(theta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\theta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,4); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(eta_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 xlabel('t [s]')
%                 ylabel('\eta_R [deg]')
%                 for j = 1:nr_steps
%                     hold on
%                     subplot(3,2,6); plot(time_ref(:,(i-1)*nr_steps+j),radtodeg(phi_R_ref(:,(i-1)*nr_steps+j)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%                     hold off
%                 end
%                 hold off
%                 hold off
%                 xlabel('t [s]')
%                 ylabel('\phi_R [deg]')
%             
%         end
%         
%     end    

end

