function Histogram(settings, pathDB)


    %Function that uses the mean and standard deviation of several
    %variables to create histograms. Besides the histograms the program
    %will generate an excel file which contains the input data of the
    %histogram and an identification of the session of wingbeats involved.

    % Concatenate the wingbeat variables of the different movie sequences
    % in large single vectors:
    

    
    for i = 1:size(pathDB.phi_body_mean,3)

        stop = find(isnan(pathDB.phi_body_mean(:,1,i))==0, 1 ,'last');
        
        if i == 1
            
            % Body
            
            phi_body = pathDB.phi_body_mean(1:stop,:,1);
            theta_body = pathDB.theta_body_mean(1:stop,:,1);
            xsi_body = pathDB.xsi_body_mean(1:stop,:,1);
            
            omegax_body = pathDB.omegax_body_mean(1:stop,:,1);
            omegay_body = pathDB.omegay_body_mean(1:stop,:,1);
            omegaz_body = pathDB.omegaz_body_mean(1:stop,:,1);
            Omega_body = pathDB.Omega_body_mean(1:stop,:,1);
            
            alfa_body = pathDB.alfa_body_mean(1:stop,:,1);
            beta_body = pathDB.beta_body_mean(1:stop,:,1);
            
            u_body = pathDB.u_body_mean(1:stop,:,1);
            v_body = pathDB.v_body_mean(1:stop,:,1);
            w_body = pathDB.w_body_mean(1:stop,:,1);
            U_body = pathDB.U_body_mean(1:stop,:,1);
            
            ax_body = pathDB.ax_body_mean(1:stop,:,1);
            ay_body = pathDB.ay_body_mean(1:stop,:,1);
            az_body = pathDB.az_body_mean(1:stop,:,1);
            a_body = pathDB.a_body_mean(1:stop,:,1);
            
            
            % Left wing
            
            phi_L_up = pathDB.phi_L_up_mean(1:stop,:,1);
            theta_L_up = pathDB.theta_L_up_mean(1:stop,:,1);
            eta_L_up = pathDB.eta_L_up_mean(1:stop,:,1);
            
            alfa_L_up = pathDB.alfa_L_up_mean(1:stop,:,1);
            beta_L_up = pathDB.beta_L_up_mean(1:stop,:,1);
            
            u_L_up = pathDB.u_L_up_mean(1:stop,:,1);
            v_L_up = pathDB.v_L_up_mean(1:stop,:,1);
            w_L_up = pathDB.w_L_up_mean(1:stop,:,1);
            U_L_up = pathDB.U_L_up_mean(1:stop,:,1);
            
            omegax_L_up = pathDB.omegax_L_up_mean(1:stop,:,1);
            omegay_L_up = pathDB.omegay_L_up_mean(1:stop,:,1);
            omegaz_L_up = pathDB.omegaz_L_up_mean(1:stop,:,1);
            Omega_L_up = pathDB.Omega_L_up_mean(1:stop,:,1);
            
            phi_L_down = pathDB.phi_L_down_mean(1:stop,:,1);
            theta_L_down = pathDB.theta_L_down_mean(1:stop,:,1);
            eta_L_down = pathDB.eta_L_down_mean(1:stop,:,1);
            
            alfa_L_down = pathDB.alfa_L_down_mean(1:stop,:,1);
            beta_L_down = pathDB.beta_L_down_mean(1:stop,:,1);
            
            u_L_down = pathDB.u_L_down_mean(1:stop,:,1);
            v_L_down = pathDB.v_L_down_mean(1:stop,:,1);
            w_L_down = pathDB.w_L_down_mean(1:stop,:,1);
            U_L_down = pathDB.U_L_down_mean(1:stop,:,1);
            
            omegax_L_down = pathDB.omegax_L_down_mean(1:stop,:,1);
            omegay_L_down = pathDB.omegay_L_down_mean(1:stop,:,1);
            omegaz_L_down = pathDB.omegaz_L_down_mean(1:stop,:,1);
            Omega_L_down = pathDB.Omega_L_down_mean(1:stop,:,1);
            
            
            % Right wing
            
            phi_R_up = pathDB.phi_R_up_mean(1:stop,:,1);
            theta_R_up = pathDB.theta_R_up_mean(1:stop,:,1);
            eta_R_up = pathDB.eta_R_up_mean(1:stop,:,1);
            
            alfa_R_up = pathDB.alfa_R_up_mean(1:stop,:,1);
            beta_R_up = pathDB.beta_R_up_mean(1:stop,:,1);
            
            u_R_up = pathDB.u_R_up_mean(1:stop,:,1);
            v_R_up = pathDB.v_R_up_mean(1:stop,:,1);
            w_R_up = pathDB.w_R_up_mean(1:stop,:,1);
            U_R_up = pathDB.U_R_up_mean(1:stop,:,1);
            
            omegax_R_up = pathDB.omegax_R_up_mean(1:stop,:,1);
            omegay_R_up = pathDB.omegay_R_up_mean(1:stop,:,1);
            omegaz_R_up = pathDB.omegaz_R_up_mean(1:stop,:,1);
            Omega_R_up = pathDB.Omega_R_up_mean(1:stop,:,1);
            
            phi_R_down = pathDB.phi_R_down_mean(1:stop,:,1);
            theta_R_down = pathDB.theta_R_down_mean(1:stop,:,1);
            eta_R_down = pathDB.eta_R_down_mean(1:stop,:,1);
            
            alfa_R_down = pathDB.alfa_R_down_mean(1:stop,:,1);
            beta_R_down = pathDB.beta_R_down_mean(1:stop,:,1);
            
            u_R_down = pathDB.u_R_down_mean(1:stop,:,1);
            v_R_down = pathDB.v_R_down_mean(1:stop,:,1);
            w_R_down = pathDB.w_R_down_mean(1:stop,:,1);
            U_R_down = pathDB.U_R_down_mean(1:stop,:,1);
            
            omegax_R_down = pathDB.omegax_R_down_mean(1:stop,:,1);
            omegay_R_down = pathDB.omegay_R_down_mean(1:stop,:,1);
            omegaz_R_down = pathDB.omegaz_R_down_mean(1:stop,:,1);
            Omega_R_down = pathDB.Omega_R_down_mean(1:stop,:,1);
            
            
            
        else
        
            % Concatenate all data into vectors
        
            phi_body = [phi_body; pathDB.phi_body_mean(1:stop,:,i)];
            theta_body = [theta_body; pathDB.theta_body_mean(1:stop,:,i)];
            xsi_body = [xsi_body; pathDB.xsi_body_mean(1:stop,:,i)];
            
            omegax_body = [omegax_body; pathDB.omegax_body_mean(1:stop,:,i)];
            omegay_body = [omegay_body; pathDB.omegay_body_mean(1:stop,:,i)];
            omegaz_body = [omegaz_body; pathDB.omegaz_body_mean(1:stop,:,i)];
            Omega_body = [Omega_body; pathDB.Omega_body_mean(1:stop,:,i)];
            
            alfa_body = [alfa_body; pathDB.alfa_body_mean(1:stop,:,i)];
            beta_body = [beta_body; pathDB.beta_body_mean(1:stop,:,i)];
            
            u_body = [u_body; pathDB.u_body_mean(1:stop,:,i)];
            v_body = [v_body; pathDB.v_body_mean(1:stop,:,i)];
            w_body = [w_body; pathDB.w_body_mean(1:stop,:,i)];
            U_body = [U_body; pathDB.U_body_mean(1:stop,:,i)];
            
            ax_body = [ax_body; pathDB.ax_body_mean(1:stop,:,i)]; 
            ay_body = [ay_body; pathDB.ay_body_mean(1:stop,:,i)];
            az_body = [az_body; pathDB.az_body_mean(1:stop,:,i)];
            a_body = [a_body; pathDB.a_body_mean(1:stop,:,i)];
            
            % Left wing
            
            phi_L_up = [phi_L_up; pathDB.phi_L_up_mean(1:stop,:,i)];
            theta_L_up = [theta_L_up; pathDB.theta_L_up_mean(1:stop,:,i)];
            eta_L_up = [eta_L_up; pathDB.eta_L_up_mean(1:stop,:,i)];
            
            alfa_L_up = [alfa_L_up; pathDB.alfa_L_up_mean(1:stop,:,i)];
            beta_L_up = [beta_L_up; pathDB.beta_L_up_mean(1:stop,:,i)];
            
            u_L_up = [u_L_up; pathDB.u_L_up_mean(1:stop,:,i)];
            v_L_up = [v_L_up; pathDB.v_L_up_mean(1:stop,:,i)];
            w_L_up = [w_L_up; pathDB.w_L_up_mean(1:stop,:,i)];
            U_L_up = [U_L_up; pathDB.U_L_up_mean(1:stop,:,i)];
            
            omegax_L_up = [omegax_L_up; pathDB.omegax_L_up_mean(1:stop,:,i)];
            omegay_L_up = [omegay_L_up; pathDB.omegay_L_up_mean(1:stop,:,i)];
            omegaz_L_up = [omegaz_L_up; pathDB.omegaz_L_up_mean(1:stop,:,i)];
            Omega_L_up = [Omega_L_up; pathDB.Omega_L_up_mean(1:stop,:,i)];
            
            phi_L_down = [phi_L_down; pathDB.phi_L_down_mean(1:stop,:,i)];
            theta_L_down = [theta_L_down; pathDB.theta_L_down_mean(1:stop,:,i)];
            eta_L_down = [eta_L_down; pathDB.eta_L_down_mean(1:stop,:,i)];
            
            alfa_L_down = [alfa_L_down; pathDB.alfa_L_down_mean(1:stop,:,i)];
            beta_L_down = [beta_L_down; pathDB.beta_L_down_mean(1:stop,:,i)];
            
            u_L_down = [u_L_down; pathDB.u_L_down_mean(1:stop,:,i)];
            v_L_down = [v_L_down; pathDB.v_L_down_mean(1:stop,:,i)];
            w_L_down = [w_L_down; pathDB.w_L_down_mean(1:stop,:,i)];
            U_L_down = [U_L_down; pathDB.U_L_down_mean(1:stop,:,i)];
            
            omegax_L_down = [omegax_L_down; pathDB.omegax_L_down_mean(1:stop,:,i)];
            omegay_L_down = [omegay_L_down; pathDB.omegay_L_down_mean(1:stop,:,i)];
            omegaz_L_down = [omegaz_L_down; pathDB.omegaz_L_down_mean(1:stop,:,i)];
            Omega_L_down = [Omega_L_down; pathDB.Omega_L_down_mean(1:stop,:,i)];
            
            % Right wing
            
            phi_R_up = [phi_R_up; pathDB.phi_R_up_mean(1:stop,:,i)];
            theta_R_up = [theta_R_up; pathDB.theta_R_up_mean(1:stop,:,i)];
            eta_R_up = [eta_R_up; pathDB.eta_R_up_mean(1:stop,:,i)];
            
            alfa_R_up = [alfa_R_up; pathDB.alfa_R_up_mean(1:stop,:,i)];
            beta_R_up = [beta_R_up; pathDB.beta_R_up_mean(1:stop,:,i)];
            
            u_R_up = [u_R_up; pathDB.u_R_up_mean(1:stop,:,i)];
            v_R_up = [v_R_up; pathDB.v_R_up_mean(1:stop,:,i)];
            w_R_up = [w_R_up; pathDB.w_R_up_mean(1:stop,:,i)];
            U_R_up = [U_R_up; pathDB.U_R_up_mean(1:stop,:,i)];
            
            omegax_R_up = [omegax_R_up; pathDB.omegax_R_up_mean(1:stop,:,i)];
            omegay_R_up = [omegay_R_up; pathDB.omegay_R_up_mean(1:stop,:,i)];
            omegaz_R_up = [omegaz_R_up; pathDB.omegaz_R_up_mean(1:stop,:,i)];
            Omega_R_up = [Omega_R_up; pathDB.Omega_R_up_mean(1:stop,:,i)];
            
            phi_R_down = [phi_R_down; pathDB.phi_R_down_mean(1:stop,:,i)];
            theta_R_down = [theta_R_down; pathDB.theta_R_down_mean(1:stop,:,i)];
            eta_R_down = [eta_R_down; pathDB.eta_R_down_mean(1:stop,:,i)];
            
            alfa_R_down = [alfa_R_down; pathDB.alfa_R_down_mean(1:stop,:,i)];
            beta_R_down = [beta_R_down; pathDB.beta_R_down_mean(1:stop,:,i)];
            
            u_R_down = [u_R_down; pathDB.u_R_down_mean(1:stop,:,i)];
            v_R_down = [v_R_down; pathDB.v_R_down_mean(1:stop,:,i)];
            w_R_down = [w_R_down; pathDB.w_R_down_mean(1:stop,:,i)];
            U_R_down = [U_R_down; pathDB.U_R_down_mean(1:stop,:,i)];
            
            omegax_R_down = [omegax_R_down; pathDB.omegax_R_down_mean(1:stop,:,i)];
            omegay_R_down = [omegay_R_down; pathDB.omegay_R_down_mean(1:stop,:,i)];
            omegaz_R_down = [omegaz_R_down; pathDB.omegaz_R_down_mean(1:stop,:,i)];
            Omega_R_down = [Omega_R_down; pathDB.Omega_R_down_mean(1:stop,:,i)];
 

        end


        

        
    end
 

%         figure()
%         hist(phi_body(:,1),50);
%         figure()
%         hist(theta_body(:,1),50);
%         figure()
%         hist(xsi_body(:,1),50);
%         
%         figure()
%         hist(omegax_body(:,1),50);
%         figure()
%         hist(omegay_body(:,1),50);
%         figure()
%         hist(omegaz_body(:,1),50);
%         figure()
%         hist(Omega_body(:,1),50);
% 
%         figure()
%         hist(u_body(:,1),50);
%         figure()
%         hist(v_body(:,1),50);
%         figure()
%         hist(w_body(:,1),50);
%         figure()
%         hist(U_body(:,1),50);
%         
%         figure()
%         hist(ax_body(:,1),50);
%         figure()
%         hist(ay_body(:,1),50);
%         figure()
%         hist(az_body(:,1),50);
%         figure()
%         hist(a_body(:,1),50);

        nr_bars = 150;


        figure()
        subplot(4,1,1); hist(radtodeg(alfa_L_up(:,1)),0.5:1:89.5)
        title('Histograms of angle of attack at the wings.')
        xlabel('Alfa left upstroke [deg]')
        ylabel('frequency')
        axis([0 90 0 200])
        subplot(4,1,2); hist(radtodeg(alfa_L_down(:,1)),-89.5:1:-0.5)
        xlabel('Alfa left downstroke [deg]')
        ylabel('frequency')
        axis([-90 0 0 200])
        subplot(4,1,3); hist(radtodeg(alfa_R_up(:,1)),0.5:1:89.5)
        xlabel('Alfa right upstroke [deg]')
        ylabel('frequency')
        axis([0 90 0 200])
        subplot(4,1,4); hist(radtodeg(alfa_R_down(:,1)),-89.5:1:-0.5)
        xlabel('Alfa right downstroke [deg]')
        ylabel('frequency')
        axis([-90 0 0 200])
        
        figure()
        subplot(4,1,1); hist(radtodeg(beta_L_up(:,1)),-49.5:0.75:49.5)
        title('Histograms of angle of sideslip on the wings.')
        xlabel('Beta left upstroke [deg]')
        ylabel('frequency')
        axis([-50 50 0 250])
        subplot(4,1,2); hist(radtodeg(beta_L_down(:,1)),-49.5:0.75:49.5)
        xlabel('Beta left downstroke [deg]')
        ylabel('frequency')
        axis([-50 50 0 250])
        subplot(4,1,3); hist(radtodeg(beta_R_up(:,1)),-49.5:0.75:49.5)
        xlabel('Beta right upstroke [deg]')
        ylabel('frequency')
        axis([-50 50 0 250])
        subplot(4,1,4); hist(radtodeg(beta_R_down(:,1)),-49.5:0.75:49.5)
        xlabel('Beta right downstroke [deg]')
        ylabel('frequency')
        axis([-50 50 0 250])
        
        figure()
        subplot(4,1,1); hist(u_L_up(:,1),-480:40:4480)
        title('Histograms wing velocity in x-direction of wing ref. frame')
        xlabel('u left upstroke [mm/s]')
        ylabel('frequency')
        axis([-500 4000 0 250])
        subplot(4,1,2); hist(u_L_down(:,1),-480:40:4480)
        xlabel('u left downstroke [mm/s]')
        ylabel('frequency')
        axis([-500 4000 0 250])
        subplot(4,1,3); hist(u_R_up(:,1),-480:40:4480)
        xlabel('u right upstroke [mm/s]')
        ylabel('frequency')
        axis([-500 4000 0 250])
        subplot(4,1,4); hist(u_R_down(:,1),-480:40:4480)
        xlabel('u right downstroke [mm/s]')
        ylabel('frequency')
        axis([-500 4000 0 250])
        
        figure()
        subplot(4,1,1); hist(v_L_up(:,1),-298:4:298)
        title('Histograms wing velocity in y-direction of wing ref. frame')
        xlabel('v left upstroke [mm/s]')
        ylabel('frequency')
        axis([-300 300 0 100])
        subplot(4,1,2); hist(v_L_down(:,1),-298:4:298)
        xlabel('v left downstroke [mm/s]')
        ylabel('frequency')
        axis([-300 300 0 100])
        subplot(4,1,3); hist(v_R_up(:,1),-298:4:298)
        xlabel('v right upstroke [mm/s]')
        ylabel('frequency')
        axis([-300 300 0 100])
        subplot(4,1,4); hist(v_R_down(:,1),-298:4:298)
        xlabel('v right downstroke [mm/s]')
        ylabel('frequency')
        axis([-300 300 0 100])
        
        figure()
        subplot(4,1,1); hist(w_L_up(:,1),1007:14:2993)
        title('Histograms wing velocity in z-direction of wing ref. frame')
        xlabel('w left upstroke [mm/s]')
        ylabel('frequency')
        axis([1000 3000 0 100])
        subplot(4,1,2); hist(w_L_down(:,1),-2993:14:-1007)
        xlabel('w left downstroke [mm/s]')
        ylabel('frequency')
        axis([-3000 -1000 0 100])
        subplot(4,1,3); hist(w_R_up(:,1),1007:14:2993)
        xlabel('w right upstroke [mm/s]')
        ylabel('frequency')
        axis([1000 3000 0 100])
        subplot(4,1,4); hist(w_R_down(:,1),-2993:14:-1007)
        xlabel('w right downstroke [mm/s]')
        ylabel('frequency')
        axis([-3000 -1000 0 100])
        
        figure()
        subplot(4,1,1); hist(U_L_up(:,1),1510:20:4490)
        title('Histograms absolute wing velocity')
        xlabel('|U| left upstroke [mm/s]')
        ylabel('frequency')
        axis([1500 4500 0 100])
        subplot(4,1,2); hist(U_L_down(:,1),1510:20:4490)
        xlabel('|U| left downstroke [mm/s]')
        ylabel('frequency')
        axis([1500 4500 0 100])
        subplot(4,1,3); hist(U_R_up(:,1),1510:20:4490)
        xlabel('|U| right upstroke [mm/s]')
        ylabel('frequency')
        axis([1500 4500 0 100])
        subplot(4,1,4); hist(U_R_down(:,1),1510:20:4490)
        xlabel('|U| right downstroke [mm/s]')
        ylabel('frequency')
        axis([1500 4500 0 100])
        

        



        

    
    
    
    
end


