function Plot_relations(settings, pathDB)


    % Program that plots relations between the averages of wing kinematic
    % parameters (alfa, velocity, beta, omega, orientation) and the
    % averages of body parameters (velocity, acceleration, omega,
    % orientation, alfa beta) during the time of the downstroke or
    % upstroke.
    
    
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
    
    % Plot relations:
    
    %----------------------------------------------------------------------
    
    % Angle of attack at left and right wing vs absolute velocity body
    
    %----------------------------------------------------------------------
    
    % |alfa_upst_L| vs |U_body|

    temp_p = polyfit(U_body(:,1),radtodeg(abs(alfa_L_up(:,1))),1);
    
    figure()
    scatter(U_body(:,1),radtodeg(abs(alfa_L_up(:,1))))
    title('Alfa upstroke left vs absolute body velocity')
    xlabel('Absolute body velocity [mm/s]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(U_body(:,1),polyval(temp_p,U_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    % |alfa_upst_R| vs |U_body|

    temp_p = polyfit(U_body(:,1),radtodeg(abs(alfa_R_up(:,1))),1);
    
    figure()
    scatter(U_body(:,1),radtodeg(abs(alfa_R_up(:,1))))
    title('Alfa upstroke right vs absolute body velocity')
    xlabel('Absolute body velocity [mm/s]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(U_body(:,1),polyval(temp_p,U_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    % |alfa_dwnst_L| vs |U_body|

    temp_p = polyfit(U_body(:,1),radtodeg(abs(alfa_L_down(:,1))),1);
    
    figure()
    scatter(U_body(:,1),radtodeg(abs(alfa_L_down(:,1))))
    title('Alfa downstroke left vs absolute body velocity')
    xlabel('Absolute body velocity [mm/s]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(U_body(:,1),polyval(temp_p,U_body(:,1)),'r')
    hold off

    
    clear temp_p
    
     % |alfa_dwnst_R| vs |U_body|

    temp_p = polyfit(U_body(:,1),radtodeg(abs(alfa_R_down(:,1))),1);
    
    figure()
    scatter(U_body(:,1),radtodeg(abs(alfa_R_down(:,1))))
    title('Alfa downstroke right vs absolute body velocity')
    xlabel('Absolute body velocity [mm/s]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(U_body(:,1),polyval(temp_p,U_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    %----------------------------------------------------------------------
    
    % Angle of attack left and right wing vs absolute body acceleration
    
    %----------------------------------------------------------------------
    
    
    % |alfa_upst_L| vs |a_body|

    temp_p = polyfit(a_body(:,1),radtodeg(abs(alfa_L_up(:,1))),1);
    
    figure()
    scatter(a_body(:,1),radtodeg(abs(alfa_L_up(:,1))))
    title(['Alfa upstroke left vs absolute body acceleration, N=' int2str(length(alfa_L_up(:,1)))])
    xlabel('Absolute body acceleration [mm/s^2]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(a_body(:,1),polyval(temp_p,a_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    % |alfa_upst_R| vs |a_body|

    temp_p = polyfit(a_body(:,1),radtodeg(abs(alfa_R_up(:,1))),1);
    
    figure()
    scatter(a_body(:,1),radtodeg(abs(alfa_R_up(:,1))))
    title(['Alfa upstroke right vs absolute body acceleration, N=' int2str(length(alfa_R_up(:,1)))])
    xlabel('Absolute body acceleration [mm/s^2]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(a_body(:,1),polyval(temp_p,a_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    % |alfa_dwnst_L| vs |a_body|

    temp_p = polyfit(a_body(:,1),radtodeg(abs(alfa_L_down(:,1))),1);
    
    figure()
    scatter(a_body(:,1),radtodeg(abs(alfa_L_down(:,1))))
    title(['Alfa downstroke left vs absolute body acceleration, N=' int2str(length(alfa_L_down(:,1)))])
    xlabel('Absolute body acceleration [mm/s^2]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(a_body(:,1),polyval(temp_p,a_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    % |alfa_dwnst_R| vs |a_body|

    temp_p = polyfit(a_body(:,1),radtodeg(abs(alfa_R_down(:,1))),1);
    
    figure()
    scatter(a_body(:,1),radtodeg(abs(alfa_R_down(:,1))))
    title(['Alfa donwstroke right vs absolute body acceleration, N=' int2str(length(alfa_R_down(:,1)))])
    xlabel('Absolute body acceleration [mm/s^2]')
    ylabel('Angle of attack [deg]')
    hold on
    plot(a_body(:,1),polyval(temp_p,a_body(:,1)),'r')
    hold off

    
    clear temp_p
    
    
    
    %----------------------------------------------------------------------
    
    % Angle of attack vs horizontal body acceleration
   
    %----------------------------------------------------------------------

    temp_p1 = polyfit(u_body(:,1),radtodeg(alfa_L_up(:,1)),1);
    temp_p2 = polyfit(u_body(:,1),radtodeg(alfa_R_up(:,1)),1);
    temp_p3 = polyfit(u_body(:,1),radtodeg(alfa_L_down(:,1)),1);
    temp_p4 = polyfit(u_body(:,1),radtodeg(alfa_R_down(:,1)),1);
    

    figure()
    subplot(2,2,1); scatter(u_body(:,1),radtodeg(alfa_L_up(:,1)))
    title('Scatter plot horizontal body velocity vs angle of attack at the wing')
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa upstroke left [deg]')
    axis([-200 600 0 90])
    hold on
    plot(u_body(:,1),polyval(temp_p2,u_body(:,1)),'r')
    hold off
    subplot(2,2,2); scatter(u_body(:,1),radtodeg(alfa_R_up(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa upstroke right [deg]')
    axis([-200 600 0 90])
    hold on
    plot(u_body(:,1),polyval(temp_p2,u_body(:,1)),'r')
    hold off
    subplot(2,2,3); scatter(u_body(:,1),-radtodeg(alfa_L_down(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa downstroke left [deg]')
    axis([-200 600 0 90])
    hold on
    plot(u_body(:,1),-polyval(temp_p3,u_body(:,1)),'r')
    hold off
    subplot(2,2,4); scatter(u_body(:,1),-radtodeg(alfa_R_down(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa downstroke right [deg]')
    axis([-200 600 0 90])
    hold on
    plot(u_body(:,1),-polyval(temp_p4,u_body(:,1)),'r')
    hold off

    clear temp_p1 temp_p2 temp_p3 temp_p4
    
    
    
        %----------------------------------------------------------------------
    
    % Angle of attack vs horizontal body acceleration
   
    %----------------------------------------------------------------------

    temp_p1 = polyfit(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_L_up(:,1)),1);
    temp_p2 = polyfit(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_R_up(:,1)),1);
    temp_p3 = polyfit(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_L_down(:,1)),1);
    temp_p4 = polyfit(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_R_down(:,1)),1);
    

    figure()
    subplot(2,2,1); scatter(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_L_up(:,1)))
    title('Scatter plot horizontal body velocity vs angle of attack at the wing')
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa upstroke left [deg]')
    axis([-200 600 0 90])
    hold on
    plot(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),polyval(temp_p1,sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1)))),'r')
    hold off
    subplot(2,2,2); scatter(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),radtodeg(alfa_R_up(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa upstroke right [deg]')
    axis([-200 600 0 90])
    hold on
    plot(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),polyval(temp_p2,sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1)))),'r')
    hold off
    subplot(2,2,3); scatter(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),-radtodeg(alfa_L_down(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa downstroke left [deg]')
    axis([-200 600 0 90])
    hold on
    plot(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),-polyval(temp_p3,sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1)))),'r')
    hold off
    subplot(2,2,4); scatter(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),-radtodeg(alfa_R_down(:,1)))
    xlabel('horizontal velocity [mm/s]')
    ylabel('alfa downstroke right [deg]')
    axis([-200 600 0 90])
    hold on
    plot(sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1))),-polyval(temp_p4,sqrt(u_body(:,1).^2+w_body(:,1).^2).*real(atan2(w_body(:,1),u_body(:,1)))),'r')
    hold off

    clear temp_p1 temp_p2 temp_p3 temp_p4
    
    
    
    
    %----------------------------------------------------------------------
    
    % Angle of attack vs horizontal body acceleration
   
    %----------------------------------------------------------------------

    temp_p1 = polyfit(ax_body(:,1),radtodeg(alfa_L_up(:,1)),1);
    temp_p2 = polyfit(ax_body(:,1),radtodeg(alfa_R_up(:,1)),1);
    temp_p3 = polyfit(ax_body(:,1),radtodeg(alfa_L_down(:,1)),1);
    temp_p4 = polyfit(ax_body(:,1),radtodeg(alfa_R_down(:,1)),1);
    

    figure()
    subplot(2,2,1); scatter(ax_body(:,1),radtodeg(alfa_L_up(:,1)))
    title('Scatter plot horizontal body acceleration vs angle of attack at the wing')
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa upstroke left [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),polyval(temp_p2,ax_body(:,1)),'r')
    hold off
    subplot(2,2,2); scatter(ax_body(:,1),radtodeg(alfa_R_up(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa upstroke right [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),polyval(temp_p2,ax_body(:,1)),'r')
    hold off
    subplot(2,2,3); scatter(ax_body(:,1),-radtodeg(alfa_L_down(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa downstroke left [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),-polyval(temp_p3,ax_body(:,1)),'r')
    hold off
    subplot(2,2,4); scatter(ax_body(:,1),-radtodeg(alfa_R_down(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa downstroke right [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),-polyval(temp_p4,ax_body(:,1)),'r')
    hold off

    clear temp_p1 temp_p2 temp_p3 temp_p4
    
    

    
    %----------------------------------------------------------------------
    
    % wing velocity vs body speed
   
    %----------------------------------------------------------------------

    temp_p1 = polyfit(ax_body(:,1),radtodeg(alfa_L_up(:,1)),1);
    temp_p2 = polyfit(ax_body(:,1),radtodeg(alfa_R_up(:,1)),1);
    temp_p3 = polyfit(ax_body(:,1),radtodeg(alfa_L_down(:,1)),1);
    temp_p4 = polyfit(ax_body(:,1),radtodeg(alfa_R_down(:,1)),1);
    

    figure()
    subplot(2,2,1); scatter(ax_body(:,1),radtodeg(alfa_L_up(:,1)))
    title('Scatter plot horizontal body acceleration vs angle of attack at the wing')
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa upstroke left [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),polyval(temp_p2,ax_body(:,1)),'r')
    hold off
    subplot(2,2,2); scatter(ax_body(:,1),radtodeg(alfa_R_up(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa upstroke right [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),polyval(temp_p2,ax_body(:,1)),'r')
    hold off
    subplot(2,2,3); scatter(ax_body(:,1),-radtodeg(alfa_L_down(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa downstroke left [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),-polyval(temp_p3,ax_body(:,1)),'r')
    hold off
    subplot(2,2,4); scatter(ax_body(:,1),-radtodeg(alfa_R_down(:,1)))
    xlabel('horizontal acceleration [mm/s^2]')
    ylabel('alfa downstroke right [deg]')
    axis([-15000 15000 0 90])
    hold on
    plot(ax_body(:,1),-polyval(temp_p4,ax_body(:,1)),'r')
    hold off

    clear temp_p1 temp_p2 temp_p3 temp_p4
    
    
    
    
end

