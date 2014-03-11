function wb_raw_filt(settings,seq_nr,qb1_f, qb2_f, qb3_f, qb4_f, wbx_f, wby_f, wbz_f, t, fig_nr1, fig_nr2,save_on_off)

    %Plots omega computed by Kalman filter and omega derived the DCM's over
    %dt.
    
    dt = t(2)-t(1);
    
    N = length(qb1_f);
    
    % Calculate the yaw, pitch and roll
    
    rn = real(atan2(2*(qb4_f.*qb1_f+qb2_f.*qb3_f),1-2*(qb1_f.^2+qb2_f.^2)));
    pn = real(asin(2*(qb4_f.*qb2_f-qb3_f.*qb1_f)));
    yn = real(atan2(2*(qb4_f.*qb3_f+qb1_f.*qb2_f),1-2*(qb2_f.^2+qb3_f.^2)));
    

    y = yn;
    r = rn;
    p = pn;
    
    om_r = zeros(1,N);
    om_p = zeros(1,N);
    om_y = zeros(1,N);

    om_r(1:N-1) = (r(2:N)-r(1:N-1))/(dt);
    om_p(1:N-1) = (p(2:N)-p(1:N-1))/(dt);
    om_y(1:N-1) = (y(2:N)-y(1:N-1))/(dt);
    
    om_g = zeros(3,N);


    for i = 1:N-1
        
        J = [1 sin(r(i))*tan(p(i)) cos(r(i))*tan(p(i)); ...
             0 cos(r(i)) -sin(r(i)); ...
             0 sin(r(i))/cos(p(i)) cos(r(i))/cos(p(i))];
        
%        R_i   = [ qb4_f(i)^2+qb1_f(i)^2-qb2_f(i)^2-qb3_f(i)^2    2*qb1_f(i)*qb2_f(i)+2*qb4_f(i)*qb3_f(i)      2*qb1_f(i)*qb3_f(i)-2*qb4_f(i)*qb2_f(i); ...
%                 2*qb1_f(i)*qb2_f(i)-2*qb4_f(i)*qb3_f(i)        qb4_f(i)^2-qb1_f(i)^2+qb2_f(i)^2-qb3_f(i)^2  2*qb2_f(i)*qb3_f(i)+2*qb4_f(i)*qb1_f(i); ...
%                 2*qb1_f(i)*qb3_f(i)+2*qb4_f(i)*qb2_f(i)        2*qb2_f(i)*qb3_f(i)-2*qb4_f(i)*qb1_f(i)      qb4_f(i)^2-qb1_f(i)^2-qb2_f(i)^2+qb3_f(i)^2];
%          
%        R_ip1   = [ qb4_f(i+1)^2+qb1_f(i+1)^2-qb2_f(i+1)^2-qb3_f(i+1)^2    2*qb1_f(i+1)*qb2_f(i+1)+2*qb4_f(i+1)*qb3_f(i+1)      2*qb1_f(i+1)*qb3_f(i+1)-2*qb4_f(i+1)*qb2_f(i+1); ...
%              2*qb1_f(i+1)*qb2_f(i+1)-2*qb4_f(i+1)*qb3_f(i+1)        qb4_f(i+1)^2-qb1_f(i+1)^2+qb2_f(i+1)^2-qb3_f(i+1)^2  2*qb2_f(i+1)*qb3_f(i+1)+2*qb4_f(i+1)*qb1_f(i+1); ...
%              2*qb1_f(i+1)*qb3_f(i+1)+2*qb4_f(i+1)*qb2_f(i+1)        2*qb2_f(i+1)*qb3_f(i+1)-2*qb4_f(i+1)*qb1_f(i+1)      qb4_f(i+1)^2-qb1_f(i+1)^2-qb2_f(i+1)^2+qb3_f(i+1)^2];
%          
%        %R_new = R_i*R_ip1';
%        
%        DR_Dt = (R_i-R_ip1)./dt;
%        
%        W_cross = R_i'*DR_Dt;
%        
%        om_g(:,i) = R_i*[(W_cross(3,2)-W_cross(2,3))/2; (-W_cross(3,1)+W_cross(1,3))/2; (W_cross(2,1)-W_cross(1,2))/2];

         
       om_g(:,i) = (J\[om_r(i); om_p(i); om_y(i)]);

    end
    
    %om_g(:,1) = om_g(:,2);
    om_g(:,N) = om_g(:,N-1);
    
    for i = 1:N

        if abs(om_g(1,i)) > 1000
            om_g(1,i) = 0;
        end
        
        if abs(om_g(2,i)) > 1000
            om_g(2,i) = 0;
        end
        
        if abs(om_g(3,i)) > 1000
            om_g(3,i) = 0;
        end
        
    end

    
    figure(fig_nr1)
    subplot(3,1,1); plot(t,om_g(1,:)',t,wbx_f,t,0)
    title('Filtered and unfiltered angular rates over time')
    xlabel('t [s]')
    ylabel('omega_x [rad/s]')
    legend('raw','filt')
    subplot(3,1,2); plot(t,om_g(2,:)',t,wby_f,t,0)
    xlabel('t [s]')
    ylabel('omega_y [rad/s]')
    legend('raw','filt')
    subplot(3,1,3); plot(t,om_g(3,:)',t,wbz_f,t,0)
    xlabel('t [s]')
    ylabel('omega_z [rad/s]')
    legend('raw','filt')
    

    figure(fig_nr2)
    subplot(3,1,1); plot(t,om_g(1,:)'-wbx_f(:),t,0)
    title('Deviation of unfiltered from filtered angular rates')
    xlabel('t [s]')
    ylabel('omega_x [rad/s]')
    legend('raw','filt')
    subplot(3,1,2); plot(t,om_g(2,:)'-wby_f(:),t,0)
    xlabel('t [s]')
    ylabel('omega_y [rad/s]')
    legend('raw','filt')
    subplot(3,1,3); plot(t,om_g(3,:)'-wbz_f(:),t,0)
    xlabel('t [s]')
    ylabel('omega_z [rad/s]')
    legend('raw','filt')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/w_body_1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/w_body/w_body_1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/w_body_2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/w_body/w_body_2_' int2str(seq_nr)], 'fig')
    
    end



end

