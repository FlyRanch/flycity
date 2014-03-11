function wR_raw_filt2(settings, seq_nr,qR1_r, qR2_r, qR3_r, qR4_r, qR1_f1, qR2_f1, qR3_f1, qR4_f1, qR1_f2, qR2_f2, qR3_f2, qR4_f2, wRx_f2, wRy_f2, wRz_f2, t, fig_nr1, fig_nr2, fig_nr3,save_on_off)

    %Plots omega computed by Kalman filter and omega derived the DCM's over
    %dt.
    
    dt = t(2)-t(1);
    
    N = length(qR1_r);
    
    % Calculate the yaw, pitch and roll
    
       
   

    
%     % Calculate the yaw, pitch and roll
%     
%     % Raw wing quaternion angular velocities:
%     
%     r_r = real(atan2(2*(qR4_r.*qR1_r+qR2_r.*qR3_r),1-2*(qR1_r.^2+qR2_r.^2)));
%     p_r = real(asin(2*(qR4_r.*qR2_r-qR3_r.*qR1_r)));
%     y_r = real(atan2(2*(qR4_r.*qR3_r+qR1_r.*qR2_r),1-2*(qR2_r.^2+qR3_r.^2)));
%     
%     om_r_r = zeros(1,N);
%     om_p_r = zeros(1,N);
%     om_y_r = zeros(1,N);
% 
%     om_r_r(2:N-1) = (r_r(3:N)-r_r(1:N-2))/(2*dt);
%     om_p_r(2:N-1) = (p_r(3:N)-p_r(1:N-2))/(2*dt);
%     om_y_r(2:N-1) = (y_r(3:N)-y_r(1:N-2))/(2*dt);
%     
%     %for i = 1:N
% 
%     om_g_r = zeros(3,N);
%     
%     for i = 1:N
%         
%        qR1_r(i) = qR1_r(i)/norm([qR1_r(i); qR2_r(i); qR3_r(i); qR4_r(i)]);
%        qR2_r(i) = qR2_r(i)/norm([qR1_r(i); qR2_r(i); qR3_r(i); qR4_r(i)]);
%        qR3_r(i) = qR3_r(i)/norm([qR1_r(i); qR2_r(i); qR3_r(i); qR4_r(i)]);
%        qR4_r(i) = qR4_r(i)/norm([qR1_r(i); qR2_r(i); qR3_r(i); qR4_r(i)]);
%         
%        J = [1 sin(r_r(i))*tan(p_r(i)) cos(r_r(i))*tan(p_r(i)); ...
%              0 cos(r_r(i)) -sin(r_r(i)); ...
%              0 sin(r_r(i))/cos(p_r(i)) cos(r_r(i))/cos(p_r(i))];
%          
% %          det(J)
% %          
% %        R   = [ qR4_r(i)^2+qR1_r(i)^2-qR2_r(i)^2-qR3_r(i)^2    2*qR1_r(i)*qR2_r(i)+2*qR4_r(i)*qR3_r(i)      2*qR1_r(i)*qR3_r(i)-2*qR4_r(i)*qR2_r(i); ...
% %                2*qR1_r(i)*qR2_r(i)-2*qR4_r(i)*qR3_r(i)        qR4_r(i)^2-qR1_r(i)^2+qR2_r(i)^2-qR3_r(i)^2  2*qR2_r(i)*qR3_r(i)+2*qR4_r(i)*qR1_r(i); ...
% %                2*qR1_r(i)*qR3_r(i)+2*qR4_r(i)*qR2_r(i)        2*qR2_r(i)*qR3_r(i)-2*qR4_r(i)*qR1_r(i)      qR4_r(i)^2-qR1_r(i)^2-qR2_r(i)^2+qR3_r(i)^2];
% 
%         
%         om_g_r(:,i) = (J\[om_r_r(i); om_p_r(i); om_y_r(i)]);
%         
%         clear J R
% 
%     end
%     
%     om_g_r(:,N) = om_g_r(:,N-1);
%     
%     for i = 1:N-1
% 
%         if abs(om_g_r(1,i)) > 10000
%             om_g_r(1,i) = 0;
%         end
%         
%         if abs(om_g_r(2,i)) > 10000
%             om_g_r(2,i) = 0;
%         end
%         
%         if abs(om_g_r(3,i)) > 10000
%             om_g_r(3,i) = 0;
%         end
%         
%     end
    
    
    % Calculate the quaternion derivative in another way:
    
    om_r_2 = zeros(3,N);
    
    for k = 2:N-1
    
    q_1 = [qR1_r(k-1); qR2_r(k-1); qR3_r(k-1); qR4_r(k-1)];
    
    q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);
    
    q_1_inv = q_1_inv./norm(q_1_inv);
    
    q_2 = [qR1_r(k+1); qR2_r(k+1); qR3_r(k+1); qR4_r(k+1)];
    
    Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
           -q_2(3) q_2(4) q_2(1) q_2(2);
           q_2(2) -q_2(1) q_2(4) q_2(3);
           -q_2(1) -q_2(2) -q_2(3) q_2(4)];    
            
    q_t = Q_2*q_1_inv;
    
    q_t = q_t./norm(q_t);
    
    omega_norm = (2/dt)*acos(q_t(4));
    
    om_r_2(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
    
    clear q_1 q_1_inv q_2 Q_2 q_t omega_norm
    
    end
    
    om_r_2(:,1) = om_r_2(:,2);
    om_r_2(:,N) = om_r_2(:,N-1);
    
    
    
    
    
%     % Body motion filtered wing quaternion angular velocities:
%     
%     r_f1 = real(atan2(2*(qR4_f1.*qR1_f1+qR2_f1.*qR3_f1),1-2*(qR1_f1.^2+qR2_f1.^2)));
%     p_f1 = real(asin(2*(qR4_f1.*qR2_f1-qR3_f1.*qR1_f1)));
%     y_f1 = real(atan2(2*(qR4_f1.*qR3_f1+qR1_f1.*qR2_f1),1-2*(qR2_f1.^2+qR3_f1.^2)));
%     
%     om_r_f1 = zeros(1,N);
%     om_p_f1 = zeros(1,N);
%     om_y_f1 = zeros(1,N);
% 
%     om_r_f1(2:N-1) = (r_f1(3:N)-r_f1(1:N-2))/(2*dt);
%     om_p_f1(2:N-1) = (p_f1(3:N)-p_f1(1:N-2))/(2*dt);
%     om_y_f1(2:N-1) = (y_f1(3:N)-y_f1(1:N-2))/(2*dt);
% 
%     om_g_f1 = zeros(3,N);
%     
%     for i = 1:N
%         
%        qR1_f1(i) = qR1_f1(i)/norm([qR1_f1(i); qR2_f1(i); qR3_f1(i); qR4_f1(i)]);
%        qR2_f1(i) = qR2_f1(i)/norm([qR1_f1(i); qR2_f1(i); qR3_f1(i); qR4_f1(i)]);
%        qR3_f1(i) = qR3_f1(i)/norm([qR1_f1(i); qR2_f1(i); qR3_f1(i); qR4_f1(i)]);
%        qR4_f1(i) = qR4_f1(i)/norm([qR1_f1(i); qR2_f1(i); qR3_f1(i); qR4_f1(i)]);
%         
%         J = [1 sin(r_f1(i))*tan(p_f1(i)) cos(r_f1(i))*tan(p_f1(i)); ...
%              0 cos(r_f1(i)) -sin(r_f1(i)); ...
%              0 sin(r_f1(i))/cos(p_f1(i)) cos(r_f1(i))/cos(p_f1(i))];
%          
% %        R   = [ qR4_f1(i)^2+qR1_f1(i)^2-qR2_f1(i)^2-qR3_f1(i)^2    2*qR1_f1(i)*qR2_f1(i)+2*qR4_f1(i)*qR3_f1(i)      2*qR1_f1(i)*qR3_f1(i)-2*qR4_f1(i)*qR2_f1(i); ...
% %                2*qR1_f1(i)*qR2_f1(i)-2*qR4_f1(i)*qR3_f1(i)        qR4_f1(i)^2-qR1_f1(i)^2+qR2_f1(i)^2-qR3_f1(i)^2  2*qR2_f1(i)*qR3_f1(i)+2*qR4_f1(i)*qR1_f1(i); ...
% %                2*qR1_f1(i)*qR3_f1(i)+2*qR4_f1(i)*qR2_f1(i)        2*qR2_f1(i)*qR3_f1(i)-2*qR4_f1(i)*qR1_f1(i)      qR4_f1(i)^2-qR1_f1(i)^2-qR2_f1(i)^2+qR3_f1(i)^2];
% 
%         
%         om_g_f1(:,i) = (J\[om_r_f1(i); om_p_f1(i); om_y_f1(i)]);
%         
%         clear J R
% 
%     end
%     
%     om_g_f1(:,N) = om_g_f1(:,N-1);
%     
%     for i = 1:N-1
% 
%         if abs(om_g_f1(1,i)) > 10000
%             om_g_f1(1,i) = 0;
%         end
%         
%         if abs(om_g_f1(2,i)) > 10000
%             om_g_f1(2,i) = 0;
%         end
%         
%         if abs(om_g_f1(3,i)) > 10000
%             om_g_f1(3,i) = 0;
%         end
%         
%     end
    
        % Calculate the quaternion derivative in another way:
    
    om_f1_2 = zeros(3,N);
    
    for k = 2:N-1
    
    q_1 = [qR1_f1(k-1); qR2_f1(k-1); qR3_f1(k-1); qR4_f1(k-1)];
    
    q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);
    
    q_1_inv = q_1_inv./norm(q_1_inv);
    
    q_2 = [qR1_f1(k+1); qR2_f1(k+1); qR3_f1(k+1); qR4_f1(k+1)];
    
    Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
           -q_2(3) q_2(4) q_2(1) q_2(2);
           q_2(2) -q_2(1) q_2(4) q_2(3);
           -q_2(1) -q_2(2) -q_2(3) q_2(4)];    
            
    q_t = Q_2*q_1_inv;
    
    q_t = q_t./norm(q_t);
    
    omega_norm = (2/dt)*acos(q_t(4));
    
    om_f1_2(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
    
    clear q_1 q_1_inv q_2 Q_2 q_t omega_norm
    
    end
    
    om_f1_2(:,1) = om_f1_2(:,2);
    om_f1_2(:,N) = om_f1_2(:,N-1);
    
    
%     % Body motion filtered wing quaternion angular velocities:
%     
%     r_f2 = real(atan2(2*(qR4_f2.*qR1_f2+qR2_f2.*qR3_f2),1-2*(qR1_f2.^2+qR2_f2.^2)));
%     p_f2 = real(asin(2*(qR4_f2.*qR2_f2-qR3_f2.*qR1_f2)));
%     y_f2 = real(atan2(2*(qR4_f2.*qR3_f2+qR1_f2.*qR2_f2),1-2*(qR2_f2.^2+qR3_f2.^2)));
%     
%     om_r_f2 = zeros(1,N);
%     om_p_f2 = zeros(1,N);
%     om_y_f2 = zeros(1,N);
% 
%     om_r_f2(2:N-1) = (r_f2(3:N)-r_f2(1:N-2))/(2*dt);
%     om_p_f2(2:N-1) = (p_f2(3:N)-p_f2(1:N-2))/(2*dt);
%     om_y_f2(2:N-1) = (y_f2(3:N)-y_f2(1:N-2))/(2*dt);
% 
%     om_g_f2 = zeros(3,N);
%     
%     for i = 1:N
%         
%        qR1_f2(i) = qR1_f2(i)/norm([qR1_f2(i); qR2_f2(i); qR3_f2(i); qR4_f2(i)]);
%        qR2_f2(i) = qR2_f2(i)/norm([qR1_f2(i); qR2_f2(i); qR3_f2(i); qR4_f2(i)]);
%        qR3_f2(i) = qR3_f2(i)/norm([qR1_f2(i); qR2_f2(i); qR3_f2(i); qR4_f2(i)]);
%        qR4_f2(i) = qR4_f2(i)/norm([qR1_f2(i); qR2_f2(i); qR3_f2(i); qR4_f2(i)]);
%         
%         J = [1 sin(r_f2(i))*tan(p_f2(i)) cos(r_f2(i))*tan(p_f2(i)); ...
%              0 cos(r_f2(i)) -sin(r_f2(i)); ...
%              0 sin(r_f2(i))/cos(p_f2(i)) cos(r_f2(i))/cos(p_f2(i))];
%          
% %        R   = [ qR4_f2(i)^2+qR1_f2(i)^2-qR2_f2(i)^2-qR3_f2(i)^2    2*qR1_f2(i)*qR2_f2(i)+2*qR4_f2(i)*qR3_f2(i)      2*qR1_f2(i)*qR3_f2(i)-2*qR4_f2(i)*qR2_f2(i); ...
% %              2*qR1_f2(i)*qR2_f2(i)-2*qR4_f2(i)*qR3_f2(i)        qR4_f2(i)^2-qR1_f2(i)^2+qR2_f2(i)^2-qR3_f2(i)^2  2*qR2_f2(i)*qR3_f2(i)+2*qR4_f2(i)*qR1_f2(i); ...
% %              2*qR1_f2(i)*qR3_f2(i)+2*qR4_f2(i)*qR2_f2(i)        2*qR2_f2(i)*qR3_f2(i)-2*qR4_f2(i)*qR1_f2(i)      qR4_f2(i)^2-qR1_f2(i)^2-qR2_f2(i)^2+qR3_f2(i)^2];
% 
%         
%         om_g_f2(:,i) = (J\[om_r_f2(i); om_p_f2(i); om_y_f2(i)]);
%         
%         clear J R
% 
%     end
%     
%     om_g_f2(:,N) = om_g_f2(:,N-1);
%     
%     for i = 1:N-1
% 
%         if abs(om_g_f2(1,i)) > 10000
%             om_g_f2(1,i) = 0;
%         end
%         
%         if abs(om_g_f2(2,i)) > 10000
%             om_g_f2(2,i) = 0;
%         end
%         
%         if abs(om_g_f2(3,i)) > 10000
%             om_g_f2(3,i) = 0;
%         end
%         
%     end
    
    % Calculate the quaternion derivative in another way:
    
    om_f2_2 = zeros(3,N);
    
    for k = 2:N-1
    
    q_1 = [qR1_f2(k-1); qR2_f2(k-1); qR3_f2(k-1); qR4_f2(k-1)];
    
    q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);
    
    q_1_inv = q_1_inv./norm(q_1_inv);
    
    q_2 = [qR1_f2(k+1); qR2_f2(k+1); qR3_f2(k+1); qR4_f2(k+1)];
    
    Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
           -q_2(3) q_2(4) q_2(1) q_2(2);
           q_2(2) -q_2(1) q_2(4) q_2(3);
           -q_2(1) -q_2(2) -q_2(3) q_2(4)];    
            
    q_t = Q_2*q_1_inv;
    
    q_t = q_t./norm(q_t);
    
    omega_norm = (2/dt)*acos(q_t(4));
    
    om_f2_2(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
    
    clear q_1 q_1_inv q_2 Q_2 q_t omega_norm
    
    end

    om_f2_2(:,1) = om_f2_2(:,2);
    om_f2_2(:,N) = om_f2_2(:,N-1);
    
    figure(fig_nr1)
    subplot(3,1,1); plot(t,om_r_2(1,:),t,om_f1_2(1,:))
    title('Filtered and body filtered angular velocities over time, right wing.')
    xlabel('t [s]')
    ylabel('omega_x [rad/s]')
    legend('raw','body-filtered')
    subplot(3,1,2); plot(t,om_r_2(2,:),t,om_f1_2(2,:))
    xlabel('t [s]')
    ylabel('omega_y [rad/s]')
    legend('raw','body-filtered')
    subplot(3,1,3); plot(t,om_r_2(3,:),t,om_f1_2(3,:))
    xlabel('t [s]')
    ylabel('omega_z [rad/s]')
    legend('raw','body-filtered')
    
    figure(fig_nr2)
    subplot(3,1,1); plot(t,om_f1_2(1,:),t,wRx_f2)
    title('Body filtered and wing filtered angular velocities over time, right wing.')
    xlabel('t [s]')
    ylabel('omega_x [rad/s]')
    legend('body-filtered','wing-filtered')
    subplot(3,1,2); plot(t,om_f1_2(2,:),t,wRy_f2)
    xlabel('t [s]')
    ylabel('omega_y [rad/s]')
    legend('body-filtered','wing-filtered')
    subplot(3,1,3); plot(t,om_f1_2(3,:),t,wRz_f2)
    xlabel('t [s]')
    ylabel('omega_z [rad/s]')
    legend('body-filtered','wing-filtered')
    
    figure(fig_nr3)
    subplot(3,1,1); plot(t,om_f2_2(1,:),t,wRx_f2)
    title('Wing filtered and filtered quaternion derived angular velocities over time, right wing.')
    xlabel('t [s]')
    ylabel('omega_x [rad/s]')
    legend('quat-derivative','wing-filtered')
    subplot(3,1,2); plot(t,om_f2_2(2,:),t,wRy_f2)
    xlabel('t [s]')
    ylabel('omega_y [rad/s]')
    legend('quat-derivative','wing-filtered')
    subplot(3,1,3); plot(t,om_f2_2(3,:),t,wRz_f2)
    xlabel('t [s]')
    ylabel('omega_z [rad/s]')
    legend('quat-derivative','wing-filtered')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/dq_r_dq_f1_R'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/w_right/dq_r_dq_f1_R_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/dq_f1_wR'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/w_right/dq_f1_wR_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/dq_f2_wR'], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(2)) '/w_right/dq_f2_wR_' int2str(seq_nr)], 'fig')
    
    end


end


