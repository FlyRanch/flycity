function Stroke_plane_pte(settings,seq_nr, phi_L, theta_L, eta_L, phi_R, theta_R, eta_R, t, fig_nr1,save_on_off)
    


  
    
    
    figure(fig_nr1)
    subplot(3,1,1); plot(t,radtodeg(phi_L),t,radtodeg(phi_R))
    title('Comparison wing kinematics left and right, stroke plane @ 55deg')
    xlabel('t [s]')
    ylabel('phi [deg]')
    legend('left','right')
    subplot(3,1,2); plot(t,radtodeg(theta_L),t,radtodeg(theta_R))
    xlabel('t [s]')
    ylabel('theta [deg]')
    legend('left','right')
    subplot(3,1,3); plot(t,radtodeg(eta_L),t,radtodeg(eta_R),t,90)
    xlabel('t [s]')
    ylabel('eta [deg]')
    legend('left','right')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/stroke_pl'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/pte_stroke_lr/stroke_pl_' int2str(seq_nr)], 'fig')
    
    end
    
end


%     %Left wingtip parameters
% 
%     Lwingtip = zeros(length(qL1),3);
%     
%     Lwt = [0; -1; 0];
%     
%     Lwingtip2 = zeros(length(qL1),3);
%     
%     Lwt2 = [0.01; -sqrt(1-0.01^2); 0];
%     
%     N = length(qL1);
%     
%     % Calculate wing_tip path
%     
%     for j = 1:N
%         
%         DCM_L_r = quat2matNEW([qL1(j) qL2(j) qL3(j) qL4(j)]);
% 
%         Lwingtip(j,:) = DCM_L_r*Lwt;
%         
%         Lwingtip2(j,:) = DCM_L_r*Lwt2;
%         
%     end
%     
%     
%     
%     
%     % Right wingtip parameters
%     
%     Rwingtip = zeros(length(qR1),3);
%     
%     Rwt = [0; 1; 0];
%     
%     Rwingtip2 = zeros(length(qR1),3);
%     
%     Rwt2 = [0.01; sqrt(1-0.01^2); 0];
%     
%     
%     % Calculate wing_tip path
%     
%     for j = 1:N
%         
%         DCM_R_r = quat2matNEW([qR1(j) qR2(j) qR3(j) qR4(j)]);
% 
%         Rwingtip(j,:) = DCM_R_r*Rwt;
%         
%         Rwingtip2(j,:) = DCM_R_r*Rwt2;
%         
%     end
%     
%     
%     
%     
% %     %Determine beginnings and ends of each stroke for the left wing.
% %     
% %     [GammaL locL] = findpeaks(Lwingtip(:,2), 'minpeakdistance' , 27);
% %       
% %     
% %     
% %     
% %     %Determine beginnings and ends of each stroke for the right wing.
% %     
% %     [GammaR locR] = findpeaks(Rwingtip(:,2), 'minpeakdistance' , 27);    
% %     
% %     
% %     % Now calculate the average wing-stroke plane for the left wing per wingbeat:
% %     
% %     Beta_avg_L = zeros(length(locL)-1,1);
% %     
% %     
% %     for k = 2:length(locL)
% %         
% %         temp_x = zeros(locL(k)-locL(k-1)-1,1);
% %         temp_z = zeros(locL(k)-locL(k-1)-1,1);
% %         
% %         j = 1;
% %         
% %         for i = (locL(k-1)):(locL(k))
% %             
% %             temp_x(j) = Lwingtip(i,1);
% %             
% %             temp_z(j) = Lwingtip(i,3);
% %             
% %             j = j+1;
% %         end
% %         
% % %         plotregression(temp_z,temp_x)
% % %         pause
% %         
% %         [r,m,b] = regression(temp_z,temp_x,'one');
% %         
% %         Beta_avg_L(k-1) =  (pi/2)-real(atan2(m,1));
% %         
% %         clear r m b temp_x temp_z
% %     end
% %     
% %     
% %     
% %     % Now calculate the average wing-stroke plane for the right wing per wingbeat:
% %     
% %     Beta_avg_R = zeros(length(locR)-1,1);
% %     
% %     
% %     for k = 2:length(locR)
% %         
% %         temp_x = zeros(locR(k)-locR(k-1)-1,1);
% %         temp_z = zeros(locR(k)-locR(k-1)-1,1);
% %         
% %         j = 1;
% %         
% %         for i = (locR(k-1)):(locR(k))
% %             
% %             temp_x(j) = Rwingtip(i,1);
% %             
% %             temp_z(j) = Rwingtip(i,3);
% %             
% %             j = j+1;
% %         end
% %         
% %         %plotregression(temp_z,temp_x)
% %         %pause
% %         
% %         [r,m,b] = regression(temp_z,temp_x,'one');
% %         
% %         Beta_avg_R(k-1) =  (pi/2)-real(atan2(m,1));
% %         
% %         clear r m b temp_x temp_z
% %     end
%     
% %     figure()
% %     plot(t(locL(1:(length(locL)-1))),radtodeg(Beta_avg_L),t(locR(1:(length(locR)-1))),radtodeg(Beta_avg_R))
%     
%     % Stroke plane angle is chosen to be constant at 55 degrees.
%     
%     beta = -(55/180)*pi;
%     
%     Rot_mat = [cos(beta) 0 -sin(beta); ...
%                    0 1 0; ...
%                    sin(beta) 0 cos(beta)];
%     
%     WTL_stroke = zeros(N,3);
%     
%     WTR_stroke = zeros(N,3);
%     
%     WTL_stroke2 = zeros(N,3);
%     
%     WTR_stroke2 = zeros(N,3);
%     
%     for j = 1:N
%         
%         % Rotate the wingtip coordinates from the body reference frame to
%         % the strokesplane coordinates.
%                
%         WTL_stroke(j,:) = Rot_mat*Lwingtip(j,:)';
%         
%         WTR_stroke(j,:) = Rot_mat*Rwingtip(j,:)';
%         
%         WTL_stroke2(j,:) = Rot_mat*Lwingtip2(j,:)';
%         
%         WTR_stroke2(j,:) = Rot_mat*Rwingtip2(j,:)';
%         
%     end
%     
% %     phi_L = real(atan2(-WTL_stroke(:,1),-WTL_stroke(:,2)));
% %     theta_L = real(atan2(-WTL_stroke(:,3),sqrt(WTL_stroke(:,1).^2+WTL_stroke(:,2).^2)));
% %     %eta_L = -real(atan2(2*(qL4.*qL1+qL2.*qL3),1-2*(qL1.^2+qL2.^2)));
% %     %eta_L = ((55/180)*pi)+real(asin(2*(qL4.*qL2-qL3.*qL1)));
% %     eta_L = -((55/180)*pi)+real(atan2(2*(qL4.*qL1+qL2.*qL3),1-2*(qL1.^2+qL2.^2)));
% %    
% %     
% %     phi_R = real(atan2(-WTR_stroke(:,1),WTR_stroke(:,2)));
% %     theta_R = real(atan2(-WTR_stroke(:,3),sqrt(WTR_stroke(:,1).^2+WTR_stroke(:,2).^2)));
% %     %eta_R = real(atan2(2*(qR4.*qR1+qR2.*qR3),1-2*(qR1.^2+qR2.^2)));
% %     %eta_R = ((55/180)*pi)+real(asin(2*(qR4.*qR2-qR3.*qR1)));
% %     eta_R = -((55/180)*pi)+real(atan2(2*(qR4.*qR1+qR2.*qR3),1-2*(qR1.^2+qR2.^2)));
% 
% 
% 
%     phi_L = -real(atan2(WTL_stroke(:,1),-WTL_stroke(:,2)));
%     theta_L = -real(asin(WTL_stroke(:,3)./sqrt(WTL_stroke(:,1).^2+WTL_stroke(:,2).^2+WTL_stroke(:,3).^2)));
%     
%     phi_L2 = -real(atan2(WTL_stroke2(:,1),-WTL_stroke2(:,2)));   
%     theta_L2 = -real(asin(WTL_stroke2(:,3)./sqrt(WTL_stroke2(:,1).^2+WTL_stroke2(:,2).^2+WTL_stroke2(:,3).^2)));
%     
%     d_phi_L = phi_L2-phi_L;
%     d_theta_L = theta_L2-theta_L;
% 
%     
%     eta_L = real(atan2(d_theta_L(:),d_phi_L(:)));
%     %eta_L = real(acos((WTL_stroke2(:,1)-WTL_stroke(:,1))/sqrt((WTL_stroke2(:,1)-WTL_stroke(:,1)).^2+(WTL_stroke2(:,2)-WTL_stroke(:,2)).^2+(WTL_stroke2(:,3)-WTL_stroke(:,3)).^2)));
%     
%     phi_R = -real(atan2(WTR_stroke(:,1),WTR_stroke(:,2)));
%     theta_R = -real(asin(WTR_stroke(:,3)./sqrt(WTR_stroke(:,1).^2+WTR_stroke(:,2).^2+WTR_stroke(:,3).^2)));
%     
%     phi_R2 = -real(atan2(WTR_stroke2(:,1),WTR_stroke2(:,2)));
%     theta_R2 = -real(asin(WTR_stroke2(:,3)./sqrt(WTR_stroke2(:,1).^2+WTR_stroke2(:,2).^2+WTR_stroke2(:,3).^2)));
%     
%     d_phi_R = phi_R2-phi_R;
%     d_theta_R = theta_R2-theta_R;
%     
%     eta_R = real(atan2(d_theta_R(:),d_phi_R(:)));
%     %eta_R = real(acos((WTR_stroke2(:,1)-WTR_stroke(:,1))/sqrt((WTR_stroke2(:,1)-WTR_stroke(:,1)).^2+(WTR_stroke2(:,2)-WTR_stroke(:,2)).^2+(WTR_stroke2(:,3)-WTR_stroke(:,3)).^2)));
    
    