n = 100;
t = [0:1/(n-1):1]';
stroke_wb_steady = fnval(stroke_steady_func,t);
pitch_wb_steady = fnval(pitch_steady_func,t);
dev_wb_steady = fnval(dev_steady_func,t);

limit = .5;
roll_limit = limit*nanstd(roll_dot_dot_mean_wb(:));
pitch_limit = limit*nanstd(pitch_dot_dot_mean_wb(:));
yaw_limit = limit*nanstd(yaw_dot_dot_mean_wb(:));
F_limit = limit*nanstd(F_mean_wb(:));

roll_dot_dot_max = max(roll_dot_dot_mean_wb);
pitch_dot_dot_max = max(pitch_dot_dot_mean_wb);
yaw_dot_dot_max = max(yaw_dot_dot_mean_wb);
F_max = max(F_mean_wb);

        t_norm = [];
        t_norm_LR = [];

        Dstroke_rollaccel_L = [];
        Dstroke_yawaccel_L = [];
        Dstroke_rollaccel_R = [];
        Dstroke_yawaccel_R = [];
        
        Dstroke_pitchaccel = [];
        Dstroke_F = [];

        Dpitch_rollaccel_L = [];
        Dpitch_yawaccel_L = [];
        Dpitch_rollaccel_R = [];
        Dpitch_yawaccel_R = [];
        
        Dpitch_pitchaccel = [];
        Dpitch_F = [];


        Ddev_rollaccel_L = [];
        Ddev_yawaccel_L = [];
        Ddev_rollaccel_R = [];
        Ddev_yawaccel_R = [];
        
        Ddev_pitchaccel = [];
        Ddev_F = [];

for seq = 1:size(n_down_start_L,2)
    counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        roll_dot_dot_mean_wb_now = roll_dot_dot_mean_wb(wb,seq);
        pitch_dot_dot_mean_wb_now = pitch_dot_dot_mean_wb(wb,seq);
        yaw_dot_dot_mean_wb_now = yaw_dot_dot_mean_wb(wb,seq);
        F_mean_wb_now = F_mean_wb(wb,seq);

        n_down_start_L_now = n_down_start_L(wb,seq);
        n_up_stop_L_now = n_down_start_L(wb+1,seq);
        n_down_start_R_now = n_down_start_R(wb,seq);
        n_up_stop_R_now = n_down_start_R(wb+1,seq);
        
        if isnan(n_down_start_L_now)==0 && isnan(n_up_stop_L_now)==0

        stroke_wb_L_now = stroke_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
        stroke_wb_R_now = stroke_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
        pitch_wb_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq));
        pitch_wb_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq));
        dev_wb_L_now = dev_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
        dev_wb_R_now = dev_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);

        t_L = [1:(length(stroke_wb_L_now)-1)/(n-1):length(stroke_wb_L_now)]';
        t_R = [1:(length(stroke_wb_R_now)-1)/(n-1):length(stroke_wb_R_now)]';

        stroke_wb_L_interp = rad2deg(interp1(stroke_wb_L_now,t_L));
        stroke_wb_R_interp = rad2deg(interp1(stroke_wb_R_now,t_R));
        pitch_wb_L_interp = rad2deg(interp1(pitch_wb_L_now,t_L));
        pitch_wb_R_interp = rad2deg(interp1(pitch_wb_R_now,t_R));
        dev_wb_L_interp = rad2deg(interp1(dev_wb_L_now,t_L));
        dev_wb_R_interp = rad2deg(interp1(dev_wb_R_now,t_R));

        Dstroke_wb_L = stroke_wb_L_interp - stroke_wb_steady;
        Dpitch_wb_L = pitch_wb_L_interp - pitch_wb_steady;
        Ddev_wb_L = dev_wb_L_interp - dev_wb_steady;
        Dstroke_wb_R = stroke_wb_R_interp - stroke_wb_steady;
        Dpitch_wb_R = pitch_wb_R_interp - pitch_wb_steady;
        Ddev_wb_R = dev_wb_R_interp - dev_wb_steady;
        
        
        
        % store data
        t_norm(end+1:end+length(t),1) = t;
        t_norm_LR(end+1:end+length(t),1) = t;
        t_norm_LR(end+1:end+length(t),1) = t;

        
        Dstroke_rollaccel_L(end+1:end+length(t),1) = Dstroke_wb_L / roll_dot_dot_mean_wb_now;
        Dstroke_rollaccel_R(end+1:end+length(t),1) = Dstroke_wb_R / roll_dot_dot_mean_wb_now;
        Ddev_rollaccel_L(end+1:end+length(t),1) = Ddev_wb_L / roll_dot_dot_mean_wb_now;
        Ddev_rollaccel_R(end+1:end+length(t),1) = Ddev_wb_R / roll_dot_dot_mean_wb_now;
        Dpitch_rollaccel_L(end+1:end+length(t),1) = Dpitch_wb_L / roll_dot_dot_mean_wb_now;
        Dpitch_rollaccel_R(end+1:end+length(t),1) = Dpitch_wb_R / roll_dot_dot_mean_wb_now;
        
        Dstroke_yawaccel_L(end+1:end+length(t),1) = Dstroke_wb_L / yaw_dot_dot_mean_wb_now;
        Dstroke_yawaccel_R(end+1:end+length(t),1) = Dstroke_wb_R / yaw_dot_dot_mean_wb_now;
        Ddev_yawaccel_L(end+1:end+length(t),1) = Ddev_wb_L / yaw_dot_dot_mean_wb_now;
        Ddev_yawaccel_R(end+1:end+length(t),1) = Ddev_wb_R / yaw_dot_dot_mean_wb_now;
        Dpitch_yawaccel_L(end+1:end+length(t),1) = Dpitch_wb_L / yaw_dot_dot_mean_wb_now;
        Dpitch_yawaccel_R(end+1:end+length(t),1) = Dpitch_wb_R / yaw_dot_dot_mean_wb_now;
        
        Dstroke_pitchaccel(end+1:end+length(t),1) = Dstroke_wb_L / pitch_dot_dot_mean_wb_now;
        Dstroke_pitchaccel(end+1:end+length(t),1) = Dstroke_wb_R / pitch_dot_dot_mean_wb_now;
        Ddev_pitchaccel(end+1:end+length(t),1) = Ddev_wb_L / pitch_dot_dot_mean_wb_now;
        Ddev_pitchaccel(end+1:end+length(t),1) = Ddev_wb_R / pitch_dot_dot_mean_wb_now;
        Dpitch_pitchaccel(end+1:end+length(t),1) = Dpitch_wb_L / pitch_dot_dot_mean_wb_now;
        Dpitch_pitchaccel(end+1:end+length(t),1) = Dpitch_wb_R / pitch_dot_dot_mean_wb_now;
        
        Dstroke_F(end+1:end+length(t),1) = Dstroke_wb_L / F_mean_wb_now;
        Dstroke_F(end+1:end+length(t),1) = Dstroke_wb_R / F_mean_wb_now;
        Ddev_F(end+1:end+length(t),1) = Ddev_wb_L / F_mean_wb_now;
        Ddev_F(end+1:end+length(t),1) = Ddev_wb_R / F_mean_wb_now;
        Dpitch_F(end+1:end+length(t),1) = Dpitch_wb_L / F_mean_wb_now;
        Dpitch_F(end+1:end+length(t),1) = Dpitch_wb_R / F_mean_wb_now;
        
        
        
%         
%         
%         
%             subplot(3,5,1)
%             plot(t,stroke_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,6)
%             plot(t,stroke_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,11)
%             plot(t,stroke_wb_L_interp-stroke_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
% 
%             subplot(3,5,2)
%             plot(t,pitch_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,7)
%             plot(t,pitch_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,12)
%             plot(t,pitch_wb_L_interp-pitch_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
% 
%             subplot(3,5,3)
%             plot(t,dev_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,8)
%             plot(t,dev_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,13)
%             plot(t,dev_wb_L_interp-dev_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
% 
%              subplot(3,5,4)
%             plot(t,aoa_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,9)
%             plot(t,aoa_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,14)
%             plot(t,aoa_wb_L_interp-aoa_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
% 
%             subplot(3,5,5)
%             plot(t,U_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,10)
%             plot(t,U_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             subplot(3,5,15)
%             plot(t,U_wb_L_interp-U_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             
% %             pause(.001)
        end
    end
end

save('WBmod_data.mat','stroke_wb_steady','pitch_wb_steady','dev_wb_steady',...
't_norm','Dstroke_rollaccel_L','Dstroke_yawaccel_L','Dstroke_rollaccel_R','Dstroke_yawaccel_R','Dstroke_pitchaccel','Dstroke_F',...
    'Dpitch_rollaccel_L','Dpitch_yawaccel_L','Dpitch_rollaccel_R','Dpitch_yawaccel_R',...
    'Dpitch_pitchaccel','Dpitch_F','Ddev_rollaccel_L','Ddev_yawaccel_L','Ddev_rollaccel_R','Ddev_yawaccel_R','Ddev_pitchaccel','Ddev_F');


% 
% calc_WBfunc_csaps    
% plot_WBfunc_csaps_LRdiff
%     
% subplot(3,5,1)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
% %     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,6)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
% %     xlabel('time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
%     grid on
% subplot(3,5,11)
% ats([0 1 -45 45])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'YTick',-90:45:90,'fontsize',8)
%     xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
%     grid on
%             
% 
% subplot(3,5,2)
% ats([0 1 0 180])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:180,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,7)
% ats([0 1 0 180])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:180,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,12)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('normalized time','fontsize',10)
% %     ylabel('Left - Right','fontsize',10)
%     grid on
%             
% subplot(3,5,3)
% ats([0 1 -30 30])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:30:180,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,8)
% ats([0 1 -30 30])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:30:180,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,13)
% ats([0 1 -45 45])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'YTick',-90:45:90,'fontsize',8)
%     xlabel('normalized time','fontsize',10)
% %     ylabel('Left - Right','fontsize',10)
%     grid on
%             
%     
% 
% subplot(3,5,4)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,9)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Right wing','fontsize',10)
%     grid on
% subplot(3,5,14)
% ats([0 1 -90 90])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('normalized time','fontsize',10)
% %     ylabel('Left - Right','fontsize',10)
%     grid on
%             
% subplot(3,5,5)
% ats([0 1 0 6])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',0:3:6,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,10)
% ats([0 1 0 6])
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
%     set(gca,'YTick',0:3:6,'fontsize',8)
% %     xlabel('time','fontsize',10)
% %     ylabel('Left wing','fontsize',10)
%     grid on
% subplot(3,5,15)
% ats([0 1 -3 3])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'YTick',-3:3:3,'fontsize',8)
%     xlabel('normalized time','fontsize',10)
% %     ylabel('Left - Right','fontsize',10)
%     grid on
%             
% mid
% hi
% 
% 
