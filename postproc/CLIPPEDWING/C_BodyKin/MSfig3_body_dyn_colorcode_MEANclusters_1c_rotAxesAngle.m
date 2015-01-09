%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,3,1)
% title('strokeplane angles','fontsize',10)
hold on
subplot(3,3,2)
% title('body euler angles','fontsize',10)
hold on
subplot(3,3,3)
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

%% MEAN clusters
for j=1:length(color_var)
    i=j;
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t-t_shift(i),drot_R_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,4)
        plot(t-t_shift(i),drot_L_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,7)
        plot(t-t_shift(i),rot_axis_pos_plot(:,i),'-','color',grey_color,'linewidth',.25)
        
        subplot(3,3,2)
        plot(t-t_shift(i),rot_dot_R_plot(:,i)/freq_steady,'-','color',grey_color,'linewidth',.25)
        subplot(3,3,5)
        plot(t-t_shift(i),rot_dot_L_plot(:,i)/freq_steady,'-','color',grey_color,'linewidth',.25)
        subplot(3,3,8)
        plot(t-t_shift(i),rot_axis_vel_plot(:,i),'-','color',grey_color,'linewidth',.25)
        
        subplot(3,3,3)
        plot(t-t_shift(i),rot_dot_dot_R_plot(:,i)/freq_steady^2,'-','color',grey_color,'linewidth',.25)
        subplot(3,3,6)
        plot(t-t_shift(i),rot_dot_dot_L_plot(:,i)/freq_steady^2,'-','color',grey_color,'linewidth',.25)
        subplot(3,3,9)
        plot(t-t_shift(i),rot_axis_accel_plot(:,i),'-','color',grey_color,'linewidth',.25)
        i
    end
end

% plot_meanclusters_bodyangles_4c
% plot_meanclusters_bodyangles_2c_rotAxesAngle
plot_meanclusters_bodyangles_1c_rotAxesAngle

% for j=1:length(subset_seqs)
%     i=subset_seqs(j);
%     size(stim_angle_vel,2) - i;
%     
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         subplot(3,3,1)
%         plot(t-t_shift(i),roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,4)
%         plot(t-t_shift(i),pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,7)
%         plot(t-t_shift(i),yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         
%         subplot(3,3,2)
%         plot(t-t_shift(i),roll_dot_plot(:,i)/freq_steady,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,5)
%         plot(t-t_shift(i),pitch_dot_plot(:,i)/freq_steady,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,8)
%         plot(t-t_shift(i),yaw_dot_plot(:,i)/freq_steady,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         
%         subplot(3,3,3)
%         plot(t-t_shift(i),roll_dot_dot_plot(:,i)/freq_steady^2,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,6)
%         plot(t-t_shift(i),pitch_dot_dot_plot(:,i)/freq_steady^2,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,9)
%         plot(t-t_shift(i),yaw_dot_dot_plot(:,i)/freq_steady^2,'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         
%     end
% end



subplot(3,3,1)
axis([t_start t_stop -30 60])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis rotation R [deg]','fontsize',10)
%    grid on 
subplot(3,3,4)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis rotation L [deg]','fontsize',10)
%    grid on 
subplot(3,3,7)
axis([t_start t_stop -90 270])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-360:90:360,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('rotation axis angle [deg]','fontsize',10)
%    grid on 


subplot(3,3,2)
axis([t_start t_stop -7.5 7.5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:2.5:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis R vel/f [-]','fontsize',10)
%    grid on 
subplot(3,3,5)
axis([t_start t_stop -7.5 7.5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:2.5:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis L vel/f [-]','fontsize',10)
%    grid on 
subplot(3,3,8)
axis([t_start t_stop -90 270])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-360:90:360,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('velocity axis [deg]','fontsize',10)
%    grid on 

subplot(3,3,3)
axis([t_start t_stop -2 2])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:1:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis R accel/f2 [-]','fontsize',10)
%    grid on 
subplot(3,3,6)
axis([t_start t_stop -2 2])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:1:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('axis L accel/f2 [-]','fontsize',10)
%    grid on 
subplot(3,3,9)
axis([t_start t_stop -90 270])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-360:90:360,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('acceleration axis [deg]','fontsize',10)
%    grid on 



