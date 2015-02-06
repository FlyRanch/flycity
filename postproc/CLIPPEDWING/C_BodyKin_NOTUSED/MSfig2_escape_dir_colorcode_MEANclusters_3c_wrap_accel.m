%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,3,1)
hold on
subplot(3,3,2)
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

Az_norm = pathDB.Fz-1;

for j=1:length(color_var)
    i=j;
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t-t_shift(i),stim_angle_vel_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,4)
        plot(t-t_shift(i),dV_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,7)
        plot(t-t_shift(i),stim_angle_accel_plot_wrap(:,i),'-','color',grey_color,'linewidth',.25)
        
%         subplot(3,3,2)
%         plot(t-t_shift(i),F_plot(:,i),'-','color',grey_color,'linewidth',.25)
%         subplot(3,3,5)
%         plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',grey_color,'linewidth',.25)
%         subplot(3,3,8)
%         plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',grey_color,'linewidth',.25)
        
        subplot(3,3,3)
        plot(t-t_shift(i),A_norm_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,6)
        plot(t-t_shift(i),Ahor_norm_plot(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,9)
        plot(t-t_shift(i),Az_norm_plot(:,i),'-','color',grey_color,'linewidth',.25)

    end
end

% 
% for j=1:length(subset_seqs)
%     i=subset_seqs(j);
%     size(stim_angle_vel,2) - i;
%     
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         subplot(3,3,1)
%         plot(t-t_shift(i),stim_angle_vel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,4)
%         plot(t-t_shift(i),dV_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,7)
%         plot(t-t_shift(i),stim_angle_accel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         
%         subplot(3,3,2)
%         plot(t-t_shift(i),F_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,5)
%         plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,8)
%         plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%     end
% end


%% plot mean clusters
% plot_meanclusters_bodydyn_5cn2c
% plot_meanclusters_bodydyn_4cn2c
% plot_meanclusters_bodydyn_4cn2c_quartiles
% plot_meanclusters_bodydyn_2c_quartiles_FdirWrap_Ahor
% plot_meanclusters_bodydyn_2c_quartiles_FdirWrap_Accel
plot_meanclusters_bodydyn_3c_quartiles_FdirWrap_Accel
% plot_meanclusters_bodydyn_3cn2c

subplot(3,3,1)
axis([t_start t_stop -90 270])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',-360:90:360,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('heading [deg]','fontsize',10) 
%    grid on 
subplot(3,3,4)
axis([t_start t_stop -.15 .3])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',-.15:.15:1,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('dV [m/s]','fontsize',10) 
%    grid on 
subplot(3,3,7)
axis([t_start t_stop -90 270])
    set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'YTick',-360:90:360,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dir F_h_o_r [deg]','fontsize',10)
%    grid on 

subplot(3,3,2)
axis([t_start t_stop 0.75 1.25])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',[0:.25:2],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F/Mg [-]','fontsize',10)
%    grid on 
subplot(3,3,5)
axis([t_start t_stop -10 10])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',[-30:5:30],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F/Mg_r_o_l_l [deg]','fontsize',10)
%    grid on 
subplot(3,3,8)
axis([t_start t_stop -10 10])
    set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'YTick',[-30:5:30],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('F/Mg_p_i_t_c_h [deg]','fontsize',10)
%    grid on 
subplot(3,3,3)
axis([t_start t_stop 0 1])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',[0:.5:10],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('a/g [-]','fontsize',10)
%    grid on 
subplot(3,3,6)
axis([t_start t_stop 0 1])
%     set(gca,'XTick',-0.06:0.03:.06) 
    set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',[0:.5:10],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('a/g horiz [-]','fontsize',10)
%    grid on 
subplot(3,3,9)
axis([t_start t_stop -.5 .5])
    set(gca,'XTick',-0.06:0.03:.06) 
%     set(gca,'XTick',-0.06:0.03:.06,'XTickLabel',[]) 
    set(gca,'YTick',[-1:.5:1],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('a/g vert [-]','fontsize',10)
%    grid on 

