%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,1,1)
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on

for j=1:length(color_var)
    i=j;
    size(stim_angle_vel,2) - i;
    
    %% cluster 3
    if abs(stim_angle_vel_mirror_pre(i))<subset_cut(3) & abs(stim_angle_vel_mirror_pre(i))>subset_cut(2) 
        if isnan(color_var(i)) == 0 && color_var(i)~=0
            subplot(3,1,1)
            plot(t-t_shift(i),stim_angle_vel_plot(:,i),'-','color',grey_color,'linewidth',.25)
            subplot(3,1,2)
            plot(t-t_shift(i),dV_plot(:,i),'-','color',grey_color,'linewidth',.25)
            subplot(3,1,3)
            plot(t-t_shift(i),stim_angle_accel_plot(:,i),'-','color',grey_color,'linewidth',.25)

    %         subplot(3,1,2)
    %         plot(t-t_shift(i),F_plot(:,i),'-','color',grey_color,'linewidth',.25)
    %         subplot(3,1,5)
    %         plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',grey_color,'linewidth',.25)
    %         subplot(3,1,8)
    %         plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',grey_color,'linewidth',.25)
        end
    end
end

%% plot mean clusters
%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% heading
subplot(3,1,1)
angle_pre = t_c3(:);
angle_post = stim_angle_vel_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_WBdp

% dV
subplot(3,1,2)
angle_pre = t_c3(:);
angle_post = dV_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_WBdp

% Fhor
subplot(3,1,3)
angle_pre = t_c3(:);
angle_post = stim_angle_accel_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_WBdplim

%%
subplot(3,1,1)
axis([t_start t_stop -180 0])
%     set(gca,'XTick',-0.04:0.02:.04) 
    set(gca,'XTick',-0.04:0.02:.04,'XTickLabel',[]) 
    set(gca,'YTick',-180:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('heading [deg]','fontsize',10) 
%    grid on 

subplot(3,1,2)
axis([t_start t_stop -.25 .5])
%     set(gca,'XTick',-0.04:0.02:.04) 
    set(gca,'XTick',-0.04:0.02:.04,'XTickLabel',[]) 
    set(gca,'YTick',-.25:.25:1,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('dV [m/s]','fontsize',10) 
%    grid on 

subplot(3,1,3)
axis([t_start t_stop -180 180])
    set(gca,'XTick',-0.04:0.02:.04) 
    set(gca,'YTick',-180:90:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dir F_h_o_r [deg]','fontsize',10)
%    grid on 
