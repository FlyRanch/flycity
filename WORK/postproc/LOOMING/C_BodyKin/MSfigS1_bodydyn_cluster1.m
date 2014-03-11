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
    
    % cluster 1
    if abs(stim_angle_vel_mirror_pre(i))<subset_cut(1) 
        if isnan(color_var(i)) == 0 && color_var(i)~=0
            subplot(3,1,1)
            plot(t-t_shift(i),roll_dot_dot_plot(:,i)/freq_steady^2,'-','color',grey_color,'linewidth',.25)
            subplot(3,1,2)
            plot(t-t_shift(i),pitch_dot_dot_plot(:,i)/freq_steady^2,'-','color',grey_color,'linewidth',.25)
            subplot(3,1,3)
            plot(t-t_shift(i),yaw_dot_dot_plot(:,i)/freq_steady^2,'-','color',grey_color,'linewidth',.25)
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
angle_pre = t_c1(:);
angle_post = roll_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_WBdp

% dV
subplot(3,1,2)
angle_pre = t_c1(:);
angle_post = pitch_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_WBdp

% Fhor
subplot(3,1,3)
angle_pre = t_c1(:);
angle_post = yaw_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_WBdp

%%
subplot(3,1,1)
axis([t_start t_stop -5 5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:2.5:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll [deg]','fontsize',10)
%    grid on 
subplot(3,1,2)
axis([t_start t_stop -5 5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-15:2.5:15,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch [deg]','fontsize',10)
%    grid on 
subplot(3,1,3)
axis([t_start t_stop -5 5])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-15:2.5:15,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('yaw [deg]','fontsize',10)
%    grid on 

