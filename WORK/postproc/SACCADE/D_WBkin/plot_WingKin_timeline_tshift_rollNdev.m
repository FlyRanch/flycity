%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,2,1)
hold on
subplot(3,2,3)
hold on
subplot(3,2,5)
hold on
subplot(3,2,2)
hold on
subplot(3,2,4)
hold on
subplot(3,2,6)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,2,1)
        plot(t-t_shift(i),F_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,3)
        plot(t-t_shift(i),F_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),F_ver_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,2,2)
        plot(t-t_shift(i),stim_angle_spn_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([t_start t_stop 0 2.5])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-1:1:3,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('F/Mg','fontsize',10) 
    grid on
subplot(3,2,3)
axis([t_start t_stop 0 2.5])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-1:1:3,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('horizontal F/Mg','fontsize',10)
    grid on
subplot(3,2,5)
axis([t_start t_stop -.5 2])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-1:1:3,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('vertical F/Mg','fontsize',10)
    grid on

subplot(3,2,2)
axis([t_start t_stop -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('SP normal','fontsize',10) 
    grid on
subplot(3,2,4)
axis([t_start t_stop -30 30])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F-SP roll','fontsize',10)
    grid on
subplot(3,2,6)
axis([t_start t_stop -30 30])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('F-SP pitch','fontsize',10)
    grid on
    
