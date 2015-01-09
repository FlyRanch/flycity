%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,3,1)
hold on
subplot(3,3,3)
hold on
subplot(3,3,5)
hold on
subplot(3,3,2)
hold on
subplot(3,3,4)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t-t_shift(i),stim_angle_spn_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,4)
        plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,7)
        plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)

        subplot(3,3,2)
        plot(t-t_shift(i),roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,5)
        plot(t-t_shift(i),yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,8)
        plot(t-t_shift(i),pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,3,3)
        plot(t-t_shift(i),roll_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,6)
        plot(t-t_shift(i),yaw_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,9)
        plot(t-t_shift(i),pitch_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([t_start t_end -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('strokeplane normal','fontsize',10) 
    grid on
subplot(3,3,4)
axis([t_start t_end -30 30])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Fsp pitch','fontsize',10)
    grid on
subplot(3,3,7)
axis([t_start t_end -30 30])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('Fsp roll','fontsize',10)
    grid on

subplot(3,3,2)
axis([t_start t_end -90 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll','fontsize',10)
    grid on
subplot(3,3,5)
axis([t_start t_end -90 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('yaw','fontsize',10)
    grid on
subplot(3,3,8)
axis([t_start t_end -45 135])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-45;0;90;135],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on


subplot(3,3,3)
axis([t_start t_end -4000 4000])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll rate','fontsize',10)
    grid on
subplot(3,3,6)
axis([t_start t_end -2000 2000])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-2000:2000:2000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('yaw rate','fontsize',10)
    grid on
subplot(3,3,9)
axis([t_start t_end -1000 3000])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-1000;0;3000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch rate','fontsize',10)
    grid on

