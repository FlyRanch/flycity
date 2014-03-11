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
        plot(t-t_shift(i),F_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,6)
        plot(t-t_shift(i),Fn_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,9)
        plot(t-t_shift(i),Ft_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([-.025 .05 -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('strokeplane normal','fontsize',10) 
    grid on
subplot(3,3,4)
axis([-.025 .05 -30 30])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Fsp pitch','fontsize',10)
    grid on
subplot(3,3,7)
axis([-.025 .05 -30 30])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:30:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('Fsp roll','fontsize',10)
    grid on

subplot(3,3,2)
axis([-.025 .05 -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll','fontsize',10)
    grid on
subplot(3,3,5)
axis([-.025 .05 -90 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('yaw','fontsize',10)
    grid on
subplot(3,3,8)
axis([-.025 .05 -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on


subplot(3,3,3)
axis([-.025 .05 0 2])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[0:1:2],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F','fontsize',10)
    grid on
subplot(3,3,6)
axis([-.025 .05 -2 2])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-2:1:2],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Fn','fontsize',10)
    grid on
subplot(3,3,9)
axis([-.025 .05 -1 2])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-2:1:2],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('Ft','fontsize',10)
    grid on

