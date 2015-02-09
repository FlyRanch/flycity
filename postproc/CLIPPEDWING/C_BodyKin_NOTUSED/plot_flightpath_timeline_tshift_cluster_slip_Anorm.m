%% plot flightpath timeline tshift clusters

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

for i=1:size(stim_angle_vel_plot,2)
    size(stim_angle_vel_plot,2) - i
    for j=2:size(stim_angle_vel_plot,1)
        if isnan(IDX_plot(j,i)) == 0 && isnan(IDX_plot(j-1,i)) == 0
            
            subplot(3,3,1)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_vel_plot(j-1,i) stim_angle_vel_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,4)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[dV_plot(j-1,i) dV_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,7)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_accel_plot(j-1,i) stim_angle_accel_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)

            subplot(3,3,2)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[roll_plot(j-1,i) roll_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,5)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[slip_plot(j-1,i) slip_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,8)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[pitch_plot(j-1,i) pitch_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            
            subplot(3,3,3)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[A_hor_plot(j-1,i) A_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,6)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[An_hor_plot(j-1,i) An_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,9)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[At_hor_plot(j-1,i) At_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j-1,i),:),'linewidth',linewidth_timelines)
        end
    end
end

subplot(3,3,1)
axis([-.025 .05 -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('heading','fontsize',10)
    grid on
subplot(3,3,4)
axis([-.025 .05 -.25 .5])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-.25:.25:1,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('dV','fontsize',10)
    grid on
subplot(3,3,7)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('A direction','fontsize',10)
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
    ylabel('slip','fontsize',10)
    grid on
subplot(3,3,8)
axis([-.025 .05 -45 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on


subplot(3,3,3)
axis([-.025 .05 0 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[0:10:20],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('A','fontsize',10)
    grid on
subplot(3,3,6)
axis([-.025 .05 -20 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-20:10:20],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('An','fontsize',10)
    grid on
subplot(3,3,9)
axis([-.025 .05 -10 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-10;0;10;20],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('At','fontsize',10)
    grid on

