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

for j=1:length(subset_seqs)
    i=subset_seqs(j);
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t-t_shift(i),Mroll(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,4)
        plot(t-t_shift(i),Mpitch(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,7)
        plot(t-t_shift(i),Myaw(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([t_start t_stop -.005 .005])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-.005:.001:.005,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('Mroll [Nm]','fontsize',10) 
%    grid on 
subplot(3,3,4)
axis([t_start t_stop -.005 .005])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-.005:.001:.005,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('Mpitch [Nm]','fontsize',10) 
%    grid on 
subplot(3,3,7)
axis([t_start t_stop -.005 .005])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-.005:.001:.005,'fontsize',8)
    xlabel('time','fontsize',10)
ylabel('Myaw [Nm]','fontsize',10) 
%    grid on 

