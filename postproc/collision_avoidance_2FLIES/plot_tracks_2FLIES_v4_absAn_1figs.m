% plot timelines for 2 FLIES
close all

fly1=1; % FLY1 head
fly2=3; % FLY2 head


subplot(3,2,1)
plot(-1,0,'^k','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',5)
hold on
plot(-1,0,'vk','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',5)
legend('Fly 1', 'Fly 2','Location','Northwest')
plot(pathDB.t,pathDB.teta(:,fly1),'-k')
% plot(pathDB.t,pathDB.teta(:,fly1),'-k')
hold on
axis([0 .15 0 20])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:10:20,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('optical angle','fontsize',18)
grid on
    
subplot(3,2,3)
plot(pathDB.t,pathDB.stim_angle_vel(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.stim_angle_vel(:,fly2),'-k')
axis([0 .15 -180 180])
set(gca,'XTick',0:.025:.125,'fontsize',12)
set(gca,'YTick',-180:180:180,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('heading','fontsize',18)
grid on

subplot(3,2,5)
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly2),'-k')
axis([0 .15 -180 180])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',-180:180:180,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('yaw','fontsize',18)
grid on

    
subplot(3,2,2)
plot(pathDB.t,pathDB.V(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.V(:,fly2),'-k')
axis([0 .15 0 .8])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:.4:1,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('V','fontsize',18)
grid on
    
subplot(3,2,4)
plot(pathDB.t,abs(pathDB.An_hor(:,fly1)),'-k')
hold on
plot(pathDB.t,abs(pathDB.An_hor(:,fly2)),'-k')
axis([0 .15 0 15])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:5:15,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('An','fontsize',18)
grid on
    
subplot(3,2,6)
plot(pathDB.t,pathDB.At_hor(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.At_hor(:,fly2),'-k')
axis([0 .15 -5 15])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',[-5;0;15],'fontsize',12)
xlabel('time','fontsize',18)
ylabel('At','fontsize',18)
grid on

    for j = 1:length(pathDB.IDX)
        length(pathDB.IDX)-j
        if isnan(pathDB.IDX(j,fly1))==0
            
            
            subplot(3,2,1)
            plot(pathDB.t(j),pathDB.teta(j),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),pathDB.teta(j),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,3)
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,5)
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)

            
            subplot(3,2,2)
            plot(pathDB.t(j),pathDB.V(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),pathDB.V(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,4)
            plot(pathDB.t(j),abs(pathDB.An_hor(j,fly1)),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),abs(pathDB.An_hor(j,fly2)),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,6)
            plot(pathDB.t(j),pathDB.At_hor(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(pathDB.t(j),pathDB.At_hor(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',5)
        end
    end

    
    saveas(gcf,['collision_response_2flies.fig'])
    saveas(gcf,['collision_response_2flies.png'])
    plot2svg('collision_response_2flies.svg')
