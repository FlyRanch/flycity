% plot timelines for 2 FLIES
close all

fly1=1; % FLY1 head
fly2=3; % FLY2 head

subplot(5,1,1)
plot(pathDB.t,pathDB.teta(:,fly1),'-k')
hold on
% plot(pathDB.t,pathDB.teta(:,fly1),'-k')
axis([0 .13 0 20])
legend('Fly 1', 'Fly 2')
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',0:10:20)
xlabel('time')
ylabel('optical angle')
grid on
    
subplot(5,1,2)
% hold off
% plot(-1,0,'^k','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',10)
% hold on
% plot(-1,0,'^k','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',10)
% legend('velocity','yaw','Location','NorthWest')
plot(pathDB.t,pathDB.stim_angle_vel(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.stim_angle_vel(:,fly2),'-k')
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly1),'-k')
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly2),'-k')
axis([0 .13 -180 180])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',-180:90:180)
xlabel('time')
ylabel('direction')
grid on

    
subplot(5,1,3)
plot(pathDB.t,pathDB.V(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.V(:,fly2),'-k')
axis([0 .13 0 1])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',0:.25:1)
xlabel('time')
ylabel('V')
grid on
    
subplot(5,1,4)
plot(pathDB.t,pathDB.An_hor(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.An_hor(:,fly2),'-k')
axis([0 .13 -10 10])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',-10:5:10)
xlabel('time')
ylabel('An')
grid on
    
subplot(5,1,5)
plot(pathDB.t,pathDB.At_hor(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.At_hor(:,fly2),'-k')
axis([0 .13 -15 15])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',-15:5:15)
xlabel('time')
ylabel('At')
grid on

    for j = 1:length(pathDB.IDX)
        length(pathDB.IDX)-j
        if isnan(pathDB.IDX(j,fly1))==0
            subplot(5,1,1)
            plot(pathDB.t(j),pathDB.teta(j),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.teta(j),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            grid on

            subplot(5,1,2)
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly1),'+k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly2),'xk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            grid on

            subplot(5,1,3)
            plot(pathDB.t(j),pathDB.V(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.V(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            grid on

            subplot(5,1,4)
            plot(pathDB.t(j),pathDB.An_hor(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.An_hor(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            grid on

            subplot(5,1,5)
            plot(pathDB.t(j),pathDB.At_hor(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.At_hor(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
            grid on
        end
    end

    
    saveas(gcf,['collision_reaponse_2flies.fig'])
    saveas(gcf,['collision_reaponse_2flies.png'])
    plot2svg

    
