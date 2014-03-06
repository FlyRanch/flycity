% plot timelines for 2 FLIES
close all

fly1=1; % FLY1 head
fly2=3; % FLY2 head

figure(1)
subplot(3,1,1)
plot(-1,0,'^k','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',10)
hold on
plot(-1,0,'vk','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',10)
legend('Fly 1', 'Fly 2','Location','Northwest')
plot(pathDB.t,pathDB.teta(:,fly1),'-k')
% plot(pathDB.t,pathDB.teta(:,fly1),'-k')
hold on
axis([0 .13 0 20])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',0:10:20)
xlabel('time')
ylabel('optical angle')
grid on
    
subplot(3,1,2)
plot(pathDB.t,pathDB.stim_angle_vel(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.stim_angle_vel(:,fly2),'-k')
axis([0 .13 -180 180])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',-180:180:180)
xlabel('time')
ylabel('heading angle')
grid on

subplot(3,1,3)
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.stim_angle_yaw(:,fly2),'-k')
axis([0 .13 -180 180])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',-180:180:180)
xlabel('time')
ylabel('yaw angle')
grid on

figure(2)    
subplot(3,1,1)
plot(pathDB.t,pathDB.V(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.V(:,fly2),'-k')
axis([0 .13 0 .6])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',0:.3:1)
xlabel('time')
ylabel('V')
grid on
    
subplot(3,1,2)
plot(pathDB.t,abs(pathDB.An_hor(:,fly1)),'-k')
hold on
plot(pathDB.t,abs(pathDB.An_hor(:,fly2)),'-k')
axis([0 .13 0 10])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',0:5:10)
xlabel('time')
ylabel('An')
grid on
    
subplot(3,1,3)
plot(pathDB.t,pathDB.At_hor(:,fly1),'-k')
hold on
plot(pathDB.t,pathDB.At_hor(:,fly2),'-k')
axis([0 .13 -5 15])
set(gca,'XTick',0:.025:.125)
set(gca,'YTick',[-5;0;15])
xlabel('time')
ylabel('At')
grid on

    for j = 1:length(pathDB.IDX)
        length(pathDB.IDX)-j
        if isnan(pathDB.IDX(j,fly1))==0
            
            figure(1)
            subplot(3,1,1)
            plot(pathDB.t(j),pathDB.teta(j),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.teta(j),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)

            subplot(3,1,2)
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.stim_angle_vel(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)

            subplot(3,1,3)
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.stim_angle_yaw(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)

            figure(2)
            subplot(3,1,1)
            plot(pathDB.t(j),pathDB.V(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.V(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)

            subplot(3,1,2)
            plot(pathDB.t(j),abs(pathDB.An_hor(j,fly1)),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),abs(pathDB.An_hor(j,fly2)),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)

            subplot(3,1,3)
            plot(pathDB.t(j),pathDB.At_hor(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly1),:),'MarkerSize',10)
            hold on
            plot(pathDB.t(j),pathDB.At_hor(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(pathDB.IDX(j,fly2),:),'MarkerSize',10)
        end
    end

    figure(1)
    saveas(gcf,['collision_response_2flies_angles.fig'])
    saveas(gcf,['collision_response_2flies_angles.png'])
    plot2svg('collision_response_2flies_angles.svg')

    figure(2)
    saveas(gcf,['collision_response_2flies_vel_accel.fig'])
    saveas(gcf,['collision_response_2flies_vel_accel.png'])
    plot2svg('collision_response_2flies_vel_accel.svg')

    
