% plot timelines for 2 FLIES
clear
clc
close all

load('flightpathDB_2flies_9clusters_Atmin_2.5n-3n2.5_response.mat')

fly1=1; % FLY1 head
fly2=3; % FLY2 head

IDX = pathDB.IDX;
cmap_k = settings.cmap_k;

t = pathDB.t;
V = pathDB.V;
At_hor = pathDB.At_hor;
An_hor = pathDB.An_hor;

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_yaw = pathDB.stim_angle_yaw;
slip = pathDB.slip;
teta = pathDB.teta;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_yaw_temp = stim_angle_yaw;
slip_temp = slip;

for i=1:size(stim_angle_vel,2)
    for j=2:size(stim_angle_vel,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip(j-1,i) = nan;
        end
    end
end

subplot(3,2,1)
plot(-1,0,'^k','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',5)
hold on
plot(-1,0,'vk','MarkerFaceColor',settings.cmap_k(5,:),'MarkerSize',5)
legend('Fly 1', 'Fly 2','Location','Northwest')
plot(t,teta(:,fly1),'-k')
% plot(t,teta(:,fly1),'-k')
hold on
axis([0 .15 0 20])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:10:20,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('optical angle','fontsize',18)
grid on
    
subplot(3,2,3)
plot(t,stim_angle_vel(:,fly1),'-k')
hold on
plot(t,stim_angle_vel(:,fly2),'-k')
axis([0 .15 -180 180])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',-180:180:180,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('heading','fontsize',18)
grid on

subplot(3,2,5)
plot(t,stim_angle_yaw(:,fly1),'-k')
hold on
plot(t,stim_angle_yaw(:,fly2),'-k')
axis([0 .15 -180 180])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',-180:180:180,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('yaw','fontsize',18)
grid on

    
subplot(3,2,2)
plot(t,V(:,fly1),'-k')
hold on
plot(t,V(:,fly2),'-k')
axis([0 .15 0 .8])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:.4:1,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('V','fontsize',18)
grid on
    
subplot(3,2,4)
plot(t,abs(An_hor(:,fly1)),'-k')
hold on
plot(t,abs(An_hor(:,fly2)),'-k')
axis([0 .15 0 15])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',0:5:15,'fontsize',12)
xlabel('time','fontsize',18)
ylabel('An','fontsize',18)
grid on
    
subplot(3,2,6)
plot(t,At_hor(:,fly1),'-k')
hold on
plot(t,At_hor(:,fly2),'-k')
axis([0 .15 -10 15])
set(gca,'XTick',0:.05:.15) 
set(gca,'YTick',[-10;0;15],'fontsize',12)
xlabel('time','fontsize',18)
ylabel('At','fontsize',18)
grid on

    for j = 1:length(IDX)
        length(IDX)-j
        if isnan(IDX(j,fly1))==0
            
            
            subplot(3,2,1)
            plot(t(j),teta(j),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),teta(j),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,3)
            plot(t(j),stim_angle_vel(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),stim_angle_vel(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,5)
            plot(t(j),stim_angle_yaw(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),stim_angle_yaw(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)

            
            subplot(3,2,2)
            plot(t(j),V(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),V(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,4)
            plot(t(j),abs(An_hor(j,fly1)),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),abs(An_hor(j,fly2)),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)

            subplot(3,2,6)
            plot(t(j),At_hor(j,fly1),'^k','MarkerFaceColor',settings.cmap_k(IDX(j,fly1),:),'MarkerSize',5)
            hold on
            plot(t(j),At_hor(j,fly2),'vk','MarkerFaceColor',settings.cmap_k(IDX(j,fly2),:),'MarkerSize',5)
        end
    end

    
    saveas(gcf,['collision_response_2flies.fig'])
    saveas(gcf,['collision_response_2flies.png'])
    plot2svg('collision_response_2flies.svg')
