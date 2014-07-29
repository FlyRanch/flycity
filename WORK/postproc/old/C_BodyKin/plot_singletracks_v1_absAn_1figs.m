% plot timelines for single tracks
clear
clc
close all

addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
load('flightpathDB_pos_qbodyEKF_9clusters_absAn_minAt_2.5n-3n2.5_response.mat')

mkdir('flightpaths')
cd('flightpaths')


IDX = pathDB.IDX;
cmap_k = settings.cmap_k;

t = pathDB.t;
V = pathDB.V;
At_hor = pathDB.At_hor;
An_hor = pathDB.An_hor;

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_yaw = pathDB.stim_angle_yaw;
slip = pathDB.slip;
teta = patternDB.teta;

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


skip = 35
for i = 156:size(V,2)
    size(V,2)-i
%     close all
    
    subplot(3,2,1)
    hold off
    plot(t,patternDB.teta(:,i),'-k')
    hold on
    axis([0 .15 0 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:90:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('optical angle','fontsize',18)
    grid on

    subplot(3,2,3)
    hold off
    plot(t,stim_angle_vel(:,i),'-k')
    hold on
    axis([0 .15 -180 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('heading','fontsize',18)
    grid on

    subplot(3,2,5)
    hold off
    plot(t,stim_angle_yaw(:,i),'-k')
    hold on
    axis([0 .15 -180 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('yaw','fontsize',18)
    grid on


    subplot(3,2,2)
    hold off
    plot(t,V(:,i),'-k')
    hold on
    axis([0 .15 0 .8])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:.4:1,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('V','fontsize',18)
    grid on

    subplot(3,2,4)
    hold off
    plot(t,abs(An_hor(:,i)),'-k')
    hold on
    axis([0 .15 0 15])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:5:15,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('An','fontsize',18)
    grid on

    subplot(3,2,6)
    hold off
    plot(t,At_hor(:,i),'-k')
    hold on
    axis([0 .15 -10 15])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',[-10;0;15],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('At','fontsize',18)
    grid on

        for j = 1:skip:length(IDX)
            length(IDX)-j;
            if isnan(IDX(j,i))==0


                subplot(3,2,1)
                plot(t(j),patternDB.teta(j,i),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)

                subplot(3,2,3)
                plot(t(j),stim_angle_vel(j,i),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)
                
                subplot(3,2,5)
                plot(t(j),stim_angle_yaw(j,i),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)

                subplot(3,2,2)
                plot(t(j),V(j,i),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)

                subplot(3,2,4)
                plot(t(j),abs(An_hor(j,i)),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)

                subplot(3,2,6)
                plot(t(j),At_hor(j,i),'ok','MarkerFaceColor',settings.cmap_k(IDX(j,i),:),'MarkerSize',5)
            end
        end


        saveas(gcf,['flightpath_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.fig'])
        saveas(gcf,['flightpath_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.png'])
        plot2svg(['flightpath_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.svg'])
end
cd ..