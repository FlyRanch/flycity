clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_9clusters_absAn_minAt_2.5n-3n2.5_response.mat')
mkdir('response_figs')
cd('response_figs')

IDX = pathDB.IDX;
cmap_k = settings.cmap_k;
cmap_360 = settings.cmap_360;

t = pathDB.t;
V = pathDB.V;
At_hor = pathDB.At_hor;
An_hor_NOabs = pathDB.An_hor;
An_hor = abs(pathDB.An_hor);

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_yaw = pathDB.stim_angle_yaw;
slip = pathDB.slip;

teta = patternDB.teta;

%% timelines NO skip NO abs(An_hor)
IDX_plot = IDX;
t_plot = t;

V_plot = V;
An_hor_plot = An_hor_NOabs;
At_hor_plot = At_hor;

stim_angle_vel_plot = stim_angle_vel;
stim_angle_yaw_plot = stim_angle_yaw;
slip_plot = slip;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_yaw_temp = stim_angle_yaw;
slip_temp = slip;

for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
    end
end

t_shift = responseDB.t_resp;
% plot_flightpath_timeline_tshift_cluster_slip
plot_flightpath_timeline_tshift_headingstart
% plot_flightpath_timeline_tshift_cluster_slip_long

subplot(3,2,4)
axis([-.05 .10 -15 15])

subplot(3,2,6)
axis([-.05 .10 -10 15])

saveas(gca,'flightpaths_headingUA_headingstart_noabs.fig')
saveas(gca,'flightpaths_headingUA_headingstart_noabs.png')
plot2svg('flightpaths_headingUA_headingstart_noabs.svg')

%% timelines NO skip INC abs(An_hor)
IDX_plot = IDX;
t_plot = t;

V_plot = V;
An_hor_plot = An_hor;
At_hor_plot = At_hor;

stim_angle_vel_plot = stim_angle_vel;
stim_angle_yaw_plot = stim_angle_yaw;
slip_plot = slip;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_yaw_temp = stim_angle_yaw;
slip_temp = slip;

for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
    end
end

t_shift = responseDB.t_resp;
% plot_flightpath_timeline_tshift_cluster_slip
plot_flightpath_timeline_tshift_headingstart
% plot_flightpath_timeline_tshift_cluster_slip_long
saveas(gca,'flightpaths_headingUA_headingstart.fig')
saveas(gca,'flightpaths_headingUA_headingstart.png')
plot2svg('flightpaths_headingUA_headingstart.svg')

%% timelines INC skip
skip = 50;

IDX_plot = IDX(1:skip:end,:);
t_plot = t(1:skip:end,:);

V_plot = V(1:skip:end,:);
An_hor_plot = An_hor(1:skip:end,:);
At_hor_plot = At_hor(1:skip:end,:);

stim_angle_vel_plot = stim_angle_vel(1:skip:end,:);
stim_angle_yaw_plot = stim_angle_yaw(1:skip:end,:);
slip_plot = slip(1:skip:end,:);

stim_angle_vel_temp = stim_angle_vel(1:skip:end,:);
stim_angle_yaw_temp = stim_angle_yaw(1:skip:end,:);
slip_temp = slip(1:skip:end,:);

for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
    end
end

% t_shift = zeros(size(V,2),1);
% plot_flightpath_timeline_tshift_cluster_slip
% saveas(gca,'flightpaths_headingUA_NOshift.fig')
% saveas(gca,'flightpaths_headingUA_NOshift.png')
% plot2svg('flightpaths_headingUA_NOshift.svg')

t_shift = responseDB.t_resp;
% plot_flightpath_timeline_tshift_cluster_slip
plot_flightpath_timeline_tshift_cluster_slip_axisfix
% plot_flightpath_timeline_tshift_cluster_slip_long
saveas(gca,'flightpaths_headingUA_tresp2.fig')
saveas(gca,'flightpaths_headingUA_tresp2.png')
plot2svg('flightpaths_headingUA_tresp2.svg')

% t_shift = responseDB.t_resp;
% % plot_flightpath_timeline_tshift_cluster_slip
% plot_flightpath_timeline_tshift_cluster_slip_long
% % plot_flightpath_timeline_tshift_cluster_slip_long
% saveas(gca,'flightpaths_headingUA_tresp_long.fig')
% saveas(gca,'flightpaths_headingUA_tresp_long.png')
% plot2svg('flightpaths_headingUA_tresp_long.svg')

% t_shift = responseDB.t_turn_start;
% plot_flightpath_timeline_tshift_cluster_slip
% saveas(gca,'flightpaths_headingUA_tturnstart.fig')
% saveas(gca,'flightpaths_headingUA_tturnstart.png')
% plot2svg('flightpaths_headingUA_tturnstart.svg')
% 
% t_shift = responseDB.t_turn_max;
% plot_flightpath_timeline_tshift_cluster_slip
% saveas(gca,'flightpaths_headingUA_tturnmax.fig')
% saveas(gca,'flightpaths_headingUA_tturnmax.png')
% plot2svg('flightpaths_headingUA_tturnmax.svg')
% 
% % slip
% t_shift = responseDB.t_resp;
% plot_sideslip_timeline
% saveas(gca,'flightpaths_slip_tresp.fig')
% saveas(gca,'flightpaths_slip_tresp.png')
% plot2svg('flightpaths_slip_tresp.svg')

%% plot histograms
close all
nx = 1000;
ny = 100;


% % no time shift
% t_shift = 0;
% t_hist = repmat(t,size(V,2),1);
% 
% plot_flightpath_timeline_histograms
% saveas(gca,'flightpath_hist_slip_notshift.fig')
% saveas(gca,'flightpath_hist_slip_notshift.png')
% plot2svg('flightpath_hist_slip_notshift.svg')
% 
% 
% % t_resp timeshift
% t_shift = responseDB.t_resp;
% t_hist = [];
% for i=1:size(V,2)    
%     t_hist = [t_hist;t-t_shift(i)];
% end
% plot_flightpath_timeline_histograms
% 
% saveas(gca,'flightpath_hist_slip_tresp.fig')
% saveas(gca,'flightpath_hist_slip_tresp.png')
% plot2svg('flightpath_hist_slip_tresp.svg')

%% hist norm
% % no time shift
% t_shift = 0;
% t_hist = repmat(t,size(V,2),1);
% plot_flightpath_timeline_histograms_norm
% 
% saveas(gca,'flightpath_hist_norm_slip_notshift.fig')
% saveas(gca,'flightpath_hist_norm_slip_notshift.png')
% plot2svg('flightpath_hist_norm_slip_notshift.svg')
% 

% t_resp timeshift
t_shift = responseDB.t_resp;
t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
plot_flightpath_timeline_histograms_norm_axisfix

saveas(gca,'flightpath_hist_norm_slip_tresp.fig')
saveas(gca,'flightpath_hist_norm_slip_tresp.png')
plot2svg('flightpath_hist_norm_slip_tresp.svg')

%% plot turn vs decel reaction time
t_Anresp = responseDB.t_turn_start;
t_Atresp = responseDB.t_decel_start;

figure
plotcolor = 'b';
plot_tAntAt

axis equal
grid on
xlabel('turn time','fontsize',18) 
ylabel('deceleration time','fontsize',18) 
set(gca,'xlim',[0 .1],'ylim',[0 .1])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2],'fontsize',12)
saveas(gca,'response_time_tturn_vs_tdecel.fig')
saveas(gca,'response_time_tturn_vs_tdecel.png')
plot2svg('response_time_tturn_vs_tdecel.svg')

%% plot turn vs accel reaction time
t_Anresp = responseDB.t_turn_start;
t_Atresp = responseDB.t_accel_start;

figure
plotcolor = 'r';
plot_tAntAt

axis equal
grid on
xlabel('turn time','fontsize',18) 
ylabel('acceleration time','fontsize',18) 
set(gca,'xlim',[0 .1],'ylim',[0 .1])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2],'fontsize',12)
saveas(gca,'response_time_tturn_vs_taccel.fig')
saveas(gca,'response_time_tturn_vs_taccel.png')
plot2svg('response_time_tturn_vs_taccel.svg')

%% plot turn vs accel and decel reaction time
t_Anresp = responseDB.t_turn_start;

figure
hold on

plotcolor = 'b';
t_Atresp = responseDB.t_decel_start;
plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plotcolor = 'r';
t_Atresp = responseDB.t_accel_start;
plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plot([0,1],[0,1],'--k')
legend('decel','accel')

axis equal
grid on
xlabel('turn response time','fontsize',18) 
ylabel('acceleration response time','fontsize',18) 
set(gca,'xlim',[0 .1],'ylim',[0 .1])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2],'fontsize',12)
saveas(gca,'response_time_tAn_vs_tAt.fig')
saveas(gca,'response_time_tAn_vs_tAt.png')
plot2svg('response_time_tAn_vs_tAt.svg')

%% calc yaw & vel turns
t = pathDB.t;
t_pre = responseDB.t_resp;
t_post = responseDB.t_turn_stop;
calc_turn_vectors

% %% turn vs heading
% figure
% plotcolor = 'r'
% angle_pre = responseDB.stim_angle_vel_pre_resp;
% turn = turn_angle_vel
% plot_angle_pre_turn_vectors
% 
% xlabel('heading pre')
% ylabel('turn angle')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'turnangle_vs_preheading.fig')
% saveas(gca,'turnangle_vs_preheading.png')
% 
% %% yaw post vs yaw pre
% figure
% plotcolor = 'b'
% angle_pre = responseDB.stim_angle_yaw_pre_resp;
% turn = turn_angle_yaw
% % plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors
% 
% xlabel('heading pre')
% ylabel('turn angle')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'yawpost_vs_yawpre.fig')
% saveas(gca,'yawpost_vs_yawpre.png')

%% turn n yaw
figure
plot(360,360,'ok','MarkerFaceColor','r','markersize',10)
hold on
plot(360,360,'ok','MarkerFaceColor','b','markersize',10)
legend('heading','yaw')

angle_pre = responseDB.stim_angle_vel_pre_resp;
turn = turn_angle_vel
plotcolor = 'r'
plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors

angle_pre = responseDB.stim_angle_yaw_pre_resp;
turn = turn_angle_yaw
plotcolor = 'b'
plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors

xlabel('heading pre','fontsize',18) 
ylabel('turn angle','fontsize',18) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',12)
saveas(gca,'heading_vs_yaw.fig')
saveas(gca,'heading_vs_yaw.png')
plot2svg('heading_vs_yaw.svg')
%% escape angle vs heading
% angle_pre = responseDB.stim_angle_vel_pre_resp;
% angle_post = responseDB.stim_angle_vel_post_resp;
% 
% figure
% plotcolor = 'k'
% plot_angle_pre_turn
% 
% xlabel('heading pre')
% ylabel('turn angle')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'escapeangle_vs_preheading.fig')
% saveas(gca,'escapeangle_vs_preheading.png')

%% slip angle effect
% figure
% plot(responseDB.slip_pre_resp,responseDB.stim_angle_vel_post_turn,'.r')
% hold on
% plot(responseDB.slip_pre_resp,responseDB.stim_angle_yaw_post_turn,'.b')
% xlabel('slip pre')
% ylabel('heading and yaw post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% figure
% plot(responseDB.slip_pre_resp,responseDB.slip_post_resp,'.r')
% xlabel('slip pre')
% ylabel('slip post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on

%% plot pre&post velocity vectors
% 
% figure
% subplot(1,2,1)
% % plot(0,0)
% % hold on
% % plot(0,0,'r')
% % legend('pre','post')
% % compass_zeroup(u_post,v_post,'r');
% % hold on
% % compass_zeroup(u_pre,v_pre,'b');
% Z = compass(u_post,v_post);
% for i=1:length(Z)
%     set(Z(i),'color',cmap(end,:))
% end
% hold on
% Y = compass(u_pre,v_pre);
% for i=1:length(Y)
%     set(Y(i),'color',cmap(1,:))
% end 
% title('horizontal flight speed vectors')
% % saveas(gca,'Vhorizontal_PreVsPost.fig')
% % saveas(gca,'Vhorizontal_PreVsPost.png')
% 
% Z = rose(deg2rad(heading_post),36);
% for i=1:length(Z)
%     set(Z(i),'color',cmap(end,:))
% end
% 
% hold on
% Y = rose(deg2rad(heading_pre),36);
% for i=1:length(Y)
%     set(Y(i),'color',cmap(1,:))
% end 
% title('heading histogram')
% % saveas(gca,'Vhorizontal_PreVsPost.fig')
% % saveas(gca,'Vhorizontal_PreVsPost.png')
% 
% %% plot pre&post vertical velocity vectors
% % figure
% % compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
% % hold on
% % compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
% subplot(1,2,2)
% Z = compass(sqrt(u_post.^2 + v_post.^2),w_post);
% for i=1:length(Z)
%     set(Z(i),'color',cmap(end,:))
% end
% hold on
% Y = compass(sqrt(u_pre.^2 + v_pre.^2),w_pre);
% for i=1:length(Y)
%     set(Y(i),'color',cmap(1,:))
% end 
% title('vertical flight speed vectors')
% saveas(gca,'V_PreVsPost.fig')
% saveas(gca,'V_PreVsPost.png')
% 
% %% plot histograms
% % figure
% subplot(3,2,1)
% dh = 360/(25)
% heading_noresp = heading(1:Nnoresp,:);
% heading_noresp = heading_noresp(isnan(heading_noresp)==0);
% hist(heading_noresp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('heading pre')
% set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[-180 180])
% % set(gca,'XTick',[-180 -90 0 90 180])
% % saveas(gca,'hist_preheading.fig')
% % saveas(gca,'hist_preheading.png')
% % plot2svg
% 
% % figure
% subplot(3,2,2)
% heading_resp = heading(Nresp:end,:);
% heading_resp = heading_resp(isnan(heading_resp)==0);
% hist(heading_resp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('heading post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% % saveas(gca,'hist_postheading.fig')
% % saveas(gca,'hist_postheading.png')
% % plot2svg
% 
% % figure
% subplot(3,2,3)
% dV = 1/25
% V_noresp = V(1:Nnoresp,:);
% V_noresp = V_noresp(isnan(V_noresp)==0);
% hist(V_noresp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('velocity pre')
% % set(gca,'ylim',[0 50000])
% 
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% % saveas(gca,'hist_preV.fig')
% % saveas(gca,'hist_preV.png')
% % plot2svg
% 
% % figure
% subplot(3,2,4)
% V_resp = V(Nresp:end,:);
% V_resp = V_resp(isnan(V_resp)==0);
% hist(V_resp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('velocity post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% % saveas(gca,'hist_postV.fig')
% % saveas(gca,'hist_postV.png')
% % plot2svg
% 
% % figure
% subplot(3,2,5)
% dA = 30/25
% A_noresp = A(1:Nnoresp,:);
% A_noresp = A_noresp(isnan(A_noresp)==0);
% hist(A_noresp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('accel pre')
% % set(gca,'ylim',[0 100000])
% 
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% % saveas(gca,'hist_preA.fig')
% % saveas(gca,'hist_preA.png')
% % plot2svg
% 
% % figure
% subplot(3,2,6)
% A_resp = A(Nresp:end,:);
% A_resp = A_resp(isnan(A_resp)==0);
% hist(A_resp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('accel post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% % saveas(gca,'hist_postA.fig')
% % saveas(gca,'hist_postA.png')
% % plot2svg
% 
% saveas(gca,'histograms_pre_post.fig')
% saveas(gca,'histograms_pre_post.png')
% % plot2svg
% 
% %% plot top view flight path vectors
% figure
% colormap(cmap)
% caxis([Nstart Nstop])
% hold on
% for j=1:size(x,2)
%     plot(x(Nstart:Nstop,j),y(Nstart:Nstop,j),'-','color',cmap(1,:))
% end
% 
% vec_col = 0;
% for i=Nstart:dN:Nstop
%     vec_col = vec_col + size(cmap,1)/Nstop;
%     for j=1:size(x,2)
% %         quiverc(x(i,:),y(i,:),u(i,:),v(i,:),'color',cmap(1,:))
%         if isnan(x(i,j)) == 0
%             quivert(x(i,j),y(i,j),u(i,j),v(i,j),i,'as',1/500,'ahr',[1 1],'nt')
%         end
%     end
%     
% end
%     
% axis equal
% colorbar
% 
% % set(gca,'xlim',[0 30])
% % set(gca,'XTick',[0 10 20 30])
% saveas(gca,'flighttracks.fig')
% saveas(gca,'flighttracks.png')
% % plot2svg

%% ON/OFF escape performance: heading, yaw, resp time, U, An, At

plot_escape_perf_ONOFF_NO165tetamax

% title('acceleration based response time')
saveas(gca,'escape performance_ONOFF.fig')
saveas(gca,'escape performance_ONOFF.png')
plot2svg('escape performance_ONOFF.svg')

%% teta_max escape performance: heading, yaw, resp time, U, An, At

plot_escape_perf_tetamax

% title('acceleration based response time')
saveas(gca,'escape performance_tetamax.fig')
saveas(gca,'escape performance_tetamax.png')
plot2svg('escape performance_tetamax.svg')

cd ..