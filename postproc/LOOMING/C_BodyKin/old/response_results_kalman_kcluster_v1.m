clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_9clusters_absAn_minAt_2.5n-3n2.5_response.mat')
mkdir('response_figs')
cd('response_figs')

t = pathDB.t;
V = pathDB.V;
At_hor = pathDB.At_hor;
An_hor = pathDB.An_hor;

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_yaw = pathDB.stim_angle_yaw;
stim_angle_vel_plot = stim_angle_vel;
stim_angle_yaw_plot = stim_angle_yaw;
for i=1:size(V,2)
    for j=2:size(V,1)
        if abs(stim_angle_vel(j,i) - stim_angle_vel(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw(j,i) - stim_angle_yaw(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
    end
end

teta = patternDB.teta;

%% plot timelines
t_shift = zeros(size(V,2),1);
plot_flightpath_timeline_tshift_ONOFF
saveas(gca,'flightpaths_headingUA_NOshift.fig')
saveas(gca,'flightpaths_headingUA_NOshift.jpg')
plot2svg('flightpaths_headingUA_NOshift.jpg')

t_shift = responseDB.t_resp;
plot_flightpath_timeline_tshift_ONOFF
saveas(gca,'flightpaths_headingUA_tresp.fig')
saveas(gca,'flightpaths_headingUA_tresp.jpg')
plot2svg('flightpaths_headingUA_tresp.jpg')

t_shift = responseDB.t_turn_start;
plot_flightpath_timeline_tshift_ONOFF
saveas(gca,'flightpaths_headingUA_tturnstart.fig')
saveas(gca,'flightpaths_headingUA_tturnstart.jpg')
plot2svg('flightpaths_headingUA_tturnstart.jpg')

t_shift = responseDB.t_turn_max;
plot_flightpath_timeline_tshift_ONOFF
saveas(gca,'flightpaths_headingUA_tturnmax.fig')
saveas(gca,'flightpaths_headingUA_tturnmax.jpg')
plot2svg('flightpaths_headingUA_tturnmax.jpg')

%% plot histograms
% no timeshift
var = stim_angle_vel(:);
t_hist = repmat(t,136,1);
plot_flightpath_timeline_histograms
pause

var = stim_angle_yaw(:);
t_hist = repmat(t,136,1);
plot_flightpath_timeline_histograms
pause

% t_resp timeshift
t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
var = stim_angle_vel(:);
plot_flightpath_timeline_histograms
pause

var = stim_angle_yaw(:);
plot_flightpath_timeline_histograms
pause

%% plot turn vs decel reaction time
t_Anresp = responseDB.t_turn_start;
t_Atresp = responseDB.t_decel_start;

figure
plot_tAntAt
% set(gca,'xlim',[0 .1],'ylim',[0 .1])
% set(gca,'XTick',[0:.05:.1])
% set(gca,'YTick',[0:.05:.1])
set(gca,'xlim',[0 .15],'ylim',[0 .15])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2])
saveas(gca,'response_time_tturn_vs_tdecel.fig')
saveas(gca,'response_time_tturn_vs_tdecel.jpg')

%% plot turn vs accel reaction time
t_Anresp = responseDB.t_turn_start;
t_Atresp = responseDB.t_accel_start;

figure
plot_tAntAt
% set(gca,'xlim',[0 .1],'ylim',[0 .1])
% set(gca,'XTick',[0:.05:.1])
% set(gca,'YTick',[0:.05:.1])
set(gca,'xlim',[0 .2],'ylim',[0 .2])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2])
saveas(gca,'response_time_tturn_vs_taccel.fig')
saveas(gca,'response_time_tturn_vs_taccel.jpg')


%% turn vs heading
angle_pre = responseDB.stim_angle_vel_pre_resp;
angle_post = responseDB.stim_angle_vel_post_turn;

figure
plotcolor = 'k'
plot_angle_pre_turn

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'turnangle_vs_preheading.fig')
saveas(gca,'turnangle_vs_preheading.jpg')

%% escape angle vs heading
angle_pre = responseDB.stim_angle_vel_pre_resp;
angle_post = responseDB.stim_angle_vel_post_resp;

figure
plotcolor = 'k'
plot_angle_pre_turn

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'escapeangle_vs_preheading.fig')
saveas(gca,'escapeangle_vs_preheading.jpg')

%% yaw post vs yaw pre
angle_pre = responseDB.stim_angle_yaw_pre_resp;
angle_post = responseDB.stim_angle_yaw_post_turn;

figure
plotcolor = 'k'
plot_angle_pre_turn_nowrap

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'yawpost_vs_yawpre.fig')
saveas(gca,'yawpost_vs_yawpre.jpg')

%% turn n yaw
figure
plot(360,360,'.r')
hold on
plot(360,360,'.b')
legend('heading','yaw')

angle_pre = responseDB.stim_angle_vel_pre_resp;
angle_post = responseDB.stim_angle_vel_post_turn;
plotcolor = 'r'
plot_angle_pre_turn

angle_pre = responseDB.stim_angle_yaw_pre_resp;
angle_post = responseDB.stim_angle_yaw_post_turn;
plotcolor = 'b'
plot_angle_pre_turn_nowrap

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'heading_vs_yaw.fig')
saveas(gca,'heading_vs_yaw.jpg')

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
% % saveas(gca,'Vhorizontal_PreVsPost.jpg')
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
% % saveas(gca,'Vhorizontal_PreVsPost.jpg')
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
% saveas(gca,'V_PreVsPost.jpg')
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
% % saveas(gca,'hist_preheading.jpg')
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
% % saveas(gca,'hist_postheading.jpg')
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
% % saveas(gca,'hist_preV.jpg')
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
% % saveas(gca,'hist_postV.jpg')
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
% % saveas(gca,'hist_preA.jpg')
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
% % saveas(gca,'hist_postA.jpg')
% % plot2svg
% 
% saveas(gca,'histograms_pre_post.fig')
% saveas(gca,'histograms_pre_post.jpg')
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
% saveas(gca,'flighttracks.jpg')
% % plot2svg

%% ON/OFF escape performance: heading, yaw, resp time, U, An, At

plot_escape_perf_ONOFF

% title('acceleration based response time')
saveas(gca,'escape performance_ONOFF.fig')
saveas(gca,'escape performance_ONOFF.jpg')
plot2svg('escape performance_ONOFF.svg')

%% teta_max escape performance: heading, yaw, resp time, U, An, At

plot_escape_perf_tetamax

% title('acceleration based response time')
saveas(gca,'escape performance_tetamax.fig')
saveas(gca,'escape performance_tetamax.jpg')
plot2svg('escape performance_tetamax.svg')

cd ..