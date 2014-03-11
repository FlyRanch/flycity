clc
clear
close all
load('flightpathDB_pos_INCq.mat')

mkdir('figs')
cd('figs')

t = pathDB.t;

x = pathDB.pos(:,:,1);
y = pathDB.pos(:,:,2);
z = pathDB.pos(:,:,3);

u = pathDB.vel(:,:,1);
v = pathDB.vel(:,:,2);
w = pathDB.vel(:,:,3);

ax = pathDB.accel(:,:,1);
ay = pathDB.accel(:,:,2);
az = pathDB.accel(:,:,3);

alpha_dot = pathDB.alpha_dot;
At = pathDB.At;
An = pathDB.An;

%% make colormap
% % white to dark
% cmap1 = [1:-(1-.75)/(length(t)-1):0.75];
% cmap2 = [1:-(1-0)/(length(t)-1):0];
% cmap3 = [1:-(1-0)/(length(t)-1):0];

% % light grey to dark
% cmap1 = [.5:-(.5-.75)/(length(t)-1):0.75];
% cmap2 = [.5:-(.5-0)/(length(t)-1):0];
% cmap3 = [.5:-(.5-0)/(length(t)-1):0];

% dark grey to dark
cmap1 = [.25:-(.25-.75)/(length(t)-1):0.75];
cmap2 = [.25:-(.25-0)/(length(t)-1):0];
cmap3 = [.25:-(.25-0)/(length(t)-1):0];

% % dark grey to dark
% cmap2 = [.25:-(.25-0)/(length(t)-1):0];
% cmap3 = [.25:-(.75-0)/(length(t)-1):0];
% cmap1 = .25* ones(size(cmap2));

% % light
% cmap2 = [1:-(1-0)/(length(t)-1):0];
% cmap3 = [1:-(1-0)/(length(t)-1):0];
% cmap1 = 1* ones(size(cmap2));

% grey to red
cmap = [cmap1; cmap2; cmap3]';

% temporal data
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;
t_start = -0.05;
t_stop = 0.2;



DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

Nnoresp=find(t==0)
Nresp=find(t==0.1)

Nstart=find(t==t_start)
Nstop=find(t==t_stop)
dN = 30;

% arena data
% arenacenter from fly positions
center = [nanmean(x(trigger_frame,:)) nanmean(y(trigger_frame,:)) nanmean(z(trigger_frame,:))];
r = 0.120; % [m] arena radius

% stimulus data
pos_stim = center - [0 r 0];
Vstim = [0;1] % direction
heading_stim = atan2(Vstim(1),Vstim(2)) *180/pi()

% reverse data with reverse pattern
for i = 1:size(x,2)
    if settings.expansion.HorPos(i) == 180
        y(:,i) = -y(:,i) + 2*center(2);
        v(:,i) = -v(:,i);
        ay(:,i) = -ay(:,i);
    end
end

% pattern direction (in fly frame of reference)
for i = 1:size(x,1)
    for j = 1:size(x,2)
        
        % 2d velocity vector
        x1 = u(i,j);
        y1 = v(i,j);
        
        % fly2stim 2d vector
        x2 = x(i,j) - pos_stim(1);
        y2 = y(i,j) - pos_stim(2);
        
        stim_angle(i,j) = atan2(x1*y2-y1*x2,x1*x2+y1*y2) *180/pi();
        
        x2=1;
        y2=0;
        dir(i,j) = -atan2(x1*y2-y1*x2,x1*x2+y1*y2) *180/pi();
        
        % remove jumps from 180deg to -180deg
        stim_angle_plot(i,j) = stim_angle(i,j);
        if i>1 && j>1 && abs(stim_angle_plot(i,j) - stim_angle_plot(i-1,j))>180
            stim_angle_plot(i,j) = nan;
        end
    end
end

side_slip = dir - pathDB.yaw_global;
side_slip(abs(side_slip)>90) = nan;

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);

% response threshold
alpha_dot_resp = 3000;
A_resp = 8;
At_resp = 3;
An_resp = 5;

% response time
t_now = t(trigger_frame:end);
for i=1:size(A,2)
    A_now = A(trigger_frame:end,i);
    At_now = At(trigger_frame:end,i);
    An_now = An(trigger_frame:end,i);
    alpha_dot_now = alpha_dot(trigger_frame:end,i);
    
    A_max(i) = max(A_now);
    t_Amax(i) = min(t_now(A_now==max(A_now)));
    tA_now = t_now(A_now>A_resp);
    if isempty(tA_now) == 1
        t_Aresp(i) = nan;
    else
        t_Aresp(i) = min(tA_now);
    end
    
    At_max(i) = max(At_now);
    t_Atmax(i) = min(t_now(At_now==max(At_now)));
    tAt_now = t_now(At_now>At_resp);
    if isempty(tAt_now) == 1
        t_Atresp(i) = nan;
    else
        t_Atresp(i) = min(tAt_now);
    end
    
    An_max(i) = max(An_now);
    t_Anmax(i) = min(t_now(An_now==max(An_now)));
    tAn_now = t_now(An_now>An_resp);
    if isempty(tAn_now) == 1
        t_Anresp(i) = nan;
    else
        t_Anresp(i) = min(tAn_now);
    end
    
    alpha_dot_max(i) = max(alpha_dot_now);
    t_alpha_dot_max(i) = min(t_now(alpha_dot_now==max(alpha_dot_now)));
    ta_now = t_now(alpha_dot_now>alpha_dot_resp);
    if isempty(ta_now) == 1
        t_alpha_dot_resp(i) = nan;
    else
        t_alpha_dot_resp(i) = min(ta_now);
    end
end


% vector at trigger (pre response)
u_pre = mean(u(trigger_frame:trigger_frame+DN_noresp,:))';
v_pre = mean(v(trigger_frame:trigger_frame+DN_noresp,:))';
w_pre = mean(w(trigger_frame:trigger_frame+DN_noresp,:))';

side_slip_pre = nanmean(side_slip(trigger_frame:trigger_frame+DN_noresp,:))';

V_pre = sqrt(u_pre.^2 + v_pre.^2 + w_pre.^2);
heading_pre = rad2deg(atan2(u_pre,v_pre));
ascend_pre = rad2deg(atan2(w_pre,sqrt(u_pre.^2 + v_pre.^2)));

% for i=1:size(u,2)
%     
%     % vector at end (post response)
%     last_frame = find(isnan(u(:,i))==0, 1, 'last' );
%     
%     z_post(i,1) = mean(z(last_frame-DN_escape:last_frame,i))';
%     
%     u_post(i,1) = mean(u(last_frame-DN_escape:last_frame,i))';
%     v_post(i,1) = mean(v(last_frame-DN_escape:last_frame,i))';
%     w_post(i,1) = mean(w(last_frame-DN_escape:last_frame,i))';
% end
% 
% V_post = sqrt(u_post.^2 + v_post.^2 + w_post.^2);
% heading_post = rad2deg(atan2(u_post,v_post));
% ascend_post = rad2deg(atan2(w_post,sqrt(u_post.^2 + v_post.^2)));
% ascend_turn = ascend_post - ascend_pre;
% 
% % turn heading < 180
% heading_turn = heading_post - heading_pre;
% heading_turn(heading_turn>180) = heading_turn(heading_turn>180) - 360;
% 
% % heading_turn2 > 0 if animal turns away from looming object, heading_turn2 < 0: towards
% heading_turn2 = abs(heading_turn);
% heading_turn2(find(abs(heading_post)) > abs(heading_pre)) = - heading_turn2((find(abs(heading_post)) > abs(heading_pre)));
% 
% % turn rate
% heading = rad2deg(atan2(u,v));
% ascend = rad2deg(atan2(w,sqrt(u.^2 + v.^2)));
% 
% turn_rate = nan(size(heading));
% for i=2:size(heading,1)-2
%  turn_rate(i,:) = heading(i+1,:) - heading(i,:);
% end
% turn_rate = turn_rate/dt;
% % 
% % turn_rate(abs(turn_rate)>8000)  = nan;
% 
% figure
% plot(0,0,'color',cmap(1,:))
% hold on
% plot(0,0,'color',cmap(end,:))
% legend('slow','fast')
% 
% for i=1:size(u,2)
%     if settings.expansion.speed(i) == 0
%         plot(t,stim_angle_plot(:,i),'color',cmap(1,:))
%     else
%         plot(t,stim_angle_plot(:,i),'color',cmap(end,:))
%     end
% end
% % grid on
% xlabel('t')
% ylabel('heading')
% axis([-0.2 0.2 -180 180])
% colormap(cmap)
% saveas(gca,'temporal_stim_angle.fig')
% saveas(gca,'temporal_stim_angle.jpg')
% 
% figure
% plot(0,0,'color',cmap(1,:))
% hold on
% plot(0,0,'color',cmap(end,:))
% legend('slow','fast')
% 
% for i=1:size(u,2)
%     if settings.expansion.speed(i) == 0
%         plot(t,V(:,i),'color',cmap(1,:))
%     else
%         plot(t,V(:,i),'color',cmap(end,:))
%     end
% end
% % grid on
% xlabel('t')
% ylabel('V')
% axis([-0.2 0.2 0 1])
% colormap(cmap)
% saveas(gca,'temporal_flightspeed.fig')
% saveas(gca,'temporal_flightspeed.jpg')
% 
% figure
% plot(0,0,'color',cmap(1,:))
% hold on
% plot(0,0,'color',cmap(end,:))
% legend('slow','fast')
% 
% for i=1:size(u,2)
%     if settings.expansion.speed(i) == 0
%         plot(t,A(:,i),'color',cmap(1,:))
%     else
%         plot(t,A(:,i),'color',cmap(end,:))
%     end
% end
% % grid on
% xlabel('t')
% ylabel('A')
% axis([-0.2 0.2 0 30])
% colormap(cmap)
% saveas(gca,'temporal_acceleration.fig')
% saveas(gca,'temporal_acceleration.jpg')
% 
% figure
% subplot(3,1,1)
% plot(0,0,'color',cmap(1,:))
% hold on
% plot(0,0,'color',cmap(end,:))
% legend('slow','fast')
% 
% for i=1:size(u,2)
%     if settings.expansion.speed(i) == 0
%         subplot(3,1,1)
%         plot(t,stim_angle_plot(:,i),'color',cmap(1,:))
%         hold on
% 
%         subplot(3,1,2)
%         plot(t,V(:,i),'color',cmap(1,:))
%         hold on
% 
%         subplot(3,1,3)
%         plot(t,A(:,i),'color',cmap(1,:))
%         hold on
%         
%     else
%         subplot(3,1,1)
%         plot(t,stim_angle_plot(:,i),'color',cmap(end,:))
%         hold on
% 
%         subplot(3,1,2)
%         plot(t,V(:,i),'color',cmap(end,:))
%         hold on
% 
%         subplot(3,1,3)
%         plot(t,A(:,i),'color',cmap(end,:))
%         hold on
%     end
% end
% subplot(3,1,1)
% % grid on
% xlabel('time (s)')
% ylabel('heading (deg)')
% axis([-0.2 0.2 -180 180])
% % set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% 
% subplot(3,1,2)
% % grid on
% xlabel('time (s)')
% ylabel('flight speed (m/s)')
% axis([-0.2 0.2 0 1])
% 
% subplot(3,1,3)
% % grid on
% xlabel('time (s)')
% ylabel('acceleration (m/s^2)')
% axis([-0.2 0.2 0 30])
% 
% saveas(gca,'flightpaths_heading_V_A.fig')
% saveas(gca,'flightpaths_heading_V_A.jpg')
% % plot2svg(flightpaths_heading_V_A,gca)
% 
% heading_pre_calc = heading_pre;
% heading_turn_calc = heading_turn;
% 
% for i=1:size(heading_turn,1)
%     if heading_pre(i) < -90 && heading_turn(i) <-90
%         heading_turn_calc(i) = heading_turn_calc(i) + 360;
%     end
% end
% 
% heading_pre_ext = heading_pre_calc;
% heading_turn_ext = heading_turn_calc;
% 
% % wrap around
% heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc -360;
% heading_turn_ext(end+1:end+length(heading_turn_calc)) = heading_turn_calc +360;
% 
% heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc +360;
% heading_turn_ext(end+1:end+length(heading_turn_calc)) = heading_turn_calc -360;
% 
% % 1st order polyfit
% headingNturn = [heading_pre_ext heading_turn_ext];
% headingNturn = sortrows(headingNturn);
% 
% [p_turn,S_turn] = polyfit(headingNturn(:,1),headingNturn(:,2),1);
% heading_sort = headingNturn(:,1);
% turn_sort = headingNturn(:,2);
% turn_sort_polyfit = polyval(p_turn,heading_sort);
% 
% p=0
% dn=20
% dm=10
% clear dturn_std_now dturn_ste_now turn_now heading_now
% for n = 1:dm:(length(heading_sort)-dn)
%     turn_sub = turn_sort(n:n+dn-1);
%     turn_polyfit_sub = turn_sort_polyfit(n:n+dn-1);
%     heading_sub = heading_sort(n:n+dn-1);
%     dturn_sub = turn_sub - turn_polyfit_sub;
% 
%     p=p+1;
%     dturn_std_now(p) = nanstd(dturn_sub);
%     dturn_ste_now(p) = nanstd(dturn_sub)./sqrt(length(dturn_sub));
%     turn_now(p) = mean(turn_sub);
%     heading_now(p) = mean(heading_sub);
% end
% 
% turn_high = turn_now + 1.96*dturn_ste_now;
% turn_low = turn_now - 1.96*dturn_ste_now;
% 
% for i=1:size(heading_turn,1)
%     if heading_pre(i) > 90 && heading_turn(i) <-180
%         heading_turn(i) = heading_turn(i) + 360;
%     end
% end
% 
% figure
% plot(0,0,'.','markersize',25,'color',cmap(1,:))
% hold on
% plot(0,0,'.','markersize',25,'color',cmap(end,:))
% legend('slow','fast')
% hold off
% hold on
% ciplot(turn_low',turn_high',heading_now',cmap(1,:))
% alpha(.25)
% plot(heading_now,turn_now,'-k')
% plot([-180,180],[180,-180],'--k')
% 
% for i=1:size(heading_pre,1)
%     if settings.expansion.speed(i) == 0
%         plot(heading_pre(i),heading_turn(i),'.','markersize',25,'color',cmap(1,:))
%         hold on
%     else
%         plot(heading_pre(i),heading_turn(i),'.','markersize',25,'color',cmap(end,:))
%         hold on
%     end
% end
% % grid on
% axis equal
% xlabel('pre-heading')
% ylabel('turn angle')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'turnangle_vs_preheading.fig')
% saveas(gca,'turnangle_vs_preheading.jpg')
% 
% 
% 
% heading_pre_calc = heading_pre;
% heading_post_calc = heading_post;
% 
% % 1st order polyfit
% preNpost = [heading_pre_calc heading_post_calc];
% preNpost = sortrows(preNpost);
% 
% [p_post,S_post] = polyfit(preNpost(:,1),preNpost(:,2),1);
% heading_sort = preNpost(:,1);
% post_sort = preNpost(:,2);
% post_sort_polyfit = polyval(p_post,heading_sort);
% 
% p=0
% dn=20
% dm=10
% clear dpost_std_now dpost_ste_now post_now heading_now
% for n = 1:dm:(length(heading_sort)-dn)
%     post_sub = post_sort(n:n+dn-1);
%     post_polyfit_sub = post_sort_polyfit(n:n+dn-1);
%     heading_sub = heading_sort(n:n+dn-1);
%     dpost_sub = post_sub - post_polyfit_sub;
% 
%     p=p+1;
%     dpost_std_now(p) = nanstd(dpost_sub);
%     dpost_ste_now(p) = nanstd(dpost_sub)./sqrt(length(dpost_sub));
%     post_now(p) = mean(post_sub);
%     heading_now(p) = mean(heading_sub);
% end
% 
% post_high = post_now + 1.96*dpost_ste_now;
% post_low = post_now - 1.96*dpost_ste_now;
% 
% figure
% plot(0,0,'color',cmap(1,:))
% hold on
% plot(0,0,'color',cmap(end,:))
% legend('slow','fast')
% hold off
% ciplot(post_low',post_high',heading_now','r')
% alpha(.5)
% hold on
% plot(heading_now,post_now,'-k')
% plot([-180,180],[0,0],'--k')
% 
% for i=1:size(heading_pre,1)
%     if settings.expansion.speed(i) == 0
%         plot(heading_pre(i),heading_post(i),'.','color',cmap(1,:))
%         hold on
%     else
%         plot(heading_pre(i),heading_post(i),'.','color',cmap(end,:))
%         hold on
%     end
% end
% % grid on
% axis equal
% xlabel('pre-heading')
% ylabel('post-heading')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'post_vs_preheading.fig')
% saveas(gca,'post_vs_preheading.jpg')
% 
% 
% 
% figure
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
% saveas(gca,'Vhorizontal_PreVsPost.fig')
% saveas(gca,'Vhorizontal_PreVsPost.jpg')
% 
% figure
% % compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
% % hold on
% % compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
% Z = compass(sqrt(u_post.^2 + v_post.^2),w_post);
% for i=1:length(Z)
%     set(Z(i),'color',cmap(end,:))
% end
% hold on
% Y = compass(sqrt(u_pre.^2 + v_pre.^2),w_pre);
% for i=1:length(Y)
%     set(Y(i),'color',cmap(1,:))
% end 
% title('flight speed vectors')
% saveas(gca,'Vvertical_PreVsPost.fig')
% saveas(gca,'Vvertical_PreVsPost.jpg')
% 
% dh = 360/(25)
% heading_noresp = heading(1:Nnoresp,:);
% heading_noresp = heading_noresp(isnan(heading_noresp)==0);
% figure
% hist(heading_noresp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% saveas(gca,'hist_preheading.fig')
% saveas(gca,'hist_preheading.jpg')
% plot2svg
% 
% heading_resp = heading(Nresp:end,:);
% heading_resp = heading_resp(isnan(heading_resp)==0);
% figure
% hist(heading_resp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% saveas(gca,'hist_postheading.fig')
% saveas(gca,'hist_postheading.jpg')
% plot2svg
% 
% V_noresp = V(1:Nnoresp,:);
% V_noresp = V_noresp(isnan(V_noresp)==0);
% figure
% hist(V_noresp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% saveas(gca,'hist_preV.fig')
% saveas(gca,'hist_preV.jpg')
% plot2svg
% 
% V_resp = V(Nresp:end,:);
% V_resp = V_resp(isnan(V_resp)==0);
% figure
% hist(V_resp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% saveas(gca,'hist_postV.fig')
% saveas(gca,'hist_postV.jpg')
% plot2svg
% 
% dA = 30/25
% A_noresp = A(1:Nnoresp,:);
% A_noresp = A_noresp(isnan(A_noresp)==0);
% figure
% hist(A_noresp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% saveas(gca,'hist_preA.fig')
% saveas(gca,'hist_preA.jpg')
% plot2svg
% 
% A_resp = A(Nresp:end,:);
% A_resp = A_resp(isnan(A_resp)==0);
% figure
% hist(A_resp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% saveas(gca,'hist_postA.fig')
% saveas(gca,'hist_postA.jpg')
% plot2svg

% plot top view flight path vectors


figure
colormap(cmap)
caxis([Nstart Nstop])
hold on
for j=1:size(x,2)
    plot(x(Nstart:Nstop,j),y(Nstart:Nstop,j),'-','color',cmap(1,:))
end

vec_col = 0;
for i=Nstart:dN:Nstop
    vec_col = vec_col + size(cmap,1)/Nstop;
    for j=1:size(x,2)
%         quiverc(x(i,:),y(i,:),u(i,:),v(i,:),'color',cmap(1,:))
        if isnan(x(i,j)) == 0
            quivert(x(i,j),y(i,j),u(i,j),v(i,j),i,'as',1/500,'ahr',[1 1],'nt')
        end
    end
    
end
    
axis equal
colorbar

% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
saveas(gca,'flighttracks.fig')
saveas(gca,'flighttracks.jpg')
plot2svg





cd ..