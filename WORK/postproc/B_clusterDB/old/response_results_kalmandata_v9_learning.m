clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_NOq.mat')
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
cmap1 = [.25:-(.25-.75)/(256-1):0.75];
cmap2 = [.25:-(.25-0)/(256-1):0];
cmap3 = [.25:-(.25-0)/(256-1):0];

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

% side_slip = dir - pathDB.yaw_global;
% side_slip(abs(side_slip)>90) = nan;

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);

% response threshold
alpha_dot_resp = 3000;
A_resp = 10;
At_resp = 10;
An_resp = 10;

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

% side_slip_pre = nanmean(side_slip(trigger_frame:trigger_frame+DN_noresp,:))';

V_pre = sqrt(u_pre.^2 + v_pre.^2 + w_pre.^2);
heading_pre = rad2deg(atan2(u_pre,v_pre));
ascend_pre = rad2deg(atan2(w_pre,sqrt(u_pre.^2 + v_pre.^2)));

for i=1:size(u,2)
    
    % vector at end (post response)
    last_frame = find(isnan(u(:,i))==0, 1, 'last' );
    
    z_post(i,1) = mean(z(last_frame-DN_escape:last_frame,i))';
    
    u_post(i,1) = mean(u(last_frame-DN_escape:last_frame,i))';
    v_post(i,1) = mean(v(last_frame-DN_escape:last_frame,i))';
    w_post(i,1) = mean(w(last_frame-DN_escape:last_frame,i))';
end

V_post = sqrt(u_post.^2 + v_post.^2 + w_post.^2);
heading_post = rad2deg(atan2(u_post,v_post));
ascend_post = rad2deg(atan2(w_post,sqrt(u_post.^2 + v_post.^2)));
ascend_turn = ascend_post - ascend_pre;

% turn heading < 180
heading_turn = heading_post - heading_pre;
heading_turn(heading_turn>180) = heading_turn(heading_turn>180) - 360;

% heading_turn2 > 0 if animal turns away from looming object, heading_turn2 < 0: towards
heading_turn2 = abs(heading_turn);
heading_turn2(find(abs(heading_post)) > abs(heading_pre)) = - heading_turn2((find(abs(heading_post)) > abs(heading_pre)));

% turn rate
heading = rad2deg(atan2(u,v));
ascend = rad2deg(atan2(w,sqrt(u.^2 + v.^2)));

turn_rate = nan(size(heading));
for i=2:size(heading,1)-2
 turn_rate(i,:) = heading(i+1,:) - heading(i,:);
end
turn_rate = turn_rate/dt;
% 
% turn_rate(abs(turn_rate)>8000)  = nan;

%% plot track data vs time
figure
subplot(3,1,1)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
plot(0,0,'c')
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t,stim_angle_plot(:,i),'r')
        subplot(3,1,2)
        hold on
        plot(t,V(:,i),'r')
        subplot(3,1,3)
        hold on
        plot(t,A(:,i),'r')
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t,stim_angle_plot(:,i),'g')
        subplot(3,1,2)
        hold on
        plot(t,V(:,i),'g')
        subplot(3,1,3)
        hold on
        plot(t,A(:,i),'g')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(3,1,1)
        hold on
        plot(t,stim_angle_plot(:,i),'b')
        subplot(3,1,2)
        hold on
        plot(t,V(:,i),'b')
        subplot(3,1,3)
        hold on
        plot(t,A(:,i),'b')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(3,1,1)
        hold on
        plot(t,stim_angle_plot(:,i),'c')
        subplot(3,1,2)
        hold on
        plot(t,V(:,i),'c')
        subplot(3,1,3)
        hold on
        plot(t,A(:,i),'c')
    end
end

subplot(3,1,1)
% grid on
xlabel('time (s)')
ylabel('heading (deg)')
axis([t_start t_stop -180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])

subplot(3,1,2)
% grid on
xlabel('time (s)')
ylabel('flight speed (m/s)')
axis([t_start t_stop 0 1])

subplot(3,1,3)
% grid on
xlabel('time (s)')
ylabel('acceleration (m/s^2)')
axis([t_start t_stop 0 30])

saveas(gca,'flightpaths_heading_V_A.fig')
saveas(gca,'flightpaths_heading_V_A.jpg')
% plot2svg(flightpaths_heading_V_A,gca)

%% plot track data vs time (COMPENSATED)
figure
subplot(3,1,1)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
plot(0,0,'c')
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),stim_angle_plot(:,i),'r')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),V(:,i),'r')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),A(:,i),'r')
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),stim_angle_plot(:,i),'g')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),V(:,i),'g')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),A(:,i),'g')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),stim_angle_plot(:,i),'b')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),V(:,i),'b')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),A(:,i),'b')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),stim_angle_plot(:,i),'c')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),V(:,i),'c')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),A(:,i),'c')
    end
end

subplot(3,1,1)
% grid on
xlabel('time (s)')
ylabel('heading (deg)')
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) -180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])

subplot(3,1,2)
% grid on
xlabel('time (s)')
ylabel('flight speed (m/s)')
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) 0 1])

subplot(3,1,3)
% grid on
xlabel('time (s)')
ylabel('acceleration (m/s^2)')
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) 0 30])

saveas(gca,'flightpaths_heading_V_A_compens.fig')
saveas(gca,'flightpaths_heading_V_A_compens.jpg')
% plot2svg(flightpaths_heading_V_A,gca)

%% plot accels vs time
figure
subplot(3,1,1)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
plot(0,0,'c')
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t,A(:,i),'r')
        subplot(3,1,2)
        hold on
        plot(t,An(:,i),'r')
        subplot(3,1,3)
        hold on
        plot(t,At(:,i),'r')
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t,A(:,i),'g')
        subplot(3,1,2)
        hold on
        plot(t,An(:,i),'g')
        subplot(3,1,3)
        hold on
        plot(t,At(:,i),'g')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(3,1,1)
        hold on
        plot(t,A(:,i),'b')
        subplot(3,1,2)
        hold on
        plot(t,An(:,i),'b')
        subplot(3,1,3)
        hold on
        plot(t,At(:,i),'b')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(3,1,1)
        hold on
        plot(t,A(:,i),'c')
        subplot(3,1,2)
        hold on
        plot(t,An(:,i),'c')
        subplot(3,1,3)
        hold on
        plot(t,At(:,i),'c')
    end
end

subplot(3,1,1)
% grid on
xlabel('time (s)')
ylabel('A')
axis([t_start t_stop -180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
axis([t_start t_stop 0 30])

subplot(3,1,2)
% grid on
xlabel('time (s)')
ylabel('An')
axis([t_start t_stop 0 30])

subplot(3,1,3)
% grid on
xlabel('time (s)')
ylabel('At')
axis([t_start t_stop -10 20])

saveas(gca,'flightpaths_As.fig')
saveas(gca,'flightpaths_As.jpg')
% plot2svg(flightpaths_As,gca)

%% plot accels vs time (COMPENSATED)
figure
subplot(3,1,1)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
plot(0,0,'c')
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),A(:,i),'r')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),An(:,i),'r')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),At(:,i),'r')
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),A(:,i),'g')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),An(:,i),'g')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),At(:,i),'g')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),A(:,i),'b')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),An(:,i),'b')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),At(:,i),'b')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(3,1,1)
        hold on
        plot(t-t_Aresp(i),A(:,i),'c')
        subplot(3,1,2)
        hold on
        plot(t-t_Aresp(i),An(:,i),'c')
        subplot(3,1,3)
        hold on
        plot(t-t_Aresp(i),At(:,i),'c')
    end
end

subplot(3,1,1)
% grid on
xlabel('time (s)')
ylabel('A')
% set(gca,'XTick',[-180 -90 0 90 180])
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) 0 30])

subplot(3,1,2)
% grid on
xlabel('time (s)')
ylabel('An')
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) 0 30])

subplot(3,1,3)
% grid on
xlabel('time (s)')
ylabel('At')
axis([t_start-nanmean(t_Aresp) t_stop-nanmean(t_Aresp) -10 20])

saveas(gca,'flightpaths_As_compens.fig')
saveas(gca,'flightpaths_As_compens.jpg')
% plot2svg(flightpaths_As,gca)

%% plot learning variables
figure
subplot(2,2,1)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
plot(0,0,'c')
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(2,2,1)
        hold on
        plot(settings.seq(i,2),V_post(i),'or')
        subplot(2,2,2)
        hold on
        plot(settings.seq(i,2),t_Aresp(i),'or')
        subplot(2,2,3)
        hold on
        plot(settings.seq(i,2),heading_post(i),'or')
        subplot(2,2,4)
        hold on
        plot(settings.seq(i,2),A_max(i),'or')
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(2,2,1)
        hold on
        plot(settings.seq(i,2),V_post(i),'og')
        subplot(2,2,2)
        hold on
        plot(settings.seq(i,2),t_Aresp(i),'og')
        subplot(2,2,3)
        hold on
        plot(settings.seq(i,2),heading_post(i),'og')
        subplot(2,2,4)
        hold on
        plot(settings.seq(i,2),A_max(i),'og')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(2,2,1)
        hold on
        plot(settings.seq(i,2),V_post(i),'ob')
        subplot(2,2,2)
        hold on
        plot(settings.seq(i,2),t_Aresp(i),'ob')
        subplot(2,2,3)
        hold on
        plot(settings.seq(i,2),heading_post(i),'ob')
        subplot(2,2,4)
        hold on
        plot(settings.seq(i,2),A_max(i),'ob')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(2,2,1)
        hold on
        plot(settings.seq(i,2),V_post(i),'oc')
        subplot(2,2,2)
        hold on
        plot(settings.seq(i,2),t_Aresp(i),'oc')
        subplot(2,2,3)
        hold on
        plot(settings.seq(i,2),heading_post(i),'oc')
        subplot(2,2,4)
        hold on
        plot(settings.seq(i,2),A_max(i),'oc')
    end
end

subplot(2,2,1)
nonan = isnan(settings.seq(:,2) + V_post);
[p,S] = polyfit(settings.seq(nonan==0,2),V_post(nonan==0),1);
plot([min(settings.seq(:,2)),max(settings.seq(:,2))],polyval(p,[min(settings.seq(:,2)),max(settings.seq(:,2))]),'k-')

subplot(2,2,2)
nonan = isnan(settings.seq(:,2) + t_Aresp');
[p,S] = polyfit(settings.seq(nonan==0,2),t_Aresp(nonan==0)',1);
plot([min(settings.seq(:,2)),max(settings.seq(:,2))],polyval(p,[min(settings.seq(:,2)),max(settings.seq(:,2))]),'k-')

subplot(2,2,3)
nonan = isnan(settings.seq(:,2) + heading_post);
[p,S] = polyfit(settings.seq(nonan==0,2),heading_post(nonan==0),1);
plot([min(settings.seq(:,2)),max(settings.seq(:,2))],polyval(p,[min(settings.seq(:,2)),max(settings.seq(:,2))]),'k-')

subplot(2,2,4)
nonan = isnan(settings.seq(:,2) + A_max');
[p,S] = polyfit(settings.seq(nonan==0,2),A_max(nonan==0)',1);
plot([min(settings.seq(:,2)),max(settings.seq(:,2))],polyval(p,[min(settings.seq(:,2)),max(settings.seq(:,2))]),'k-')

subplot(2,2,1)
% grid on
xlabel('trigger number')
ylabel('response speed')
% axis([t_start t_stop -180 180])
% % set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])

subplot(2,2,2)
% grid on
xlabel('trigger number')
ylabel('response time')
% axis([t_start t_stop 0 1])

subplot(2,2,3)
% grid on
xlabel('trigger number')
ylabel('response heading')
% axis([t_start t_stop 0 30])

subplot(2,2,4)
% grid on
xlabel('trigger number')
ylabel('response A max')
% axis([t_start t_stop 0 30])

saveas(gca,'learning metrics.fig')
saveas(gca,'learning metrics.jpg')
% plot2svg(flightpaths_As,gca)

%% plot linear vs transverse reaction time

% 1st order polyfit
t_AtNAn = [t_Atresp' t_Anresp'];
t_AtNAn = sortrows(t_AtNAn);

nonan = isnan(t_AtNAn(:,1)) + isnan(t_AtNAn(:,2));
t_AtNAn = t_AtNAn(nonan==0,:);

[p_t,S_t] = polyfit(t_AtNAn(:,1),t_AtNAn(:,2),1);
t_At_sort = t_AtNAn(:,1);
t_An_sort = t_AtNAn(:,2);
t_An_sort_polyfit = polyval(p_t,t_At_sort);

p=0
dn=20
dm=10
clear dt_An_std_now dt_An_ste_now t_An_now t_At_now
for n = 1:dm:(length(t_At_sort)-dn)
    t_An_sub = t_An_sort(n:n+dn-1);
    t_An_polyfit_sub = t_An_sort_polyfit(n:n+dn-1);
    t_At_sub = t_At_sort(n:n+dn-1);
    dt_An_sub = t_An_sub - t_An_polyfit_sub;

    p=p+1;
    dt_An_std_now(p) = nanstd(dt_An_sub);
    dt_An_ste_now(p) = nanstd(dt_An_sub)./sqrt(length(dt_An_sub));
    t_An_now(p) = nanmean(t_An_sub);
    t_At_now(p) = nanmean(t_At_sub);
end

t_An_high = t_An_now + 1.96*dt_An_ste_now;
t_An_low = t_An_now - 1.96*dt_An_ste_now;

figure
plot(0,0,'.r','markersize',25)
hold on
plot(0,0,'.g','markersize',25)
plot(0,0,'.b','markersize',25)
plot(0,0,'.c','markersize',25)
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

% hold off
ciplot(t_An_low',t_An_high',t_At_now',cmap(1,:))
hold on
alpha(.25)
plot(t_At_now,t_An_now,'-k')
plot([0,max(t_Atresp)],[0,max(t_Atresp)],'--k')

for i=1:size(t_Anresp,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        plot(t_Atresp(i),t_Anresp(i),'.r','markersize',25)
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        plot(t_Atresp(i),t_Anresp(i),'.g','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        plot(t_Atresp(i),t_Anresp(i),'.b','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        plot(t_Atresp(i),t_Anresp(i),'.c','markersize',25)
    end
end
% grid on
axis equal
xlabel('At')
ylabel('An')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'At_vs_An_response_time.fig')
saveas(gca,'At_vs_An_response_time.jpg')


%% turn vs heading
heading_pre_calc = heading_pre;
heading_turn_calc = heading_turn;

for i=1:size(heading_turn,1)
    if heading_pre(i) < -90 && heading_turn(i) <-90
        heading_turn_calc(i) = heading_turn_calc(i) + 360;
    end
end

heading_pre_ext = heading_pre_calc;
heading_turn_ext = heading_turn_calc;

% wrap around
heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc -360;
heading_turn_ext(end+1:end+length(heading_turn_calc)) = heading_turn_calc +360;

heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc +360;
heading_turn_ext(end+1:end+length(heading_turn_calc)) = heading_turn_calc -360;

% 1st order polyfit
headingNturn = [heading_pre_ext heading_turn_ext];
headingNturn = sortrows(headingNturn);

[p_turn,S_turn] = polyfit(headingNturn(:,1),headingNturn(:,2),1);
heading_sort = headingNturn(:,1);
turn_sort = headingNturn(:,2);
turn_sort_polyfit = polyval(p_turn,heading_sort);

p=0
dn=20
dm=10
clear dturn_std_now dturn_ste_now turn_now heading_now
for n = 1:dm:(length(heading_sort)-dn)
    turn_sub = turn_sort(n:n+dn-1);
    turn_polyfit_sub = turn_sort_polyfit(n:n+dn-1);
    heading_sub = heading_sort(n:n+dn-1);
    dturn_sub = turn_sub - turn_polyfit_sub;

    p=p+1;
    dturn_std_now(p) = nanstd(dturn_sub);
    dturn_ste_now(p) = nanstd(dturn_sub)./sqrt(length(dturn_sub));
    turn_now(p) = mean(turn_sub);
    heading_now(p) = mean(heading_sub);
end

turn_high = turn_now + 1.96*dturn_ste_now;
turn_low = turn_now - 1.96*dturn_ste_now;

for i=1:size(heading_turn,1)
    if heading_pre(i) > 90 && heading_turn(i) <-180
        heading_turn(i) = heading_turn(i) + 360;
    end
end

figure
plot(0,0,'.r','markersize',25)
hold on
plot(0,0,'.g','markersize',25)
plot(0,0,'.b','markersize',25)
plot(0,0,'.c','markersize',25)
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

% hold off
ciplot(turn_low',turn_high',heading_now',cmap(1,:))
hold on
alpha(.25)
plot(heading_now,turn_now,'-k')
plot([-180,180],[180,-180],'--k')

for i=1:size(heading_pre,1)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        plot(heading_pre(i),heading_turn(i),'.r','markersize',25)
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        plot(heading_pre(i),heading_turn(i),'.g','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        plot(heading_pre(i),heading_turn(i),'.b','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        plot(heading_pre(i),heading_turn(i),'.c','markersize',25)
    end
end
% grid on
axis equal
xlabel('pre-heading')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'turnangle_vs_preheading.fig')
saveas(gca,'turnangle_vs_preheading.jpg')

%% plot heading pre vs post

heading_pre_calc = heading_pre;
heading_post_calc = heading_post;

heading_pre_ext = heading_pre_calc;
heading_post_ext = heading_post_calc;

% wrap around
heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc -360;
heading_post_ext(end+1:end+length(heading_post_calc)) = heading_post_calc;

heading_pre_ext(end+1:end+length(heading_pre_calc)) = heading_pre_calc +360;
heading_post_ext(end+1:end+length(heading_post_calc)) = heading_post_calc;

% 1st order polyfit
preNpost = [heading_pre_calc heading_post_calc];
preNpost = sortrows(preNpost);

[p_post,S_post] = polyfit(preNpost(:,1),preNpost(:,2),1);
heading_sort = preNpost(:,1);
post_sort = preNpost(:,2);
post_sort_polyfit = polyval(p_post,heading_sort);

p=0
dn=20
dm=10
clear dpost_std_now dpost_ste_now post_now heading_now
for n = 1:dm:(length(heading_sort)-dn)
    post_sub = post_sort(n:n+dn-1);
    post_polyfit_sub = post_sort_polyfit(n:n+dn-1);
    heading_sub = heading_sort(n:n+dn-1);
    dpost_sub = post_sub - post_polyfit_sub;

    p=p+1;
    dpost_std_now(p) = nanstd(dpost_sub);
    dpost_ste_now(p) = nanstd(dpost_sub)./sqrt(length(dpost_sub));
    post_now(p) = mean(post_sub);
    heading_now(p) = mean(heading_sub);
end

post_high = post_now + 1.96*dpost_ste_now;
post_low = post_now - 1.96*dpost_ste_now;

figure
plot(0,0,'.r','markersize',25)
hold on
plot(0,0,'.g','markersize',25)
plot(0,0,'.b','markersize',25)
plot(0,0,'.c','markersize',25)
legend('stepwise slow 165deg','stepwise fast 165deg','continuous fast 165deg','continuous fast 64deg')

% hold off
ciplot(post_low',post_high',heading_now',cmap(1,:))
hold on
alpha(.25)
plot(heading_now,post_now,'-k')
plot([-180,180],[0,0],'--k')

for i=1:size(heading_pre,1)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        plot(heading_pre(i),heading_post(i),'.r','markersize',25)
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        plot(heading_pre(i),heading_post(i),'.g','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        plot(heading_pre(i),heading_post(i),'.b','markersize',25)
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        plot(heading_pre(i),heading_post(i),'.c','markersize',25)
    end
end

% grid on
axis equal
xlabel('pre-heading')
ylabel('post-heading')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])
saveas(gca,'post_vs_preheading.fig')
saveas(gca,'post_vs_preheading.jpg')


%% plot pre&post velocity vectors

figure
subplot(1,2,1)
% plot(0,0)
% hold on
% plot(0,0,'r')
% legend('pre','post')
% compass_zeroup(u_post,v_post,'r');
% hold on
% compass_zeroup(u_pre,v_pre,'b');
Z = compass(u_post,v_post);
for i=1:length(Z)
    set(Z(i),'color',cmap(end,:))
end
hold on
Y = compass(u_pre,v_pre);
for i=1:length(Y)
    set(Y(i),'color',cmap(1,:))
end 
title('horizontal flight speed vectors')
% saveas(gca,'Vhorizontal_PreVsPost.fig')
% saveas(gca,'Vhorizontal_PreVsPost.jpg')

%% plot pre&post vertical velocity vectors
% figure
% compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
% hold on
% compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
subplot(1,2,2)
Z = compass(sqrt(u_post.^2 + v_post.^2),w_post);
for i=1:length(Z)
    set(Z(i),'color',cmap(end,:))
end
hold on
Y = compass(sqrt(u_pre.^2 + v_pre.^2),w_pre);
for i=1:length(Y)
    set(Y(i),'color',cmap(1,:))
end 
title('vertical flight speed vectors')
saveas(gca,'V_PreVsPost.fig')
saveas(gca,'V_PreVsPost.jpg')

%% plot histograms
% figure
subplot(3,2,1)
dh = 360/(25)
heading_noresp = heading(1:Nnoresp,:);
heading_noresp = heading_noresp(isnan(heading_noresp)==0);
hist(heading_noresp,-180+dh/2:dh:180-dh/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('heading pre')
set(gca,'ylim',[0 10000])

set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% saveas(gca,'hist_preheading.fig')
% saveas(gca,'hist_preheading.jpg')
% plot2svg

% figure
subplot(3,2,2)
heading_resp = heading(Nresp:end,:);
heading_resp = heading_resp(isnan(heading_resp)==0);
hist(heading_resp,-180+dh/2:dh:180-dh/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('heading post')
% set(gca,'ylim',[0 10000])

set(gca,'xlim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
% saveas(gca,'hist_postheading.fig')
% saveas(gca,'hist_postheading.jpg')
% plot2svg

% figure
subplot(3,2,3)
dV = 1/25
V_noresp = V(1:Nnoresp,:);
V_noresp = V_noresp(isnan(V_noresp)==0);
hist(V_noresp,0+dV/2:dV:1-dV/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('velocity pre')
% set(gca,'ylim',[0 50000])

set(gca,'xlim',[0 1])
set(gca,'XTick',[0 .5 1])
% saveas(gca,'hist_preV.fig')
% saveas(gca,'hist_preV.jpg')
% plot2svg

% figure
subplot(3,2,4)
V_resp = V(Nresp:end,:);
V_resp = V_resp(isnan(V_resp)==0);
hist(V_resp,0+dV/2:dV:1-dV/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('velocity post')
% set(gca,'ylim',[0 10000])

set(gca,'xlim',[0 1])
set(gca,'XTick',[0 .5 1])
% saveas(gca,'hist_postV.fig')
% saveas(gca,'hist_postV.jpg')
% plot2svg

% figure
subplot(3,2,5)
dA = 30/25
A_noresp = A(1:Nnoresp,:);
A_noresp = A_noresp(isnan(A_noresp)==0);
hist(A_noresp,0+dA/2:dA:30-dA/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('accel pre')
% set(gca,'ylim',[0 100000])

set(gca,'xlim',[0 30])
set(gca,'XTick',[0 10 20 30])
% saveas(gca,'hist_preA.fig')
% saveas(gca,'hist_preA.jpg')
% plot2svg

% figure
subplot(3,2,6)
A_resp = A(Nresp:end,:);
A_resp = A_resp(isnan(A_resp)==0);
hist(A_resp,0+dA/2:dA:30-dA/2)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',cmap(1,:),'EdgeColor','w')

title('accel post')
% set(gca,'ylim',[0 10000])

set(gca,'xlim',[0 30])
set(gca,'XTick',[0 10 20 30])
% saveas(gca,'hist_postA.fig')
% saveas(gca,'hist_postA.jpg')
% plot2svg

saveas(gca,'histograms_pre_post.fig')
saveas(gca,'histograms_pre_post.jpg')
% plot2svg

%% plot top view flight path vectors
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

%% escape performance: heading, U, A, time
k=0;
l=0;
m=0;
n=0;

figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on
for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(2,2,1)
        plot(.5,heading_post(i),'ko','MarkerSize',5)
        subplot(2,2,2)
        plot(.5,V_post(i),'ko','MarkerSize',5)
        subplot(2,2,3)
        plot(.5,A_max(i),'ko','MarkerSize',5)
        subplot(2,2,4)
        plot(.5,t_Aresp(i),'ko','MarkerSize',5)
        
        k = k+1;
        H1(k) = heading_post(i);
        V1(k) = V_post(i);
        A1(k) = A_max(i);
        T1(k) = t_Aresp(i);
        
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(2,2,1)
        plot(1,heading_post(i),'ko','MarkerSize',5)
        subplot(2,2,2)
        plot(1,V_post(i),'ko','MarkerSize',5)
        subplot(2,2,3)
        plot(1,A_max(i),'ko','MarkerSize',5)
        subplot(2,2,4)
        plot(1,t_Aresp(i),'ko','MarkerSize',5)
        
        l = l+1;
        H2(l) = heading_post(i);
        V2(l) = V_post(i);
        A2(l) = A_max(i);
        T2(l) = t_Aresp(i);
        
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 165
        subplot(2,2,1)
        plot(1.5,heading_post(i),'ko','MarkerSize',5)
        subplot(2,2,2)
        plot(1.5,V_post(i),'ko','MarkerSize',5)
        subplot(2,2,3)
        plot(1.5,A_max(i),'ko','MarkerSize',5)
        subplot(2,2,4)
        plot(1.5,t_Aresp(i),'ko','MarkerSize',5)
        
        m = m+1;
        H3(m) = heading_post(i);
        V3(m) = V_post(i);
        A3(m) = A_max(i);
        T3(m) = t_Aresp(i);
        
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1 && settings.expansion.maxangle(i) == 64
        subplot(2,2,1)
        plot(2,heading_post(i),'ko','MarkerSize',5)
        subplot(2,2,2)
        plot(2,V_post(i),'ko','MarkerSize',5)
        subplot(2,2,3)
        plot(2,A_max(i),'ko','MarkerSize',5)
        subplot(2,2,4)
        plot(2,t_Aresp(i),'ko','MarkerSize',5)
        
        n = n+1;
        H4(n) = heading_post(i);
        V4(n) = V_post(i);
        A4(n) = A_max(i);
        T4(n) = t_Aresp(i);
    end
end

subplot(2,2,1)
errorbar(.5,nanmean(H1),nanstd(H1),'dr','linewidth',2,'MarkerSize',10)
errorbar(1,nanmean(H2),nanstd(H2),'dg','linewidth',2,'MarkerSize',10)
errorbar(1.5,nanmean(H3),nanstd(H3),'db','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(H4),nanstd(H4),'dc','linewidth',2,'MarkerSize',10)

subplot(2,2,2)
errorbar(.5,nanmean(V1),nanstd(V1),'dr','linewidth',2,'MarkerSize',10)
errorbar(1,nanmean(V2),nanstd(V2),'db','linewidth',2,'MarkerSize',10)
errorbar(1.5,nanmean(V3),nanstd(V3),'dg','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(V4),nanstd(V4),'dc','linewidth',2,'MarkerSize',10)

subplot(2,2,3)
errorbar(.5,nanmean(A1),nanstd(A1),'dr','linewidth',2,'MarkerSize',10)
errorbar(1,nanmean(A2),nanstd(A2),'db','linewidth',2,'MarkerSize',10)
errorbar(1.5,nanmean(A3),nanstd(A3),'dg','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(A4),nanstd(A4),'dc','linewidth',2,'MarkerSize',10)

subplot(2,2,4)
errorbar(.5,nanmean(T1),nanstd(T1),'dr','linewidth',2,'MarkerSize',10)
errorbar(1,nanmean(T2),nanstd(T2),'db','linewidth',2,'MarkerSize',10)
errorbar(1.5,nanmean(T3),nanstd(T3),'dg','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(T4),nanstd(T4),'dc','linewidth',2,'MarkerSize',10)

subplot(2,2,1)
ylabel('escape heading')
set(gca,'XTick',[.5:.5:2])
set(gca,'XTickLabel',{'step slow 165deg','step fast 165deg','cont fast 165deg','cont fast 64deg'})
set(gca,'YTick',[-90:30:90])
axis([0 2.5 -90 90])

subplot(2,2,2)
ylabel('escape speed')
set(gca,'XTick',[.5:.5:2])
set(gca,'XTickLabel',{'step slow 165deg','step fast 165deg','cont fast 165deg','cont fast 64deg'})
set(gca,'YTick',[0:.5:1])
axis([0 2.5 0 1])

subplot(2,2,3)
ylabel('accel max')
set(gca,'XTick',[.5:.5:2])
set(gca,'XTickLabel',{'step slow 165deg','step fast 165deg','cont fast 165deg','cont fast 64deg'})
set(gca,'YTick',[0:5:30])
axis([0 2.5 0 30])

subplot(2,2,4)
ylabel('response time')
set(gca,'XTick',[.5:.5:2])
set(gca,'XTickLabel',{'step slow 165deg','step fast 165deg','cont fast 165deg','cont fast 64deg'})
set(gca,'YTick',[0:.05:.15])
axis([0 2.5 0 .15])

% title('acceleration based response time')
saveas(gca,'escape performance.fig')
saveas(gca,'escape performance.jpg')
% plot2svg

cd ..