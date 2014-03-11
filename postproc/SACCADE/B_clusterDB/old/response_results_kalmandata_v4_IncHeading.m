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

% temporal data
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

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

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')

for i=1:size(u,2)
    if settings.expansion.speed(i) == 0
        plot(t,stim_angle_plot(:,i))
    elseif settings.expansion.HorPos(i) == 0
        plot(t,stim_angle_plot(:,i),'r')
    elseif settings.expansion.HorPos(i) == 180
        plot(t,stim_angle_plot(:,i),'g')
    end
end
grid on
xlabel('t')
ylabel('heading')
axis([-0.2 0.2 -180 180])

saveas(gca,'temporal_stim_angle.fig')
saveas(gca,'temporal_stim_angle.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('left slip','right slip','slip too high')

for i=1:size(u,2)
    if side_slip_pre(i) < 0
        plot(t,stim_angle_plot(:,i))
    elseif side_slip_pre(i) > 0
        plot(t,stim_angle_plot(:,i),'r')
    else
        plot(t,stim_angle_plot(:,i),'g')
    end
end
grid on
xlabel('t')
ylabel('heading')
axis([-0.2 0.2 -180 180])

saveas(gca,'temporal_stim_angle_sideslip.fig')
saveas(gca,'temporal_stim_angle_sideslip.jpg')

figure
subplot(3,1,1)
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')

for i=1:size(u,2)
    if settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i))

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i))

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i))
        
    elseif settings.expansion.HorPos(i) == 0
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i),'r')

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i),'r')

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i),'r')
        
    elseif settings.expansion.HorPos(i) == 180
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i),'g')

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i),'g')

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i),'g')
    end
end
subplot(3,1,1)
grid on
xlabel('t')
ylabel('V')
axis([0 0.25 0 1])

subplot(3,1,2)
grid on
xlabel('t')
ylabel('A')
axis([0 0.25 0 30])

subplot(3,1,3)
grid on
xlabel('t')
ylabel('angular velocity')
axis([0 0.25 0 8000])

saveas(gca,'flightpaths_turn_V_A.fig')
saveas(gca,'flightpaths_turn_V_A.jpg')

figure
subplot(3,1,1)
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
plot(0,0,'c')
legend('stepwise slow','stepwise fast','continuous slow','continuous fast')

for i=1:size(u,2)
    if settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i))

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i))

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i))
        
    elseif settings.expansion.stepwise(i) == 1 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i),'r')

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i),'r')

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i),'r')
        
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 0
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i),'g')

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i),'g')

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i),'g')
        
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 1
        subplot(3,1,1)
        hold on
        plot(t(:),V(:,i),'c')

        subplot(3,1,2)
        hold on
        plot(t(:),A(:,i),'c')

        subplot(3,1,3)
        hold on
        plot(t(:),alpha_dot(:,i),'c')
    end
end
subplot(3,1,1)
grid on
xlabel('t')
ylabel('V')
axis([0 0.25 0 1])

subplot(3,1,2)
grid on
xlabel('t')
ylabel('A')
axis([0 0.25 0 30])

subplot(3,1,3)
grid on
xlabel('t')
ylabel('angular velocity')
axis([0 0.25 0 8000])

saveas(gca,'flightpaths_turn_V_A_stepwise.fig')
saveas(gca,'flightpaths_turn_V_A_stepwise.jpg')

figure
% subplot(2,2,1)
% plot(0,0)
% hold on
% plot(0,0,'r')
% plot(0,0,'g')
% plot(0,0,'c')
% legend('A','At','An')

for i=1:size(u,2)
    if settings.expansion.speed(i) == 0
        subplot(2,2,1)
        hold on
        plot(t(:),A(:,i))

        subplot(2,2,2)
        hold on
        plot(t(:),At(:,i))

        subplot(2,2,3)
        hold on
        plot(t(:),An(:,i))
        
        subplot(2,2,4)
        hold on
        plot(t(:),alpha_dot(:,i))
        
    elseif settings.expansion.HorPos(i) == 0
        subplot(2,2,1)
        hold on
        plot(t(:),A(:,i),'r')

        subplot(2,2,2)
        hold on
        plot(t(:),At(:,i),'r')

        subplot(2,2,3)
        hold on
        plot(t(:),An(:,i),'r')

        subplot(2,2,4)
        hold on
        plot(t(:),alpha_dot(:,i),'r')

    elseif settings.expansion.HorPos(i) == 180
        subplot(2,2,1)
        hold on
        plot(t(:),A(:,i),'g')

        subplot(2,2,2)
        hold on
        plot(t(:),At(:,i),'g')

        subplot(2,2,3)
        hold on
        plot(t(:),An(:,i),'g')

        subplot(2,2,4)
        hold on
        plot(t(:),alpha_dot(:,i),'g')
    end
end
subplot(2,2,1)
grid on
xlabel('t')
ylabel('A')
axis([0 0.25 0 30])

subplot(2,2,2)
grid on
xlabel('t')
ylabel('At')
axis([0 0.25 -10 20])

subplot(2,2,3)
grid on
xlabel('t')
ylabel('An')
axis([0 0.25 0 30])

subplot(2,2,4)
grid on
xlabel('t')
ylabel('angular velocity')
axis([0 0.25 0 8000])

saveas(gca,'flightpaths_accel.fig')
saveas(gca,'flightpaths_accel.jpg')

subplot(2,2,1)
plot([0 1],[A_resp A_resp],'--k')

subplot(2,2,2)
plot([0 1],[At_resp At_resp],'--k')

subplot(2,2,3)
plot([0 1],[An_resp An_resp],'--k')

subplot(2,2,4)
plot([0 1],[alpha_dot_resp alpha_dot_resp],'--k')

saveas(gca,'flightpaths_accel_incResponse.fig')
saveas(gca,'flightpaths_accel_incResponse.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(heading_pre,1)
    if settings.expansion.speed(i) == 0
        plot(side_slip_pre(i),heading_post(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(side_slip_pre(i),heading_post(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(side_slip_pre(i),heading_post(i),'.g')
    end
end
plot([-180,180],[0,0],'--k')
grid on
axis equal
xlabel('sideslip pre')
ylabel('heading post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'side_slip_pre_vs_heading_post.fig')
saveas(gca,'side_slip_pre_vs_heading_post.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(heading_pre,1)
    if settings.expansion.speed(i) == 0
        plot(side_slip_pre(i),heading_turn(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(side_slip_pre(i),heading_turn(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(side_slip_pre(i),heading_turn(i),'.g')
    end
end
plot([-180,180],[0,0],'--k')
grid on
axis equal
xlabel('sideslip pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'side_slip_pre_vs_turn.fig')
saveas(gca,'side_slip_pre_vs_turn.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(heading_pre,1)
    if settings.expansion.speed(i) == 0
        plot(heading_pre(i),heading_post(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(heading_pre(i),heading_post(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(heading_pre(i),heading_post(i),'.g')
    end
end
plot([-180,180],[0,0],'--k')
grid on
axis equal
xlabel('heading pre')
ylabel('heading post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'heading_pre_vs_post.fig')
saveas(gca,'heading_pre_vs_post.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('left slip','right slip','slip too high')
for i=1:size(heading_pre,1)
    if side_slip_pre(i) < 0
        plot(heading_pre(i),heading_post(i),'.')
    elseif side_slip_pre(i) > 0
        plot(heading_pre(i),heading_post(i),'.r')
    else
        plot(heading_pre(i),heading_post(i),'.g')
    end
end
plot([-180,180],[0,0],'--k')
grid on
axis equal
xlabel('heading pre')
ylabel('heading post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'heading_pre_vs_post_sideslip.fig')
saveas(gca,'heading_pre_vs_post_sideslip.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(heading_pre,1)
    if settings.expansion.speed(i) == 0
        plot(heading_pre(i),heading_turn(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(heading_pre(i),heading_turn(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(heading_pre(i),heading_turn(i),'.g')
    end
end
plot([-180,180],[180,-180],'--k')
grid on
axis equal
xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle_vs_preheading.fig')
saveas(gca,'turnangle_vs_preheading.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
legend('stepwise','continuous')
for i=1:size(heading_pre,1)
    if settings.expansion.stepwise(i) == 1
        plot(heading_pre(i),heading_turn(i),'.')
    elseif settings.expansion.stepwise(i) == 0
        plot(heading_pre(i),heading_turn(i),'.r')
    end
end
plot([-180,180],[180,-180],'--k')
grid on
axis equal
xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle_vs_preheading_stepwise.fig')
saveas(gca,'turnangle_vs_preheading_stepwise.jpg')


figure
plot(0,0)
hold on
plot(0,0,'g')
legend('left slip','right slip','slip too high')
for i=1:size(heading_pre,1)
    if side_slip_pre(i) < 0
        plot(heading_pre(i),heading_turn(i),'.')
    elseif side_slip_pre(i) > 0
        plot(heading_pre(i),heading_turn(i),'.r')
    else
        plot(heading_pre(i),heading_turn(i),'.g')
    end
end
plot([-180,180],[180,-180],'--k')
grid on
axis equal
xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle_vs_preheading_sideslip.fig')
saveas(gca,'turnangle_vs_preheading_sideslip.jpg')


figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(ascend_pre,1)
    if settings.expansion.speed(i) == 0
        plot(ascend_pre(i),ascend_post(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(ascend_pre(i),ascend_post(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(ascend_pre(i),ascend_post(i),'.g')
    end
end
plot([-180,180],[-180,180],'--k')
grid on
axis equal
xlabel('ascend pre')
ylabel('ascend post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'ascend_pre_vs_post.fig')
saveas(gca,'ascend_pre_vs_post.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(ascend_pre,1)
    if settings.expansion.speed(i) == 0
        plot(ascend_pre(i),ascend_turn(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(ascend_pre(i),ascend_turn(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(ascend_pre(i),ascend_turn(i),'.g')
    end
end
plot([-180,180],[180,-180],'--k')
grid on
axis equal
xlabel('ascend pre')
ylabel('ascend change')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'ascendchange_vs_ascend_pre.fig')
saveas(gca,'ascendchange_vs_ascend_pre.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(w_pre,1)
    if settings.expansion.speed(i) == 0
        plot(w_pre(i),w_post(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(w_pre(i),w_post(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(w_pre(i),w_post(i),'.g')
    end
end
grid on
axis equal
xlabel('Vz_p_r_e')
ylabel('Vz_p_o_s_t')
saveas(gca,'Vz_preVSpost.fig')
saveas(gca,'Vz_preVSpost.jpg')

figure
plot(0,0)
hold on
plot(0,0,'r')
plot(0,0,'g')
legend('slow','fast','reverse')
for i=1:size(V_pre,1)
    if settings.expansion.speed(i) == 0
        plot(V_pre(i),V_post(i),'.')
    elseif settings.expansion.HorPos(i) == 0
        plot(V_pre(i),V_post(i),'.r')
    elseif settings.expansion.HorPos(i) == 180
        plot(V_pre(i),V_post(i),'.g')
    end
end
grid on
axis equal
xlabel('V_p_r_e')
ylabel('V_p_o_s_t')
saveas(gca,'Vpre_vs_Vpost.fig')
saveas(gca,'Vpre_vs_Vpost.jpg')

figure
% plot(0,0)
% hold on
% plot(0,0,'r')
% legend('pre','post')
compass_zeroup(u_post,v_post,'r')
hold on
compass(u_pre,v_pre,'b')
title('horizontal flight speed vectors')
saveas(gca,'Vhorizontal_PreVsPost.fig')
saveas(gca,'Vhorizontal_PreVsPost.jpg')

figure
compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
hold on
compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
title('flight speed vectors')
saveas(gca,'Vvertical_PreVsPost.fig')
saveas(gca,'Vvertical_PreVsPost.jpg')

%
figure
hold on
errorbar(-53,mean(ascend_post(settings.expansion.VerPos==-1)),std(ascend_post(settings.expansion.VerPos==-1)),'db','linewidth',2,'MarkerSize',6,'MarkerFaceColor','b')
errorbar(0,mean(ascend_post(settings.expansion.VerPos==0)),std(ascend_post(settings.expansion.VerPos==0)),'dr','linewidth',2,'MarkerSize',6,'MarkerFaceColor','r')
errorbar(53,mean(ascend_post(settings.expansion.VerPos==1)),std(ascend_post(settings.expansion.VerPos==1)),'dg','linewidth',2,'MarkerSize',6,'MarkerFaceColor','g')
plot(53*settings.expansion.VerPos,ascend_post,'ok')

grid on
xlabel('vertical looming angle')
ylabel('ascend angle')
% title('adcend angle versus vertical looming angle')
set(gca,'xlim',[-60 60])
saveas(gca,'adcend angle versus vertical looming angle.fig')
saveas(gca,'adcend angle versus vertical looming angle.jpg')

figure
hold on
errorbar(-53,mean(z_post(settings.expansion.VerPos==-1)),std(z_post(settings.expansion.VerPos==-1)),'db','linewidth',2,'MarkerSize',6,'MarkerFaceColor','b')
errorbar(0,mean(z_post(settings.expansion.VerPos==0)),std(z_post(settings.expansion.VerPos==0)),'dr','linewidth',2,'MarkerSize',6,'MarkerFaceColor','r')
errorbar(53,mean(z_post(settings.expansion.VerPos==1)),std(z_post(settings.expansion.VerPos==1)),'dg','linewidth',2,'MarkerSize',6,'MarkerFaceColor','g')
plot(53*settings.expansion.VerPos,z_post,'ok')

grid on
xlabel('vertical looming angle')
ylabel('escape height')
% title('escape height versus vertical looming angle')
set(gca,'xlim',[-60 60])
saveas(gca,'escape height versus vertical looming angle.fig')
saveas(gca,'escape height versus vertical looming angle.jpg')

% 
figure
hold on
errorbar(1,mean(V_post(settings.expansion.speed==1)),std(V_post(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,mean(V_post(settings.expansion.speed==0)),std(V_post(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,V_post,'ko')

grid on
ylabel('escape speed')

set(gca,'XTick',[0 1])
set(gca,'XTickLabel',{'slow','fast'})

% title('escape speed versus looming speed')
set(gca,'xlim',[-.5 1.5])
saveas(gca,'escape speed versus looming speed.fig')
saveas(gca,'escape speed versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Aresp(settings.expansion.speed==1)),nanstd(t_Aresp(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_Aresp(settings.expansion.speed==0)),nanstd(t_Aresp(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_Aresp,'ko')

grid on
ylabel('response time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('acceleration based response time')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'response time (A) versus looming speed.fig')
saveas(gca,'response time (A) versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Amax(settings.expansion.speed==1)),nanstd(t_Amax(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_Amax(settings.expansion.speed==0)),nanstd(t_Amax(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_Amax,'ko')

grid on
ylabel('Amax time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('time of maximum acceleration')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'time of maximum acceleration versus looming speed.fig')
saveas(gca,'time of maximum acceleration versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Atresp(settings.expansion.speed==1)),nanstd(t_Atresp(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_Atresp(settings.expansion.speed==0)),nanstd(t_Atresp(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_Atresp,'ko')

grid on
ylabel('response time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('At based response time')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'response time (At) versus looming speed.fig')
saveas(gca,'response time (At) versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Atmax(settings.expansion.speed==1)),nanstd(t_Atmax(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_Atmax(settings.expansion.speed==0)),nanstd(t_Atmax(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_Atmax,'ko')

grid on
ylabel('Atmax time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('time of maximum At')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'time of maximum At versus looming speed.fig')
saveas(gca,'time of maximum At versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Anresp(settings.expansion.speed==1)),nanstd(t_Anresp(settings.expansion.speed==1)),'b','linewidth',2)
errorbar(0,nanmean(t_Anresp(settings.expansion.speed==0)),nanstd(t_Anresp(settings.expansion.speed==0)),'r','linewidth',2)

plot(1,nanmean(t_Anresp(settings.expansion.speed==1)),'db','MarkerSize',10,'MarkerFaceColor','b')
plot(0,nanmean(t_Anresp(settings.expansion.speed==0)),'dr','MarkerSize',10,'MarkerFaceColor','r')

plot(settings.expansion.speed,t_Anresp,'ko')

grid on
ylabel('response time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('An based response time')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'response time (An) versus looming speed.fig')
saveas(gca,'response time (An) versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_Anmax(settings.expansion.speed==1)),nanstd(t_Anmax(settings.expansion.speed==1)),'b','linewidth',2)
errorbar(0,nanmean(t_Anmax(settings.expansion.speed==0)),nanstd(t_Anmax(settings.expansion.speed==0)),'r','linewidth',2)

plot(1,nanmean(t_Anmax(settings.expansion.speed==1)),'db','MarkerSize',10,'MarkerFaceColor','b')
plot(0,nanmean(t_Anmax(settings.expansion.speed==0)),'dr','MarkerSize',10,'MarkerFaceColor','r')

plot(settings.expansion.speed,t_Anmax,'ko')

grid on
ylabel('Anmax time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('time of maximum An')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'time of maximum An versus looming speed.fig')
saveas(gca,'time of maximum An versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_alpha_dot_resp(settings.expansion.speed==1)),nanstd(t_alpha_dot_resp(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_alpha_dot_resp(settings.expansion.speed==0)),nanstd(t_alpha_dot_resp(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_alpha_dot_resp,'ko')

grid on
ylabel('response time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('angular velocity based response time')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'response time (alpha_dot) versus looming speed.fig')
saveas(gca,'response time (alpha_dot) versus looming speed.jpg')

figure
hold on
errorbar(1,nanmean(t_alpha_dot_max(settings.expansion.speed==1)),nanstd(t_alpha_dot_max(settings.expansion.speed==1)),'db','linewidth',2,'MarkerSize',10,'MarkerFaceColor','b')
errorbar(0,nanmean(t_alpha_dot_max(settings.expansion.speed==0)),nanstd(t_alpha_dot_max(settings.expansion.speed==0)),'dr','linewidth',2,'MarkerSize',10,'MarkerFaceColor','r')
plot(settings.expansion.speed,t_alpha_dot_max,'ko')

grid on
ylabel('angular velocity time')

set(gca,'XTick',[0 1])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'slow','fast'})

title('time of maximum angular velocity')
set(gca,'xlim',[-.5 1.5])
set(gca,'ylim',[0 .25]) 
saveas(gca,'time of maximum angular velocity versus looming speed.fig')
saveas(gca,'time of maximum angular velocity versus looming speed.jpg')

figure
hold on
plot(1,t_Aresp,'ko')
plot(2,t_Atresp,'ko')
plot(3,t_Anresp,'ko')
plot(4,t_alpha_dot_resp,'ko')

errorbar(1,nanmean(t_Aresp),nanstd(t_Aresp),'db','linewidth',2,'MarkerSize',6,'MarkerFaceColor','b')
errorbar(2,nanmean(t_Atresp),nanstd(t_Atresp),'dr','linewidth',2,'MarkerSize',6,'MarkerFaceColor','r')
errorbar(3,nanmean(t_Anresp),nanstd(t_Anresp),'dg','linewidth',2,'MarkerSize',6,'MarkerFaceColor','g')
errorbar(4,nanmean(t_alpha_dot_resp),nanstd(t_alpha_dot_resp),'dc','linewidth',2,'MarkerSize',6,'MarkerFaceColor','c')

grid on
ylabel('response time')

set(gca,'XTick',[1:4])
set(gca,'YTick',[0:.05:.3])
set(gca,'XTickLabel',{'A','At','An','angular velocity'})

% title('acceleration based response times')
set(gca,'xlim',[0 5])
set(gca,'ylim',[0 .25])
saveas(gca,'acceleration based response times.fig')
saveas(gca,'acceleration based response times.jpg')

figure
% subplot(2,2,1)
% plot(0,0)
% hold on
% plot(0,0,'r')
% plot(0,0,'g')
% plot(0,0,'c')
% legend('slow','fast','reverse')
% 
for i=1:size(u,2)
        subplot(2,2,1)
        log_hist(A(isnan(A)==0),100)
        subplot(2,2,2)
        log_hist(At(isnan(At)==0),100)
        subplot(2,2,3)
        log_hist(An(isnan(An)==0),100)
        subplot(2,2,4)
        log_hist(alpha_dot(isnan(alpha_dot)==0),100)
end
subplot(2,2,1)
grid on
xlabel('A')
ylabel('N')
% axis([0 0.25 0 30])

subplot(2,2,2)
grid on
xlabel('At')
ylabel('N')
% axis([0 0.25 -10 20])

subplot(2,2,3)
grid on
xlabel('An')
ylabel('N')
% axis([0 0.25 0 30])

subplot(2,2,4)
grid on
xlabel('angular velocity')
ylabel('N')
% axis([0 0.25 0 8000])

saveas(gca,'hist_accel.fig')
saveas(gca,'hist_accel.jpg')

figure
% subplot(2,2,1)
% plot(0,0)
% hold on
% plot(0,0,'r')
% plot(0,0,'g')
% plot(0,0,'c')
% legend('slow','fast','reverse')
set(gca,'xscale','log') 
% 
for i=1:size(u,2)
        subplot(2,2,1)
        log_hist(A(isnan(A(:))==0),100)
        subplot(2,2,2)
        log_hist(At(isnan(At(:))==0),100)
        subplot(2,2,3)
        log_hist(An(isnan(An(:))==0),100)
        subplot(2,2,4)
        log_hist(alpha_dot(isnan(alpha_dot(:))==0),100)
end
subplot(2,2,1)
grid on
xlabel('A')
ylabel('N')
% axis([0 0.25 0 30])

subplot(2,2,2)
grid on
xlabel('At')
ylabel('N')
% axis([0 0.25 -10 20])

subplot(2,2,3)
grid on
xlabel('An')
ylabel('N')
% axis([0 0.25 0 30])

subplot(2,2,4)
grid on
xlabel('angular velocity')
ylabel('N')
% axis([0 0.25 0 8000])

saveas(gca,'hist_accel_fromtrig.fig')
saveas(gca,'hist_accel_fromtrig.jpg')

cd ..