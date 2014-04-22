% clc
% clear
close all
% load('flightpathDB_kalman.mat')

% stimulus direction
Vstim = [0;1]
heading_stim = atan2(Vstim(1),Vstim(2)) *180/pi()

trigger_frame = find(t == min(abs(t)));

dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

% vector at trigger (pre response)
u_pre = mean(u(trigger_frame:trigger_frame+DN_noresp,:))';
v_pre = mean(v(trigger_frame:trigger_frame+DN_noresp,:))';
w_pre = mean(w(trigger_frame:trigger_frame+DN_noresp,:))';

V_pre = sqrt(u_pre.^2 + v_pre.^2 + w_pre.^2);
heading_pre = rad2deg(atan2(u_pre,v_pre));
ascend_pre = rad2deg(atan2(w_pre,sqrt(u_pre.^2 + v_pre.^2)));

for i=1:size(u,2)
    
    % vector at end (post response)
    last_frame = find(isnan(u(:,i))==0, 1, 'last' );
    
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

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);

% turn rate
heading = rad2deg(atan2(u,v));
ascend = rad2deg(atan2(w,sqrt(u.^2 + v.^2)));

turn_rate = nan(size(heading));
for i=2:size(heading,1)-2
 turn_rate(i,:) = heading(i+1,:) - heading(i,:);
end
turn_rate = turn_rate/dt;

turn_rate(abs(turn_rate)>10000)  = nan;

figure
for i=1:size(u,2)
    
    subplot(3,1,1)
    hold on
    plot(t,turn_rate(:,i))
    
    subplot(3,1,2)
    hold on
    plot(t,V(:,i))
    
    subplot(3,1,3)
    hold on
    plot(t,A(:,i))
    
end
subplot(3,1,1)
grid on
xlabel('t')
ylabel('turn rate')

subplot(3,1,2)
grid on
xlabel('t')
ylabel('V')

subplot(3,1,3)
grid on
xlabel('t')
ylabel('A')

saveas(gca,'flightpaths_turn_V_A.fig')
saveas(gca,'flightpaths_turn_V_A.jpg')

figure
plot(heading_pre,heading_post,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('heading pre')
ylabel('heading post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'heading_pre_vs_post.fig')
saveas(gca,'heading_pre_vs_post.jpg')

figure
plot(heading_pre,heading_turn,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle_vs_preheading.fig')
saveas(gca,'turnangle_vs_preheading.jpg')

figure
plot(ascend_pre,ascend_post,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('ascend pre')
ylabel('ascend post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'ascend_pre_vs_post.fig')
saveas(gca,'ascend_pre_vs_post.jpg')

figure
plot(ascend_pre,ascend_turn,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('ascend pre')
ylabel('ascend change')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'ascendchange_vs_ascend_pre.fig')
saveas(gca,'ascendchange_vs_ascend_pre.jpg')

figure
plot(heading_pre,heading_turn2,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('heading pre')
ylabel('turn angle (away)')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle2_vs_preheading.fig')
saveas(gca,'turnangle2_vs_preheading.jpg')

figure
plot(w_pre,w_post,'*')
hold on
grid on
axis equal
xlabel('Vz pre')
ylabel('Vz post')
saveas(gca,'Vz_preVSpost.fig')
saveas(gca,'Vz_preVSpost.jpg')

figure
plot(V_pre,V_post,'*')
hold on
grid on
axis equal
xlabel('V pre')
ylabel('V post')
saveas(gca,'Vpre_vs_Vpost.fig')
saveas(gca,'Vpre_vs_Vpost.jpg')

figure
compass_zeroup(u_post,v_post,'r')
hold on
compass(u_pre,v_pre,'b')
title('horizontal flight speed')
saveas(gca,'Vhorizontal_PreVsPost.fig')
saveas(gca,'Vhorizontal_PreVsPost.jpg')

figure
compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
hold on
compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
title('vertical flight speed')
saveas(gca,'Vvertical_PreVsPost.fig')
saveas(gca,'Vvertical_PreVsPost.jpg')

