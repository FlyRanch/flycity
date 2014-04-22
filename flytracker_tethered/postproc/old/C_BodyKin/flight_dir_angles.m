clc
clear
close all
load('flightpathDB_10steps_all.mat')

% stimulus direction
Vstim = [0;1]
direction_stim = atan2(Vstim(1),Vstim(2)) *180/pi()

trigger_frame = find(DBt == min(abs(DBt)));

dt = DBt(2) - DBt(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

% vector at trigger (pre response)
dt_pre = (DBt(trigger_frame+DN_noresp,:)-DBt(trigger_frame,:))';
dx_pre = (DBx(trigger_frame+DN_noresp,:)-DBx(trigger_frame,:))';
dy_pre = (DBy(trigger_frame+DN_noresp,:)-DBy(trigger_frame,:))';
dz_pre = (DBz(trigger_frame+DN_noresp,:)-DBz(trigger_frame,:))';

vx_pre = dx_pre / dt_pre/1000;
vy_pre = dy_pre / dt_pre/1000;
vz_pre = dz_pre / dt_pre/1000;

V_pre = sqrt(vx_pre.^2 + vy_pre.^2 + vz_pre.^2);
direction_pre = rad2deg(atan2(dx_pre,dy_pre));
ascend_pre = rad2deg(atan2(dz_pre,sqrt(dx_pre.^2 + dy_pre.^2)));

for i=1:size(DBx,2)
    
%     direction_pre(i,1) = atan2(dx_pre(i),dy_pre(i)) *180/pi();
    
    % vector at end (post response)
    last_frame = find(isnan(DBx(:,i))==0, 1, 'last' );
    
    dt_post = DBt(last_frame)-DBt(last_frame-DN_escape);
    dx_post(i,1) = DBx(last_frame,i)-DBx(last_frame-DN_escape,i);
    dy_post(i,1) = DBy(last_frame,i)-DBy(last_frame-DN_escape,i);
    dz_post(i,1) = DBz(last_frame,i)-DBz(last_frame-DN_escape,i);
end

vx_post = dx_post / dt_post/1000;
vy_post = dy_post / dt_post/1000;
vz_post = dz_post / dt_post/1000;

V_post = sqrt(vx_post.^2 + vy_post.^2 + vz_post.^2);
direction_post = rad2deg(atan2(dx_post,dy_post));
ascend_post = rad2deg(atan2(dz_post,sqrt(dx_post.^2 + dy_post.^2)));
ascend_turn = ascend_post - ascend_pre;

% turn direction < 180
direction_turn = direction_post - direction_pre;
direction_turn(direction_turn>180) = direction_turn(direction_turn>180) - 360;

% direction_turn2 > 0 if animal turns away from looming object, direction_turn2 < 0: towards
direction_turn2 = abs(direction_turn);
direction_turn2(find(abs(direction_post)) > abs(direction_pre)) = - direction_turn2((find(abs(direction_post)) > abs(direction_pre)));

figure
plot(direction_pre,direction_post,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('direction pre')
ylabel('direction post')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'direction_pre_vs_post.fig')
saveas(gca,'direction_pre_vs_post.jpg')

figure
plot(direction_pre,direction_turn,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('direction pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle_vs_predirection.fig')
saveas(gca,'turnangle_vs_predirection.jpg')

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
plot(direction_pre,direction_turn2,'*')
hold on
plot([-180,180],[-180,180],'-k')
grid on
axis equal
xlabel('direction pre')
ylabel('turn angle (away)')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
saveas(gca,'turnangle2_vs_predirection.fig')
saveas(gca,'turnangle2_vs_predirection.jpg')

figure
plot(vz_pre,vz_post,'*')
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
compass_zeroup(vx_post,vy_post,'r')
hold on
compass(vx_pre,vy_pre,'b')
title('horizontal flight speed')
saveas(gca,'Vhorizontal_PreVsPost.fig')
saveas(gca,'Vhorizontal_PreVsPost.jpg')

figure
compass_zeroright(sqrt(vx_post.^2 + vy_post.^2),vz_post,'r')
hold on
compass(sqrt(vx_pre.^2 + vy_pre.^2),vz_pre,'b')
title('vertical flight speed')
saveas(gca,'Vvertical_PreVsPost.fig')
saveas(gca,'Vvertical_PreVsPost.jpg')

