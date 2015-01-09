% kmean flighttracks

clc
clear
close all

loadname = ('flightpathDB_pos_qbodyEKF_INCroll.mat');
if exist(loadname) ~= 2
    loadname = ('flightpathDB_pos_qbodyEKF_NOroll.mat');
end

load(loadname)

%% EM based manual clustering
k=9;
At_hor_thresh =1.4;
At_hor_thresh_min =-1.7;
An_hor_thresh =1.9;

% A_hor_thresh = 2.7;
A_hor_thresh = 2.75; % new kalman filter

A_thresh = 3; % NOT FROM EM

%% constants
var_name = dir('flyVar*');
var_name = var_name.name;
load(var_name);
strokeplane_angle = strokeplane_WBkin;

% % strokeplane_angle = -55; % [deg]
% % strokeplane_angle = -45; % [deg]
% strokeplane_angle = -47.5; % [deg]
% g = 9.81; % [m/s2]

%% rotation axis angle
% % rot_axis_angle = 45; % [deg]
% rot_axis_angle = 30; % [deg]
% rot_axis_angle = 30.1; % [deg] !!! average accel axis angle during SACCADE

%% average torque axis angle during SACCADE
% 
% % HI inertia inc added mass sphere, -.025<tpre<-.005, .03<tpost<-.06
% rot_axis_angle = 25.5; % [deg] !!!
% 
% % LOW inertia inc added mass disc, -.025<tpre<-.005, .03<tpost<-.06
% rot_axis_angle = 23.2; % [deg] !!!
% 
% % MEDIUM inertia mean of added mass sphere/disc, -.025<tpre<-.005, .03<tpost<-.06
% rot_axis_angle = 24.5; % [deg] !!!

%% responseDB info

% % start at 150 frames after trigger
% nstart = 150;
% 
% trigger_frame = find(pathDB.t == min(abs(pathDB.t)));
% start_frame = trigger_frame + nstart;

% start at frame 1
start_frame = 1;
start_frame = 145;
start_frame = 100;

savename = ([loadname(1:end-4),'_strokeplane',num2str(-strokeplane_angle),'deg_rotaxis',num2str(rot_axis_angle),'deg_startframe',num2str(start_frame),'.mat']);
% savename = ([loadname(1:end-4),'_2clusters_Ahor',num2str(A_hor_thresh),'mps2_strokeplane',num2str(-strokeplane_angle),'deg_startframe',num2str(start_frame),'.mat']);

%% plotting clusters & timelines?
toplot_clusterdistr =1;
toplot_clusterdistr =0;
skip = 10

toplot=1;
toplot=0;
if toplot == 1
    figdir = ['accel_plots_',num2str(At_hor_thresh),'n',num2str(At_hor_thresh_min),'n',num2str(An_hor_thresh),'mps2_startframe',num2str(start_frame)];
%     figdir = ['accel_plots_Ahor',num2str(A_hor_thresh),'mps2_startframe',num2str(start_frame)];
    mkdir(figdir)
end

%% cluster cmap
cmap_k = [0 0 0;  0 0 0; 0 0 0;... % black
         0 0 1; .5 .5 .5; 1 0 0;... % blue ; grey; red;
         0 1 0;  1 1 0; 1 .5  0];   % green; yellow; orange

% circular cmap
cmap_360r = [zeros(45,1); [0:1/(45-1):1]'; ones(3*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(45,1)];
cmap_360g = [ones(2*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]'];
cmap_360b = [[1:-1/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]';ones(3*45,1)];

cmap_360 = [cmap_360r cmap_360g cmap_360b];

%%
t = pathDB.t;
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);

x = pathDB.pos(:,:,1);
y = pathDB.pos(:,:,2);
z = pathDB.pos(:,:,3);

u = pathDB.vel(:,:,1);
v = pathDB.vel(:,:,2);
w = pathDB.vel(:,:,3);

ax = pathDB.accel(:,:,1);
ay = pathDB.accel(:,:,2);
az = pathDB.accel(:,:,3);

yaw_global = pathDB.yaw_global;
pitch_global = pathDB.pitch_global;
roll_global = pathDB.roll_global;

yaw = pathDB.yaw;
pitch = pathDB.pitch;
roll = pathDB.roll;

roll_dot_body  = radtodeg(pathDB.omega(:,:,1));
pitch_dot_body = radtodeg(pathDB.omega(:,:,2));
yaw_dot_body   = radtodeg(pathDB.omega(:,:,3));

%% REMOVE FIRST FRAMES (start_frame)

for i = 1:size(x,2)
    
    x_now = x(:,i);
    start_frame_now = start_frame + find(isnan(x_now)==0, 1 )-1;
    
    x(1:start_frame_now,i)=nan;
    y(1:start_frame_now,i)=nan;
    z(1:start_frame_now,i)=nan;
    
    u(1:start_frame_now,i)=nan;
    v(1:start_frame_now,i)=nan;
    w(1:start_frame_now,i)=nan;
    
    ax(1:start_frame_now,i)=nan;
    ay(1:start_frame_now,i)=nan;
    az(1:start_frame_now,i)=nan;
    
    yaw_global(1:start_frame_now,i)=nan;
    pitch_global(1:start_frame_now,i)=nan;
    roll_global(1:start_frame_now,i)=nan;
    
    yaw(1:start_frame_now,i)=nan;
    pitch(1:start_frame_now,i)=nan;
    roll(1:start_frame_now,i)=nan;
    
    yaw_dot_body(1:start_frame_now,i)=nan;
    pitch_dot_body(1:start_frame_now,i)=nan;
    roll_dot_body(1:start_frame_now,i)=nan;
end
    
%% calc attitude in stroke pane coord system
roll_dot_sp = nan(size(roll_dot_body));
pitch_dot_sp = nan(size(roll_dot_body));
yaw_dot_sp = nan(size(roll_dot_body));
for i = 1:size(roll_dot_body,1)
    for j = 1:size(roll_dot_body,2)
        if isnan(roll_dot_body(i,j))==0
            [roll_dot_sp(i,j) pitch_dot_sp(i,j) yaw_dot_sp(i,j)] = body2strokeplane(roll_dot_body(i,j),pitch_dot_body(i,j),yaw_dot_body(i,j),strokeplane_angle);
        end
    end
end

% body orientation
qbody = pathDB.qbody;

% % yawrollpitch angles
% [roll2,pitch2,yaw2] = qbody2angles_manual_yawrollpitch(qbody);
% 
% % adjust roll & pitch < -45
% roll2 = -roll2;
% yaw2(pitch2<-15) = -yaw2(pitch2<-15);
% pitch2(pitch2<-15) = pitch2(pitch2<-15) + 180;

% temporal data
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);

% arena data
% arenacenter from fly positions
center = [nanmean(x(trigger_frame,:)) nanmean(y(trigger_frame,:)) nanmean(z(trigger_frame,:))];
% center = [nanmean(x(1,:)) nanmean(y(1,:)) nanmean(z(1,:))];
r_arena = 0.120; % [m] arena radius

% stimulus data
pos_stim = center - [0 r_arena 0];
Vstim = [0;1] % direction
heading_stim = atan2(Vstim(1),Vstim(2)) *180/pi()

% reverse data with reverse pattern
% for i = 1:size(x,2)
%     if settings.expansion.HorPos(i) == 180
%         y(:,i) = -y(:,i) + 2*center(2);
%         v(:,i) = -v(:,i);
%         ay(:,i) = -ay(:,i);
%         yaw(:,i) = -yaw(:,i);
%         roll(:,i) = -roll(:,i);
%         yaw_global(:,i) = -yaw_global(:,i);
%         roll_global(:,i) = -roll_global(:,i);
%         yaw_dot_body(:,i) = -yaw_dot_body(:,i);
%         roll_dot_body(:,i) = -roll_dot_body(:,i);
%         yaw_dot_sp(:,i) = -yaw_dot_sp(:,i);
%         roll_dot_sp(:,i) = -roll_dot_sp(:,i);
% %         yaw2(:,i) = -yaw2(:,i);
% %         roll2(:,i) = -roll2(:,i);
%     end
% end

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
V_hor = sqrt(u.^2 + v.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);
A_hor = sqrt(ax.^2 + ay.^2);

%% OMEGA LnR rot axes from roll_dot & pitch_dot
% rot_dot_sp_L = (roll_dot_sp - pitch_dot_sp) * cosd(rot_axis_angle);
% rot_dot_sp_R = (roll_dot_sp + pitch_dot_sp) * cosd(rot_axis_angle);
% 
% rot_dot_body_L = (roll_dot_body - pitch_dot_body) * cosd(rot_axis_angle);
% rot_dot_body_R = (roll_dot_body + pitch_dot_body) * cosd(rot_axis_angle);

rot_dot_sp_L = roll_dot_sp*sind(rot_axis_angle) - pitch_dot_sp*cosd(rot_axis_angle);
rot_dot_sp_R = roll_dot_sp*cosd(rot_axis_angle) + pitch_dot_sp*sind(rot_axis_angle);

rot_dot_body_L = roll_dot_body*sind(rot_axis_angle) - pitch_dot_body*cosd(rot_axis_angle);
rot_dot_body_R = roll_dot_body*cosd(rot_axis_angle) + pitch_dot_body*sind(rot_axis_angle);

%% force production NORMALIZED
% wl = 3 % [mm] hydei wing length
% M = 1e-6 * .1078 * wl^(3.008) % Hydei mass
% Fg = g*M;
% Fx = ax*M;
% Fy = ay*M;
% Fz = az*M + Fg;

Fx = ax/g;
Fy = ay/g;
Fz = az/g + 1;

F = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
F_hor = sqrt(Fx.^2 + Fy.^2);

% pattern direction (in fly frame of reference)
stim_angle_vel = nan(size(x));
stim_angle_accel = nan(size(x));
stim_angle_yaw = nan(size(x));
stim_angle_F = nan(size(x));
stim_angle_spn = nan(size(x));

angle_vel = nan(size(x));
angle_accel = nan(size(x));
angle_yaw = nan(size(x));
angle_F = nan(size(x));
angle_spn = nan(size(x));

accel_angle_hor_vel = nan(size(x));
accel_angle_hor_body = nan(size(x));

F_angle_hor_vel = nan(size(x));
F_angle_hor_body = nan(size(x));

Fb_pitch = nan(size(x));
Fb_roll = nan(size(x));

Fsp_pitch = nan(size(x));
Fsp_roll = nan(size(x));

omega_turn_body = nan(size(x));
omega_turn_axis_body = nan(size(x));

omega_turn_sp = nan(size(x));
omega_turn_axis_sp = nan(size(x));

droll_body = nan(size(x));
dpitch_body = nan(size(x));
dyaw_body = nan(size(x));

droll_sp = nan(size(x));
dpitch_sp = nan(size(x));
dyaw_sp = nan(size(x));

drot_sp_L = nan(size(rot_dot_sp_L));
drot_sp_R = nan(size(rot_dot_sp_R));

drot_body_L = nan(size(rot_dot_body_L));
drot_body_R = nan(size(rot_dot_body_R));

slip = nan(size(x));
% slip2 = nan(size(x));

slip_aero = nan(size(x));
pitch_aero = nan(size(x));

slip_aero_hor = nan(size(x));
pitch_aero_hor = nan(size(x));

slip_body = nan(size(x));
pitch_body = nan(size(x));

An = nan(size(x));
At = nan(size(x));
An_hor = nan(size(x));
At_hor = nan(size(x));

Fn = nan(size(x));
Ft = nan(size(x));
Fn_hor = nan(size(x));
Ft_hor = nan(size(x));

for j = 1:size(x,2)
    size(x,2)-j
    for i = 1:size(x,1)
        if isnan(x(i,j))==0
            
            % omega data
            omega_turn_body(i,j) = norm([pitch_dot_body(i,j),yaw_dot_body(i,j)]);
            omega_turn_axis_body(i,j) = radtodeg(atan2(yaw_dot_body(i,j),pitch_dot_body(i,j)));
            
            omega_turn_sp(i,j) = norm([pitch_dot_sp(i,j),yaw_dot_sp(i,j)]);
            omega_turn_axis_sp(i,j) = radtodeg(atan2(yaw_dot_sp(i,j),pitch_dot_sp(i,j)));
            
%             if i > trigger_frame
                % previous d's
                if isnan(droll_body(i-1,j))==0
                    droll_body_prev = droll_body(i-1,j);
                    dpitch_body_prev = dpitch_body(i-1,j);
                    dyaw_body_prev = dyaw_body(i-1,j);
                    
                    droll_sp_prev = droll_sp(i-1,j);
                    dpitch_sp_prev = dpitch_sp(i-1,j);
                    dyaw_sp_prev = dyaw_sp(i-1,j);
                    
                    drot_sp_L_prev = drot_sp_L(i-1,j);
                    drot_sp_R_prev = drot_sp_R(i-1,j);
                    
                    drot_body_L_prev = drot_body_L(i-1,j);
                    drot_body_R_prev = drot_body_R(i-1,j);
                else
                    droll_body_prev = 0;
                    dpitch_body_prev = 0;
                    dyaw_body_prev = 0;
                    
                    droll_sp_prev = 0;
                    dpitch_sp_prev = 0;
                    dyaw_sp_prev = 0;
                    
                    drot_sp_L_prev = 0;
                    drot_sp_R_prev = 0;
                    
                    drot_body_L_prev = 0;
                    drot_body_R_prev = 0;
                end
                % integrate
                droll_body(i,j) = droll_body_prev + roll_dot_body(i,j)*dt;
                dpitch_body(i,j) = dpitch_body_prev + pitch_dot_body(i,j)*dt;
                dyaw_body(i,j) = dyaw_body_prev + yaw_dot_body(i,j)*dt;

                droll_sp(i,j) = droll_sp_prev + roll_dot_sp(i,j)*dt;
                dpitch_sp(i,j) = dpitch_sp_prev + pitch_dot_sp(i,j)*dt;
                dyaw_sp(i,j) = dyaw_sp_prev + yaw_dot_sp(i,j)*dt;
                
                drot_sp_L(i,j) = drot_sp_L_prev + rot_dot_sp_L(i,j)*dt;
                drot_sp_R(i,j) = drot_sp_R_prev + rot_dot_sp_R(i,j)*dt;
                
                drot_body_L(i,j) = drot_body_L_prev + rot_dot_body_L(i,j)*dt;
                drot_body_R(i,j) = drot_body_R_prev + rot_dot_body_R(i,j)*dt;
%             end
            
            % force & vel vectors in body coord system (UNMIRROR for calc with q)
            if settings.expansion.HorPos(j) == 180
                [x_Vb y_Vb z_Vb] = world2body_q(u(i,j),-v(i,j),w(i,j),qbody(i,j,:));
                [x_Vb_hor y_Vb_hor z_Vb_hor] = world2body_q(u(i,j),-v(i,j),0,qbody(i,j,:));
                [x_Fb y_Fb z_Fb] = world2body_q(Fx(i,j),-Fy(i,j),Fz(i,j),qbody(i,j,:));
                
                y_Vb = -y_Vb;
                y_Vb_hor = -y_Vb_hor;
                y_Fb = -y_Fb;
                y_Fsp = -y_Fsp;
            else
                [x_Vb y_Vb z_Vb] = world2body_q(u(i,j),v(i,j),w(i,j),qbody(i,j,:));
                [x_Vb_hor y_Vb_hor z_Vb_hor] = world2body_q(u(i,j),v(i,j),0,qbody(i,j,:));
                [x_Fb y_Fb z_Fb] = world2body_q(Fx(i,j),Fy(i,j),Fz(i,j),qbody(i,j,:));
            end
            Vb = norm([x_Vb y_Vb z_Vb]);
            Vb_hor = norm([x_Vb y_Vb z_Vb]);
            
            % strokeplane normal in world coordinates
            [x_spn y_spn z_spn] = strokeplane_normal2world(strokeplane_angle,qbody(i,j,:));
            
            % force vector in SP coord system
            [x_Fsp y_Fsp z_Fsp] = body2strokeplane(x_Fb,y_Fb,z_Fb,strokeplane_angle);
            
            % slip & pitch in aero frame of reference
            pitch_aero(i,j) = atan2(z_Vb,x_Vb) * 180/pi();
            slip_aero(i,j) =  asind(y_Vb/Vb);

            % slip & pitch in HORIZONTAL aero frame of reference
            pitch_aero_hor(i,j) = atan2(z_Vb_hor,x_Vb_hor) * 180/pi();
            slip_aero_hor(i,j) =  asind(y_Vb_hor/Vb_hor);

            % force vector pitch & yaw
            Fb_pitch(i,j) =          atan2(-z_Fb,x_Fb) *180/pi();
            Fb_roll(i,j)  =          atan2(y_Fb,-z_Fb) *180/pi();
            Fsp_pitch(i,j) =         atan2(x_Fsp,-z_Fsp) *180/pi();
            Fsp_roll(i,j)  =         atan2(y_Fsp,-z_Fsp) *180/pi();
            
            % distance from center
            r(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2 + (z(i,j) - center(3))^2);
            r_hor(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2);
            
            % fly2stim 2d vector
            x_stim = x(i,j) - pos_stim(1);
            y_stim = y(i,j) - pos_stim(2);
            
            % 2d yaw vector
            x_yaw = cosd(yaw(i,j));
            y_yaw = sind(yaw(i,j));
%             x_yaw2 = cosd(yaw2(i,j));
%             y_yaw2 = sind(yaw2(i,j));
            
            % 2d pitch vector
            x_pitch = cosd(pitch(i,j));
            z_pitch = sind(pitch(i,j));
            
            % velocity vectors
            x_vel = u(i,j);
            y_vel = v(i,j);
            z_vel = w(i,j);
            
            xy_vel = sqrt(x_vel.^2 + y_vel.^2);
            
            vel_now(1:3,1) = [x_vel y_vel z_vel];
            vel_hor_now(1:2,1) = [x_vel y_vel];
            
            % acceleration vectors
            x_a = ax(i,j);
            y_a = ay(i,j);
            z_a = ay(i,j);
            
            accel_now(1:3,1) = [x_a y_a z_a];
            accel_hor_now(1:2,1) = [x_a y_a];
            
            % F vectors
            x_F = Fx(i,j);
            y_F = Fy(i,j);
            z_F = Fy(i,j);
            
            F_now(1:3,1) = [x_F y_F z_F];
            F_hor_now(1:2,1) = [x_F y_F];
                        
            stim_angle_vel(i,j)     = atan2(x_vel*y_stim-y_vel*x_stim,x_vel*x_stim+y_vel*y_stim) *180/pi();
            stim_angle_accel(i,j)   = atan2(x_a*y_stim-y_a*x_stim,x_a*x_stim+y_a*y_stim) *180/pi();
            stim_angle_F(i,j)       = atan2(x_F*y_stim-y_F*x_stim,x_F*x_stim+y_F*y_stim) *180/pi();
            stim_angle_yaw(i,j)     = atan2(x_yaw*y_stim-y_yaw*x_stim,x_yaw*x_stim+y_yaw*y_stim) *180/pi();
            stim_angle_spn(i,j)     = atan2(x_spn*y_stim-y_spn*x_stim,x_spn*x_stim+y_spn*y_stim) *180/pi();
            
            angle_vel(i,j)     = atan2(x_vel,y_vel) *180/pi();
            angle_accel(i,j)   = atan2(x_a,y_a) *180/pi();
            angle_F(i,j)       = atan2(x_F,y_F) *180/pi();
            angle_yaw(i,j)     = atan2(x_yaw,y_yaw) *180/pi();
            angle_spn(i,j)     = atan2(x_spn,y_spn) *180/pi();
            
            slip(i,j) =             atan2(x_yaw*y_vel-y_yaw*x_vel,x_yaw*x_vel+y_yaw*y_vel) *180/pi();
            x_slip = cosd(slip(i,j));
            y_slip = sind(slip(i,j));
            
%             slip2(i,j) =             atan2(x_yaw2*y_vel-y_yaw2*x_vel,x_yaw2*x_vel+y_yaw2*y_vel) *180/pi();
%             pitch_vel(i,j) =        atan2(x_pitch*z_vel-z_pitch*xy_vel,x_pitch*xy_vel+z_pitch*z_vel) *180/pi();

            % rotate pitch & slip around x-axis with roll angle
            Xx = 1;
            Xy = tand(slip(i,j));
            Xz = sqrt(1+Xy^2) * tand(pitch_global(i,j));
            
            Xsp = [Xx;Xy;Xz];
            Xsp_b = rollrotx(Xsp,roll(i,j));
            
            slip_body(i,j) = atand(Xsp_b(2));
            pitch_body(i,j) = -atand(Xsp_b(3) / sqrt(1+Xsp_b(2)^2));

            
            accel_angle_hor_vel(i,j)    = atan2(x_a*y_vel-y_a*x_vel,x_a*x_vel+y_a*y_vel) *180/pi();
            accel_angle_hor_body(i,j)   = atan2(x_a*y_yaw-y_a*x_yaw,x_a*x_yaw+y_a*y_yaw) *180/pi();
            
            At(i,j) = dot(vel_now,accel_now) / norm(vel_now);
            An(i,j) = sqrt(norm(accel_now)^2 - At(i,j)^2);

            At_hor(i,j) = dot(vel_hor_now,accel_hor_now) / norm(vel_hor_now);
            An_hor(i,j) = sign(accel_angle_hor_vel(i,j))*sqrt(norm(accel_hor_now)^2 - At_hor(i,j)^2);

            F_angle_hor_vel(i,j)    =  atan2(x_F*y_vel-y_F*x_vel,x_F*x_vel+y_F*y_vel) *180/pi();
            F_angle_hor_body(i,j)   = atan2(x_F*y_yaw-y_F*x_yaw,x_F*x_yaw+y_F*y_yaw) *180/pi();
            
            Ft(i,j) = dot(vel_now,F_now) / norm(vel_now);
            Fn(i,j) = sqrt(norm(F_now)^2 - Ft(i,j)^2);

            Ft_hor(i,j) = dot(vel_hor_now,F_hor_now) / norm(vel_hor_now);
            Fn_hor(i,j) = sign(F_angle_hor_vel(i,j))*sqrt(norm(F_hor_now)^2 - Ft_hor(i,j)^2);
            
%             hold off
%             quiver(x_vel,y_vel,'b')
%             hold on
%             quiver(x_a,y_a,'b')
%             quiver(x_yaw,y_yaw,'r')
%             quiver(x_stim,y_stim,'g')
%             quiver(x_yaxis,y_yaxis,'k')
%             quiver(x_pitch,z_pitch,'k')
%             quiver(xy_vel,z_vel,'b')
%             axis equal
%             legend('vel','accel','yaw','stim','y-axis')

    %         % remove jumps from 180deg to -180deg
    %         stim_angle_plot(i,j) = stim_angle(i,j);
    %         if i>1 && j>1 && abs(stim_angle_plot(i,j) - stim_angle_plot(i-1,j))>180
    %             stim_angle_plot(i,j) = nan;
    %         end
        end
    end
end

% kcluster db
tk=t;
for i=1:size(V,2)-1
    tk = [tk;t];
end
t = tk;

stim_angle_velk = stim_angle_vel(:);
stim_angle_accelk = stim_angle_accel(:);
r_hork = r_hor(:);
Vk = V(:);

Ak=A(:);
Ank=An(:);
Atk=At(:);
A_hork=A_hor(:);
At_hork=At_hor(:);
% An_hork=An_hor(:);
An_hork=abs(An_hor(:));
% alpha_dot_hork = alpha_dot_hor(:);

%% cluster variables
X_accel = [tk An_hork At_hork];
% X_accel = [tk An_hork At_hork stim_angle_velk stim_angle_accelk Vk];

% %% k-mean cluster
% 
% % inc seeds
% k=9;
% seeds = [0.1 0 0; 0.1 10 0; 0.1 -5 0; 0.1 0 10; 0.1 0 -5; 0.1 10 10; 0.1 10 -5; 0.1 -10 10; 0.1 -10 -5];
% IDXk = kmeans(X_accel,k,'start', seeds);
% 
% % no seeds
% k = 7
% IDXk = kmeans(X_accel,k);

%% Manual cluster An & At HOR
for i=1:length(An_hork)
    if isnan(An_hork(i)) == 1
        IDXk(i,1) = nan;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) < At_hor_thresh_min
        IDXk(i,1) = 1;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 3;
    elseif An_hork(i) < -An_hor_thresh
        IDXk(i,1) = 2;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) < At_hor_thresh_min
        IDXk(i,1) = 4;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 6;
    elseif abs(An_hork(i)) < An_hor_thresh
        IDXk(i,1) = 5;
    elseif An_hork(i) > An_hor_thresh && At_hork(i) < At_hor_thresh_min
        IDXk(i,1) = 7;
    elseif An_hork(i) > An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 9;
    elseif An_hork(i) > An_hor_thresh
        IDXk(i,1) = 8;
    end
end

%% Manual cluster A
% for i=1:length(Ak)
%     if isnan(Ak(i)) == 1
%         IDXk(i,1) = nan;
%     elseif Ak(i) < A_thresh
%         IDXk(i,1) = 5;
%     else
%         IDXk(i,1) = 9;
%     end
% end

%% Manual cluster A hor
% for i=1:length(A_hork)
%     if isnan(A_hork(i)) == 1
%         IDXk(i,1) = nan;
%     elseif A_hork(i) < A_hor_thresh
%         IDXk(i,1) = 5;
%     else
%         IDXk(i,1) = 9;
%     end
% end

%% subsets
for i = 1:size(A,2)
    IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
end     

%% plot data
if toplot_clusterdistr ==1
    figure
%     skip = 100
    for i = 1:skip:size(IDXk)
        length(IDXk)-i
        hold on
        if isnan(IDXk(i)) == 0
    %         plot3(X_accel(i,1),X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
            plot(X_accel(i,2),X_accel(i,3),'ok','markerfacecolor',cmap_k(IDXk(i),:),'markersize',5)
    %         plot(X_accel(i,1),X_accel(i,2),'ok','markerfacecolor',cmap_k(IDXk(i),:),'markersize',5)
        end
    end
    grid on
    axis equal
    xlabel('An','fontsize',18)
    ylabel('At','fontsize',18)
    set(gca,'xlim',[0 20],'ylim',[-20 20])
    set(gca,'XTick',[-20:10:20],'fontsize',12)
    set(gca,'YTick',[-20:10:20],'fontsize',12)

    saveas(gca,['kmeanbased_manualclusters_absAn_skip',num2str(skip),'.fig'])
    saveas(gca,['kmeanbased_manualclusters_absAn_skip',num2str(skip),'.png'])
    plot2svg(['kmeanbased_manualclusters_absAn_skip',num2str(skip),'.svg'])
end

pathDB.stim_angle_vel = stim_angle_vel;
pathDB.stim_angle_accel = stim_angle_accel;
pathDB.stim_angle_F = stim_angle_F;
pathDB.stim_angle_yaw = stim_angle_yaw;
pathDB.stim_angle_spn = stim_angle_spn;

pathDB.angle_vel = angle_vel;
pathDB.angle_accel = angle_accel;
pathDB.angle_F = angle_F;
pathDB.angle_yaw = angle_yaw;
pathDB.angle_spn = angle_spn;

pathDB.accel_angle_hor_vel = accel_angle_hor_vel;
pathDB.accel_angle_hor_body = accel_angle_hor_body;

pathDB.F_angle_hor_vel = F_angle_hor_vel;
pathDB.F_angle_hor_body = F_angle_hor_body;

pathDB.roll_global = pathDB.roll;

pathDB.slip = slip;
pathDB.pitch = pitch;
pathDB.roll = roll;

pathDB.yaw_dot_body = yaw_dot_body;
pathDB.pitch_dot_body = pitch_dot_body;
pathDB.roll_dot_body = roll_dot_body;

pathDB.dyaw_body = dyaw_body;
pathDB.dpitch_body = dpitch_body;
pathDB.droll_body = droll_body;

pathDB.yaw_dot_sp = yaw_dot_sp;
pathDB.pitch_dot_sp = pitch_dot_sp;
pathDB.roll_dot_sp = roll_dot_sp;

pathDB.dyaw_sp = dyaw_sp;
pathDB.dpitch_sp = dpitch_sp;
pathDB.droll_sp = droll_sp;

pathDB.drot_body_L = drot_body_L;
pathDB.drot_body_R = drot_body_R;

pathDB.drot_sp_L = drot_sp_L;
pathDB.drot_sp_R = drot_sp_R;

pathDB.rot_dot_body_L = rot_dot_body_L;
pathDB.rot_dot_body_R = rot_dot_body_R;

pathDB.rot_dot_sp_L = rot_dot_sp_L;
pathDB.rot_dot_sp_R = rot_dot_sp_R;

pathDB.omega_turn_body = omega_turn_body;
pathDB.omega_turn_axis_body = omega_turn_axis_body;

pathDB.omega_turn_sp = omega_turn_sp;
pathDB.omega_turn_axis_sp = omega_turn_axis_sp;

% pathDB.slip2 = slip2;
% pathDB.pitch2 = pitch2;
% pathDB.yaw2 = yaw2;
% pathDB.roll2 = roll2;

pathDB.slip_aero = slip_aero;
pathDB.pitch_aero = pitch_aero;

pathDB.slip_aero_hor = slip_aero_hor;
pathDB.pitch_aero_hor = pitch_aero_hor;

pathDB.slip_body = slip_body;
pathDB.pitch_body = pitch_body;

pathDB.V = V;
pathDB.V_hor = V_hor;

pathDB.A = A;
pathDB.At = At;
pathDB.An = An;

pathDB.A_hor = A_hor;
pathDB.At_hor = At_hor;
pathDB.An_hor = An_hor;

pathDB.F = F;
pathDB.Ft = Ft;
pathDB.Fn = Fn;

pathDB.F_hor = F_hor;
pathDB.Ft_hor = Ft_hor;
pathDB.Fn_hor = Fn_hor;

pathDB.Fx = Fx;
pathDB.Fy = Fy;
pathDB.Fz = Fz;

pathDB.Fsp_pitch = Fsp_pitch;
pathDB.Fsp_roll = Fsp_roll;

pathDB.Fb_pitch = Fb_pitch;
pathDB.Fb_roll = Fb_roll;

pathDB.IDX = IDX;

settings.g=g;
settings.strokeplane_angle=strokeplane_angle;
settings.rot_axis_angle=rot_axis_angle;

settings.cmap_k=cmap_k;
% settings.cmap_k_abs=cmap_k_abs;
settings.cmap_360=cmap_360;
settings.kmclust_k=k;
settings.kmclust_At_hor_thresh=At_hor_thresh;
settings.kmclust_At_hor_thresh_min=At_hor_thresh_min;
settings.kmclust_An_hor_thresh=An_hor_thresh;

%% run responseDB file
batch_calc_responseDB_subfunc_SACCADE

save(savename,'pathDB','patternDB','settings','responseDB')

