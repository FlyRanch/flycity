% kmean flighttracks

clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')

loadname = ('flightpathDB_pos_qbodyEKF_INCroll');
% loadname = ('flightpathDB_pos_qbodyEKF_NOroll');
savename = ([loadname,'_9clusters_']);
load([loadname,'.mat'])

plotting =0

% EM based manual clustering
k=9;
At_hor_thresh =2.5;
At_hor_thresh_min =-3.3;
An_hor_thresh =3.0;

t = pathDB.t;
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t = t;

x = pathDB.pos(:,:,1);
y = pathDB.pos(:,:,2);
z = pathDB.pos(:,:,3);

u = pathDB.vel(:,:,1);
v = pathDB.vel(:,:,2);
w = pathDB.vel(:,:,3);

ax = pathDB.accel(:,:,1);
ay = pathDB.accel(:,:,2);
az = pathDB.accel(:,:,3);

% body orientation
yaw = pathDB.yaw_global;
pitch = pathDB.pitch_global;
roll = pathDB.roll;

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
for i = 1:size(x,2)
    if settings.expansion.HorPos(i) == 180
        y(:,i) = -y(:,i) + 2*center(2);
        v(:,i) = -v(:,i);
        ay(:,i) = -ay(:,i);
        yaw(:,i) = -yaw(:,i);
        roll(:,i) = -roll(:,i);
    end
end

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
V_hor = sqrt(u.^2 + v.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);
A_hor = sqrt(ax.^2 + ay.^2);

%% force production NORMALIZED
g = 9.81; % [m/s2]
strokeplane_angle = -55; % [deg]
% wl = 3 % [mm] hydei wing length
% M = 1e-6 * .1078 * wl^(3.008) % Hydei mass
% Fg = g*M;
% Fx = ax*M;
% Fy = ay*M;
% Fz = az*M - Fg;

Fx = ax/g;
Fy = ay/g;
Fz = az/g - 1;

F = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
F_hor = sqrt(Fx.^2 + Fy.^2);

% forces in body coord system
Fbx = nan(size(x));
Fby = nan(size(x));
Fbz = nan(size(x));

% forces in strokeplane coord system
Fspx = nan(size(x));
Fspy = nan(size(x));
Fspz = nan(size(x));

% pattern direction (in fly frame of reference)
stim_angle_vel = nan(size(x));
stim_angle_accel = nan(size(x));
stim_angle_yaw = nan(size(x));

accel_angle_hor_vel = nan(size(x));
accel_angle_hor_body = nan(size(x));

F_angle_hor_vel = nan(size(x));
F_angle_hor_body = nan(size(x));

heading = nan(size(x));
slip = nan(size(x));
pitch_vel = nan(size(x));

An = nan(size(x));
At = nan(size(x));
An_hor = nan(size(x));
At_hor = nan(size(x));

for j = 1:size(x,2)
    size(x,2)-j
    for i = 1:size(x,1)
        if isnan(x(i,j))==0
            
            % force vector in body & strokeplane coord system
            [Fbx(i,j) Fby(i,j) Fbz(i,j)] = world2body_euler(Fx(i,j),Fy(i,j),Fz(i,j),roll(i,j),pitch(i,j),yaw(i,j));
            [Fspx(i,j) Fspy(i,j) Fspz(i,j)] = body2strokeplane(Fbx(i,j),Fby(i,j),Fbz(i,j),strokeplane_angle);
            
            Fsp_pitch(i,j) = atand(-Fspx(i,j) / Fspz(i,j));
            Fsp_roll(i,j) = atand(-Fspy(i,j) / Fspz(i,j));
            
            % distance from center
            r(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2 + (z(i,j) - center(3))^2);
            r_hor(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2);
            
            % fly2stim 2d vector
            x_stim = x(i,j) - pos_stim(1);
            y_stim = y(i,j) - pos_stim(2);
            
            % 2d yaw vector
            x_yaw = cosd(yaw(i,j));
            y_yaw = sind(yaw(i,j));
            
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
            
            % vector in y-direction
            x_yaxis=0;
            y_yaxis=1;
            
            stim_angle_vel(i,j) =   atan2(x_vel*y_stim-y_vel*x_stim,x_vel*x_stim+y_vel*y_stim) *180/pi();
            stim_angle_accel(i,j) = atan2(x_a*y_stim-y_a*x_stim,x_a*x_stim+y_a*y_stim) *180/pi();
            stim_angle_F(i,j) = atan2(x_F*y_stim-y_F*x_stim,x_F*x_stim+y_F*y_stim) *180/pi();
            stim_angle_yaw(i,j) =   atan2(x_yaw*y_stim-y_yaw*x_stim,x_yaw*x_stim+y_yaw*y_stim) *180/pi();
            
            heading(i,j) =          atan2(x_vel*y_yaxis-y_vel*x_yaxis,x_vel*x_yaxis+y_vel*y_yaxis) *180/pi();
            slip(i,j) =             atan2(x_yaw*y_vel-y_yaw*x_vel,x_yaw*x_vel+y_yaw*y_vel) *180/pi();
            pitch_vel(i,j) =        atan2(x_pitch*z_vel-z_pitch*xy_vel,x_pitch*xy_vel+z_pitch*z_vel) *180/pi();
            
            accel_angle_hor_vel(i,j) =  atan2(x_a*y_vel-y_a*x_vel,x_a*x_vel+y_a*y_vel) *180/pi();
            accel_angle_hor_body(i,j) = atan2(x_a*y_yaw-y_a*x_yaw,x_a*x_yaw+y_a*y_yaw) *180/pi();
            
            At(i,j) = dot(vel_now,accel_now) / norm(vel_now);
            An(i,j) = sqrt(norm(accel_now)^2 - At(i,j)^2);

            At_hor(i,j) = dot(vel_hor_now,accel_hor_now) / norm(vel_hor_now);
            An_hor(i,j) = sign(accel_angle_hor_vel(i,j))*sqrt(norm(accel_hor_now)^2 - At_hor(i,j)^2);

            F_angle_hor_vel(i,j) =  atan2(x_F*y_vel-y_F*x_vel,x_F*x_vel+y_F*y_vel) *180/pi();
            F_angle_hor_body(i,j) = atan2(x_F*y_yaw-y_F*x_yaw,x_F*x_yaw+y_F*y_yaw) *180/pi();
            
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

% %% cluster variables
% X_accel = [tk An_hork At_hork];
% % X_accel = [tk An_hork At_hork stim_angle_velk stim_angle_accelk Vk];
% 
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

%% Manual cluster
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

% subsets
for i = 1:size(A,2)
    IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
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
     

%% plot data
if plotting ==1
    figure
    skip = 100
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

pathDB.accel_angle_hor_vel = accel_angle_hor_vel;
pathDB.accel_angle_hor_body = accel_angle_hor_body;

pathDB.F_angle_hor_vel = F_angle_hor_vel;
pathDB.F_angle_hor_body = F_angle_hor_body;

pathDB.heading = heading;
pathDB.slip = slip;
pathDB.pitch = pitch;
pathDB.pitch_vel = pitch_vel;
pathDB.roll_global = pathDB.roll;
pathDB.roll = roll;

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

pathDB.Fbx = Fbx;
pathDB.Fby = Fby;
pathDB.Fbz = Fbz;

pathDB.Fspx = Fspx;
pathDB.Fspy = Fspy;
pathDB.Fspz = Fspz;

pathDB.Fsp_pitch = Fsp_pitch;
pathDB.Fsp_roll = Fsp_roll;

pathDB.IDX = IDX;

settings.g=g;
settings.strokeplane_angle=strokeplane_angle;

settings.cmap_k=cmap_k;
% settings.cmap_k_abs=cmap_k_abs;
settings.cmap_360=cmap_360;
settings.kmclust_k=k;
settings.kmclust_At_hor_thresh=At_hor_thresh;
settings.kmclust_At_hor_thresh_min=At_hor_thresh_min;
settings.kmclust_An_hor_thresh=An_hor_thresh;

save([savename,num2str(At_hor_thresh),'n',num2str(At_hor_thresh_min),'n',num2str(An_hor_thresh),'.mat'],...
    'pathDB','patternDB','settings')

