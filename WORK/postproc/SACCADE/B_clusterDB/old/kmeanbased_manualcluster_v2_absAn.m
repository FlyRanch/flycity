% kmean flighttracks

clc
clear
close all

% addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
% load('flightpathDB_pos_qbodyEKF.mat')
load('flightpathDB_pos_qbodyEKF_posneg.mat')

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

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);
A_hor = sqrt(ax.^2 + ay.^2);

% A
An = abs(pathDB.An);
At = pathDB.At;
An_hor = abs(pathDB.An_hor);
At_hor = pathDB.At_hor;
alpha_dot_hor = abs(pathDB.alpha_dot_hor);

% body orientation
yaw = pathDB.yaw_global;
pitch = pathDB.pitch_global;

% temporal data
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;
t_start = -0.05;
t_stop = 0.2;

% arena data
% arenacenter from fly positions
% center = [nanmean(x(1,:)) nanmean(y(1,:)) nanmean(z(1,:))];
center = [nanmean(x(trigger_frame,:)) nanmean(y(trigger_frame,:)) nanmean(z(trigger_frame,:))];
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
    end
end


% pattern direction (in fly frame of reference)
for j = 1:size(x,2)
    size(x,2)-j
    for i = 1:size(x,1)
        if isnan(x(i,j))==0
            % distance from center
            r(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2 + (z(i,j) - center(3))^2);
            r_hor(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2);

            
            % fly2stim 2d vector
            x_stim = x(i,j) - pos_stim(1);
            y_stim = y(i,j) - pos_stim(2);
            
            % 2d yaw vector
            x_yaw = cosd(yaw(i,j));
            y_yaw = sind(yaw(i,j));
            
            % 2d velocity vector
            x_vel = u(i,j);
            y_vel = v(i,j);
            
            % vector in y-direction
            x_yaxis=0;
            y_yaxis=1;
            
            stim_angle_vel(i,j) =   atan2(x_vel*y_stim-y_vel*x_stim,x_vel*x_stim+y_vel*y_stim) *180/pi();
            stim_angle_yaw(i,j) =   atan2(x_yaw*y_stim-y_yaw*x_stim,x_yaw*x_stim+y_yaw*y_stim) *180/pi();
            heading(i,j) =          atan2(x_vel*y_yaxis-y_vel*x_yaxis,x_vel*x_yaxis+y_vel*y_yaxis) *180/pi();
            slip(i,j) =             atan2(x_yaw*y_vel-y_yaw*x_vel,x_yaw*x_vel+y_yaw*y_vel) *180/pi();
            
%             hold off
%             quiver(x_vel,y_vel,'b')
%             hold on
%             quiver(x_yaw,y_yaw,'r')
%             quiver(x_stim,y_stim,'g')
%             quiver(x_yaxis,y_yaxis,'k')
%             axis equal
%             legend('vel','yaw','stim','y-axis')



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
r_hork = r_hor(:);
Vk = V(:);

Ak=A(:);
Ank=An(:);
Atk=At(:);
A_hork=A_hor(:);
At_hork=At_hor(:);
An_hork=An_hor(:);

alpha_dot_hork = alpha_dot_hor(:);

%% cluster variables
X_accel = [tk An_hork At_hork];
% X_accel = [tk An_hork At_hork stim_angle_velk Vk];

%% k-mean cluster

% % inc seeds
% k=9;
% seeds = [0.1 0 0; 0.1 10 0; 0.1 -5 0; 0.1 0 10; 0.1 0 -5; 0.1 10 10; 0.1 10 -5; 0.1 -10 10; 0.1 -10 -5];
% IDXk = kmeans(X_accel,k,'start', seeds);
% 
% % no seeds
% k = 7
% IDXk = kmeans(X_accel,k);

% manual clustering
k=9;
At_hor_thresh =2.5;
An_hor_thresh =2.5;

for i=1:length(An_hork)
    if isnan(An_hork(i)) == 1
        IDXk(i,1) = nan;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) < -At_hor_thresh
        IDXk(i,1) = 1;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 3;
    elseif An_hork(i) < -An_hor_thresh
        IDXk(i,1) = 2;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) < -At_hor_thresh
        IDXk(i,1) = 4;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 6;
    elseif abs(An_hork(i)) < An_hor_thresh
        IDXk(i,1) = 5;
    elseif An_hork(i) > An_hor_thresh && At_hork(i) < -At_hor_thresh
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

%% plot data
% cmap = cmap(1:end/k:end,:);

% cluster cmap
cmap_k = [1 1 1;  1  1  1; 1 1 1;... % black
         0 0 1; .5 .5 .5; 1 0 0;... % blue ; grey; red;
         0 1 0;  1 1 0; 1 .5  0];   % green; yellow; orange
     
% cmap_k = [.5 0 1;  0  0  1; 0 1 1;... % dark purple; blue; light blue
%          1 0 1; .5 .5 .5; 0 1 0;... % light purple; grey; green;
%          1 0 0;  1 .5  0; 1 1 0];   % red; orange; yellow
%      
% cmap_k = [0 1 0;  1  1  0; 1 .5 0;... % green; yellow; orange
%           0 1 1; .5 .5 .5; 1  0 0;... % cyan; grey; red;
%           0 0 1; .5  0  1; 1  0 1];   % blue; purple; magenta
% 
% cmap_k_abs = [0 1 0;  1  1  0; 1 .5 0;... % green; yellow; orange
%               0 0 1; .5 .5 .5; 1  0 0]; % blue; grey; red;
% 
% cmap_k = [0 0 1; .5 .5 .5; 1  0 0;... % blue; grey; red;
%           0 1 0;  1  1  0; 1 .5 0]; % green; yellow; orange
         

% circular cmap
cmap_360r = [zeros(45,1); [0:1/(45-1):1]'; ones(3*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(45,1)];
cmap_360g = [ones(2*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]'];
cmap_360b = [[1:-1/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]';ones(3*45,1)];

cmap_360 = [cmap_360r cmap_360g cmap_360b];
     

%% state clusters
figure
for i = 1:100:size(IDXk)
    length(IDXk)-i
    hold on
    if isnan(IDXk(i)) == 0
%         plot3(X_accel(i,1),X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
        plot(X_accel(i,2),X_accel(i,3),'.','color',cmap_k(IDXk(i),:))
    end
end
grid on
xlabel('An hor')
ylabel('At hor')

pathDB.stim_angle_vel = stim_angle_vel;
pathDB.stim_angle_yaw = stim_angle_yaw;
pathDB.heading = heading;
pathDB.slip = slip;
pathDB.V = V;
pathDB.A = A;
pathDB.A_hor = A_hor;
pathDB.IDX = IDX;

settings.cmap_k=cmap_k;
% settings.cmap_k_abs=cmap_k_abs;
settings.cmap_360=cmap_360;
settings.kmclust_k=k;
settings.kmclust_At_hor_thresh=At_hor_thresh;
settings.kmclust_An_hor_thresh=An_hor_thresh;

save(['flightpathDB_pos_qbodyEKF_9clusters_absAn_',num2str(At_hor_thresh),'n',num2str(An_hor_thresh),'.mat'],'pathDB','settings')

