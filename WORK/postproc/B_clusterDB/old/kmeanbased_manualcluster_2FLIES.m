% kmean flighttracks

clc
clear
close all

% addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos.mat')

t = pathDB.t;
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

% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);
A_hor = sqrt(ax.^2 + ay.^2);

% A
An = pathDB.An;
At = pathDB.At;
An_hor = pathDB.An_hor;
At_hor = pathDB.At_hor;
alpha_dot_hor = pathDB.alpha_dot_hor;

% kcluster db
tk=t;
for i=1:size(V,2)-1
    tk = [tk;t];
end
t = tk;

% headingk = stim_angle(:);
% r_hork = r_hor(:);
% Vk = V(:);
% 
Ak=A(:);
Ank=An(:);
Atk=At(:);
A_hork=A_hor(:);
At_hork=At_hor(:);
An_hork=An_hor(:);

alpha_dot_hork = alpha_dot_hor(:);

% angle of fly2 in fly1 view, and vice versa
for i = 1:size(x,1)
    %% Fly1
    
    % stimulus direction: horizontal head to head vector
    x_stim = x(i,1) - x(i,3);
    y_stim = y(i,1) - y(i,3);
    
    % yaw angle: horizontal fly vector (head to tail)
    x_yaw = x(i,1) - x(i,2); 
    y_yaw = y(i,1) - y(i,2);

    % 2d velocity vector
    x_vel = u(i,1);
    y_vel = v(i,1);

    % vector in y-direction
    x_yaxis=0;
    y_yaxis=1;
    
    stim_angle_vel(i,1:2) =   atan2(x_vel*y_stim-y_vel*x_stim,x_vel*x_stim+y_vel*y_stim) *180/pi();
    stim_angle_yaw(i,1:2) =   atan2(x_yaw*y_stim-y_yaw*x_stim,x_yaw*x_stim+y_yaw*y_stim) *180/pi();
    heading(i,1:2) =          atan2(x_vel*y_yaxis-y_vel*x_yaxis,x_vel*x_yaxis+y_vel*y_yaxis) *180/pi();
    slip(i,1:2) =             atan2(x_yaw*y_vel-y_yaw*x_vel,x_yaw*x_vel+y_yaw*y_vel) *180/pi();
            
    
%     % remove jumps from 180deg to -180deg
%     stim_angle_plot(i,1:2) = stim_angle(i,1:2);
%     if i>1 && abs(stim_angle_plot(i,1) - stim_angle_plot(i-1,1))>180
%         stim_angle_plot(i,1:2) = nan;
%     end
    
    %% Fly2
    
    % stimulus direction: horizontal head to head vector
    x_stim = x(i,3) - x(i,1);
    y_stim = y(i,3) - y(i,1);
    
    % yaw angle: horizontal fly vector (head to tail)
    x_yaw = x(i,3) - x(i,4); 
    y_yaw = y(i,3) - y(i,4);

    % 2d velocity vector
    x_vel = u(i,3);
    y_vel = v(i,3);

    % vector in y-direction
    x_yaxis=0;
    y_yaxis=1;
    
    stim_angle_vel(i,3:4) =   atan2(x_vel*y_stim-y_vel*x_stim,x_vel*x_stim+y_vel*y_stim) *180/pi();
    stim_angle_yaw(i,3:4) =   atan2(x_yaw*y_stim-y_yaw*x_stim,x_yaw*x_stim+y_yaw*y_stim) *180/pi();
    heading(i,3:4) =          atan2(x_vel*y_yaxis-y_vel*x_yaxis,x_vel*x_yaxis+y_vel*y_yaxis) *180/pi();
    slip(i,3:4) =             atan2(x_yaw*y_vel-y_yaw*x_vel,x_yaw*x_vel+y_yaw*y_vel) *180/pi();

%     % remove jumps from 180deg to -180deg
%     stim_angle_plot(i,3:4) = stim_angle(i,3:4);
%     if i>1 && abs(stim_angle_plot(i,1) - stim_angle_plot(i-1,1))>180
%         stim_angle_plot(i,3:4) = nan;
%     end
end

%% cluster variables
X_accel = [tk An_hork At_hork];
% X_accel = [tk headingk r_hork Vk A_hork An_hork At_hork alpha_dot_hork];
% X_accel = [tk headingk r_hork Vk A_hork abs(An_hork) At_hork abs(alpha_dot_hork)];
% X_accel = [tk Ank Atk];
% X_accel = [tk Ank Atk];

%% k-mean cluster

% % inc seeds
% k=9;
% seeds = [0.1 0 0; 0.1 10 0; 0.1 -5 0; 0.1 0 10; 0.1 0 -5; 0.1 10 10; 0.1 10 -5; 0.1 -10 10; 0.1 -10 -5];
% IDXk = kmeans(X_accel,k,'start', seeds);
% 
% % no seeds
% k = 6
% IDXk = kmeans(X_accel,k);

% manual clustering
k=9;
At_hor_thresh =3;
An_hor_thresh =3.5;

for i=1:length(An_hork)
    if isnan(An_hork(i)) == 1
        IDXk(i,1) = nan;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) < -At_hor_thresh
        IDXk(i,1) = 1;
    elseif An_hork(i) < -An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
        IDXk(i,1) = 2;
    elseif An_hork(i) < -An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 3;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) < -At_hor_thresh
        IDXk(i,1) = 4;
    elseif abs(An_hork(i)) < An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
        IDXk(i,1) = 5;
    elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 6;
    elseif An_hork(i) > An_hor_thresh && At_hork(i) < -At_hor_thresh
        IDXk(i,1) = 7;
    elseif An_hork(i) > An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
        IDXk(i,1) = 8;
    elseif An_hork(i) > An_hor_thresh && At_hork(i) > At_hor_thresh
        IDXk(i,1) = 9;
    end
end

% subsets
for i = 1:size(A,2)
    IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
end

%% plot data
% cmap = cmap(1:end/k:end,:);

% cluster cmap
cmap_k = [.5 0 1;  0  0  1; 0 1 1;... % dark purple; blue; light blue
         1 0 1; .5 .5 .5; 0 1 0;... % light purple; grey; green;
         1 0 0;  1 .5  0; 1 1 0];   % red; orange; yellow
     
cmap_k = [0 1  0;  1  1  0; 1 .5  0;... % green; yellow; orange
         0 1 1; .5 .5 .5; 1 0 0;... % cyan; grey; red;
         0 0 1; .5  0 1; 1 0 1];   % blue; purple; magenta

cmap_k_abs = [0 1 0;  1  1  0; 1 .5 0;... % green; yellow; orange
              0 0 1; .5 .5 .5; 1  0 0]; % blue; grey; red;

% circular cmap
cmap_360r = [zeros(45,1); [0:1/(45-1):1]'; ones(3*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(45,1)];
cmap_360g = [ones(2*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]'];
cmap_360b = [[1:-1/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]';ones(3*45,1)];

cmap_360 = [cmap_360r cmap_360g cmap_360b];

%% state clusters
figure
for i = 1:size(IDXk)
    size(IDXk)-i
    hold on
    if isnan(IDXk(i)) == 0
%             plot3(X_accel(i,1),X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
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
settings.cmap_k_abs=cmap_k_abs;
settings.cmap_360=cmap_360;
settings.kmclust_k=k;
settings.kmclust_At_hor_thresh=At_hor_thresh;
settings.kmclust_An_hor_thresh=An_hor_thresh;

save('flightpathDB_2flies_9clusters','pathDB','settings')

