% kmean flighttracks

clc
clear
close all

% addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_posneg.mat')

t = pathDB.t;
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t = t(trigger_frame:end,:);

x = pathDB.pos(trigger_frame:end,:,1);
y = pathDB.pos(trigger_frame:end,:,2);
z = pathDB.pos(trigger_frame:end,:,3);

u = pathDB.vel(trigger_frame:end,:,1);
v = pathDB.vel(trigger_frame:end,:,2);
w = pathDB.vel(trigger_frame:end,:,3);

ax = pathDB.accel(trigger_frame:end,:,1);
ay = pathDB.accel(trigger_frame:end,:,2);
az = pathDB.accel(trigger_frame:end,:,3);


% vel & accel
V = sqrt(u.^2 + v.^2 + w.^2);
A = sqrt(ax.^2 + ay.^2 + az.^2);
A_hor = sqrt(ax.^2 + ay.^2);

% A
An = pathDB.An(trigger_frame:end,:);
At = pathDB.At(trigger_frame:end,:);
An_hor = pathDB.An_hor(trigger_frame:end,:);
At_hor = pathDB.At_hor(trigger_frame:end,:);
alpha_dot_hor = pathDB.alpha_dot_hor(trigger_frame:end,:);


% temporal data
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;
t_start = -0.05;
t_stop = 0.2;

% arena data
% arenacenter from fly positions
center = [nanmean(x(1,:)) nanmean(y(1,:)) nanmean(z(1,:))];
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
for i = 1:size(x,1)
    for j = 1:size(x,2)
        size(x,1) - i
        % distance from center
        r(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2 + (z(i,j) - center(3))^2);
        r_hor(i,j) = sqrt((x(i,j) - center(1))^2 + (y(i,j) - center(2))^2);
        
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

% max accels & normalized
for i=1:size(A,2)
    A_hor_norm(:,i) = A_hor(:,i)./(max(A_hor(:,i))-min(A_hor(:,i)));
    At_hor_norm(:,i) = At_hor(:,i)./(max(At_hor(:,i))-min(At_hor(:,i)));
    An_hor_norm(:,i) = An_hor(:,i)./(max(An_hor(:,i))-min(An_hor(:,i)));
    
    A_hor_now = A_hor(:,i);
    At_hor_now = At_hor(:,i);
    An_hor_now = An_hor(:,i);
    alpha_dot_hor_now = alpha_dot_hor(:,i);
    
    A_hor_max(i) = max(A_hor_now);
    tA_hor_max(i) = min(t(A_hor_now==max(A_hor_now)));
    if tA_hor_max(i) == t(1);
        A_hor_max(i) = nan;
        tA_hor_max(i) = nan;
    elseif tA_hor_max(i) == t(end);
        A_hor_max(i) = nan;
        tA_hor_max(i) = nan;
    end
    
    At_hor_max(i) = max(At_hor_now);
    tAt_hor_max(i) = min(t(At_hor_now==max(At_hor_now)));
    if tAt_hor_max(i) == t(1);
        At_hor_max(i) = nan;
        tAt_hor_max(i) = nan;
    elseif tAt_hor_max(i) == t(end);
        At_hor_max(i) = nan;
        tAt_hor_max(i) = nan;
    end
    
    An_hor_max(i) = max(An_hor_now);
    tAn_hor_max(i) = min(t(An_hor_now==max(An_hor_now)));
    if tAn_hor_max(i) == t(1);
        An_hor_max(i) = nan;
        tAn_hor_max(i) = nan;
    elseif tAn_hor_max(i) == t(end);
        An_hor_max(i) = nan;
        tAn_hor_max(i) = nan;
    end
    
    alpha_dot_hor_max(i) = max(alpha_dot_hor_now);
    talpha_dot_hor_max(i) = min(t(alpha_dot_hor_now==max(alpha_dot_hor_now)));
    if talpha_dot_hor_max(i) == t(1);
        alpha_dot_hor_max(i) = nan;
        talpha_dot_hor_max(i) = nan;
    elseif talpha_dot_hor_max(i) == t(end);
        alpha_dot_hor_max(i) = nan;
        talpha_dot_hor_max(i) = nan;
    end
end

% kcluster db
tk=t;
for i=1:size(V,2)-1
    tk = [tk;t];
end
t = tk;

headingk = stim_angle(:);
r_hork = r_hor(:);
Vk = V(:);

Ak=A(:);
Ank=An(:);
Atk=At(:);
A_hork=A_hor(:);
At_hork=At_hor(:);
An_hork=An_hor(:);
A_hor_nork=A_hor_norm(:);
At_hor_nork=At_hor_norm(:);
An_hor_nork=An_hor_norm(:);

alpha_dot_hork = alpha_dot_hor(:);

%% cluster variables
% X_accel = [tk An_hor_nork At_hor_nork];
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
% no seeds
k = 9
IDXk = kmeans(X_accel,k);

% % manual clustering
% k=9;
% At_hor_thresh =3;
% An_hor_thresh =3.5;
% 
% for i=1:length(An_hork)
%     if isnan(An_hork(i)) == 1
%         IDXk(i,1) = nan;
%     elseif An_hork(i) < -An_hor_thresh && At_hork(i) < -At_hor_thresh
%         IDXk(i,1) = 1;
%     elseif An_hork(i) < -An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
%         IDXk(i,1) = 2;
%     elseif An_hork(i) < -An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 3;
%     elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) < -At_hor_thresh
%         IDXk(i,1) = 4;
%     elseif abs(An_hork(i)) < An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
%         IDXk(i,1) = 5;
%     elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 6;
%     elseif An_hork(i) > An_hor_thresh && At_hork(i) < -At_hor_thresh
%         IDXk(i,1) = 7;
%     elseif An_hork(i) > An_hor_thresh && abs(At_hork(i)) < At_hor_thresh
%         IDXk(i,1) = 8;
%     elseif An_hork(i) > An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 9;
%     end
% end

% subsets
for i = 1:size(A,2)
    IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
end

%% plot data
cmap = colormap('jet');
cmap = cmap(1:end/k:end,:);

cmap = [.5 0 1;  0  0  1; 0 1 1;...
         1 0 1; .5 .5 .5; 0 1 0;...
         1 0 0;  1 .5  0; 1 1 0];
%      
% cmap = [ 0 0 1; 0 1 1; 0  1 0;...
%         .5 0 1; 0 0 0; 1  1 0;...
%          1 0 1; 1 0 0; 1 .5 0];
     
% close all
figure
for i = 1:100:size(IDXk)
    i
    hold on
    if isnan(IDXk(i)) == 0
%             plot3(X_accel(i,1),X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
        plot(X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
    end
end
grid on
xlabel('An hor')
ylabel('At hor')
%     xlabel('t')
%     ylabel('An hor')
%     zlabel('At hor')

    
% flightpaths
figure
for j = 1:size(A,2)
    
            subplot(3,1,1)
            hold off

            subplot(3,1,2)
            hold off

            subplot(3,1,3)
            hold off

            
    for i = 1:5:size(A,1)
        if isnan(IDX(i,j)) == 0
            
            subplot(3,1,1)
            plot(t(i),stim_angle(i,j),'.','color',cmap(IDX(i,j),:))
            ylabel('heading')
            hold on

            subplot(3,1,2)
            plot(t(i),An_hor(i,j),'.','color',cmap(IDX(i,j),:))
            ylabel('An')
            hold on
            
            subplot(3,1,3)
            plot(t(i),At_hor(i,j),'.','color',cmap(IDX(i,j),:))
            ylabel('At')
            hold on
        end
    end
    pause
end


    
for j = 1:size(A,2)
    
            subplot(3,1,1)
            hold off
            subplot(3,1,2)
            hold off
            subplot(3,1,3)
            hold off
            
    for i = 1:5:size(A,1)
        if isnan(IDX(i,j)) == 0
            
            subplot(3,1,1)
            plot(t(i),stim_angle(i,j),'.','color',cmap(IDX(i,j),:))
            hold on

            subplot(3,1,2)
            plot(t(i),V(i,j),'.','color',cmap(IDX(i,j),:))
            hold on
            
            subplot(3,1,3)
            plot(t(i),alpha_dot_hor(i,j),'.','color',cmap(IDX(i,j),:))
            hold on
        end
    end
    pause
end


pathDB.IDX = IDX;
settings.kmclust_k=cmap;
settings.kmclust_k=k;
settings.kmclust_At_hor_thresh=At_hor_thresh;
settings.kmclust_An_hor_thresh=An_hor_thresh;

save('flightpathDB_pos_qbodyEKF_9clusters','pathDB','settings')

