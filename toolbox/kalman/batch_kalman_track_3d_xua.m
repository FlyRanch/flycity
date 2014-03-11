% Kalmar filter 3D flight track
% State = (x y z xdot ydot zdot xdotdot ydotdot zdotdot). We only observe (x y z).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
% 
% clc
% clear
% close all
% 
% addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/Kalman')
% addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/KPMstats')
% addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/KPMtools')
% 
% load('flightpathDB_noskip.mat')
% % load('flightpathDB_skip10.mat')

savefile = 'flightpathDB_kalman.mat';

ss = 9; % state size
os = 3; % observation size
F = [1 0 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0; 0 0 1 0 0 1 0 0 0; 0 0 0 1 0 0 1 0 0; 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
H = [1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];

%% filter settings

if length(DBt)<1000
    % skip 10
    filter.x = 0.1;
    filter.u = 0.01;
    filter.a = 0.001;
    filter.R = 10;
else
    % no skip
    filter.x = 0.0001;
    filter.u = 0.00001;
    filter.a = 0.000001;
    filter.R = 10000;
end

Q = eye(ss);
Q(1:3,:) = filter.x*Q(1:3,:);
Q(4:6,:) = filter.u*Q(4:6,:);
Q(7:9,:) = filter.a*Q(7:9,:);
% Q = 0.1*eye(ss);
R = filter.R*eye(os);

%% batch loop
dt = DBt(2)-DBt(1);
t = DBt;

x = nan(size(DBx));
y = nan(size(DBx));
z = nan(size(DBx));
u = nan(size(DBx));
v = nan(size(DBx));
w = nan(size(DBx));
ax = nan(size(DBx));
ay = nan(size(DBx));
az = nan(size(DBx));

for i=1:size(DBseq,1)
    

    % measured values
    start = find(isnan(DBx(:,i))==0, 1 );
    stop = find(isnan(DBx(:,i))==0, 1, 'last' );
    
    t_obs = DBt(start:stop);
    Y_obs = ([DBx(start:stop,i) DBy(start:stop,i) DBz(start:stop,i)]./1000)';

    % initial state
    init_pos = Y_obs(:,1);
    
    if length(DBt)<1000
        init_dx = Y_obs(:,2)-Y_obs(:,1);
    else
        init_dx = (Y_obs(:,11)-Y_obs(:,1))/10;
    end
    
    init_a = [0 0 0]';
    
    initx = [init_pos;init_dx;init_a];
    initV = eye(ss); % initial state covariance...

    % initx = [10 10 1 0]';
    % initV = 10*eye(ss);

%     [xfilt, Vfilt, VVfilt, loglik] = kalman_filter(Y_obs, F, H, Q, R, initx, initV);
    [xsmooth, Vsmooth] = kalman_smoother(Y_obs, F, H, Q, R, initx, initV);

    x(start:stop,i) = xsmooth(1,:);
    y(start:stop,i) = xsmooth(2,:);
    z(start:stop,i) = xsmooth(3,:);
    
    u(start:stop,i) = xsmooth(4,:)./dt;
    v(start:stop,i) = xsmooth(5,:)./dt;
    w(start:stop,i) = xsmooth(6,:)./dt;
    
    ax(start:stop,i) = xsmooth(7,:)./dt./dt;
    ay(start:stop,i) = xsmooth(8,:)./dt./dt;
    az(start:stop,i) = xsmooth(9,:)./dt./dt;

    figure(1)
    plot3(Y_obs(1,:), Y_obs(2,:), Y_obs(3,:), 'bo');
    hold on
%     plot3(xfilt(1,:), xfilt(2,:), xfilt(3,:), 'g+:');
    plot3(xsmooth(1,:), xsmooth(2,:), xsmooth(3,:), 'rx:');
    hold off
    grid on
    axis equal
%     legend('observed', 'filtered', 'smoothed')
    legend('observed', 'smoothed')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['3D path ' num2str(DBseq(i,1)) 'S' num2str(DBseq(i,2))])
%     pause

    figure(2)
    title(['kalman results ' num2str(DBseq(i,1)) 'S' num2str(DBseq(i,2))])
    subplot(3,1,1)
    hold off
    plot(t_obs,xsmooth(1,:),'b')
    hold on
    plot(t_obs,xsmooth(2,:),'r')
    plot(t_obs,xsmooth(3,:),'g')
    grid on
    legend('x', 'y', 'z','Location','NorthWest')
    xlabel('t')
    ylabel('x')

    subplot(3,1,2)
    hold off
    plot(t_obs,xsmooth(4,:)./dt,'b')
    hold on
    plot(t_obs,xsmooth(5,:)./dt,'r')
    plot(t_obs,xsmooth(6,:)./dt,'g')
    plot(t_obs,sqrt(xsmooth(4,:).^2+xsmooth(5,:).^2+xsmooth(6,:).^2)./dt,'k')
    grid on
    legend('u', 'v', 'w','U','Location','NorthWest')
    xlabel('t')
    ylabel('u')

    subplot(3,1,3)
    hold off
    plot(t_obs,xsmooth(7,:)./dt./dt,'b')
    hold on
    plot(t_obs,xsmooth(8,:)./dt./dt,'r')
    plot(t_obs,xsmooth(9,:)./dt./dt,'g')
    plot(t_obs,sqrt(xsmooth(7,:).^2+xsmooth(8,:).^2+xsmooth(9,:).^2)./dt./dt,'k')
    grid on
    legend('ax', 'ay', 'az','A','Location','NorthWest')
    xlabel('t')
    ylabel('a')
    
%     pause
end

t = DBt;
x_obs = DBx;
y_obs = DBy;
z_obs = DBz;

save(savefile,'DBseq','t','x_obs','y_obs','z_obs','x','y','z','u','v','w','ax','ay','az','filter')