% Kalmar filter 3D flight track
% State = (x y z xdot ydot zdot xdotdot ydotdot zdotdot). We only observe (x y z).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
% 

function [pos,vel,accel] = batch_kalman_pos(pos_obs,t,settings,plotting,saving)



%% filter settings
filter = settings.kalman_pos;

Q = eye(filter.ss);
Q(1:3,:) = filter.x*Q(1:3,:);
Q(4:6,:) = filter.u*Q(4:6,:);
Q(7:9,:) = filter.a*Q(7:9,:);
% Q = 0.1*eye(filter.ss);
R = filter.R*eye(filter.os);

%% batch loop
dt = t(2)-t(1);

x_obs = pos_obs(:,:,1);
y_obs = pos_obs(:,:,2);
z_obs = pos_obs(:,:,3);

x = nan(size(x_obs));
y = nan(size(x_obs));
z = nan(size(x_obs));
u = nan(size(x_obs));
v = nan(size(x_obs));
w = nan(size(x_obs));
ax = nan(size(x_obs));
ay = nan(size(x_obs));
az = nan(size(x_obs));

%seq loop
for i=1:size(x_obs,2)
    counter = size(x_obs,2) - i

    % measured values
    start = find(isnan(x_obs(:,i))==0, 1 );
    stop = find(isnan(x_obs(:,i))==0, 1, 'last' );
    
    t_now = t(start:stop);
    Y = ([x_obs(start:stop,i) y_obs(start:stop,i) z_obs(start:stop,i)]./1000)';

    
    %% 1st kalman smoother
    % initial state
    init_pos = Y(:,1);
    init_dx = (Y(:,11)-Y(:,1))/10;
    init_a = [0 0 0]';
    
    initx = [init_pos;init_dx;init_a];
    initV = eye(filter.ss); % initial state covariance...
    % initV = 10*eye(filter.ss);

%     [xfilt, Vfilt, VVfilt, loglik] = kalman_filter(Y, filter.F, filter.H, Q, R, initx, initV);
    [xsmooth, Vsmooth] = kalman_smoother(Y, filter.F, filter.H, Q, R, initx, initV);
    
    %% 2nd kalman smoother
    % initial state
    init_pos = xsmooth(1:3,1);
    init_dx = xsmooth(4:6,1);
    init_a = xsmooth(7:9,1);
    
    initx = [init_pos;init_dx;init_a];
    initV = eye(filter.ss); % initial state covariance...
    % initV = 10*eye(filter.ss);

%     [xfilt, Vfilt, VVfilt, loglik] = kalman_filter(Y, filter.F, filter.H, Q, R, initx, initV);
    [xsmooth, Vsmooth] = kalman_smoother(Y, filter.F, filter.H, Q, R, initx, initV);
    

%     % johan's smoother
%     [ xsmooth ] = XYZ_body( Y*1000, dt );
%     xsmooth = xsmooth/1000;
    %%

    x(start:stop,i) = xsmooth(1,:);
    y(start:stop,i) = xsmooth(2,:);
    z(start:stop,i) = xsmooth(3,:);
    
    u(start:stop,i) = xsmooth(4,:)./dt;
    v(start:stop,i) = xsmooth(5,:)./dt;
    w(start:stop,i) = xsmooth(6,:)./dt;
    
    ax(start:stop,i) = xsmooth(7,:)./dt./dt;
    ay(start:stop,i) = xsmooth(8,:)./dt./dt;
    az(start:stop,i) = xsmooth(9,:)./dt./dt;

    if plotting == 1
        figure(1)
        plot3(Y(1,:), Y(2,:), Y(3,:), 'b.');
        hold on
    %     plot3(xfilt(1,:), xfilt(2,:), xfilt(3,:), 'g+:');
        plot3(xsmooth(1,:), xsmooth(2,:), xsmooth(3,:), '.r');
        hold off
        grid on
        axis equal
    %     legend('observed', 'filtered', 'smoothed')
        legend('observed', 'smoothed')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
        if saving == 1
            title(['3D path ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2))])
        %     saveas(gca,['3D path ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.fig'])
            saveas(gca,['3D path ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.jpg'])
        end
    %     pause

        figure(2)
        subplot(3,1,1)
        hold off
        plot(t_now,xsmooth(1,:),'b')
        hold on
        plot(t_now,xsmooth(2,:),'r')
        plot(t_now,xsmooth(3,:),'g')
        grid on
        legend('x', 'y', 'z','Location','NorthWest')
        xlabel('t')
        ylabel('x')

        subplot(3,1,2)
        hold off
        plot(t_now,xsmooth(4,:)./dt,'b')
        hold on
        plot(t_now,xsmooth(5,:)./dt,'r')
        plot(t_now,xsmooth(6,:)./dt,'g')
        plot(t_now,sqrt(xsmooth(4,:).^2+xsmooth(5,:).^2+xsmooth(6,:).^2)./dt,'k')
        grid on
        legend('u', 'v', 'w','U','Location','NorthWest')
        xlabel('t')
        ylabel('u')

        subplot(3,1,3)
        hold off
        plot(t_now,xsmooth(7,:)./dt./dt,'b')
        hold on
        plot(t_now,xsmooth(8,:)./dt./dt,'r')
        plot(t_now,xsmooth(9,:)./dt./dt,'g')
        plot(t_now,sqrt(xsmooth(7,:).^2+xsmooth(8,:).^2+xsmooth(9,:).^2)./dt./dt,'k')
        grid on
        legend('ax', 'ay', 'az','A','Location','NorthWest')
        xlabel('t')
        ylabel('a')
        
        if saving ==1
            title(['kalman results ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2))])
        %     saveas(gca,['xua ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.fig'])
            saveas(gca,['xua ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.jpg'])
        end
    %     pause
    end
end

pos = nan(size(pos_obs));
vel = nan(size(pos_obs));
accel = nan(size(pos_obs));

pos(:,:,1) = x;
pos(:,:,2) = y;
pos(:,:,3) = z;

vel(:,:,1) = u;
vel(:,:,2) = v;
vel(:,:,3) = w;

accel(:,:,1) = ax;
accel(:,:,2) = ay;
accel(:,:,3) = az;
