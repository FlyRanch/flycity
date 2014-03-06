% Kalmar filter 2D flight track: only x & y.
% State = (x y xdot ydot). We only observe (x y).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)

clc
clear
close all

addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/flytracker/kalman/KalmanAll/KPMtools')

load('flightpathDB_10steps.mat')
savefile = 'flightpathDB_kalman.mat';

ss = 9; % state size
os = 3; % observation size
F = [1 0 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0; 0 0 1 0 0 1 0 0 0; 0 0 0 1 0 0 1 0 0; 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
H = [1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];

% filter settings
filter_x = 0.1;
filter_u = 0.01;
filter_a = 0.001;

Q = eye(ss);
Q(1:3,:) = filter_x*Q(1:3,:);
Q(4:6,:) = filter_u*Q(4:6,:);
Q(7:9,:) = filter_a*Q(7:9,:);

% Q = 0.1*eye(ss);

R = 10*eye(os);

% measurement values
t = DBt(isnan(DBx)==0);
y = ([DBx(isnan(DBx)==0) DBy(isnan(DBx)==0) DBz(isnan(DBx)==0)]./1000)';

% initial state
init_pos = y(:,1);
init_dx = y(:,2)-y(:,1);
initx = [init_pos;init_dx; [0 0 0]'];
initV = eye(ss); % initial state covariance...

% initx = [10 10 1 0]';
% initV = 10*eye(ss);

[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(y, F, H, Q, R, initx, initV);
[xsmooth, Vsmooth] = kalman_smoother(y, F, H, Q, R, initx, initV);

dt = DBt(2)-DBt(1);
x = xsmooth(1:3,:)';
u = (xsmooth(4:6,:)./dt)';
a = (xsmooth(7:9,:)/dt/dt)';

X = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
U = sqrt(u(:,1).^2+u(:,2).^2+u(:,3).^2);
A = sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2);

figure
plot3(y(1,:), y(2,:), y(3,:), 'bo');
hold on
plot3(xfilt(1,:), xfilt(2,:), xfilt(3,:), 'rx:');
plot3(xsmooth(1,:), xsmooth(2,:), xsmooth(3,:), 'g+:');
hold off
grid on
legend('observed', 'filtered', 'smoothed')
xlabel('x')
ylabel('y')

figure
subplot(3,1,1)
hold on
plot(t,x(:,1),'b')
plot(t,x(:,2),'r')
plot(t,x(:,3),'g')
% plot(t,X,'c')
grid on
legend('x', 'y', 'z')
xlabel('t')
ylabel('x')

subplot(3,1,2)
hold on
plot(t,u(:,1),'b')
plot(t,u(:,2),'r')
plot(t,u(:,3),'g')
plot(t,U,'c')
grid on
legend('u', 'v', 'w','U')
xlabel('t')
ylabel('u')

subplot(3,1,3)
hold on
plot(t,a(:,1),'b')
plot(t,a(:,2),'r')
plot(t,a(:,3),'g')
plot(t,A,'c')
grid on
legend('ax', 'ay', 'az','A')
xlabel('t')
ylabel('a')

save(savefile,'DBseq','t','x_obs','y_obs','z_obs','x','y','z','u','v','w','ax','ay','az','filter')
