% Kalmar filter 2D flight track: only x & y.
% State = (x y xdot ydot). We only observe (x y).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)

clc
clear

load('flightpathDB_10steps.mat')

ss = 4; % state size
os = 2; % observation size
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
H = [1 0 0 0; 0 1 0 0];
Q = 0.1*eye(ss);
R = 1*eye(os);

% measurement values
y = [DBx(isnan(DBx)==0) DBy(isnan(DBx)==0)]';

% initial state
init_pos = y(:,1);
init_dx = y(:,2)-y(:,1);
initx = [init_pos;init_dx];
initV = eye(ss); % initial state covariance...
% initx = [10 10 1 0]';
% initV = 10*eye(ss);

[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(y, F, H, Q, R, initx, initV);
[xsmooth, Vsmooth] = kalman_smoother(y, F, H, Q, R, initx, initV);

% dfilt = x([1 2],:) - xfilt([1 2],:);
% mse_filt = sqrt(sum(sum(dfilt.^2)))
% 
% dsmooth = x([1 2],:) - xsmooth([1 2],:);
% mse_smooth = sqrt(sum(sum(dsmooth.^2)))


plot(y(1,:), y(2,:), 'bo');
hold on
plot(xfilt(1,:), xfilt(2,:), 'rx:');
plot(xsmooth(1,:), xsmooth(2,:), 'g+:');
hold off
legend('observed', 'filtered', 'smoothed')
xlabel('x')
ylabel('y')

