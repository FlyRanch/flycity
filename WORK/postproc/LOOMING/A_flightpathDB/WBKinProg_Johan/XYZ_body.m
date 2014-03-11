function [ X_smooth ] = XYZ_body( Y, dt )

KF_func = @KF;

% Initialization

Q =  [dt^2 0 0 0 0 0 0 0 0; ...
     0 dt^2 0 0 0 0 0 0 0; ...
     0 0 dt^2 0 0 0 0 0 0; ...
     0 0 0 1 0 0 0 0 0; ...
     0 0 0 0 1 0 0 0 0; ...
     0 0 0 0 0 1 0 0 0; ...
     0 0 0 0 0 0 1/dt^2 0 0; ...
     0 0 0 0 0 0 0 1/dt^2 0; ...
     0 0 0 0 0 0 0 0 1/dt^2];


R = eye(3);

 
P1 = eye(9);

x1 = [Y(:,1); (Y(:,20)-Y(:,1))/(20*dt); zeros(3,1)];

% Run KF for first time:

[xs1, Ps1] = KF_func( Y, Q, R, x1, P1, dt );

x1 = xs1(:,1);

P1 = Ps1(:,:,1);

[xs2, Ps2] = KF_func( Y, Q, R, x1, P1, dt );

X_smooth(1:3,:) = xs2(1:3,:);
X_smooth(4:6,:) = xs2(4:6,:);
X_smooth(7:9,:) = xs2(7:9,:);

clear xs1 xs2