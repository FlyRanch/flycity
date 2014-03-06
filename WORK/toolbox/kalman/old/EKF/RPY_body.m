function [ X_smooth ] = RPY_body( Y, dt )

EKF_func = @EKF2;       %Extended Kalman Filter



if size(Y,1) == 4

d = 100;

a = 0;
b = 0;
c = 0.003;

Q = [a 0 0 0 0 0 0 0 0 0; ...
     0 a 0 0 0 0 0 0 0 0; ...
     0 0 a 0 0 0 0 0 0 0; ...
     0 0 0 b 0 0 0 0 0 0; ...
     0 0 0 0 b 0 0 0 0 0; ...
     0 0 0 0 0 b 0 0 0 0; ...
     0 0 0 0 0 0 c 0 0 0; ...
     0 0 0 0 0 0 0 c 0 0; ...
     0 0 0 0 0 0 0 0 c 0; ...
     0 0 0 0 0 0 0 0 0 c];

R = [1 0 0 0; ...
     0 1 0 0; ...
     0 0 1 0; ...
     0 0 0 1];

H = [0 0 0 0 0 0 1 0 0 0; ...
     0 0 0 0 0 0 0 1 0 0; ...
     0 0 0 0 0 0 0 0 1 0; ...
     0 0 0 0 0 0 0 0 0 1];


x1 = [zeros(6,1); Y(:,1)];

P1 = eye(10);

% Call Extended Kalman Filter for smoothened data for the first time

[ xs1, Ps1 ] = EKF_func( Y, Q, R, H, x1, P1, dt );

x1 = xs1(:,1);

P1 = Ps1(:,:,1);

% Call Extended Kalman Filter for smoothened data for the second time

[ xs2, Ps2 ] = EKF_func( Y, Q, R, H, x1, P1, dt );



elseif size(Y,1) == 7
    


a = 0.0001/dt;
b = 0.0001;
c = 0.0001;
d = 1;



Q = [a 0 0 0 0 0 0 0 0 0; ...
     0 a 0 0 0 0 0 0 0 0; ...
     0 0 a 0 0 0 0 0 0 0; ...
     0 0 0 b 0 0 0 0 0 0; ...
     0 0 0 0 c 0 0 0 0 0; ...
     0 0 0 0 0 b 0 0 0 0; ...
     0 0 0 0 0 0 d 0 0 0; ...
     0 0 0 0 0 0 0 d 0 0; ...
     0 0 0 0 0 0 0 0 d 0; ...
     0 0 0 0 0 0 0 0 0 d];

R = eye(7);

H = [0 0 0 1 0 0 0 0 0 0; ...
     0 0 0 0 1 0 0 0 0 0; ...
     0 0 0 0 0 1 0 0 0 0; ...
     0 0 0 0 0 0 1 0 0 0; ...
     0 0 0 0 0 0 0 1 0 0; ...
     0 0 0 0 0 0 0 0 1 0; ...
     0 0 0 0 0 0 0 0 0 1];


x1 = [zeros(3,1); Y(:,1)];

P1 = eye(10);

% Call Extended Kalman Filter for smoothened data for the first time

[ xs1, Ps1 ] = EKF_func( Y, Q, R, H, x1, P1, dt );

x1 = xs1(:,1);

P1 = Ps1(:,:,1);

% Call Extended Kalman Filter for smoothened data for the second time

[ xs2, Ps2 ] = EKF_func( Y, Q, R, H, x1, P1, dt );


end

%Return smoothened data
X_smooth(1,:) = xs2(4,:);
X_smooth(2,:) = xs2(5,:);
X_smooth(3,:) = xs2(6,:);
X_smooth(4,:) = xs2(7,:);
X_smooth(5,:) = xs2(8,:);
X_smooth(6,:) = xs2(9,:);
X_smooth(7,:) = xs2(10,:);

end




