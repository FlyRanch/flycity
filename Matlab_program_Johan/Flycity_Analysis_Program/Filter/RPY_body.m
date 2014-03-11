function [ X_smooth ] = RPY_body( Y, dt, Q_vec1, Q_vec2 )

EKF_func = @EKF;       %Extended Kalman Filter



if size(Y,1) == 4

Q = [Q_vec1(1) 0 0 0 0 0 0 0 0 0; ...
     0 Q_vec1(2) 0 0 0 0 0 0 0 0; ...
     0 0 Q_vec1(3) 0 0 0 0 0 0 0; ...
     0 0 0 Q_vec1(4) 0 0 0 0 0 0; ...
     0 0 0 0 Q_vec1(5) 0 0 0 0 0; ...
     0 0 0 0 0 Q_vec1(6) 0 0 0 0; ...
     0 0 0 0 0 0 Q_vec1(7) 0 0 0; ...
     0 0 0 0 0 0 0 Q_vec1(8) 0 0; ...
     0 0 0 0 0 0 0 0 Q_vec1(9) 0; ...
     0 0 0 0 0 0 0 0 0 Q_vec1(10)];


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


Q = [Q_vec2(1) 0 0 0 0 0 0 0 0 0; ...
     0 Q_vec2(2) 0 0 0 0 0 0 0 0; ...
     0 0 Q_vec2(3) 0 0 0 0 0 0 0; ...
     0 0 0 Q_vec2(4) 0 0 0 0 0 0; ...
     0 0 0 0 Q_vec2(5) 0 0 0 0 0; ...
     0 0 0 0 0 Q_vec2(6) 0 0 0 0; ...
     0 0 0 0 0 0 Q_vec2(7) 0 0 0; ...
     0 0 0 0 0 0 0 Q_vec2(8) 0 0; ...
     0 0 0 0 0 0 0 0 Q_vec2(9) 0; ...
     0 0 0 0 0 0 0 0 0 Q_vec2(10)];
 
 R = [1 0 0 0 0 0 0; ...
      0 1 0 0 0 0 0; ...
      0 0 1 0 0 0 0; ...
      0 0 0 1 0 0 0; ...
      0 0 0 0 1 0 0; ...
      0 0 0 0 0 1 0; ...
      0 0 0 0 0 0 1];

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

% figure()
% subplot(3,1,1); plot(X_smooth(1,:));
% subplot(3,1,2); plot(X_smooth(2,:));
% subplot(3,1,3); plot(X_smooth(3,:));
% 
% figure()
% subplot(4,1,1); plot(X_smooth(4,:));
% subplot(4,1,2); plot(X_smooth(5,:));
% subplot(4,1,3); plot(X_smooth(6,:));
% subplot(4,1,4); plot(X_smooth(7,:));
% 
% pause

end




