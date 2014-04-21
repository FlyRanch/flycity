function [ x_smooth, P_smooth ] = EKF2( Y, Q, R, H, x1, P1, dt )


%Initialization

PHI_func = @EKF_PHI;

F_func = @EKF_F5;

M = size(Q,1);

N = size(Y,2);

O = size(R,1);

x = x1;

P = P1;

% Create matrices to save x_min, p_min, x_plus and p_plus.

x_min = zeros(M,N);
p_min = zeros(M,M,N);
x_plus = zeros(M,N);
p_plus = zeros(M,M,N);


% Forward Kalman Filter step.

for k=1:N
    
      
    x = PHI_func(x,dt)*x;
    
    x(M-3:M) = x(M-3:M)./norm(x(M-3:M));
    
    P = PHI_func(x,dt)*P*PHI_func(x,dt)'+Q;
    
    x_min(:,k) = x;
    
    p_min(:,:,k) = P;
    
    K = P*H'/(H*P*H'+R);
        
    x = x + K*(Y(:,k)-H*x);
    
    P = (eye(M)-K*H)*P;
    
    x_plus(:,k) = x;
    
    p_plus(:,:,k) = P;
    
end





% Backward Rauch-Tung-Striebel smoother step.

xs = zeros(M,N);

Ps = zeros(M,M,N);
    
    
for k=N:-1:1
    
    if k == N
               
        xs(:,k) = x_plus(:,k);
        
        Ps(:,:,k) = p_plus(:,:,k);
        
    else
        
        P_pred = PHI_func(x_plus(:,k),dt)*p_plus(:,:,k)*PHI_func(x_plus(:,k),dt)' + Q;
        C = p_plus(:,:,k)*PHI_func(x_plus(:,k),dt)';
        
        D = C / P_pred;

%         xs(:,k) = x_plus(:,k) + D*(xs(:,k+1) - x_min(:,k+1));

        xt = x_plus(:,k) + D*(xs(:,k+1) - x_min(:,k+1));
        
        xs(1:6,k) = xt(1:6);
        
        xs(7:10,k) = xt(7:10)./norm(xt(7:10));
        
        clear xt
        
        Ps(:,:,k) = p_plus(:,:,k) + D*(Ps(:,:,k+1)-p_min(:,:,k+1))*D';
        
    end
    
end

%     % Now transform omega_x, omega_y, omega_z from global coordinates back to
%     % body coordinates:
% 
%     for k = 1:N
% 
%          q_1(1) = xs(7,k);
%          q_1(2) = xs(8,k);
%          q_1(3) = xs(9,k);
%          q_1(4) = xs(10,k);
% 
%          omega_world = [ xs(4,k); xs(5,k); xs(6,k)];
% 
%          R   = [ q_1(4)^2+q_1(1)^2-q_1(2)^2-q_1(3)^2    2*q_1(1)*q_1(2)+2*q_1(4)*q_1(3)      2*q_1(1)*q_1(3)-2*q_1(4)*q_1(2); ...
%                  2*q_1(1)*q_1(2)-2*q_1(4)*q_1(3)        q_1(4)^2-q_1(1)^2+q_1(2)^2-q_1(3)^2  2*q_1(2)*q_1(3)+2*q_1(4)*q_1(1); ...
%                  2*q_1(1)*q_1(3)+2*q_1(4)*q_1(2)        2*q_1(2)*q_1(3)-2*q_1(4)*q_1(1)      q_1(4)^2-q_1(1)^2-q_1(2)^2+q_1(3)^2];  
% 
%          omega_body = R*omega_world;
% 
%          xs(4,k) = omega_body(1);
%          xs(5,k) = omega_body(2);
%          xs(6,k) = omega_body(3);
% 
%          clear omega_world omega_body q_1 R
% 
%     end

     
     x_smooth = xs;
     
     P_smooth = Ps;


end

