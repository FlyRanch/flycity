function [ x_smooth, P_smooth ] = EKF( Y, Q, R, H, x1, P1, dt )


%Initialization

PHI_func = @EKF_PHI;

N = size(Y,2);

x = x1;

P = P1;

% Create matrices to save x_min, p_min, x_plus and p_plus.

x_min = zeros(7,N);
p_min = zeros(7,7,N);
x_plus = zeros(7,N);
p_plus = zeros(7,7,N);

% Forward Kalman Filter step.

for k=1:N
    
    x = PHI_func(x,dt)*x;
    
    P = PHI_func(x,dt)*P*PHI_func(x,dt)'+Q;
    
    x_min(:,k) = x;
    
    p_min(:,:,k) = P;
    
    K = P*H'/(H*P*H'+R);
    
    x = x + K*(Y(:,k)-H*x);
    
    P = (eye(7)-K*H)*P;
    
    x_plus(:,k) = x;
    
    p_plus(:,:,k) = P;
    
end


% Backward Rauch-Tung-Striebel smoother step.

xs = zeros(7,N);

Ps = zeros(7,7,N);

for k=N:-1:1
    
    if k == N
               
        xs(:,k) = x_plus(:,k);
        
        Ps(:,:,k) = p_plus(:,:,k);
        
    else
        
        A = p_plus(:,:,k)*(PHI_func(x_plus(:,k),dt)')/(p_min(:,:,k+1));
        
        xs(:,k) = x_plus(:,k) + A*(xs(:,k+1)-x_min(:,k+1));
        
        Ps(:,:,k) = p_plus(:,:,k)+A*(Ps(:,:,k+1)-p_min(:,:,k+1))*A';
        
    end
    
end

x_smooth = xs;

P_smooth = Ps;


