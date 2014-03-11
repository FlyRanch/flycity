function [ x_smooth, P_smooth ] = EKF( Y, Q, R, H, x1, P1, dt )


%Initialization

PHI_func = @EKF_PHI;

F_func = @EKF_F;

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
     
for k=1:N
    
    xs(7:10,k) = xs(7:10,k)./norm(xs(7:10,k));
    
end

     x_smooth = xs;
     
     P_smooth = Ps;


end

