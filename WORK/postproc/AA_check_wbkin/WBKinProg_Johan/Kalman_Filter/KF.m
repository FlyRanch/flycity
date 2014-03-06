function [ xsmooth, Psmooth ] = KF( Y, Q, R, x1, P1, dt )

%Initialization

lx = size(Q,1);


F = [1 0 0 dt 0 0 0.5*dt^2 0 0; ...
    0 1 0 0 dt 0 0 0.5*dt^2 0; ...
    0 0 1 0 0 dt 0 0 0.5*dt^2; ...
    0 0 0 1 0 0 dt 0 0; ...
    0 0 0 0 1 0 0 dt 0; ...
    0 0 0 0 0 1 0 0 dt; ...
    0 0 0 0 0 0 1 0 0; ...
    0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 1];

H = [eye(3) zeros(3,6)];

N = size(Y,2);

x = x1;

P = P1;

% Create matrices to save x_min, p_min, x_plus and p_plus.

x_min = zeros(lx,N);
p_min = zeros(lx,lx,N);
x_plus = zeros(lx,N);
p_plus = zeros(lx,lx,N);


% Forward Kalman Filter step.

for k=1:N
    
    x = F*x;
    
    P = F*P*F'+Q;
    
    x_min(:,k) = x;
    
    p_min(:,:,k) = P;
    
    K = P*H'/(H*P*H'+R);
        
    x = x + K*(Y(:,k)-H*x);
    
    P = (eye(lx)-K*H)*P;
    
    x_plus(:,k) = x;
    
    p_plus(:,:,k) = P;
    
end

% Backward Rauch-Tung-Striebel smoother step.

xs = zeros(lx,N);
Ps = zeros(lx,lx,N);

for k=N:-1:1
    
    if k == N
        
        xs(:,k) = x_plus(:,k);
        
        Ps(:,:,k) = p_plus(:,:,k);
                
    else
        A = p_plus(:,:,k)*(F)'/(p_min(:,:,k+1));
        
        xs(:,k) = x_plus(:,k) + A*(xs(:,k+1)-x_min(:,k+1));
        
        Ps(:,:,k) = p_plus(:,:,k)+A*(Ps(:,:,k+1)-p_min(:,:,k+1))*A';
                
    end
    
end

% Return filtered and smoothened data + new initial conditions for x1 and
% P1

xsmooth = xs;

Psmooth = Ps;

