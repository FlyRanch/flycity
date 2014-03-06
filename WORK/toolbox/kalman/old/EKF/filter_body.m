function [q1,q2,q3,q4,omega_1,omega_2,omega_3] = filter_body(q1_raw,q2_raw,q3_raw,q4_raw,dt)

savefile = 'pathDB2.mat';

% Program that contains a Kalman filter which filters the body flight
% trajectory and the body orientation. The output will be the filtered body
% position plus the filtered body orientation, captured in a quaternion.

N = length(q1_raw);

% batch loop wings
q1 = zeros(1,N);
q2 = zeros(1,N);
q3 = zeros(1,N);
q4 = zeros(1,N);
omega_1 = zeros(1,N);
omega_2 = zeros(1,N);
omega_3 = zeros(1,N);
% 
% for j=1:N
    

    
    omega_on_off = 0;
      
    [Y1] = MakeY(q1_raw, q2_raw, q3_raw, q4_raw, dt, omega_on_off);

    %Run Extended Kalman filter + smoother 
    [X_smooth1] = RPY_body(Y1,dt);
    
    % Generate a new set of measurements by calculating omega from the
    % smoothened quaternions
    
    omega_on_off = 1;
      
    [Y2] = MakeY(X_smooth1(4,:), X_smooth1(5,:), X_smooth1(6,:), X_smooth1(7,:), dt, omega_on_off);

    %Run Extended Kalman filter + smoother 
    [X_smooth2] = RPY_body(Y2,dt);
    

    %Save smooth data (Note that the obtained omega's have to be converted
    %to the body frame of reference)
    omega_1 = X_smooth2(1,:);
    omega_2 = X_smooth2(2,:);
    omega_3 = X_smooth2(3,:);
    q1= X_smooth2(4,:);
    q2= X_smooth2(5,:);
    q3= X_smooth2(6,:);
    q4= X_smooth2(7,:);

   
%     clear Y1 Y2 start stop X_smooth1 X_smooth2 X_s Y1 q1_raw q2_raw q3_raw q4_raw q1_filt q2_filt q3_filt q4_filt
%     
% end


