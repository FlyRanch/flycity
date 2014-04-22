function [temp] = filter_body_nosave(settings,pathDB)

savefile = 'pathDB2.mat';

% Program that contains a Kalman filter which filters the body flight
% trajectory and the body orientation. The output will be the filtered body
% position plus the filtered body orientation, captured in a quaternion.

dt = pathDB.t(2)-pathDB.t(1);
t = pathDB.t;

% batch loop body
x = nan(size(pathDB.x));
y = nan(size(pathDB.x));
z = nan(size(pathDB.x));
u = nan(size(pathDB.x));
v = nan(size(pathDB.x));
w = nan(size(pathDB.x));
ax = nan(size(pathDB.x));
ay = nan(size(pathDB.x));
az = nan(size(pathDB.x));

for i=1:size(pathDB.x,2)
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.x(:,i))==0, 1 );
    stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
    
    % select measurements from NaN values
    Y = ([pathDB.x(start:stop,i) pathDB.y(start:stop,i) pathDB.z(start:stop,i)])';
    
    %Run Linear Kalman filter + smoother (with or without ALS Q and R optimization)
    [X_smooth] = XYZ_body(Y,dt);
    
    
    %Save smooth data
    x(start:stop,i) = X_smooth(1,:);
    y(start:stop,i) = X_smooth(2,:);
    z(start:stop,i) = X_smooth(3,:);
    
    u(start:stop,i) = X_smooth(4,:);
    v(start:stop,i) = X_smooth(5,:);
    w(start:stop,i) = X_smooth(6,:);
    
    ax(start:stop,i) = X_smooth(7,:);
    ay(start:stop,i) = X_smooth(8,:);
    az(start:stop,i) = X_smooth(9,:);
       
    clear Y start stop X_smooth
end



% batch loop wings
q1 = nan(size(pathDB.x));
q2 = nan(size(pathDB.x));
q3 = nan(size(pathDB.x));
q4 = nan(size(pathDB.x));
omega_1 = nan(size(pathDB.x));
omega_2 = nan(size(pathDB.x));
omega_3 = nan(size(pathDB.x));

for j=1:size(pathDB.x,2)
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.x(:,j))==0, 1 );
    stop = find(isnan(pathDB.x(:,j))==0, 1, 'last' );
   
    q1_raw = pathDB.qb1(start:stop,j);
    q2_raw = pathDB.qb2(start:stop,j);
    q3_raw = pathDB.qb3(start:stop,j);
    q4_raw = pathDB.qb4(start:stop,j);
    
    omega_on_off = 0;
      
    [Y1] = MakeY(q1_raw, q2_raw, q3_raw, q4_raw, start, stop, dt, omega_on_off);

    %Run Extended Kalman filter + smoother 
    [X_smooth1] = RPY_body(Y1,dt);
    
    % Generate a new set of measurements by calculating omega from the
    % smoothened quaternions
    
    omega_on_off = 1;
      
    [Y2] = MakeY(X_smooth1(4,:), X_smooth1(5,:), X_smooth1(6,:), X_smooth1(7,:), start, stop, dt, omega_on_off);

    %Run Extended Kalman filter + smoother 
    [X_smooth2] = RPY_body(Y2,dt);
    

    %Save smooth data (Note that the obtained omega's have to be converted
    %to the body frame of reference)
    omega_1(start:stop,j) = X_smooth2(1,:);
    omega_2(start:stop,j) = X_smooth2(2,:);
    omega_3(start:stop,j) = X_smooth2(3,:);
    q1(start:stop,j) = X_smooth2(4,:);
    q2(start:stop,j) = X_smooth2(5,:);
    q3(start:stop,j) = X_smooth2(6,:);
    q4(start:stop,j) = X_smooth2(7,:);

    
%     wb_raw_filt(X_smooth1(4,:), X_smooth1(5,:), X_smooth1(6,:), X_smooth1(7,:), X_smooth1(1,:), X_smooth1(2,:), X_smooth1(3,:), pathDB.t(start:stop), 1, 2)
%     pause
   
    clear Y1 Y2 start stop X_smooth1 X_smooth2 X_s Y1 q1_raw q2_raw q3_raw q4_raw q1_filt q2_filt q3_filt q4_filt
    
end

        temp.x = x;
        temp.y = y;
        temp.z = z;
        
        temp.u = u;
        temp.v = v;
        temp.w = w;

        temp.ax = ax;
        temp.ay = ay;
        temp.az = az;

        temp.q1 = q1;
        temp.q2 = q2;
        temp.q3 = q3;
        temp.q4 = q4;
        
        temp.omega_1 = omega_1;
        temp.omega_2 = omega_2;
        temp.omega_3 = omega_3;

% save(savefile,'t','x','y','z','u','v','w','ax','ay','az','q1','q2','q3','q4','omega_1','omega_2','omega_3')

