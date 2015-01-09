function [ output_args ] = filter_wing(settings, pathDB)


    savefile = 'pathDB4.mat';

    % Program that contains a Kalman filter which filters the measured wing
    % quaternions of the left and right wing. The output will be the filtered
    % wing quaternions.

    dt = pathDB.t(2)-pathDB.t(1);
    t = pathDB.t;

    % batch loop body
    qL1_filt2 = nan(size(pathDB.x));
    qL2_filt2 = nan(size(pathDB.x));
    qL3_filt2 = nan(size(pathDB.x));
    qL4_filt2 = nan(size(pathDB.x));
    
    qR1_filt2 = nan(size(pathDB.x));
    qR2_filt2 = nan(size(pathDB.x));
    qR3_filt2 = nan(size(pathDB.x));
    qR4_filt2 = nan(size(pathDB.x));
    
    omega1_L = nan(size(pathDB.x));
    omega2_L = nan(size(pathDB.x));
    omega3_L = nan(size(pathDB.x));
    
    omega1_R = nan(size(pathDB.x));
    omega2_R = nan(size(pathDB.x));
    omega3_R =  nan(size(pathDB.x));
    
    
    
    for i=1:size(pathDB.x,2)
        
        % start and stop point for the measurements
        start = find(isnan(pathDB.x(:,i))==0, 1 );
        stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
        
        % Remove NaN values from input
        qL1 = pathDB.qL1_filt1(start:stop,i);
        qL2 = pathDB.qL2_filt1(start:stop,i);
        qL3 = pathDB.qL3_filt1(start:stop,i);
        qL4 = pathDB.qL4_filt1(start:stop,i);
               
        qR1 = pathDB.qR1_filt1(start:stop,i);
        qR2 = pathDB.qR2_filt1(start:stop,i);
        qR3 = pathDB.qR3_filt1(start:stop,i);
        qR4 = pathDB.qR4_filt1(start:stop,i);

        
        % Filter quaternion data such that there will be no jumps between
        % positive and negative quaternions.
        
        omega_on_off = 0;
        
        [YL1, YR1] = MakeYL_YR( qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4, start, stop, omega_on_off, dt);
               
        %Run Extended Kalman filter + smoother for the left wing.
        
        [xs_L1] = RPY_wing(YL1,dt);
        
        %Run Extended Kalman filter + smoother for the right wing.
        
        [xs_R1] = RPY_wing(YR1,dt);
        
        
        % Now use the smoothened values of xs_L1 and xs_R1 to include the
        % omega values in the filter:
        
        omega_on_off = 1;
        
        [YL2, YR2] = MakeYL_YR( xs_L1(4,:), xs_L1(5,:), xs_L1(6,:), xs_L1(7,:), xs_R1(4,:), xs_R1(5,:), xs_R1(6,:), xs_R1(7,:), start, stop, omega_on_off, dt);
        
        %Run Extended Kalman filter + smoother for the left wing.
        
        [xs_L2] = RPY_wing(YL2,dt);
        
        %Run Extended Kalman filter + smoother for the right wing.
        
        [xs_R2] = RPY_wing(YR2,dt);

        
        
        
        % Store filtered results
        omega1_L(start:stop,i) = xs_L2(1,:);
        omega2_L(start:stop,i) = xs_L2(2,:);
        omega3_L(start:stop,i) = xs_L2(3,:);
        qL1_filt2(start:stop,i) = xs_L2(4,:);
        qL2_filt2(start:stop,i) = xs_L2(5,:);
        qL3_filt2(start:stop,i) = xs_L2(6,:);
        qL4_filt2(start:stop,i) = xs_L2(7,:);
        
        omega1_R(start:stop,i) = xs_R2(1,:);
        omega2_R(start:stop,i) = xs_R2(2,:);
        omega3_R(start:stop,i) = xs_R2(3,:);
        qR1_filt2(start:stop,i) = xs_R2(4,:);
        qR2_filt2(start:stop,i) = xs_R2(5,:);
        qR3_filt2(start:stop,i) = xs_R2(6,:);
        qR4_filt2(start:stop,i) = xs_R2(7,:);
        
        clear YL YR start stop xs_L xs_R qL1 qL2 qL3 qL4 qR1 qR2 qR3 qR4
    end

        
   save(savefile,'omega1_L', 'omega2_L', 'omega3_L', 'qL1_filt2', 'qL2_filt2', 'qL3_filt2', 'qL4_filt2', 'omega1_R', 'omega2_R', 'omega3_R', 'qR1_filt2', 'qR2_filt2', 'qR3_filt2', 'qR4_filt2') 

end

