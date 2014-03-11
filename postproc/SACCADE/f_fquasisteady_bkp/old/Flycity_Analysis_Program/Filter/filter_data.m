function filter_data( settings, pathDB )

    % Use a Linear Kalman Filter to filter the body position data and an
    % Extended Kalman Filter to filter body and wing orientation.
    
    savefile = 'pathDB2.mat';
    
    nr_frames = settings.frame_end;
    nr_seq = settings.nr_of_seq;
    
    
    xyz =       nan(nr_frames,3,nr_seq);
    uvw =       nan(nr_frames,3,nr_seq);
    a_xyz =     nan(nr_frames,3,nr_seq);
    
    qB =        nan(nr_frames,4,nr_seq);
    wB =        nan(nr_frames,3,nr_seq);
    
    qL =        nan(nr_frames,4,nr_seq);
    wL =        nan(nr_frames,3,nr_seq);
    
    qR =        nan(nr_frames,4,nr_seq);
    wR =        nan(nr_frames,3,nr_seq);
    
    % raw data:
    
    dt = pathDB.dt;
    
    xyz_raw =   pathDB.raw.xyz;
    qB_raw =    pathDB.raw.qB;
    qL_raw =    pathDB.raw.qL;
    qR_raw =    pathDB.raw.qR;
    
    % Obtain filter setttings:
    
    Q_xyz    = settings.filter_set.xyz;
    Q_rpy_b1 = settings.filter_set.qB1;
    Q_rpy_b2 = settings.filter_set.qB2;
    Q_rpy_w1 = settings.filter_set.qW1;
    Q_rpy_w2 = settings.filter_set.qW2;
    
    % Filter body position data:
    
    for i = 1:nr_seq
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
                
        Y = (xyz_raw(start:stop,:,i))';
        
        [Xs] = XYZ_body(Y,dt,Q_xyz);
        
        xyz(start:stop,:,i) = [Xs(1,:); Xs(2,:); Xs(3,:)]';
        uvw(start:stop,:,i) = [Xs(4,:); Xs(5,:); Xs(6,:)]';
        a_xyz(start:stop,:,i) = [Xs(7,:); Xs(8,:); Xs(9,:)]';
        
    end
    
    clear Xs Y start stop
    
    % Filter body orientation data:
    
    for i = 1:nr_seq
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
                
        N = stop-start+1;
        
        % Filter round 1:
        
        omega_on_off = 0;
        
        q1_raw = qB_raw(start:stop,1,i);
        q2_raw = qB_raw(start:stop,2,i);
        q3_raw = qB_raw(start:stop,3,i);
        q4_raw = qB_raw(start:stop,4,i);
      
        [Y1] = MakeY(q1_raw, q2_raw, q3_raw, q4_raw, N, dt, omega_on_off);
        
        [Xs1] = RPY_body(Y1,dt,Q_rpy_b1,Q_rpy_b2);
        
        % Filter round 2:
        
        omega_on_off = 1;
        
        [Y2] = MakeY(Xs1(4,:), Xs1(5,:), Xs1(6,:), Xs1(7,:), N, dt, omega_on_off);
        
        [Xs2] = RPY_body(Y2,dt,Q_rpy_b1,Q_rpy_b2);
        
        qB(start:stop,:,i) = [Xs2(4,:); Xs2(5,:); Xs2(6,:); Xs2(7,:)]';
        wB(start:stop,:,i) = [Xs2(1,:); Xs2(2,:); Xs2(3,:)]';
        
    end
    
    clear Xs1 Y1 Xs2 Y2 start stop q1_raw q2_raw q3_raw q4_raw
    
    
    % Filter wing orientation data:
    
    for i = 1:nr_seq
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
                
        N = stop-start+1;
        
        % Filter round 1:
        
        omega_on_off = 0;
        
        qL1_raw = qL_raw(start:stop,1,i);
        qL2_raw = qL_raw(start:stop,2,i);
        qL3_raw = qL_raw(start:stop,3,i);
        qL4_raw = qL_raw(start:stop,4,i);
        
        qR1_raw = qR_raw(start:stop,1,i);
        qR2_raw = qR_raw(start:stop,2,i);
        qR3_raw = qR_raw(start:stop,3,i);
        qR4_raw = qR_raw(start:stop,4,i);
        
        [YL1, YR1] = MakeYL_YR( qL1_raw, qL2_raw, qL3_raw, qL4_raw, qR1_raw, qR2_raw, qR3_raw, qR4_raw, N, omega_on_off, dt);
               
        %Run Extended Kalman filter + smoother for the left wing.
        
        [xs_L1] = RPY_wing(YL1,dt,Q_rpy_w1,Q_rpy_w2);
        
        %Run Extended Kalman filter + smoother for the right wing.
        
        [xs_R1] = RPY_wing(YR1,dt,Q_rpy_w1,Q_rpy_w2);
        
        % Filter round 2:
        
        omega_on_off = 1;
        
        [YL2, YR2] = MakeYL_YR( xs_L1(4,:), xs_L1(5,:), xs_L1(6,:), xs_L1(7,:), xs_R1(4,:), xs_R1(5,:), xs_R1(6,:), xs_R1(7,:), N, omega_on_off, dt);
               
        %Run Extended Kalman filter + smoother for the left wing.
        
        [xs_L2] = RPY_wing(YL2,dt,Q_rpy_w1,Q_rpy_w2);
        
        %Run Extended Kalman filter + smoother for the right wing.
        
        [xs_R2] = RPY_wing(YR2,dt,Q_rpy_w1,Q_rpy_w2);
        
        qL(start:stop,:,i) = [xs_L2(4,:); xs_L2(5,:); xs_L2(6,:); xs_L2(7,:)]';
        wL(start:stop,:,i) = [xs_L2(1,:); xs_L2(2,:); xs_L2(3,:)]';
        
        qR(start:stop,:,i) = [xs_R2(4,:); xs_R2(5,:); xs_R2(6,:); xs_R2(7,:)]';
        wR(start:stop,:,i) = [xs_R2(1,:); xs_R2(2,:); xs_R2(3,:)]';
        
    end
    
    save(savefile,'xyz','uvw','a_xyz','qB','wB','qL','wL','qR','wR')

end

