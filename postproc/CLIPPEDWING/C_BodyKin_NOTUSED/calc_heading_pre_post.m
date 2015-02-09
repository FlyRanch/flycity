
for i = 1:length(stim_angle_vel_post_turn)
    %% heading
    if isnan(stim_angle_vel_pre_resp(i)) == 0
        heading_pre(i,1) = stim_angle_vel_pre_resp(i);
    elseif isnan(stim_angle_vel_pre_turn(i)) == 0
        heading_pre(i,1) = stim_angle_vel_pre_turn(i);
    else
        heading_temp = stim_angle_vel(:,i);
        heading_pre(i,1) = heading_temp(t==0);
        if isnan(heading_pre(i)) == 1
            heading_temp = heading_temp(isnan(heading_temp)==0);
            heading_pre(i,1) = heading_temp(1);
        end
    end
    
    if isnan(stim_angle_vel_post_resp(i)) == 0
        heading_post(i,1) = stim_angle_vel_post_resp(i);
    elseif isnan(stim_angle_vel_post_turn(i)) == 0
        heading_post(i,1) = stim_angle_vel_post_turn(i);
    else
        heading_temp = stim_angle_vel(:,i);
        heading_temp = heading_temp(isnan(heading_temp)==0);
        heading_post(i,1) = heading_temp(end);
    end
    
    %% A direction
    % pre
    if isnan(stim_angle_accel_pre_resp(i)) == 0
        Adir_pre(i,1) = stim_angle_accel_pre_resp(i);
    elseif isnan(stim_angle_accel_pre_turn(i)) == 0
        Adir_pre(i,1) = stim_angle_accel_pre_turn(i);
%     elseif isnan(stim_angle_accel_pre_accel(i)) == 0
%         Adir_pre(i,1) = stim_angle_accel_pre_accel(i);
    else
        Adir_pre(i,1)=nan;
    end
    
    % post
    if isnan(stim_angle_accel_post_resp(i)) == 0
        Adir_post(i,1) = stim_angle_accel_post_resp(i);
    elseif isnan(stim_angle_accel_post_turn(i)) == 0
        Adir_post(i,1) = stim_angle_accel_post_turn(i);
%     elseif isnan(stim_angle_accel_post_accel(i)) == 0
%         Adir_post(i,1) = stim_angle_accel_post_accel(i);
    else
        Adir_temp = stim_angle_accel(:,i);
        Adir_temp = Adir_temp(isnan(Adir_temp)==0);
        Adir_post(i,1) = Adir_temp(end);
    end
    
    % mean
    if isnan(stim_angle_accel_resp_mean(i)) == 0
        Adir_mean(i,1) = stim_angle_accel_resp_mean(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(stim_angle_accel(n_pre(i):n_post(i),i));
        Adir_mean(i,1) = mu;
    elseif isnan(n_pre(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(stim_angle_accel(n_pre(i):end,i));
        Adir_mean(i,1) = mu;
    else
        Adir_mean(i,1) = nan;
    end
    
    % std
    if isnan(stim_angle_accel_resp_std(i)) == 0
        Adir_std(i,1) = stim_angle_accel_resp_std(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        
        [s s0] = circ_std_deg_nonan(stim_angle_accel(n_pre(i):n_post(i),i));
        Adir_std(i,1) = s0;
    elseif isnan(n_pre(i)) == 0
        
        [s s0] = circ_std_deg_nonan(stim_angle_accel(n_pre(i):end,i));
        Adir_std(i,1) = s0;
    else
        Adir_std(i,1) = nan;
    end
    
    % min
    if isnan(stim_angle_accel_resp_min(i)) == 0
        Adir_min(i,1) = stim_angle_accel_resp_min(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        Adir_min(i,1) = min(stim_angle_accel(n_pre(i):n_post(i),i));
    elseif isnan(n_pre(i)) == 0
        Adir_min(i,1) = min(stim_angle_accel(n_pre(i):end,i));
    else
        Adir_min(i,1) = nan;
    end
    
    % max
    if isnan(stim_angle_accel_resp_max(i)) == 0
        Adir_max(i,1) = stim_angle_accel_resp_max(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        Adir_max(i,1) = max(stim_angle_accel(n_pre(i):n_post(i),i));
    elseif isnan(n_pre(i)) == 0
        Adir_max(i,1) = max(stim_angle_accel(n_pre(i):end,i));
    else
        Adir_max(i,1) = nan;
    end

    %% A norm
    % pre
    if isnan(A_hor_pre_resp(i)) == 0
        A_pre(i,1) = A_hor_pre_resp(i);
    elseif isnan(A_hor_pre_turn(i)) == 0
        A_pre(i,1) = A_hor_pre_turn(i);
    elseif isnan(A_hor_pre_accel(i)) == 0
        A_pre(i,1) = A_hor_pre_accel(i);
    else
        A_pre(i,1)=nan;
    end
    
    % post
    if isnan(A_hor_post_resp(i)) == 0
        A_post(i,1) = A_hor_post_resp(i);
    elseif isnan(A_hor_post_turn(i)) == 0
        A_post(i,1) = A_hor_post_turn(i);
    elseif isnan(A_hor_post_accel(i)) == 0
        A_post(i,1) = A_hor_post_accel(i);
    else
        A_temp = A_hor(:,i);
        A_temp = A_temp(isnan(A_temp)==0);
        A_post(i,1) = A_temp(end);
    end
    

    %% A direction relative to body
    % pre
    if isnan(accel_angle_hor_body_pre_resp(i)) == 0
        Adir_body_pre(i,1) = accel_angle_hor_body_pre_resp(i);
    elseif isnan(accel_angle_hor_body_pre_turn(i)) == 0
        Adir_body_pre(i,1) = accel_angle_hor_body_pre_turn(i);
    elseif isnan(accel_angle_hor_body_pre_accel(i)) == 0
        Adir_body_pre(i,1) = accel_angle_hor_body_pre_accel(i);
    else
        Adir_body_pre(i,1) = nan;
    end
    
    % post
    if isnan(accel_angle_hor_body_post_resp(i)) == 0
        Adir_body_post(i,1) = accel_angle_hor_body_post_resp(i);
    elseif isnan(accel_angle_hor_body_post_turn(i)) == 0
        Adir_body_post(i,1) = accel_angle_hor_body_post_turn(i);
    elseif isnan(accel_angle_hor_body_post_accel(i)) == 0
        Adir_body_post(i,1) = accel_angle_hor_body_post_accel(i);
    else
        Adir_body_temp = accel_angle_hor_body(:,i);
        Adir_body_temp = Adir_body_temp(isnan(Adir_body_temp)==0);
        Adir_body_post(i,1) = Adir_body_temp(end);
    end
    
    % mean
    if isnan(accel_angle_hor_body_resp_mean(i)) == 0
        Adir_body_mean(i,1) = accel_angle_hor_body_resp_mean(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(accel_angle_hor_body(n_pre(i):n_post(i),i));
        Adir_body_mean(i,1) = mu;
    elseif isnan(n_pre(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(accel_angle_hor_body(n_pre(i):end,i));
        Adir_body_mean(i,1) = mu;
    else
        Adir_body_mean(i,1) = nan;
    end
    
    %% A direction relative to V_hor
    % pre
    if isnan(accel_angle_hor_vel_pre_resp(i)) == 0
        Adir_vel_pre(i,1) = accel_angle_hor_vel_pre_resp(i);
    elseif isnan(accel_angle_hor_vel_pre_turn(i)) == 0
        Adir_vel_pre(i,1) = accel_angle_hor_vel_pre_turn(i);
    elseif isnan(accel_angle_hor_vel_pre_accel(i)) == 0
        Adir_vel_pre(i,1) = accel_angle_hor_vel_pre_accel(i);
    else
        Adir_vel_pre(i,1) = nan;
    end
    
    % post
    if isnan(accel_angle_hor_vel_post_resp(i)) == 0
        Adir_vel_post(i,1) = accel_angle_hor_vel_post_resp(i);
    elseif isnan(accel_angle_hor_vel_post_turn(i)) == 0
        Adir_vel_post(i,1) = accel_angle_hor_vel_post_turn(i);
    elseif isnan(accel_angle_hor_vel_post_accel(i)) == 0
        Adir_vel_post(i,1) = accel_angle_hor_vel_post_accel(i);
    else
        Adir_vel_temp = accel_angle_hor_vel(:,i);
        Adir_vel_temp = Adir_vel_temp(isnan(Adir_vel_temp)==0);
        Adir_vel_post(i,1) = Adir_vel_temp(end);
    end
    
    % mean
    if isnan(accel_angle_hor_vel_resp_mean(i)) == 0
        Adir_vel_mean(i,1) = accel_angle_hor_vel_resp_mean(i);
    elseif isnan(n_pre(i)) == 0 && isnan(n_post(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(accel_angle_hor_vel(n_pre(i):n_post(i),i));
        Adir_vel_mean(i,1) = mu;
    elseif isnan(n_pre(i)) == 0
        
        [mu ul ll] = circ_mean_deg_nonan(accel_angle_hor_vel(n_pre(i):end,i));
        Adir_vel_mean(i,1) = mu;
    else
        Adir_vel_mean(i,1) = nan;
    end
    
    
    %% yaw
    if isnan(stim_angle_yaw_pre_resp(i)) == 0
        yaw_pre(i,1) = stim_angle_yaw_pre_resp(i);
    elseif isnan(stim_angle_yaw_pre_turn(i)) == 0
        yaw_pre(i,1) = stim_angle_yaw_pre_turn(i);
    else
        yaw_temp = stim_angle_yaw(:,i);
        yaw_pre(i,1) = yaw_temp(t==0);
        if isnan(yaw_pre(i)) == 1
            yaw_temp = yaw_temp(isnan(yaw_temp)==0);
            yaw_pre(i,1) = yaw_temp(1);
        end
    end
    
    if isnan(stim_angle_yaw_post_turn(i)) == 0
        yaw_post(i,1) = stim_angle_yaw_post_turn(i);
    elseif isnan(stim_angle_yaw_post_resp(i)) == 0
        yaw_post(i,1) = stim_angle_yaw_post_resp(i);
    else
        yaw_temp = stim_angle_yaw(:,i);
        yaw_temp = yaw_temp(isnan(yaw_temp)==0);
        yaw_post(i,1) = yaw_temp(end);
    end
    
    %% slip
    if isnan(slip_pre_resp(i)) == 0
        slip_pre(i,1) = slip_pre_resp(i);
    elseif isnan(slip_pre_turn(i)) == 0
        slip_pre(i,1) = slip_pre_turn(i);
    else
        slip_temp = slip(:,i);
        slip_pre(i,1) = slip_temp(t==0);
        if isnan(slip_pre(i)) == 1
            slip_temp = slip_temp(isnan(slip_temp)==0);
            slip_pre(i,1) = slip_temp(1);
        end
    end
    
    if isnan(slip_post_turn(i)) == 0
        slip_post(i,1) = slip_post_turn(i);
    elseif isnan(slip_post_resp(i)) == 0
        slip_post(i,1) = slip_post_resp(i);
    else
        slip_temp = slip(:,i);
        slip_temp = slip_temp(isnan(slip_temp)==0);
        slip_post(i,1) = slip_temp(end);
    end
    
    %% pitch
    if isnan(pitch_pre_resp(i)) == 0
        pitch_pre(i,1) = pitch_pre_resp(i);
    elseif isnan(pitch_pre_turn(i)) == 0
        pitch_pre(i,1) = pitch_pre_turn(i);
    else
        pitch_temp = pitch(:,i);
        pitch_pre(i,1) = pitch_temp(t==0);
        if isnan(pitch_pre(i)) == 1
            pitch_temp = pitch_temp(isnan(pitch_temp)==0);
            pitch_pre(i,1) = pitch_temp(1);
        end
    end
    
    if isnan(pitch_post_turn(i)) == 0
        pitch_post(i,1) = pitch_post_turn(i);
    elseif isnan(pitch_post_resp(i)) == 0
        pitch_post(i,1) = pitch_post_resp(i);
    else
        pitch_temp = pitch(:,i);
        pitch_temp = pitch_temp(isnan(pitch_temp)==0);
        pitch_post(i,1) = pitch_temp(end);
    end
    
    %% roll
    if isnan(roll_pre_resp(i)) == 0
        roll_pre(i,1) = roll_pre_resp(i);
    elseif isnan(roll_pre_turn(i)) == 0
        roll_pre(i,1) = roll_pre_turn(i);
    else
        roll_temp = roll(:,i);
        roll_pre(i,1) = roll_temp(t==0);
        if isnan(roll_pre(i)) == 1
            roll_temp = roll_temp(isnan(roll_temp)==0);
            roll_pre(i,1) = roll_temp(1);
        end
    end
    
    if isnan(roll_post_turn(i)) == 0
        roll_post(i,1) = roll_post_turn(i);
    elseif isnan(roll_post_resp(i)) == 0
        roll_post(i,1) = roll_post_resp(i);
    else
        roll_temp = roll(:,i);
        roll_temp = roll_temp(isnan(roll_temp)==0);
        roll_post(i,1) = roll_temp(end);
    end
    
    %% V
    if isnan(V_pre_resp(i)) == 0
        V_pre(i,1) = V_pre_resp(i);
    elseif isnan(V_pre_accel(i)) == 0
        V_pre(i,1) = V_pre_accel(i);
    else
        V_temp = V(:,i);
        V_pre(i,1) = V_temp(t==0);
        if isnan(V_pre(i)) == 1
            V_temp = V_temp(isnan(V_temp)==0);
            V_pre(i,1) = V_temp(1);
        end
    end
    
    if isnan(V_post_resp(i)) == 0
        V_post(i,1) = V_post_resp(i);
    elseif isnan(V_post_accel(i)) == 0
        V_post(i,1) = V_post_accel(i);
    else
        V_temp = V(:,i);
        V_temp = V_temp(isnan(V_temp)==0);
        V_post(i,1) = V_temp(end);
    end
end
