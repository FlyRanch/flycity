
    %% pre
    if sum(isnan(f_wb_seq_pre(:,wb_now))==0) >= wb_min4mean_clusters    % only mean if Nwb>wb_min4mean

        % mean
        % wb kin
        f_wb_mean_now = nanmean(f_wb_seq_pre(:,wb_now)')';
        
        stroke_wb_L_mean_now = nanmean(stroke_wb_L_seq_bins_pre(:,:,wb_now)')';
        stroke_wb_R_mean_now = nanmean(stroke_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        dev_wb_L_mean_now = nanmean(dev_wb_L_seq_bins_pre(:,:,wb_now)')';
        dev_wb_R_mean_now = nanmean(dev_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        pitch_wb_L_mean_now = nanmean(pitch_wb_L_seq_bins_pre(:,:,wb_now)')';
        pitch_wb_R_mean_now = nanmean(pitch_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        % body kin
        roll_mean_wb_mean_now = nanmean(roll_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_mean_wb_mean_now = nanmean(roll_dot_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_dot_mean_wb_mean_now = nanmean(roll_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        pitch_mean_wb_mean_now = nanmean(pitch_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_mean_wb_mean_now = nanmean(pitch_dot_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_dot_mean_wb_mean_now = nanmean(pitch_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        yaw_mean_wb_mean_now = nanmean(yaw_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_mean_wb_mean_now = nanmean(yaw_dot_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_dot_mean_wb_mean_now = nanmean(yaw_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
%         drot_L_mean_wb_mean_now = nanmean(drot_L_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_L_mean_wb_mean_now = nanmean(rot_dot_L_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_dot_L_mean_wb_mean_now = nanmean(rot_dot_dot_L_mean_wb_seq_pre(:,wb_now)')';
%         
%         drot_R_mean_wb_mean_now = nanmean(drot_R_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_R_mean_wb_mean_now = nanmean(rot_dot_R_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_dot_R_mean_wb_mean_now = nanmean(rot_dot_dot_R_mean_wb_seq_pre(:,wb_now)')';
%         
%         % torque
%         Mroll_mean_wb_accel_mean_now = nanmean(Mroll_mean_wb_accel_seq_pre(:,wb_now)')';
%         Mroll_mean_wb_damp_mean_now = nanmean(Mroll_mean_wb_damp_seq_pre(:,wb_now)')';
%         Mroll_mean_wb_mean_now = nanmean(Mroll_mean_wb_seq_pre(:,wb_now)')';
%         
%         Mpitch_mean_wb_accel_mean_now = nanmean(Mpitch_mean_wb_accel_seq_pre(:,wb_now)')';
%         Mpitch_mean_wb_damp_mean_now = nanmean(Mpitch_mean_wb_damp_seq_pre(:,wb_now)')';
%         Mpitch_mean_wb_mean_now = nanmean(Mpitch_mean_wb_seq_pre(:,wb_now)')';
%         
%         Myaw_mean_wb_accel_mean_now = nanmean(Myaw_mean_wb_accel_seq_pre(:,wb_now)')';
%         Myaw_mean_wb_damp_mean_now = nanmean(Myaw_mean_wb_damp_seq_pre(:,wb_now)')';
%         Myaw_mean_wb_mean_now = nanmean(Myaw_mean_wb_seq_pre(:,wb_now)')';
%         
%         M_L_mean_wb_accel_mean_now = nanmean(M_L_mean_wb_accel_seq_pre(:,wb_now)')';
%         M_L_mean_wb_damp_mean_now = nanmean(M_L_mean_wb_damp_seq_pre(:,wb_now)')';
%         M_L_mean_wb_mean_now = nanmean(M_L_mean_wb_seq_pre(:,wb_now)')';
%         
%         M_R_mean_wb_accel_mean_now = nanmean(M_R_mean_wb_accel_seq_pre(:,wb_now)')';
%         M_R_mean_wb_damp_mean_now = nanmean(M_R_mean_wb_damp_seq_pre(:,wb_now)')';
%         M_R_mean_wb_mean_now = nanmean(M_R_mean_wb_seq_pre(:,wb_now)')';
        
        %% ste
        n_now = sum(isnan(f_wb_seq_pre(:,wb_now))==0);
        f_wb_ste_now = nanstd(f_wb_seq_pre(:,wb_now)')'/sqrt(n_now);

        stroke_wb_L_ste_now = nanstd(stroke_wb_L_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        stroke_wb_R_ste_now = nanstd(stroke_wb_R_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        
        dev_wb_L_ste_now = nanstd(dev_wb_L_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        dev_wb_R_ste_now = nanstd(dev_wb_R_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        
        pitch_wb_L_ste_now = nanstd(pitch_wb_L_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        pitch_wb_R_ste_now = nanstd(pitch_wb_R_seq_bins_pre(:,:,wb_now)')'/sqrt(n_now);
        
        % body kin
        roll_mean_wb_ste_now = nanstd(roll_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        roll_dot_mean_wb_ste_now = nanstd(roll_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        roll_dot_dot_mean_wb_ste_now = nanstd(roll_dot_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        
        pitch_mean_wb_ste_now = nanstd(pitch_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        pitch_dot_mean_wb_ste_now = nanstd(pitch_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        pitch_dot_dot_mean_wb_ste_now = nanstd(pitch_dot_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        
        yaw_mean_wb_ste_now = nanstd(yaw_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        yaw_dot_mean_wb_ste_now = nanstd(yaw_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        yaw_dot_dot_mean_wb_ste_now = nanstd(yaw_dot_dot_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        
%         drot_L_mean_wb_ste_now = nanstd(drot_L_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         rot_dot_L_mean_wb_ste_now = nanstd(rot_dot_L_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         rot_dot_dot_L_mean_wb_ste_now = nanstd(rot_dot_dot_L_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         drot_R_mean_wb_ste_now = nanstd(drot_R_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         rot_dot_R_mean_wb_ste_now = nanstd(rot_dot_R_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         rot_dot_dot_R_mean_wb_ste_now = nanstd(rot_dot_dot_R_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         % torque
%         Mroll_mean_wb_accel_ste_now = nanstd(Mroll_mean_wb_accel_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Mroll_mean_wb_damp_ste_now = nanstd(Mroll_mean_wb_damp_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Mroll_mean_wb_ste_now = nanstd(Mroll_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         Mpitch_mean_wb_accel_ste_now = nanstd(Mpitch_mean_wb_accel_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Mpitch_mean_wb_damp_ste_now = nanstd(Mpitch_mean_wb_damp_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Mpitch_mean_wb_ste_now = nanstd(Mpitch_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         Myaw_mean_wb_accel_ste_now = nanstd(Myaw_mean_wb_accel_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Myaw_mean_wb_damp_ste_now = nanstd(Myaw_mean_wb_damp_seq_pre(:,wb_now)')'/sqrt(n_now);
%         Myaw_mean_wb_ste_now = nanstd(Myaw_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         M_L_mean_wb_accel_ste_now = nanstd(M_L_mean_wb_accel_seq_pre(:,wb_now)')'/sqrt(n_now);
%         M_L_mean_wb_damp_ste_now = nanstd(M_L_mean_wb_damp_seq_pre(:,wb_now)')'/sqrt(n_now);
%         M_L_mean_wb_ste_now = nanstd(M_L_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
%         
%         M_R_mean_wb_accel_ste_now = nanstd(M_R_mean_wb_accel_seq_pre(:,wb_now)')'/sqrt(n_now);
%         M_R_mean_wb_damp_ste_now = nanstd(M_R_mean_wb_damp_seq_pre(:,wb_now)')'/sqrt(n_now);
%         M_R_mean_wb_ste_now = nanstd(M_R_mean_wb_seq_pre(:,wb_now)')'/sqrt(n_now);
        
        %% all
        % wb kin
        f_wb_now = (f_wb_seq_pre(:,wb_now)')';
        
        stroke_wb_L_now = (stroke_wb_L_seq_bins_pre(:,:,wb_now)')';
        stroke_wb_R_now = (stroke_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        dev_wb_L_now = (dev_wb_L_seq_bins_pre(:,:,wb_now)')';
        dev_wb_R_now = (dev_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        pitch_wb_L_now = (pitch_wb_L_seq_bins_pre(:,:,wb_now)')';
        pitch_wb_R_now = (pitch_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        % body kin
        roll_mean_wb_now = (roll_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_mean_wb_now = (roll_dot_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_dot_mean_wb_now = (roll_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        pitch_mean_wb_now = (pitch_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_mean_wb_now = (pitch_dot_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_dot_mean_wb_now = (pitch_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        yaw_mean_wb_now = (yaw_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_mean_wb_now = (yaw_dot_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_dot_mean_wb_now = (yaw_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
%         drot_L_mean_wb_now = (drot_L_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_L_mean_wb_now = (rot_dot_L_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_dot_L_mean_wb_now = (rot_dot_dot_L_mean_wb_seq_pre(:,wb_now)')';
%         
%         drot_R_mean_wb_now = (drot_R_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_R_mean_wb_now = (rot_dot_R_mean_wb_seq_pre(:,wb_now)')';
%         rot_dot_dot_R_mean_wb_now = (rot_dot_dot_R_mean_wb_seq_pre(:,wb_now)')';
%         
%         % torque
%         Mroll_mean_wb_accel_now = (Mroll_mean_wb_accel_seq_pre(:,wb_now)')';
%         Mroll_mean_wb_damp_now = (Mroll_mean_wb_damp_seq_pre(:,wb_now)')';
%         Mroll_mean_wb_now = (Mroll_mean_wb_seq_pre(:,wb_now)')';
%         
%         Mpitch_mean_wb_accel_now = (Mpitch_mean_wb_accel_seq_pre(:,wb_now)')';
%         Mpitch_mean_wb_damp_now = (Mpitch_mean_wb_damp_seq_pre(:,wb_now)')';
%         Mpitch_mean_wb_now = (Mpitch_mean_wb_seq_pre(:,wb_now)')';
%         
%         Myaw_mean_wb_accel_now = (Myaw_mean_wb_accel_seq_pre(:,wb_now)')';
%         Myaw_mean_wb_damp_now = (Myaw_mean_wb_damp_seq_pre(:,wb_now)')';
%         Myaw_mean_wb_now = (Myaw_mean_wb_seq_pre(:,wb_now)')';
%         
%         M_L_mean_wb_accel_now = (M_L_mean_wb_accel_seq_pre(:,wb_now)')';
%         M_L_mean_wb_damp_now = (M_L_mean_wb_damp_seq_pre(:,wb_now)')';
%         M_L_mean_wb_now = (M_L_mean_wb_seq_pre(:,wb_now)')';
%         
%         M_R_mean_wb_accel_now = (M_R_mean_wb_accel_seq_pre(:,wb_now)')';
%         M_R_mean_wb_damp_now = (M_R_mean_wb_damp_seq_pre(:,wb_now)')';
%         M_R_mean_wb_now = (M_R_mean_wb_seq_pre(:,wb_now)')';
        
        %% add to list
        
        %% mean
        % wb kin
        f_wb_seq_pre_mean = [f_wb_seq_pre_mean f_wb_mean_now];
        
        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean stroke_wb_L_mean_now(1:end-1)];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean stroke_wb_R_mean_now(1:end-1)];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean dev_wb_L_mean_now(1:end-1)];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean dev_wb_R_mean_now(1:end-1)];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean pitch_wb_L_mean_now(1:end-1)];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean pitch_wb_R_mean_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean = [roll_mean_wb_seq_pre_mean roll_mean_wb_mean_now];
        roll_dot_mean_wb_seq_pre_mean = [roll_dot_mean_wb_seq_pre_mean roll_dot_mean_wb_mean_now];
        roll_dot_dot_mean_wb_seq_pre_mean = [roll_dot_dot_mean_wb_seq_pre_mean roll_dot_dot_mean_wb_mean_now];
        
        pitch_mean_wb_seq_pre_mean = [pitch_mean_wb_seq_pre_mean pitch_mean_wb_mean_now];
        pitch_dot_mean_wb_seq_pre_mean = [pitch_dot_mean_wb_seq_pre_mean pitch_dot_mean_wb_mean_now];
        pitch_dot_dot_mean_wb_seq_pre_mean = [pitch_dot_dot_mean_wb_seq_pre_mean pitch_dot_dot_mean_wb_mean_now];
        
        yaw_mean_wb_seq_pre_mean = [yaw_mean_wb_seq_pre_mean yaw_mean_wb_mean_now];
        yaw_dot_mean_wb_seq_pre_mean = [yaw_dot_mean_wb_seq_pre_mean yaw_dot_mean_wb_mean_now];
        yaw_dot_dot_mean_wb_seq_pre_mean = [yaw_dot_dot_mean_wb_seq_pre_mean yaw_dot_dot_mean_wb_mean_now];
        
%         drot_L_mean_wb_seq_pre_mean = [drot_L_mean_wb_seq_pre_mean drot_L_mean_wb_mean_now];
%         rot_dot_L_mean_wb_seq_pre_mean = [rot_dot_L_mean_wb_seq_pre_mean rot_dot_L_mean_wb_mean_now];
%         rot_dot_dot_L_mean_wb_seq_pre_mean = [rot_dot_dot_L_mean_wb_seq_pre_mean rot_dot_dot_L_mean_wb_mean_now];
%         
%         drot_R_mean_wb_seq_pre_mean = [drot_R_mean_wb_seq_pre_mean drot_R_mean_wb_mean_now];
%         rot_dot_R_mean_wb_seq_pre_mean = [rot_dot_R_mean_wb_seq_pre_mean rot_dot_R_mean_wb_mean_now];
%         rot_dot_dot_R_mean_wb_seq_pre_mean = [rot_dot_dot_R_mean_wb_seq_pre_mean rot_dot_dot_R_mean_wb_mean_now];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_mean = [Mroll_mean_wb_accel_seq_pre_mean Mroll_mean_wb_accel_mean_now];
%         Mroll_mean_wb_damp_seq_pre_mean = [Mroll_mean_wb_damp_seq_pre_mean Mroll_mean_wb_damp_mean_now];
%         Mroll_mean_wb_seq_pre_mean = [Mroll_mean_wb_seq_pre_mean Mroll_mean_wb_mean_now];
% 
%         Mpitch_mean_wb_accel_seq_pre_mean = [Mpitch_mean_wb_accel_seq_pre_mean Mpitch_mean_wb_accel_mean_now];
%         Mpitch_mean_wb_damp_seq_pre_mean = [Mpitch_mean_wb_damp_seq_pre_mean Mpitch_mean_wb_damp_mean_now];
%         Mpitch_mean_wb_seq_pre_mean = [Mpitch_mean_wb_seq_pre_mean Mpitch_mean_wb_mean_now];
% 
%         Myaw_mean_wb_accel_seq_pre_mean = [Myaw_mean_wb_accel_seq_pre_mean Myaw_mean_wb_accel_mean_now];
%         Myaw_mean_wb_damp_seq_pre_mean = [Myaw_mean_wb_damp_seq_pre_mean Myaw_mean_wb_damp_mean_now];
%         Myaw_mean_wb_seq_pre_mean = [Myaw_mean_wb_seq_pre_mean Myaw_mean_wb_mean_now];
% 
%         M_L_mean_wb_accel_seq_pre_mean = [M_L_mean_wb_accel_seq_pre_mean M_L_mean_wb_accel_mean_now];
%         M_L_mean_wb_damp_seq_pre_mean = [M_L_mean_wb_damp_seq_pre_mean M_L_mean_wb_damp_mean_now];
%         M_L_mean_wb_seq_pre_mean = [M_L_mean_wb_seq_pre_mean M_L_mean_wb_mean_now];
% 
%         M_R_mean_wb_accel_seq_pre_mean = [M_R_mean_wb_accel_seq_pre_mean M_R_mean_wb_accel_mean_now];
%         M_R_mean_wb_damp_seq_pre_mean = [M_R_mean_wb_damp_seq_pre_mean M_R_mean_wb_damp_mean_now];
%         M_R_mean_wb_seq_pre_mean = [M_R_mean_wb_seq_pre_mean M_R_mean_wb_mean_now];

        %% ste
        % wb kin
        f_wb_seq_pre_ste = [f_wb_seq_pre_ste f_wb_ste_now];
        
        stroke_wb_L_seq_bins_pre_ste = [stroke_wb_L_seq_bins_pre_ste stroke_wb_L_ste_now(1:end-1)];
        stroke_wb_R_seq_bins_pre_ste = [stroke_wb_R_seq_bins_pre_ste stroke_wb_R_ste_now(1:end-1)];
        
        dev_wb_L_seq_bins_pre_ste = [dev_wb_L_seq_bins_pre_ste dev_wb_L_ste_now(1:end-1)];
        dev_wb_R_seq_bins_pre_ste = [dev_wb_R_seq_bins_pre_ste dev_wb_R_ste_now(1:end-1)];
        
        pitch_wb_L_seq_bins_pre_ste = [pitch_wb_L_seq_bins_pre_ste pitch_wb_L_ste_now(1:end-1)];
        pitch_wb_R_seq_bins_pre_ste = [pitch_wb_R_seq_bins_pre_ste pitch_wb_R_ste_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_pre_ste = [roll_mean_wb_seq_pre_ste roll_mean_wb_ste_now];
        roll_dot_mean_wb_seq_pre_ste = [roll_dot_mean_wb_seq_pre_ste roll_dot_mean_wb_ste_now];
        roll_dot_dot_mean_wb_seq_pre_ste = [roll_dot_dot_mean_wb_seq_pre_ste roll_dot_dot_mean_wb_ste_now];
        
        pitch_mean_wb_seq_pre_ste = [pitch_mean_wb_seq_pre_ste pitch_mean_wb_ste_now];
        pitch_dot_mean_wb_seq_pre_ste = [pitch_dot_mean_wb_seq_pre_ste pitch_dot_mean_wb_ste_now];
        pitch_dot_dot_mean_wb_seq_pre_ste = [pitch_dot_dot_mean_wb_seq_pre_ste pitch_dot_dot_mean_wb_ste_now];
        
        yaw_mean_wb_seq_pre_ste = [yaw_mean_wb_seq_pre_ste yaw_mean_wb_ste_now];
        yaw_dot_mean_wb_seq_pre_ste = [yaw_dot_mean_wb_seq_pre_ste yaw_dot_mean_wb_ste_now];
        yaw_dot_dot_mean_wb_seq_pre_ste = [yaw_dot_dot_mean_wb_seq_pre_ste yaw_dot_dot_mean_wb_ste_now];
        
%         drot_L_mean_wb_seq_pre_ste = [drot_L_mean_wb_seq_pre_ste drot_L_mean_wb_ste_now];
%         rot_dot_L_mean_wb_seq_pre_ste = [rot_dot_L_mean_wb_seq_pre_ste rot_dot_L_mean_wb_ste_now];
%         rot_dot_dot_L_mean_wb_seq_pre_ste = [rot_dot_dot_L_mean_wb_seq_pre_ste rot_dot_dot_L_mean_wb_ste_now];
%         
%         drot_R_mean_wb_seq_pre_ste = [drot_R_mean_wb_seq_pre_ste drot_R_mean_wb_ste_now];
%         rot_dot_R_mean_wb_seq_pre_ste = [rot_dot_R_mean_wb_seq_pre_ste rot_dot_R_mean_wb_ste_now];
%         rot_dot_dot_R_mean_wb_seq_pre_ste = [rot_dot_dot_R_mean_wb_seq_pre_ste rot_dot_dot_R_mean_wb_ste_now];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_ste = [Mroll_mean_wb_accel_seq_pre_ste Mroll_mean_wb_accel_ste_now];
%         Mroll_mean_wb_damp_seq_pre_ste = [Mroll_mean_wb_damp_seq_pre_ste Mroll_mean_wb_damp_ste_now];
%         Mroll_mean_wb_seq_pre_ste = [Mroll_mean_wb_seq_pre_ste Mroll_mean_wb_ste_now];
% 
%         Mpitch_mean_wb_accel_seq_pre_ste = [Mpitch_mean_wb_accel_seq_pre_ste Mpitch_mean_wb_accel_ste_now];
%         Mpitch_mean_wb_damp_seq_pre_ste = [Mpitch_mean_wb_damp_seq_pre_ste Mpitch_mean_wb_damp_ste_now];
%         Mpitch_mean_wb_seq_pre_ste = [Mpitch_mean_wb_seq_pre_ste Mpitch_mean_wb_ste_now];
% 
%         Myaw_mean_wb_accel_seq_pre_ste = [Myaw_mean_wb_accel_seq_pre_ste Myaw_mean_wb_accel_ste_now];
%         Myaw_mean_wb_damp_seq_pre_ste = [Myaw_mean_wb_damp_seq_pre_ste Myaw_mean_wb_damp_ste_now];
%         Myaw_mean_wb_seq_pre_ste = [Myaw_mean_wb_seq_pre_ste Myaw_mean_wb_ste_now];
% 
%         M_L_mean_wb_accel_seq_pre_ste = [M_L_mean_wb_accel_seq_pre_ste M_L_mean_wb_accel_ste_now];
%         M_L_mean_wb_damp_seq_pre_ste = [M_L_mean_wb_damp_seq_pre_ste M_L_mean_wb_damp_ste_now];
%         M_L_mean_wb_seq_pre_ste = [M_L_mean_wb_seq_pre_ste M_L_mean_wb_ste_now];
% 
%         M_R_mean_wb_accel_seq_pre_ste = [M_R_mean_wb_accel_seq_pre_ste M_R_mean_wb_accel_ste_now];
%         M_R_mean_wb_damp_seq_pre_ste = [M_R_mean_wb_damp_seq_pre_ste M_R_mean_wb_damp_ste_now];
%         M_R_mean_wb_seq_pre_ste = [M_R_mean_wb_seq_pre_ste M_R_mean_wb_ste_now];

        %% mean ALL
        % wb kin
        f_wb_seq_pre_mean_all = [f_wb_seq_pre_mean_all; f_wb_mean_now];
        
        stroke_wb_L_seq_bins_pre_mean_all = [stroke_wb_L_seq_bins_pre_mean_all; stroke_wb_L_mean_now(1:end-1)];
        stroke_wb_R_seq_bins_pre_mean_all = [stroke_wb_R_seq_bins_pre_mean_all; stroke_wb_R_mean_now(1:end-1)];
        
        dev_wb_L_seq_bins_pre_mean_all = [dev_wb_L_seq_bins_pre_mean_all; dev_wb_L_mean_now(1:end-1)];
        dev_wb_R_seq_bins_pre_mean_all = [dev_wb_R_seq_bins_pre_mean_all; dev_wb_R_mean_now(1:end-1)];
        
        pitch_wb_L_seq_bins_pre_mean_all = [pitch_wb_L_seq_bins_pre_mean_all; pitch_wb_L_mean_now(1:end-1)];
        pitch_wb_R_seq_bins_pre_mean_all = [pitch_wb_R_seq_bins_pre_mean_all; pitch_wb_R_mean_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean_all = [roll_mean_wb_seq_pre_mean_all; roll_mean_wb_mean_now];
        roll_dot_mean_wb_seq_pre_mean_all = [roll_dot_mean_wb_seq_pre_mean_all; roll_dot_mean_wb_mean_now];
        roll_dot_dot_mean_wb_seq_pre_mean_all = [roll_dot_dot_mean_wb_seq_pre_mean_all; roll_dot_dot_mean_wb_mean_now];
        
        pitch_mean_wb_seq_pre_mean_all = [pitch_mean_wb_seq_pre_mean_all; pitch_mean_wb_mean_now];
        pitch_dot_mean_wb_seq_pre_mean_all = [pitch_dot_mean_wb_seq_pre_mean_all; pitch_dot_mean_wb_mean_now];
        pitch_dot_dot_mean_wb_seq_pre_mean_all = [pitch_dot_dot_mean_wb_seq_pre_mean_all; pitch_dot_dot_mean_wb_mean_now];
        
        yaw_mean_wb_seq_pre_mean_all = [yaw_mean_wb_seq_pre_mean_all; yaw_mean_wb_mean_now];
        yaw_dot_mean_wb_seq_pre_mean_all = [yaw_dot_mean_wb_seq_pre_mean_all; yaw_dot_mean_wb_mean_now];
        yaw_dot_dot_mean_wb_seq_pre_mean_all = [yaw_dot_dot_mean_wb_seq_pre_mean_all; yaw_dot_dot_mean_wb_mean_now];
        
%         drot_L_mean_wb_seq_pre_mean_all = [drot_L_mean_wb_seq_pre_mean_all; drot_L_mean_wb_mean_now];
%         rot_dot_L_mean_wb_seq_pre_mean_all = [rot_dot_L_mean_wb_seq_pre_mean_all; rot_dot_L_mean_wb_mean_now];
%         rot_dot_dot_L_mean_wb_seq_pre_mean_all = [rot_dot_dot_L_mean_wb_seq_pre_mean_all; rot_dot_dot_L_mean_wb_mean_now];
%         
%         drot_R_mean_wb_seq_pre_mean_all = [drot_R_mean_wb_seq_pre_mean_all; drot_R_mean_wb_mean_now];
%         rot_dot_R_mean_wb_seq_pre_mean_all = [rot_dot_R_mean_wb_seq_pre_mean_all; rot_dot_R_mean_wb_mean_now];
%         rot_dot_dot_R_mean_wb_seq_pre_mean_all = [rot_dot_dot_R_mean_wb_seq_pre_mean_all; rot_dot_dot_R_mean_wb_mean_now];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_mean_all = [Mroll_mean_wb_accel_seq_pre_mean_all; Mroll_mean_wb_accel_mean_now];
%         Mroll_mean_wb_damp_seq_pre_mean_all = [Mroll_mean_wb_damp_seq_pre_mean_all; Mroll_mean_wb_damp_mean_now];
%         Mroll_mean_wb_seq_pre_mean_all = [Mroll_mean_wb_seq_pre_mean_all; Mroll_mean_wb_mean_now];
% 
%         Mpitch_mean_wb_accel_seq_pre_mean_all = [Mpitch_mean_wb_accel_seq_pre_mean_all; Mpitch_mean_wb_accel_mean_now];
%         Mpitch_mean_wb_damp_seq_pre_mean_all = [Mpitch_mean_wb_damp_seq_pre_mean_all; Mpitch_mean_wb_damp_mean_now];
%         Mpitch_mean_wb_seq_pre_mean_all = [Mpitch_mean_wb_seq_pre_mean_all; Mpitch_mean_wb_mean_now];
% 
%         Myaw_mean_wb_accel_seq_pre_mean_all = [Myaw_mean_wb_accel_seq_pre_mean_all; Myaw_mean_wb_accel_mean_now];
%         Myaw_mean_wb_damp_seq_pre_mean_all = [Myaw_mean_wb_damp_seq_pre_mean_all; Myaw_mean_wb_damp_mean_now];
%         Myaw_mean_wb_seq_pre_mean_all = [Myaw_mean_wb_seq_pre_mean_all; Myaw_mean_wb_mean_now];
% 
%         M_L_mean_wb_accel_seq_pre_mean_all = [M_L_mean_wb_accel_seq_pre_mean_all; M_L_mean_wb_accel_mean_now];
%         M_L_mean_wb_damp_seq_pre_mean_all = [M_L_mean_wb_damp_seq_pre_mean_all; M_L_mean_wb_damp_mean_now];
%         M_L_mean_wb_seq_pre_mean_all = [M_L_mean_wb_seq_pre_mean_all; M_L_mean_wb_mean_now];
% 
%         M_R_mean_wb_accel_seq_pre_mean_all = [M_R_mean_wb_accel_seq_pre_mean_all; M_R_mean_wb_accel_mean_now];
%         M_R_mean_wb_damp_seq_pre_mean_all = [M_R_mean_wb_damp_seq_pre_mean_all; M_R_mean_wb_damp_mean_now];
%         M_R_mean_wb_seq_pre_mean_all = [M_R_mean_wb_seq_pre_mean_all; M_R_mean_wb_mean_now];

        %% ste ALL
        % wb kin
        f_wb_seq_pre_ste_all = [f_wb_seq_pre_ste_all; f_wb_ste_now];
        
        stroke_wb_L_seq_bins_pre_ste_all = [stroke_wb_L_seq_bins_pre_ste_all; stroke_wb_L_ste_now(1:end-1)];
        stroke_wb_R_seq_bins_pre_ste_all = [stroke_wb_R_seq_bins_pre_ste_all; stroke_wb_R_ste_now(1:end-1)];
        
        dev_wb_L_seq_bins_pre_ste_all = [dev_wb_L_seq_bins_pre_ste_all; dev_wb_L_ste_now(1:end-1)];
        dev_wb_R_seq_bins_pre_ste_all = [dev_wb_R_seq_bins_pre_ste_all; dev_wb_R_ste_now(1:end-1)];
        
        pitch_wb_L_seq_bins_pre_ste_all = [pitch_wb_L_seq_bins_pre_ste_all; pitch_wb_L_ste_now(1:end-1)];
        pitch_wb_R_seq_bins_pre_ste_all = [pitch_wb_R_seq_bins_pre_ste_all; pitch_wb_R_ste_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_pre_ste_all = [roll_mean_wb_seq_pre_ste_all; roll_mean_wb_ste_now];
        roll_dot_mean_wb_seq_pre_ste_all = [roll_dot_mean_wb_seq_pre_ste_all; roll_dot_mean_wb_ste_now];
        roll_dot_dot_mean_wb_seq_pre_ste_all = [roll_dot_dot_mean_wb_seq_pre_ste_all; roll_dot_dot_mean_wb_ste_now];
        
        pitch_mean_wb_seq_pre_ste_all = [pitch_mean_wb_seq_pre_ste_all; pitch_mean_wb_ste_now];
        pitch_dot_mean_wb_seq_pre_ste_all = [pitch_dot_mean_wb_seq_pre_ste_all; pitch_dot_mean_wb_ste_now];
        pitch_dot_dot_mean_wb_seq_pre_ste_all = [pitch_dot_dot_mean_wb_seq_pre_ste_all; pitch_dot_dot_mean_wb_ste_now];
        
        yaw_mean_wb_seq_pre_ste_all = [yaw_mean_wb_seq_pre_ste_all; yaw_mean_wb_ste_now];
        yaw_dot_mean_wb_seq_pre_ste_all = [yaw_dot_mean_wb_seq_pre_ste_all; yaw_dot_mean_wb_ste_now];
        yaw_dot_dot_mean_wb_seq_pre_ste_all = [yaw_dot_dot_mean_wb_seq_pre_ste_all; yaw_dot_dot_mean_wb_ste_now];
        
%         drot_L_mean_wb_seq_pre_ste_all = [drot_L_mean_wb_seq_pre_ste_all; drot_L_mean_wb_ste_now];
%         rot_dot_L_mean_wb_seq_pre_ste_all = [rot_dot_L_mean_wb_seq_pre_ste_all; rot_dot_L_mean_wb_ste_now];
%         rot_dot_dot_L_mean_wb_seq_pre_ste_all = [rot_dot_dot_L_mean_wb_seq_pre_ste_all; rot_dot_dot_L_mean_wb_ste_now];
%         
%         drot_R_mean_wb_seq_pre_ste_all = [drot_R_mean_wb_seq_pre_ste_all; drot_R_mean_wb_ste_now];
%         rot_dot_R_mean_wb_seq_pre_ste_all = [rot_dot_R_mean_wb_seq_pre_ste_all; rot_dot_R_mean_wb_ste_now];
%         rot_dot_dot_R_mean_wb_seq_pre_ste_all = [rot_dot_dot_R_mean_wb_seq_pre_ste_all; rot_dot_dot_R_mean_wb_ste_now];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_ste_all = [Mroll_mean_wb_accel_seq_pre_ste_all; Mroll_mean_wb_accel_ste_now];
%         Mroll_mean_wb_damp_seq_pre_ste_all = [Mroll_mean_wb_damp_seq_pre_ste_all; Mroll_mean_wb_damp_ste_now];
%         Mroll_mean_wb_seq_pre_ste_all = [Mroll_mean_wb_seq_pre_ste_all; Mroll_mean_wb_ste_now];
% 
%         Mpitch_mean_wb_accel_seq_pre_ste_all = [Mpitch_mean_wb_accel_seq_pre_ste_all; Mpitch_mean_wb_accel_ste_now];
%         Mpitch_mean_wb_damp_seq_pre_ste_all = [Mpitch_mean_wb_damp_seq_pre_ste_all; Mpitch_mean_wb_damp_ste_now];
%         Mpitch_mean_wb_seq_pre_ste_all = [Mpitch_mean_wb_seq_pre_ste_all; Mpitch_mean_wb_ste_now];
% 
%         Myaw_mean_wb_accel_seq_pre_ste_all = [Myaw_mean_wb_accel_seq_pre_ste_all; Myaw_mean_wb_accel_ste_now];
%         Myaw_mean_wb_damp_seq_pre_ste_all = [Myaw_mean_wb_damp_seq_pre_ste_all; Myaw_mean_wb_damp_ste_now];
%         Myaw_mean_wb_seq_pre_ste_all = [Myaw_mean_wb_seq_pre_ste_all; Myaw_mean_wb_ste_now];
% 
%         M_L_mean_wb_accel_seq_pre_ste_all = [M_L_mean_wb_accel_seq_pre_ste_all; M_L_mean_wb_accel_ste_now];
%         M_L_mean_wb_damp_seq_pre_ste_all = [M_L_mean_wb_damp_seq_pre_ste_all; M_L_mean_wb_damp_ste_now];
%         M_L_mean_wb_seq_pre_ste_all = [M_L_mean_wb_seq_pre_ste_all; M_L_mean_wb_ste_now];
% 
%         M_R_mean_wb_accel_seq_pre_ste_all = [M_R_mean_wb_accel_seq_pre_ste_all; M_R_mean_wb_accel_ste_now];
%         M_R_mean_wb_damp_seq_pre_ste_all = [M_R_mean_wb_damp_seq_pre_ste_all; M_R_mean_wb_damp_ste_now];
%         M_R_mean_wb_seq_pre_ste_all = [M_R_mean_wb_seq_pre_ste_all; M_R_mean_wb_ste_now];

        %% all
        % wb kin
        f_wb_seq_pre_all = [f_wb_seq_pre_all f_wb_now];
        
        stroke_wb_L_seq_bins_pre_all = [stroke_wb_L_seq_bins_pre_all; stroke_wb_L_now(1:end-1,:)];
        stroke_wb_R_seq_bins_pre_all = [stroke_wb_R_seq_bins_pre_all; stroke_wb_R_now(1:end-1,:)];
        
        dev_wb_L_seq_bins_pre_all = [dev_wb_L_seq_bins_pre_all; dev_wb_L_now(1:end-1,:)];
        dev_wb_R_seq_bins_pre_all = [dev_wb_R_seq_bins_pre_all; dev_wb_R_now(1:end-1,:)];
        
        pitch_wb_L_seq_bins_pre_all = [pitch_wb_L_seq_bins_pre_all; pitch_wb_L_now(1:end-1,:)];
        pitch_wb_R_seq_bins_pre_all = [pitch_wb_R_seq_bins_pre_all; pitch_wb_R_now(1:end-1,:)];
        
        % body kin
        roll_mean_wb_seq_pre_all = [roll_mean_wb_seq_pre_all roll_mean_wb_now];
        roll_dot_mean_wb_seq_pre_all = [roll_dot_mean_wb_seq_pre_all roll_dot_mean_wb_now];
        roll_dot_dot_mean_wb_seq_pre_all = [roll_dot_dot_mean_wb_seq_pre_all roll_dot_dot_mean_wb_now];
        
        pitch_mean_wb_seq_pre_all = [pitch_mean_wb_seq_pre_all pitch_mean_wb_now];
        pitch_dot_mean_wb_seq_pre_all = [pitch_dot_mean_wb_seq_pre_all pitch_dot_mean_wb_now];
        pitch_dot_dot_mean_wb_seq_pre_all = [pitch_dot_dot_mean_wb_seq_pre_all pitch_dot_dot_mean_wb_now];
        
        yaw_mean_wb_seq_pre_all = [yaw_mean_wb_seq_pre_all yaw_mean_wb_now];
        yaw_dot_mean_wb_seq_pre_all = [yaw_dot_mean_wb_seq_pre_all yaw_dot_mean_wb_now];
        yaw_dot_dot_mean_wb_seq_pre_all = [yaw_dot_dot_mean_wb_seq_pre_all yaw_dot_dot_mean_wb_now];
        
%         drot_L_mean_wb_seq_pre_all = [drot_L_mean_wb_seq_pre_all drot_L_mean_wb_now];
%         rot_dot_L_mean_wb_seq_pre_all = [rot_dot_L_mean_wb_seq_pre_all rot_dot_L_mean_wb_now];
%         rot_dot_dot_L_mean_wb_seq_pre_all = [rot_dot_dot_L_mean_wb_seq_pre_all rot_dot_dot_L_mean_wb_now];
%         
%         drot_R_mean_wb_seq_pre_all = [drot_R_mean_wb_seq_pre_all drot_R_mean_wb_now];
%         rot_dot_R_mean_wb_seq_pre_all = [rot_dot_R_mean_wb_seq_pre_all rot_dot_R_mean_wb_now];
%         rot_dot_dot_R_mean_wb_seq_pre_all = [rot_dot_dot_R_mean_wb_seq_pre_all rot_dot_dot_R_mean_wb_now];
% 
%         % torque
%         Mroll_mean_wb_accel_seq_pre_all = [Mroll_mean_wb_accel_seq_pre_all Mroll_mean_wb_accel_now];
%         Mroll_mean_wb_damp_seq_pre_all = [Mroll_mean_wb_damp_seq_pre_all Mroll_mean_wb_damp_now];
%         Mroll_mean_wb_seq_pre_all = [Mroll_mean_wb_seq_pre_all Mroll_mean_wb_now];
% 
%         Mpitch_mean_wb_accel_seq_pre_all = [Mpitch_mean_wb_accel_seq_pre_all Mpitch_mean_wb_accel_now];
%         Mpitch_mean_wb_damp_seq_pre_all = [Mpitch_mean_wb_damp_seq_pre_all Mpitch_mean_wb_damp_now];
%         Mpitch_mean_wb_seq_pre_all = [Mpitch_mean_wb_seq_pre_all Mpitch_mean_wb_now];
% 
%         Myaw_mean_wb_accel_seq_pre_all = [Myaw_mean_wb_accel_seq_pre_all Myaw_mean_wb_accel_now];
%         Myaw_mean_wb_damp_seq_pre_all = [Myaw_mean_wb_damp_seq_pre_all Myaw_mean_wb_damp_now];
%         Myaw_mean_wb_seq_pre_all = [Myaw_mean_wb_seq_pre_all Myaw_mean_wb_now];
% 
%         M_L_mean_wb_accel_seq_pre_all = [M_L_mean_wb_accel_seq_pre_all M_L_mean_wb_accel_now];
%         M_L_mean_wb_damp_seq_pre_all = [M_L_mean_wb_damp_seq_pre_all M_L_mean_wb_damp_now];
%         M_L_mean_wb_seq_pre_all = [M_L_mean_wb_seq_pre_all M_L_mean_wb_now];
% 
%         M_R_mean_wb_accel_seq_pre_all = [M_R_mean_wb_accel_seq_pre_all M_R_mean_wb_accel_now];
%         M_R_mean_wb_damp_seq_pre_all = [M_R_mean_wb_damp_seq_pre_all M_R_mean_wb_damp_now];
%         M_R_mean_wb_seq_pre_all = [M_R_mean_wb_seq_pre_all M_R_mean_wb_now];

    else
        %% mean
        % wb kin
        f_wb_seq_pre_mean = [f_wb_seq_pre_mean nan];

        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean nan(n_bins-1,1)];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean = [roll_mean_wb_seq_pre_mean nan];
        roll_dot_mean_wb_seq_pre_mean = [roll_dot_mean_wb_seq_pre_mean nan];
        roll_dot_dot_mean_wb_seq_pre_mean = [roll_dot_dot_mean_wb_seq_pre_mean nan];
        
        pitch_mean_wb_seq_pre_mean = [pitch_mean_wb_seq_pre_mean nan];
        pitch_dot_mean_wb_seq_pre_mean = [pitch_dot_mean_wb_seq_pre_mean nan];
        pitch_dot_dot_mean_wb_seq_pre_mean = [pitch_dot_dot_mean_wb_seq_pre_mean nan];
        
        yaw_mean_wb_seq_pre_mean = [yaw_mean_wb_seq_pre_mean nan];
        yaw_dot_mean_wb_seq_pre_mean = [yaw_dot_mean_wb_seq_pre_mean nan];
        yaw_dot_dot_mean_wb_seq_pre_mean = [yaw_dot_dot_mean_wb_seq_pre_mean nan];
        
%         drot_L_mean_wb_seq_pre_mean = [drot_L_mean_wb_seq_pre_mean nan];
%         rot_dot_L_mean_wb_seq_pre_mean = [rot_dot_L_mean_wb_seq_pre_mean nan];
%         rot_dot_dot_L_mean_wb_seq_pre_mean = [rot_dot_dot_L_mean_wb_seq_pre_mean nan];
%         
%         drot_R_mean_wb_seq_pre_mean = [drot_R_mean_wb_seq_pre_mean nan];
%         rot_dot_R_mean_wb_seq_pre_mean = [rot_dot_R_mean_wb_seq_pre_mean nan];
%         rot_dot_dot_R_mean_wb_seq_pre_mean = [rot_dot_dot_R_mean_wb_seq_pre_mean nan];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_mean = [Mroll_mean_wb_accel_seq_pre_mean nan];
%         Mroll_mean_wb_damp_seq_pre_mean = [Mroll_mean_wb_damp_seq_pre_mean nan];
%         Mroll_mean_wb_seq_pre_mean = [Mroll_mean_wb_seq_pre_mean nan];
% 
%         Mpitch_mean_wb_accel_seq_pre_mean = [Mpitch_mean_wb_accel_seq_pre_mean nan];
%         Mpitch_mean_wb_damp_seq_pre_mean = [Mpitch_mean_wb_damp_seq_pre_mean nan];
%         Mpitch_mean_wb_seq_pre_mean = [Mpitch_mean_wb_seq_pre_mean nan];
% 
%         Myaw_mean_wb_accel_seq_pre_mean = [Myaw_mean_wb_accel_seq_pre_mean nan];
%         Myaw_mean_wb_damp_seq_pre_mean = [Myaw_mean_wb_damp_seq_pre_mean nan];
%         Myaw_mean_wb_seq_pre_mean = [Myaw_mean_wb_seq_pre_mean nan];
% 
%         M_L_mean_wb_accel_seq_pre_mean = [M_L_mean_wb_accel_seq_pre_mean nan];
%         M_L_mean_wb_damp_seq_pre_mean = [M_L_mean_wb_damp_seq_pre_mean nan];
%         M_L_mean_wb_seq_pre_mean = [M_L_mean_wb_seq_pre_mean nan];
% 
%         M_R_mean_wb_accel_seq_pre_mean = [M_R_mean_wb_accel_seq_pre_mean nan];
%         M_R_mean_wb_damp_seq_pre_mean = [M_R_mean_wb_damp_seq_pre_mean nan];
%         M_R_mean_wb_seq_pre_mean = [M_R_mean_wb_seq_pre_mean nan];

        %% ste
        % wb kin
        f_wb_seq_pre_ste = [f_wb_seq_pre_ste nan];

        stroke_wb_L_seq_bins_pre_ste = [stroke_wb_L_seq_bins_pre_ste nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_pre_ste = [stroke_wb_R_seq_bins_pre_ste nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_pre_ste = [dev_wb_L_seq_bins_pre_ste nan(n_bins-1,1)];
        dev_wb_R_seq_bins_pre_ste = [dev_wb_R_seq_bins_pre_ste nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_pre_ste = [pitch_wb_L_seq_bins_pre_ste nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_pre_ste = [pitch_wb_R_seq_bins_pre_ste nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_pre_ste = [roll_mean_wb_seq_pre_ste nan];
        roll_dot_mean_wb_seq_pre_ste = [roll_dot_mean_wb_seq_pre_ste nan];
        roll_dot_dot_mean_wb_seq_pre_ste = [roll_dot_dot_mean_wb_seq_pre_ste nan];
        
        pitch_mean_wb_seq_pre_ste = [pitch_mean_wb_seq_pre_ste nan];
        pitch_dot_mean_wb_seq_pre_ste = [pitch_dot_mean_wb_seq_pre_ste nan];
        pitch_dot_dot_mean_wb_seq_pre_ste = [pitch_dot_dot_mean_wb_seq_pre_ste nan];
        
        yaw_mean_wb_seq_pre_ste = [yaw_mean_wb_seq_pre_ste nan];
        yaw_dot_mean_wb_seq_pre_ste = [yaw_dot_mean_wb_seq_pre_ste nan];
        yaw_dot_dot_mean_wb_seq_pre_ste = [yaw_dot_dot_mean_wb_seq_pre_ste nan];
        
%         drot_L_mean_wb_seq_pre_ste = [drot_L_mean_wb_seq_pre_ste nan];
%         rot_dot_L_mean_wb_seq_pre_ste = [rot_dot_L_mean_wb_seq_pre_ste nan];
%         rot_dot_dot_L_mean_wb_seq_pre_ste = [rot_dot_dot_L_mean_wb_seq_pre_ste nan];
%         
%         drot_R_mean_wb_seq_pre_ste = [drot_R_mean_wb_seq_pre_ste nan];
%         rot_dot_R_mean_wb_seq_pre_ste = [rot_dot_R_mean_wb_seq_pre_ste nan];
%         rot_dot_dot_R_mean_wb_seq_pre_ste = [rot_dot_dot_R_mean_wb_seq_pre_ste nan];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_ste = [Mroll_mean_wb_accel_seq_pre_ste nan];
%         Mroll_mean_wb_damp_seq_pre_ste = [Mroll_mean_wb_damp_seq_pre_ste nan];
%         Mroll_mean_wb_seq_pre_ste = [Mroll_mean_wb_seq_pre_ste nan];
% 
%         Mpitch_mean_wb_accel_seq_pre_ste = [Mpitch_mean_wb_accel_seq_pre_ste nan];
%         Mpitch_mean_wb_damp_seq_pre_ste = [Mpitch_mean_wb_damp_seq_pre_ste nan];
%         Mpitch_mean_wb_seq_pre_ste = [Mpitch_mean_wb_seq_pre_ste nan];
% 
%         Myaw_mean_wb_accel_seq_pre_ste = [Myaw_mean_wb_accel_seq_pre_ste nan];
%         Myaw_mean_wb_damp_seq_pre_ste = [Myaw_mean_wb_damp_seq_pre_ste nan];
%         Myaw_mean_wb_seq_pre_ste = [Myaw_mean_wb_seq_pre_ste nan];
% 
%         M_L_mean_wb_accel_seq_pre_ste = [M_L_mean_wb_accel_seq_pre_ste nan];
%         M_L_mean_wb_damp_seq_pre_ste = [M_L_mean_wb_damp_seq_pre_ste nan];
%         M_L_mean_wb_seq_pre_ste = [M_L_mean_wb_seq_pre_ste nan];
% 
%         M_R_mean_wb_accel_seq_pre_ste = [M_R_mean_wb_accel_seq_pre_ste nan];
%         M_R_mean_wb_damp_seq_pre_ste = [M_R_mean_wb_damp_seq_pre_ste nan];
%         M_R_mean_wb_seq_pre_ste = [M_R_mean_wb_seq_pre_ste nan];

        %% mean ALL
        % wb kin
        f_wb_seq_pre_mean_all = [f_wb_seq_pre_mean_all; nan];

        stroke_wb_L_seq_bins_pre_mean_all = [stroke_wb_L_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_pre_mean_all = [stroke_wb_R_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_pre_mean_all = [dev_wb_L_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        dev_wb_R_seq_bins_pre_mean_all = [dev_wb_R_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_pre_mean_all = [pitch_wb_L_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_pre_mean_all = [pitch_wb_R_seq_bins_pre_mean_all; nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean_all = [roll_mean_wb_seq_pre_mean_all; nan];
        roll_dot_mean_wb_seq_pre_mean_all = [roll_dot_mean_wb_seq_pre_mean_all; nan];
        roll_dot_dot_mean_wb_seq_pre_mean_all = [roll_dot_dot_mean_wb_seq_pre_mean_all; nan];
        
        pitch_mean_wb_seq_pre_mean_all = [pitch_mean_wb_seq_pre_mean_all; nan];
        pitch_dot_mean_wb_seq_pre_mean_all = [pitch_dot_mean_wb_seq_pre_mean_all; nan];
        pitch_dot_dot_mean_wb_seq_pre_mean_all = [pitch_dot_dot_mean_wb_seq_pre_mean_all; nan];
        
        yaw_mean_wb_seq_pre_mean_all = [yaw_mean_wb_seq_pre_mean_all; nan];
        yaw_dot_mean_wb_seq_pre_mean_all = [yaw_dot_mean_wb_seq_pre_mean_all; nan];
        yaw_dot_dot_mean_wb_seq_pre_mean_all = [yaw_dot_dot_mean_wb_seq_pre_mean_all; nan];
        
%         drot_L_mean_wb_seq_pre_mean_all = [drot_L_mean_wb_seq_pre_mean_all; nan];
%         rot_dot_L_mean_wb_seq_pre_mean_all = [rot_dot_L_mean_wb_seq_pre_mean_all; nan];
%         rot_dot_dot_L_mean_wb_seq_pre_mean_all = [rot_dot_dot_L_mean_wb_seq_pre_mean_all; nan];
%         
%         drot_R_mean_wb_seq_pre_mean_all = [drot_R_mean_wb_seq_pre_mean_all; nan];
%         rot_dot_R_mean_wb_seq_pre_mean_all = [rot_dot_R_mean_wb_seq_pre_mean_all; nan];
%         rot_dot_dot_R_mean_wb_seq_pre_mean_all = [rot_dot_dot_R_mean_wb_seq_pre_mean_all; nan];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_mean_all = [Mroll_mean_wb_accel_seq_pre_mean_all; nan];
%         Mroll_mean_wb_damp_seq_pre_mean_all = [Mroll_mean_wb_damp_seq_pre_mean_all; nan];
%         Mroll_mean_wb_seq_pre_mean_all = [Mroll_mean_wb_seq_pre_mean_all; nan];
% 
%         Mpitch_mean_wb_accel_seq_pre_mean_all = [Mpitch_mean_wb_accel_seq_pre_mean_all; nan];
%         Mpitch_mean_wb_damp_seq_pre_mean_all = [Mpitch_mean_wb_damp_seq_pre_mean_all; nan];
%         Mpitch_mean_wb_seq_pre_mean_all = [Mpitch_mean_wb_seq_pre_mean_all; nan];
% 
%         Myaw_mean_wb_accel_seq_pre_mean_all = [Myaw_mean_wb_accel_seq_pre_mean_all; nan];
%         Myaw_mean_wb_damp_seq_pre_mean_all = [Myaw_mean_wb_damp_seq_pre_mean_all; nan];
%         Myaw_mean_wb_seq_pre_mean_all = [Myaw_mean_wb_seq_pre_mean_all; nan];
% 
%         M_L_mean_wb_accel_seq_pre_mean_all = [M_L_mean_wb_accel_seq_pre_mean_all; nan];
%         M_L_mean_wb_damp_seq_pre_mean_all = [M_L_mean_wb_damp_seq_pre_mean_all; nan];
%         M_L_mean_wb_seq_pre_mean_all = [M_L_mean_wb_seq_pre_mean_all; nan];
% 
%         M_R_mean_wb_accel_seq_pre_mean_all = [M_R_mean_wb_accel_seq_pre_mean_all; nan];
%         M_R_mean_wb_damp_seq_pre_mean_all = [M_R_mean_wb_damp_seq_pre_mean_all; nan];
%         M_R_mean_wb_seq_pre_mean_all = [M_R_mean_wb_seq_pre_mean_all; nan];
% 
        %% ste ALL
        % wb kin
        f_wb_seq_pre_ste_all = [f_wb_seq_pre_ste_all; nan];

        stroke_wb_L_seq_bins_pre_ste_all = [stroke_wb_L_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_pre_ste_all = [stroke_wb_R_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_pre_ste_all = [dev_wb_L_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        dev_wb_R_seq_bins_pre_ste_all = [dev_wb_R_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_pre_ste_all = [pitch_wb_L_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_pre_ste_all = [pitch_wb_R_seq_bins_pre_ste_all; nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_pre_ste_all = [roll_mean_wb_seq_pre_ste_all; nan];
        roll_dot_mean_wb_seq_pre_ste_all = [roll_dot_mean_wb_seq_pre_ste_all; nan];
        roll_dot_dot_mean_wb_seq_pre_ste_all = [roll_dot_dot_mean_wb_seq_pre_ste_all; nan];
        
        pitch_mean_wb_seq_pre_ste_all = [pitch_mean_wb_seq_pre_ste_all; nan];
        pitch_dot_mean_wb_seq_pre_ste_all = [pitch_dot_mean_wb_seq_pre_ste_all; nan];
        pitch_dot_dot_mean_wb_seq_pre_ste_all = [pitch_dot_dot_mean_wb_seq_pre_ste_all; nan];
        
        yaw_mean_wb_seq_pre_ste_all = [yaw_mean_wb_seq_pre_ste_all; nan];
        yaw_dot_mean_wb_seq_pre_ste_all = [yaw_dot_mean_wb_seq_pre_ste_all; nan];
        yaw_dot_dot_mean_wb_seq_pre_ste_all = [yaw_dot_dot_mean_wb_seq_pre_ste_all; nan];
        
%         drot_L_mean_wb_seq_pre_ste_all = [drot_L_mean_wb_seq_pre_ste_all; nan];
%         rot_dot_L_mean_wb_seq_pre_ste_all = [rot_dot_L_mean_wb_seq_pre_ste_all; nan];
%         rot_dot_dot_L_mean_wb_seq_pre_ste_all = [rot_dot_dot_L_mean_wb_seq_pre_ste_all; nan];
%         
%         drot_R_mean_wb_seq_pre_ste_all = [drot_R_mean_wb_seq_pre_ste_all; nan];
%         rot_dot_R_mean_wb_seq_pre_ste_all = [rot_dot_R_mean_wb_seq_pre_ste_all; nan];
%         rot_dot_dot_R_mean_wb_seq_pre_ste_all = [rot_dot_dot_R_mean_wb_seq_pre_ste_all; nan];
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre_ste_all = [Mroll_mean_wb_accel_seq_pre_ste_all; nan];
%         Mroll_mean_wb_damp_seq_pre_ste_all = [Mroll_mean_wb_damp_seq_pre_ste_all; nan];
%         Mroll_mean_wb_seq_pre_ste_all = [Mroll_mean_wb_seq_pre_ste_all; nan];
% 
%         Mpitch_mean_wb_accel_seq_pre_ste_all = [Mpitch_mean_wb_accel_seq_pre_ste_all; nan];
%         Mpitch_mean_wb_damp_seq_pre_ste_all = [Mpitch_mean_wb_damp_seq_pre_ste_all; nan];
%         Mpitch_mean_wb_seq_pre_ste_all = [Mpitch_mean_wb_seq_pre_ste_all; nan];
% 
%         Myaw_mean_wb_accel_seq_pre_ste_all = [Myaw_mean_wb_accel_seq_pre_ste_all; nan];
%         Myaw_mean_wb_damp_seq_pre_ste_all = [Myaw_mean_wb_damp_seq_pre_ste_all; nan];
%         Myaw_mean_wb_seq_pre_ste_all = [Myaw_mean_wb_seq_pre_ste_all; nan];
% 
%         M_L_mean_wb_accel_seq_pre_ste_all = [M_L_mean_wb_accel_seq_pre_ste_all; nan];
%         M_L_mean_wb_damp_seq_pre_ste_all = [M_L_mean_wb_damp_seq_pre_ste_all; nan];
%         M_L_mean_wb_seq_pre_ste_all = [M_L_mean_wb_seq_pre_ste_all; nan];
% 
%         M_R_mean_wb_accel_seq_pre_ste_all = [M_R_mean_wb_accel_seq_pre_ste_all; nan];
%         M_R_mean_wb_damp_seq_pre_ste_all = [M_R_mean_wb_damp_seq_pre_ste_all; nan];
%         M_R_mean_wb_seq_pre_ste_all = [M_R_mean_wb_seq_pre_ste_all; nan];

        %% all
        % wb kin
        f_wb_seq_pre_all = [f_wb_seq_pre_all nan(max(seq_nr),1)];

        stroke_wb_L_seq_bins_pre_all = [stroke_wb_L_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        stroke_wb_R_seq_bins_pre_all = [stroke_wb_R_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        
        dev_wb_L_seq_bins_pre_all = [dev_wb_L_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        dev_wb_R_seq_bins_pre_all = [dev_wb_R_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        
        pitch_wb_L_seq_bins_pre_all = [pitch_wb_L_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        pitch_wb_R_seq_bins_pre_all = [pitch_wb_R_seq_bins_pre_all; nan(n_bins-1,max(seq_nr))];
        
        % body kin
        roll_mean_wb_seq_pre_all = [roll_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        roll_dot_mean_wb_seq_pre_all = [roll_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        roll_dot_dot_mean_wb_seq_pre_all = [roll_dot_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        
        pitch_mean_wb_seq_pre_all = [pitch_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        pitch_dot_mean_wb_seq_pre_all = [pitch_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        pitch_dot_dot_mean_wb_seq_pre_all = [pitch_dot_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        
        yaw_mean_wb_seq_pre_all = [yaw_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        yaw_dot_mean_wb_seq_pre_all = [yaw_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        yaw_dot_dot_mean_wb_seq_pre_all = [yaw_dot_dot_mean_wb_seq_pre_all nan(max(seq_nr),1)];
        
%         drot_L_mean_wb_seq_pre_all = [drot_L_mean_wb_seq_pre_all nan(max(seq_nr),1)];
%         rot_dot_L_mean_wb_seq_pre_all = [rot_dot_L_mean_wb_seq_pre_all nan(max(seq_nr),1)];
%         rot_dot_dot_L_mean_wb_seq_pre_all = [rot_dot_dot_L_mean_wb_seq_pre_all nan(max(seq_nr),1)];
%         
%         drot_R_mean_wb_seq_pre_all = [drot_R_mean_wb_seq_pre_all nan(max(seq_nr),1)];
%         rot_dot_R_mean_wb_seq_pre_all = [rot_dot_R_mean_wb_seq_pre_all nan(max(seq_nr),1)];
%         rot_dot_dot_R_mean_wb_seq_pre_all = [rot_dot_dot_R_mean_wb_seq_pre_all nan(max(seq_nr),1)];
% 
%         % torque
%         Mroll_mean_wb_accel_seq_pre_all = [Mroll_mean_wb_accel_seq_pre_all nan(max(seq_nr),1)];
%         Mroll_mean_wb_damp_seq_pre_all = [Mroll_mean_wb_damp_seq_pre_all nan(max(seq_nr),1)];
%         Mroll_mean_wb_seq_pre_all = [Mroll_mean_wb_seq_pre_all nan(max(seq_nr),1)];
% 
%         Mpitch_mean_wb_accel_seq_pre_all = [Mpitch_mean_wb_accel_seq_pre_all nan(max(seq_nr),1)];
%         Mpitch_mean_wb_damp_seq_pre_all = [Mpitch_mean_wb_damp_seq_pre_all nan(max(seq_nr),1)];
%         Mpitch_mean_wb_seq_pre_all = [Mpitch_mean_wb_seq_pre_all nan(max(seq_nr),1)];
% 
%         Myaw_mean_wb_accel_seq_pre_all = [Myaw_mean_wb_accel_seq_pre_all nan(max(seq_nr),1)];
%         Myaw_mean_wb_damp_seq_pre_all = [Myaw_mean_wb_damp_seq_pre_all nan(max(seq_nr),1)];
%         Myaw_mean_wb_seq_pre_all = [Myaw_mean_wb_seq_pre_all nan(max(seq_nr),1)];
% 
%         M_L_mean_wb_accel_seq_pre_all = [M_L_mean_wb_accel_seq_pre_all nan(max(seq_nr),1)];
%         M_L_mean_wb_damp_seq_pre_all = [M_L_mean_wb_damp_seq_pre_all nan(max(seq_nr),1)];
%         M_L_mean_wb_seq_pre_all = [M_L_mean_wb_seq_pre_all nan(max(seq_nr),1)];
% 
%         M_R_mean_wb_accel_seq_pre_all = [M_R_mean_wb_accel_seq_pre_all nan(max(seq_nr),1)];
%         M_R_mean_wb_damp_seq_pre_all = [M_R_mean_wb_damp_seq_pre_all nan(max(seq_nr),1)];
%         M_R_mean_wb_seq_pre_all = [M_R_mean_wb_seq_pre_all nan(max(seq_nr),1)];

    end

