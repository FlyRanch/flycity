

rot_axis_pos_mirror = atan2d(pitch,roll_mirror);
rot_axis_vel_mirror = atan2d(pitch_dot,roll_dot_mirror);
rot_axis_accel_mirror = atan2d(pitch_dot_dot,roll_dot_dot_mirror);

Maccel_axis = atan2d(Mpitch_accel,Mroll_accel);
Mdamp_axis = atan2d(Mpitch_damp,Mroll_damp);
M_axis = atan2d(Mpitch,Mroll);

% no angles < -90
rot_axis_pos_plot = rot_axis_pos_mirror;
rot_axis_pos_plot(rot_axis_pos_plot<-90) = rot_axis_pos_plot(rot_axis_pos_plot<-90)+360;

rot_axis_vel_plot = rot_axis_vel_mirror;
rot_axis_vel_plot(rot_axis_vel_plot<-90) = rot_axis_vel_plot(rot_axis_vel_plot<-90)+360;

rot_axis_accel_plot = rot_axis_accel_mirror;
rot_axis_accel_plot(rot_axis_accel_plot<-90) = rot_axis_accel_plot(rot_axis_accel_plot<-90)+360;

Maccel_axis_plot = Maccel_axis;
Maccel_axis_plot(Maccel_axis_plot<-90) = Maccel_axis_plot(Maccel_axis_plot<-90)+360;

Mdamp_axis_plot = Mdamp_axis;
Mdamp_axis_plot(Mdamp_axis_plot<-90) = Mdamp_axis_plot(Mdamp_axis_plot<-90)+360;

M_axis_plot = M_axis;
M_axis_plot(M_axis_plot<-90) = M_axis_plot(M_axis_plot<-90)+360;

% remove jumps from plots (+/-180deg)
rot_axis_pos_temp = rot_axis_pos_plot;
rot_axis_vel_temp = rot_axis_vel_plot;
rot_axis_accel_temp = rot_axis_accel_plot;

Maccel_axis_temp = Maccel_axis_plot;
Mdamp_axis_temp = Mdamp_axis_plot;
M_axis_temp = M_axis_plot;

for i=1:size(rot_axis_vel_plot,2)
    for j=2:size(rot_axis_vel_plot,1)
        if abs(rot_axis_pos_temp(j,i) - rot_axis_pos_temp(j-1,i)) > 80
            rot_axis_pos_plot(j-1,i) = nan;
        end
        
        if abs(rot_axis_vel_temp(j,i) - rot_axis_vel_temp(j-1,i)) > 80
            rot_axis_vel_plot(j-1,i) = nan;
        end
        
        if abs(rot_axis_accel_temp(j,i) - rot_axis_accel_temp(j-1,i)) > 80
            rot_axis_accel_plot(j-1,i) = nan;
        end
        
        if abs(Maccel_axis_temp(j,i) - Maccel_axis_temp(j-1,i)) > 80
            Maccel_axis_plot(j-1,i) = nan;
        end
        
        if abs(Mdamp_axis_temp(j,i) - Mdamp_axis_temp(j-1,i)) > 80
            Mdamp_axis_plot(j-1,i) = nan;
        end
        
        if abs(M_axis_temp(j,i) - M_axis_temp(j-1,i)) > 80
            M_axis_plot(j-1,i) = nan;
        end
    end
end

for i=1:size(rot_axis_vel_plot,2)
    
    t_now = t-t_shift(i);
    
    %% yaw angles
    dyaw_now = dyaw_mirror(:,i);
    
    dyaw_rot1 = dyaw_now(t_now>t_rot1_min & t_now<t_rot1_max);
    dyaw_rot1 = dyaw_rot1(isnan(dyaw_rot1)==0);
    dyaw_rot1_mean(i,1) = nanmean(dyaw_rot1);
    
    dyaw_rot2 = dyaw_now(t_now>t_rot2_min & t_now<t_rot2_max);
    dyaw_rot2 = dyaw_rot2(isnan(dyaw_rot2)==0);
    dyaw_rot2_mean(i,1) = nanmean(dyaw_rot2);
    
    
    %% pitch & roll accels
    pitch_dot_dot_now = pitch_dot_dot(:,i);    
    roll_dot_dot_now = roll_dot_dot_mirror(:,i);    
%     rot_axis_accel_now = rot_axis_accel_mirror(:,i);    
    
    % rot1
    pitch_dot_dot_rot1 = pitch_dot_dot_now(t_now>t_rot1_min & t_now<t_rot1_max);
    roll_dot_dot_rot1 = roll_dot_dot_now(t_now>t_rot1_min & t_now<t_rot1_max);
%     rot_axis_accel_rot1 = rot_axis_accel_now(t_now>t_rot1_min & t_now<t_rot1_max);
    
    pitch_dot_dot_rot1 = pitch_dot_dot_rot1(isnan(pitch_dot_dot_rot1)==0);
    roll_dot_dot_rot1 = roll_dot_dot_rot1(isnan(roll_dot_dot_rot1)==0);
%     rot_axis_accel_rot1 = rot_axis_accel_rot1(isnan(rot_axis_accel_rot1)==0);
    
    pitch_dot_dot_rot1_mean(i,1) = nanmean(pitch_dot_dot_rot1);
    roll_dot_dot_rot1_mean(i,1) = nanmean(roll_dot_dot_rot1);
%     rot_axis_accel_rot1_mean(i,1) = nanmean(rot_axis_accel_rot1);
    
    % rot2
    pitch_dot_dot_rot2 = pitch_dot_dot_now(t_now>t_rot2_min & t_now<t_rot2_max);
    roll_dot_dot_rot2 = roll_dot_dot_now(t_now>t_rot2_min & t_now<t_rot2_max);
%     rot_axis_accel_rot2 = rot_axis_accel_now(t_now>t_rot2_min & t_now<t_rot2_max);
    
    pitch_dot_dot_rot2 = pitch_dot_dot_rot2(isnan(pitch_dot_dot_rot2)==0);
    roll_dot_dot_rot2 = roll_dot_dot_rot2(isnan(roll_dot_dot_rot2)==0);
%     rot_axis_accel_rot2 = rot_axis_accel_rot2(isnan(rot_axis_accel_rot2)==0);
    
    pitch_dot_dot_rot2_mean(i,1) = nanmean(pitch_dot_dot_rot2);
    roll_dot_dot_rot2_mean(i,1) = nanmean(roll_dot_dot_rot2);
%     rot_axis_accel_rot2_mean(i,1) = nanmean(rot_axis_accel_rot2);


    %% TORQUES NORM
    Mpitch_norm_now = Mpitch_norm(:,i);    
    Mroll_norm_now = Mroll_norm(:,i);    
    Myaw_norm_now = Myaw_norm(:,i);    
    M_L_norm_now = M_L_norm(:,i);    
    M_R_norm_now = M_R_norm(:,i);    
    M_axis1_norm_now = M_axis1_norm(:,i);    
    M_axis1normal_norm_now = M_axis1normal_norm(:,i);    
    M_axis2_norm_now = M_axis2_norm(:,i);    
    M_axis2normal_norm_now = M_axis2normal_norm(:,i);    
    
    % t_start < t_stop
    Mpitch_norm_startstop = Mpitch_norm_now(t_now>t_start & t_now<t_stop);
    Mroll_norm_startstop = Mroll_norm_now(t_now>t_start & t_now<t_stop);
    Myaw_norm_startstop = Myaw_norm_now(t_now>t_start & t_now<t_stop);
    M_L_norm_startstop = M_L_norm_now(t_now>t_start & t_now<t_stop);
    M_R_norm_startstop = M_R_norm_now(t_now>t_start & t_now<t_stop);
    M_axis1_norm_startstop = M_axis1_norm_now(t_now>t_start & t_now<t_stop);
    M_axis1normal_norm_startstop = M_axis1normal_norm_now(t_now>t_start & t_now<t_stop);
    M_axis2_norm_startstop = M_axis2_norm_now(t_now>t_start & t_now<t_stop);
    M_axis2normal_norm_startstop = M_axis2normal_norm_now(t_now>t_start & t_now<t_stop);
    
    % 0 < t_stop
    Mpitch_norm_0stop = Mpitch_norm_now(t_now>0 & t_now<t_stop);
    Mroll_norm_0stop = Mroll_norm_now(t_now>0 & t_now<t_stop);
    Myaw_norm_0stop = Myaw_norm_now(t_now>0 & t_now<t_stop);
    M_L_norm_0stop = M_L_norm_now(t_now>0 & t_now<t_stop);
    M_R_norm_0stop = M_R_norm_now(t_now>0 & t_now<t_stop);
    M_axis1_norm_0stop = M_axis1_norm_now(t_now>0 & t_now<t_stop);
    M_axis1normal_norm_0stop = M_axis1normal_norm_now(t_now>0 & t_now<t_stop);
    M_axis2_norm_0stop = M_axis2_norm_now(t_now>0 & t_now<t_stop);
    M_axis2normal_norm_0stop = M_axis2normal_norm_now(t_now>0 & t_now<t_stop);
    
    % MIN&MAX&MEAN TORQUES
    Mroll_norm_max_startstop(i,1) = max(Mroll_norm_startstop);
    Mpitch_norm_max_startstop(i,1) = max(Mpitch_norm_startstop);
    Myaw_norm_max_startstop(i,1) = max(Myaw_norm_startstop);
    M_L_norm_max_startstop(i,1) = max(M_L_norm_startstop);
    M_R_norm_max_startstop(i,1) = max(M_R_norm_startstop);
    M_axis1_norm_max_startstop(i,1) = max(M_axis1_norm_startstop);
    M_axis1normal_norm_max_startstop(i,1) = max(M_axis1normal_norm_startstop);
    M_axis2_norm_max_startstop(i,1) = max(M_axis2_norm_startstop);
    M_axis2normal_norm_max_startstop(i,1) = max(M_axis2normal_norm_startstop);
    
    Mroll_norm_min_startstop(i,1) = min(Mroll_norm_startstop);
    Mpitch_norm_min_startstop(i,1) = min(Mpitch_norm_startstop);
    Myaw_norm_min_startstop(i,1) = min(Myaw_norm_startstop);
    M_L_norm_min_startstop(i,1) = min(M_L_norm_startstop);
    M_R_norm_min_startstop(i,1) = min(M_R_norm_startstop);
    M_axis1_norm_min_startstop(i,1) = min(M_axis1_norm_startstop);
    M_axis1normal_norm_min_startstop(i,1) = min(M_axis1normal_norm_startstop);
    M_axis2_norm_min_startstop(i,1) = min(M_axis2_norm_startstop);
    M_axis2normal_norm_min_startstop(i,1) = min(M_axis2normal_norm_startstop);
    
    Mroll_norm_min_0stop(i,1) = min(Mroll_norm_0stop);
    Mpitch_norm_min_0stop(i,1) = min(Mpitch_norm_0stop);
    Myaw_norm_min_0stop(i,1) = min(Myaw_norm_0stop);
    M_L_norm_min_0stop(i,1) = min(M_L_norm_0stop);
    M_R_norm_min_0stop(i,1) = min(M_R_norm_0stop);
    M_axis1_norm_min_0stop(i,1) = min(M_axis1_norm_0stop);
    M_axis1normal_norm_min_0stop(i,1) = min(M_axis1normal_norm_0stop);
    M_axis2_norm_min_0stop(i,1) = min(M_axis2_norm_0stop);
    M_axis2normal_norm_min_0stop(i,1) = min(M_axis2normal_norm_0stop);
    
    Mroll_norm_mean_startstop(i,1) = nanmean(Mroll_norm_startstop);
    Mpitch_norm_mean_startstop(i,1) = nanmean(Mpitch_norm_startstop);
    Myaw_norm_mean_startstop(i,1) = nanmean(Myaw_norm_startstop);
    M_L_norm_mean_startstop(i,1) = nanmean(M_L_norm_startstop);
    M_R_norm_mean_startstop(i,1) = nanmean(M_R_norm_startstop);
    M_axis1_norm_mean_startstop(i,1) = nanmean(M_axis1_norm_startstop);
    M_axis1normal_norm_mean_startstop(i,1) = nanmean(M_axis1normal_norm_startstop);
    M_axis2_norm_mean_startstop(i,1) = nanmean(M_axis2_norm_startstop);
    M_axis2normal_norm_mean_startstop(i,1) = nanmean(M_axis2normal_norm_startstop);
    
    Mroll_norm_mean_0stop(i,1) = nanmean(Mroll_norm_0stop);
    Mpitch_norm_mean_0stop(i,1) = nanmean(Mpitch_norm_0stop);
    Myaw_norm_mean_0stop(i,1) = nanmean(Myaw_norm_0stop);
    M_L_norm_mean_0stop(i,1) = nanmean(M_L_norm_0stop);
    M_R_norm_mean_0stop(i,1) = nanmean(M_R_norm_0stop);
    M_axis1_norm_mean_0stop(i,1) = nanmean(M_axis1_norm_0stop);
    M_axis1normal_norm_mean_0stop(i,1) = nanmean(M_axis1normal_norm_0stop);
    M_axis2_norm_mean_0stop(i,1) = nanmean(M_axis2_norm_0stop);
    M_axis2normal_norm_mean_0stop(i,1) = nanmean(M_axis2normal_norm_0stop);
    
    %% TORQUES
    Mpitch_now = Mpitch(:,i);    
    Mroll_now = Mroll(:,i);    
    Myaw_now = Myaw(:,i);    
    M_L_now = M_L(:,i);    
    M_R_now = M_R(:,i);    
    M_axis1_now = M_axis1(:,i);    
    M_axis1normal_now = M_axis1normal(:,i);    
    M_axis2_now = M_axis2(:,i);    
    M_axis2normal_now = M_axis2normal(:,i);    
    
    % t_start < t_stop
    Mpitch_startstop = Mpitch_now(t_now>t_start & t_now<t_stop);
    Mroll_startstop = Mroll_now(t_now>t_start & t_now<t_stop);
    Myaw_startstop = Myaw_now(t_now>t_start & t_now<t_stop);
    M_L_startstop = M_L_now(t_now>t_start & t_now<t_stop);
    M_R_startstop = M_R_now(t_now>t_start & t_now<t_stop);
    M_axis1_startstop = M_axis1_now(t_now>t_start & t_now<t_stop);
    M_axis1normal_startstop = M_axis1normal_now(t_now>t_start & t_now<t_stop);
    M_axis2_startstop = M_axis2_now(t_now>t_start & t_now<t_stop);
    M_axis2normal_startstop = M_axis2normal_now(t_now>t_start & t_now<t_stop);
    
    % 0 < t_stop
    Mpitch_0stop = Mpitch_now(t_now>0 & t_now<t_stop);
    Mroll_0stop = Mroll_now(t_now>0 & t_now<t_stop);
    Myaw_0stop = Myaw_now(t_now>0 & t_now<t_stop);
    M_L_0stop = M_L_now(t_now>0 & t_now<t_stop);
    M_R_0stop = M_R_now(t_now>0 & t_now<t_stop);
    M_axis1_0stop = M_axis1_now(t_now>0 & t_now<t_stop);
    M_axis1normal_0stop = M_axis1normal_now(t_now>0 & t_now<t_stop);
    M_axis2_0stop = M_axis2_now(t_now>0 & t_now<t_stop);
    M_axis2normal_0stop = M_axis2normal_now(t_now>0 & t_now<t_stop);
    
    % MIN&MAX&MEAN TORQUES
    Mroll_max_startstop(i,1) = max(Mroll_startstop);
    Mpitch_max_startstop(i,1) = max(Mpitch_startstop);
    Myaw_max_startstop(i,1) = max(Myaw_startstop);
    M_L_max_startstop(i,1) = max(M_L_startstop);
    M_R_max_startstop(i,1) = max(M_R_startstop);
    M_axis1_max_startstop(i,1) = max(M_axis1_startstop);
    M_axis1normal_max_startstop(i,1) = max(M_axis1normal_startstop);
    M_axis2_max_startstop(i,1) = max(M_axis2_startstop);
    M_axis2normal_max_startstop(i,1) = max(M_axis2normal_startstop);
    
    Mroll_min_startstop(i,1) = min(Mroll_startstop);
    Mpitch_min_startstop(i,1) = min(Mpitch_startstop);
    Myaw_min_startstop(i,1) = min(Myaw_startstop);
    M_L_min_startstop(i,1) = min(M_L_startstop);
    M_R_min_startstop(i,1) = min(M_R_startstop);
    M_axis1_min_startstop(i,1) = min(M_axis1_startstop);
    M_axis1normal_min_startstop(i,1) = min(M_axis1normal_startstop);
    M_axis2_min_startstop(i,1) = min(M_axis2_startstop);
    M_axis2normal_min_startstop(i,1) = min(M_axis2normal_startstop);
    
    Mroll_min_0stop(i,1) = min(Mroll_0stop);
    Mpitch_min_0stop(i,1) = min(Mpitch_0stop);
    Myaw_min_0stop(i,1) = min(Myaw_0stop);
    M_L_min_0stop(i,1) = min(M_L_0stop);
    M_R_min_0stop(i,1) = min(M_R_0stop);
    M_axis1_min_0stop(i,1) = min(M_axis1_0stop);
    M_axis1normal_min_0stop(i,1) = min(M_axis1normal_0stop);
    M_axis2_min_0stop(i,1) = min(M_axis2_0stop);
    M_axis2normal_min_0stop(i,1) = min(M_axis2normal_0stop);
    
    Mroll_mean_startstop(i,1) = nanmean(Mroll_startstop);
    Mpitch_mean_startstop(i,1) = nanmean(Mpitch_startstop);
    Myaw_mean_startstop(i,1) = nanmean(Myaw_startstop);
    M_L_mean_startstop(i,1) = nanmean(M_L_startstop);
    M_R_mean_startstop(i,1) = nanmean(M_R_startstop);
    M_axis1_mean_startstop(i,1) = nanmean(M_axis1_startstop);
    M_axis1normal_mean_startstop(i,1) = nanmean(M_axis1normal_startstop);    
    M_axis2_mean_startstop(i,1) = nanmean(M_axis2_startstop);
    M_axis2normal_mean_startstop(i,1) = nanmean(M_axis2normal_startstop);
    
    Mroll_mean_0stop(i,1) = nanmean(Mroll_0stop);
    Mpitch_mean_0stop(i,1) = nanmean(Mpitch_0stop);
    Myaw_mean_0stop(i,1) = nanmean(Myaw_0stop);
    M_L_mean_0stop(i,1) = nanmean(M_L_0stop);
    M_R_mean_0stop(i,1) = nanmean(M_R_0stop);
    M_axis1_mean_0stop(i,1) = nanmean(M_axis1_0stop);
    M_axis1normal_mean_0stop(i,1) = nanmean(M_axis1normal_0stop);
    M_axis2_mean_0stop(i,1) = nanmean(M_axis2_0stop);
    M_axis2normal_mean_0stop(i,1) = nanmean(M_axis2normal_0stop);
    
    %% TORQUE AXES
    % rot1
    Mpitch_rot1 = Mpitch_now(t_now>t_rot1_min & t_now<t_rot1_max);
    Mroll_rot1 = Mroll_now(t_now>t_rot1_min & t_now<t_rot1_max);
    
    Mpitch_rot1 = Mpitch_rot1(isnan(Mpitch_rot1)==0);
    Mroll_rot1 = Mroll_rot1(isnan(Mroll_rot1)==0);
    
    Mpitch_rot1_mean(i,1) = nanmean(Mpitch_rot1);
    Mroll_rot1_mean(i,1) = nanmean(Mroll_rot1);
    
    % rot2
    Mpitch_rot2 = Mpitch_now(t_now>t_rot2_min & t_now<t_rot2_max);
    Mroll_rot2 = Mroll_now(t_now>t_rot2_min & t_now<t_rot2_max);
    
    Mpitch_rot2 = Mpitch_rot2(isnan(Mpitch_rot2)==0);
    Mroll_rot2 = Mroll_rot2(isnan(Mroll_rot2)==0);
    
    Mpitch_rot2_mean(i,1) = nanmean(Mpitch_rot2);
    Mroll_rot2_mean(i,1) = nanmean(Mroll_rot2);
    
    %% Torque max & min
    
end

% yaw turns 
Dyaw_rot1rot2_mean = dyaw_rot2_mean - dyaw_rot1_mean;
Dyaw_rot1rot2_mean_mean = circ_mean_deg_nonan(Dyaw_rot1rot2_mean);

% mean rot axes angle
rot_axis_accel_rot1_mean = atan2d(pitch_dot_dot_rot1_mean,roll_dot_dot_rot1_mean);
rot_axis_accel_rot1_mean_mean = circ_mean_deg_nonan(rot_axis_accel_rot1_mean);

rot_axis_accel_rot2_mean = atan2d(pitch_dot_dot_rot2_mean,roll_dot_dot_rot2_mean);
rot_axis_accel_rot2_mean_mean = circ_mean_deg_nonan(rot_axis_accel_rot2_mean);

Drot_axis_accel_rot1rot2_mean_mean = rot_axis_accel_rot1_mean_mean - rot_axis_accel_rot2_mean_mean;
rot_axis_accel_rot1rot2_mean_mean = circ_mean_deg_nonan([rot_axis_accel_rot1_mean;rot_axis_accel_rot2_mean+180]);

% mean torque axes angle
Maxis_rot1_mean = atan2d(Mpitch_rot1_mean,Mroll_rot1_mean);
Maxis_rot1_mean_mean = circ_mean_deg_nonan(Maxis_rot1_mean);

Maxis_rot2_mean = atan2d(Mpitch_rot2_mean,Mroll_rot2_mean);
Maxis_rot2_mean_mean = circ_mean_deg_nonan(Maxis_rot2_mean);

DMaxis_rot1rot2_mean_mean = Maxis_rot1_mean_mean - Maxis_rot2_mean_mean
Maxis_rot1rot2_mean_mean = circ_mean_deg_nonan([Maxis_rot1_mean;Maxis_rot2_mean+180])


