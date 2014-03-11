

rot_axis_pos_mirror = atan2d(pitch,roll_mirror);
rot_axis_vel_mirror = atan2d(pitch_dot,roll_dot_mirror);
rot_axis_accel_mirror = atan2d(pitch_dot_dot,roll_dot_dot_mirror);

% no angles < -90
rot_axis_pos_plot = rot_axis_pos_mirror;
rot_axis_pos_plot(rot_axis_pos_plot<-90) = rot_axis_pos_plot(rot_axis_pos_plot<-90)+360;

rot_axis_vel_plot = rot_axis_vel_mirror;
rot_axis_vel_plot(rot_axis_vel_plot<-90) = rot_axis_vel_plot(rot_axis_vel_plot<-90)+360;

rot_axis_accel_plot = rot_axis_accel_mirror;
rot_axis_accel_plot(rot_axis_accel_plot<-90) = rot_axis_accel_plot(rot_axis_accel_plot<-90)+360;

% remove jumps from plots (+/-180deg)
rot_axis_pos_temp = rot_axis_pos_plot;
rot_axis_vel_temp = rot_axis_vel_plot;
rot_axis_accel_temp = rot_axis_accel_plot;

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
    
end

% yaw turns 
Dyaw_rot1rot2_mean = dyaw_rot2_mean - dyaw_rot1_mean;
Dyaw_rot1rot2_mean_mean = circ_mean_deg_nonan(Dyaw_rot1rot2_mean);

% mean rot axes angle
rot_axis_accel_rot1_mean = atan2d(pitch_dot_dot_rot1_mean,roll_dot_dot_rot1_mean);
rot_axis_accel_rot1_mean_mean = circ_mean_deg_nonan(rot_axis_accel_rot1_mean);

% mean rot axes angle
rot_axis_accel_rot2_mean = atan2d(pitch_dot_dot_rot2_mean,roll_dot_dot_rot2_mean);
rot_axis_accel_rot2_mean_mean = circ_mean_deg_nonan(rot_axis_accel_rot2_mean);

Drot_axis_accel_rot1rot2_mean_mean = rot_axis_accel_rot1_mean_mean - rot_axis_accel_rot2_mean_mean;
rot_axis_accel_rot1rot2_mean_mean = circ_mean_deg_nonan([rot_axis_accel_rot1_mean;rot_axis_accel_rot2_mean+180]);


