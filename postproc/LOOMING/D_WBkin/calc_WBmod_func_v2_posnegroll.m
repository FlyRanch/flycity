n = 100;
t = [0:1/(n-1):1]';
stroke_wb_steady = fnval(stroke_steady_func,t);
pitch_wb_steady = fnval(pitch_steady_func,t);
dev_wb_steady = fnval(dev_steady_func,t);

limit = .5;
roll_dot_dot_limit = limit*nanstd(roll_dot_dot_mean_wb(:));
pitch_dot_dot_limit = limit*nanstd(pitch_dot_dot_mean_wb(:));
yaw_dot_dot_limit = limit*nanstd(yaw_dot_dot_mean_wb(:));
F_limit = limit*nanstd(F_mean_wb(:));

roll_dot_dot_max = max(abs(roll_dot_dot_mean_wb(:)));
pitch_dot_dot_max = max(abs(pitch_dot_dot_mean_wb(:)));
yaw_dot_dot_max = max(abs(yaw_dot_dot_mean_wb(:)));
F_max = max(abs(F_mean_wb(:)));

        t_roll = [];
        t_roll_LR = [];

        t_yaw = [];
        t_yaw_LR = [];

        t_pitch = [];
        t_pitch_LR = [];

        t_F = [];
        t_F_LR = [];

        Dstroke_rollaccel_L = [];
        Dstroke_rollaccel_R = [];
        Dpitch_rollaccel_L = [];
        Dpitch_rollaccel_R = [];
        Ddev_rollaccel_L = [];
        Ddev_rollaccel_R = [];
        
        Ddstroke_rollaccel = [];
        Ddpitch_rollaccel = [];
        Dddev_rollaccel = [];
        
        Dstroke_yawaccel_L = [];
        Dstroke_yawaccel_R = [];
        Dpitch_yawaccel_L = [];
        Dpitch_yawaccel_R = [];
        Ddev_yawaccel_L = [];
        Ddev_yawaccel_R = [];
        
        Ddstroke_yawaccel = [];
        Ddpitch_yawaccel = [];
        Dddev_yawaccel = [];
        
        Dstroke_pitchaccel = [];
        Dpitch_pitchaccel = [];
        Ddev_pitchaccel = [];
        
        Dstroke_F = [];
        Dpitch_F = [];
        Ddev_F = [];

for seq = 1:size(n_down_start_L,2)
    counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        roll_dot_dot_mean_wb_now = roll_dot_dot_mean_wb(wb,seq);
        pitch_dot_dot_mean_wb_now = pitch_dot_dot_mean_wb(wb,seq);
        yaw_dot_dot_mean_wb_now = yaw_dot_dot_mean_wb(wb,seq);
        F_mean_wb_now = F_mean_wb(wb,seq);

        n_down_start_L_now = n_down_start_L(wb,seq);
        n_up_stop_L_now = n_down_start_L(wb+1,seq);
        n_down_start_R_now = n_down_start_R(wb,seq);
        n_up_stop_R_now = n_down_start_R(wb+1,seq);
        
        if isnan(n_down_start_L_now)==0 && isnan(n_up_stop_L_now)==0

        stroke_wb_L_now = stroke_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
        stroke_wb_R_now = stroke_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
        pitch_wb_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq));
        pitch_wb_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq));
        dev_wb_L_now = dev_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
        dev_wb_R_now = dev_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);

        t_L = [1:(length(stroke_wb_L_now)-1)/(n-1):length(stroke_wb_L_now)]';
        t_R = [1:(length(stroke_wb_R_now)-1)/(n-1):length(stroke_wb_R_now)]';

        stroke_wb_L_interp = rad2deg(interp1(stroke_wb_L_now,t_L));
        stroke_wb_R_interp = rad2deg(interp1(stroke_wb_R_now,t_R));
        pitch_wb_L_interp = rad2deg(interp1(pitch_wb_L_now,t_L));
        pitch_wb_R_interp = rad2deg(interp1(pitch_wb_R_now,t_R));
        dev_wb_L_interp = rad2deg(interp1(dev_wb_L_now,t_L));
        dev_wb_R_interp = rad2deg(interp1(dev_wb_R_now,t_R));

        Dstroke_wb_L = stroke_wb_L_interp - stroke_wb_steady;
        Dpitch_wb_L = pitch_wb_L_interp - pitch_wb_steady;
        Ddev_wb_L = dev_wb_L_interp - dev_wb_steady;
        Dstroke_wb_R = stroke_wb_R_interp - stroke_wb_steady;
        Dpitch_wb_R = pitch_wb_R_interp - pitch_wb_steady;
        Ddev_wb_R = dev_wb_R_interp - dev_wb_steady;
        
        
        
        % store data
        if abs(roll_dot_dot_mean_wb_now) > abs(roll_dot_dot_limit)
            if roll_dot_dot_mean_wb_now > 0
                t_roll(end+1:end+length(t),1) = t;
                t_roll_LR(end+1:end+length(t),1) = t;
                t_roll_LR(end+1:end+length(t),1) = t;

                Dstroke_rollaccel_L(end+1:end+length(t),1) = Dstroke_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dstroke_rollaccel_R(end+1:end+length(t),1) = Dstroke_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddev_rollaccel_L(end+1:end+length(t),1) = Ddev_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddev_rollaccel_R(end+1:end+length(t),1) = Ddev_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dpitch_rollaccel_L(end+1:end+length(t),1) = Dpitch_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dpitch_rollaccel_R(end+1:end+length(t),1) = Dpitch_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;

                Ddstroke_rollaccel(end+1:end+length(t),1) = (Dstroke_wb_L-Dstroke_wb_R) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddpitch_rollaccel(end+1:end+length(t),1) = (Dpitch_wb_L-Dpitch_wb_R) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dddev_rollaccel(end+1:end+length(t),1) = (Ddev_wb_L-Ddev_wb_R) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
            else
                t_roll(end+1:end+length(t),1) = t;
                t_roll_LR(end+1:end+length(t),1) = t;
                t_roll_LR(end+1:end+length(t),1) = t;

                Dstroke_rollaccel_L(end+1:end+length(t),1) = Dstroke_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dstroke_rollaccel_R(end+1:end+length(t),1) = Dstroke_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddev_rollaccel_L(end+1:end+length(t),1) = Ddev_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddev_rollaccel_R(end+1:end+length(t),1) = Ddev_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dpitch_rollaccel_L(end+1:end+length(t),1) = Dpitch_wb_R / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dpitch_rollaccel_R(end+1:end+length(t),1) = Dpitch_wb_L / roll_dot_dot_mean_wb_now * roll_dot_dot_max;

                Ddstroke_rollaccel(end+1:end+length(t),1) = (Dstroke_wb_R-Dstroke_wb_L) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Ddpitch_rollaccel(end+1:end+length(t),1) = (Dpitch_wb_R-Dpitch_wb_L) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
                Dddev_rollaccel(end+1:end+length(t),1) = (Ddev_wb_R-Ddev_wb_L) / roll_dot_dot_mean_wb_now * roll_dot_dot_max;
            end
        end
        
        if abs(yaw_dot_dot_mean_wb_now) > abs(yaw_dot_dot_limit)
            t_yaw(end+1:end+length(t),1) = t;
            t_yaw_LR(end+1:end+length(t),1) = t;
            t_yaw_LR(end+1:end+length(t),1) = t;

            Dstroke_yawaccel_L(end+1:end+length(t),1) = Dstroke_wb_L / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Dstroke_yawaccel_R(end+1:end+length(t),1) = Dstroke_wb_R / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Ddev_yawaccel_L(end+1:end+length(t),1) = Ddev_wb_L / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Ddev_yawaccel_R(end+1:end+length(t),1) = Ddev_wb_R / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Dpitch_yawaccel_L(end+1:end+length(t),1) = Dpitch_wb_L / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Dpitch_yawaccel_R(end+1:end+length(t),1) = Dpitch_wb_R / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            
            Ddstroke_yawaccel(end+1:end+length(t),1) = (Dstroke_wb_L-Dstroke_wb_R) / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Ddpitch_yawaccel(end+1:end+length(t),1) = (Dpitch_wb_L-Dpitch_wb_R) / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
            Dddev_yawaccel(end+1:end+length(t),1) = (Ddev_wb_L-Ddev_wb_R) / yaw_dot_dot_mean_wb_now * yaw_dot_dot_max;
        end
        
        if abs(pitch_dot_dot_mean_wb_now) > abs(pitch_dot_dot_limit)
            t_pitch(end+1:end+length(t),1) = t;
            t_pitch_LR(end+1:end+length(t),1) = t;
            t_pitch_LR(end+1:end+length(t),1) = t;

            Dstroke_pitchaccel(end+1:end+length(t),1) = Dstroke_wb_L / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
            Dstroke_pitchaccel(end+1:end+length(t),1) = Dstroke_wb_R / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
            Ddev_pitchaccel(end+1:end+length(t),1) = Ddev_wb_L / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
            Ddev_pitchaccel(end+1:end+length(t),1) = Ddev_wb_R / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
            Dpitch_pitchaccel(end+1:end+length(t),1) = Dpitch_wb_L / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
            Dpitch_pitchaccel(end+1:end+length(t),1) = Dpitch_wb_R / pitch_dot_dot_mean_wb_now * pitch_dot_dot_max;
        end
        
        if abs(F_mean_wb_now) > abs(F_limit) + 1
            t_F(end+1:end+length(t),1) = t;
            t_F_LR(end+1:end+length(t),1) = t;
            t_F_LR(end+1:end+length(t),1) = t;

            Dstroke_F(end+1:end+length(t),1) = Dstroke_wb_L / F_mean_wb_now * F_max;
            Dstroke_F(end+1:end+length(t),1) = Dstroke_wb_R / F_mean_wb_now * F_max;
            Ddev_F(end+1:end+length(t),1) = Ddev_wb_L / F_mean_wb_now * F_max;
            Ddev_F(end+1:end+length(t),1) = Ddev_wb_R / F_mean_wb_now * F_max;
            Dpitch_F(end+1:end+length(t),1) = Dpitch_wb_L / F_mean_wb_now * F_max;
            Dpitch_F(end+1:end+length(t),1) = Dpitch_wb_R / F_mean_wb_now * F_max;
        end
        
        
        end
    end
end

save('WBmod_data.mat','stroke_wb_steady','pitch_wb_steady','dev_wb_steady',...
    'roll_dot_dot_max','yaw_dot_dot_max','pitch_dot_dot_max','F_max',...
    't_roll','t_roll_LR','t_yaw','t_yaw_LR','t_pitch','t_pitch_LR','t_F','t_F_LR',...
	'Dstroke_rollaccel_L','Dstroke_yawaccel_L','Dstroke_rollaccel_R','Dstroke_yawaccel_R','Dstroke_pitchaccel','Dstroke_F',...
    'Dpitch_rollaccel_L','Dpitch_yawaccel_L','Dpitch_rollaccel_R','Dpitch_yawaccel_R','Dpitch_pitchaccel','Dpitch_F',...
    'Ddev_rollaccel_L','Ddev_yawaccel_L','Ddev_rollaccel_R','Ddev_yawaccel_R','Ddev_pitchaccel','Ddev_F');

%% plot heatmaps

% bins
nx=50;
ny=100;
binx = 0:1/(nx-1):1;

%% roll accel
figure
plot_WBmod_rollaccel

%% yaw accel
figure
plot_WBmod_yawaccel

%% pitch accel & force
figure
plot_WBmod_pitchaccelNforce

