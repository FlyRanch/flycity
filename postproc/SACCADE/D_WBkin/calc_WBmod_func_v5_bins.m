t = [0:1/(n-1):1]';

stroke_wb_steady = stroke_wb_steady_bins_meanCIstd(:,1);
pitch_wb_steady = pitch_wb_steady_bins_meanCIstd(:,1);
dev_wb_steady = dev_wb_steady_bins_meanCIstd(:,1);

%     limit = .5;
%     rollaccel_limit = limit*nanstd(roll_dot_dot_mean_wb(:));
%     pitchaccel_limit = limit*nanstd(pitch_dot_dot_mean_wb(:));
%     yawaccel_limit = limit*nanstd(yaw_dot_dot_mean_wb(:));
%     Fenhance_limit = limit*nanstd(F_mean_wb(:));

rollaccel_3std = 3*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_3std = 3*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_3std = 3*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_3std = 3*nanstd(F_mean_wb(:));

        t_roll = [];
        t_yaw = [];
        t_pitch = [];
        t_F = [];

        Dstroke_rollaccel_up = [];
        Dstroke_rollaccel_down = [];
        Dpitch_rollaccel_up = [];
        Dpitch_rollaccel_down = [];
        Ddev_rollaccel_up = [];
        Ddev_rollaccel_down = [];
        
        Ddstroke_rollaccel = [];
        Ddpitch_rollaccel = [];
        Dddev_rollaccel = [];
        
        Dstroke_yawaccel_right = [];
        Dstroke_yawaccel_left = [];
        Dpitch_yawaccel_right = [];
        Dpitch_yawaccel_left = [];
        Ddev_yawaccel_right = [];
        Ddev_yawaccel_left = [];
        
        Ddstroke_yawaccel = [];
        Ddpitch_yawaccel = [];
        Dddev_yawaccel = [];
        
        Dstroke_pitchaccel = [];
        Dpitch_pitchaccel = [];
        Ddev_pitchaccel = [];
        
        
        Dstroke_rollaccel_up_active = [];
        Dstroke_rollaccel_down_active = [];
        Dpitch_rollaccel_up_active = [];
        Dpitch_rollaccel_down_active = [];
        Ddev_rollaccel_up_active = [];
        Ddev_rollaccel_down_active = [];
        
        Ddstroke_rollaccel_active = [];
        Ddpitch_rollaccel_active = [];
        Dddev_rollaccel_active = [];
        
        Dstroke_yawaccel_right_active = [];
        Dstroke_yawaccel_left_active = [];
        Dpitch_yawaccel_right_active = [];
        Dpitch_yawaccel_left_active = [];
        Ddev_yawaccel_right_active = [];
        Ddev_yawaccel_left_active = [];
        
        Ddstroke_yawaccel_active = [];
        Ddpitch_yawaccel_active = [];
        Dddev_yawaccel_active = [];
        
        Dstroke_pitchaccel_active = [];
        Dpitch_pitchaccel_active = [];
        Ddev_pitchaccel_active = [];
        
        
        
        Dstroke_F = [];
        Dpitch_F = [];
        Ddev_F = [];
        
        
        

% for seq = 1:size(n_down_start_L,2)
%     counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        rollaccel_mean_wb_now = roll_dot_dot_mean_wb(wb);
        pitchaccel_mean_wb_now = pitch_dot_dot_mean_wb(wb);
        yawaccel_mean_wb_now = yaw_dot_dot_mean_wb(wb);
        Fenhance_mean_wb_now = F_mean_wb(wb);

        rollvel_mean_wb_now = roll_dot_mean_wb(wb);
        pitchvel_mean_wb_now = pitch_dot_mean_wb(wb);
        yawvel_mean_wb_now = yaw_dot_mean_wb(wb);
        
        % active roll
        rollaccel_mean_wb_now_active = rollaccel_mean_wb_now + k_damp_roll*rollvel_mean_wb_now;
        pitchaccel_mean_wb_now_active = pitchaccel_mean_wb_now + k_damp_pitch*pitchvel_mean_wb_now;
        yawaccel_mean_wb_now_active = yawaccel_mean_wb_now + k_damp_yaw*yawvel_mean_wb_now;

        n_down_start_L_now = n_down_start_L(wb);
        n_up_stop_L_now = n_down_start_L(wb+1);
        n_down_start_R_now = n_down_start_R(wb);
        n_up_stop_R_now = n_down_start_R(wb+1);
        
        if isnan(n_down_start_L_now)==0 && isnan(n_up_stop_L_now)==0

        stroke_wb_L_now = stroke_L_mirror(n_down_start_L_now:n_up_stop_L_now);
        stroke_wb_R_now = stroke_R_mirror(n_down_start_R_now:n_up_stop_R_now);
        pitch_wb_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_up_stop_L_now));
        pitch_wb_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_up_stop_R_now));
        dev_wb_L_now = dev_L_mirror(n_down_start_L_now:n_up_stop_L_now);
        dev_wb_R_now = dev_R_mirror(n_down_start_R_now:n_up_stop_R_now);

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
        if abs(rollaccel_mean_wb_now) > abs(rollaccel_limit)
            if rollaccel_mean_wb_now > 0
                t_roll(:,end+1) = t;

                Dstroke_rollaccel_up(:,end+1) = Dstroke_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dstroke_rollaccel_down(:,end+1) = Dstroke_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddev_rollaccel_up(:,end+1) = Ddev_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddev_rollaccel_down(:,end+1) = Ddev_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dpitch_rollaccel_up(:,end+1) = Dpitch_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dpitch_rollaccel_down(:,end+1) = Dpitch_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;

                Ddstroke_rollaccel(:,end+1) = (Dstroke_wb_L-Dstroke_wb_R) / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddpitch_rollaccel(:,end+1) = (Dpitch_wb_L-Dpitch_wb_R) / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dddev_rollaccel(:,end+1) = (Ddev_wb_L-Ddev_wb_R) / abs(rollaccel_mean_wb_now) * rollaccel_3std;

                % active roll
                Dstroke_rollaccel_up_active(:,end+1) = Dstroke_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dstroke_rollaccel_down_active(:,end+1) = Dstroke_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddev_rollaccel_up_active(:,end+1) = Ddev_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddev_rollaccel_down_active(:,end+1) = Ddev_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dpitch_rollaccel_up_active(:,end+1) = Dpitch_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dpitch_rollaccel_down_active(:,end+1) = Dpitch_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;

                Ddstroke_rollaccel_active(:,end+1) = (Dstroke_wb_L-Dstroke_wb_R) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddpitch_rollaccel_active(:,end+1) = (Dpitch_wb_L-Dpitch_wb_R) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dddev_rollaccel_active(:,end+1) = (Ddev_wb_L-Ddev_wb_R) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                
%                 hold off
%                 plot(Ddev_rollaccel_up(:,end),'r')
%                 hold on
%                 plot(Ddev_rollaccel_down(:,end),'m')
%                 plot(Ddev_rollaccel_up_active(:,end),'b')
%                 plot(Ddev_rollaccel_down_active(:,end),'c')
%                 pause
                
            else
                t_roll(:,end+1) = t;

                Dstroke_rollaccel_up(:,end+1) = Dstroke_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dstroke_rollaccel_down(:,end+1) = Dstroke_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddev_rollaccel_up(:,end+1) = Ddev_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddev_rollaccel_down(:,end+1) = Ddev_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dpitch_rollaccel_up(:,end+1) = Dpitch_wb_R / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dpitch_rollaccel_down(:,end+1) = Dpitch_wb_L / abs(rollaccel_mean_wb_now) * rollaccel_3std;

                Ddstroke_rollaccel(:,end+1) = (Dstroke_wb_R-Dstroke_wb_L) / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Ddpitch_rollaccel(:,end+1) = (Dpitch_wb_R-Dpitch_wb_L) / abs(rollaccel_mean_wb_now) * rollaccel_3std;
                Dddev_rollaccel(:,end+1) = (Ddev_wb_R-Ddev_wb_L) / abs(rollaccel_mean_wb_now) * rollaccel_3std;

                % active roll
                Dstroke_rollaccel_up_active(:,end+1) = Dstroke_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dstroke_rollaccel_down_active(:,end+1) = Dstroke_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddev_rollaccel_up_active(:,end+1) = Ddev_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddev_rollaccel_down_active(:,end+1) = Ddev_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dpitch_rollaccel_up_active(:,end+1) = Dpitch_wb_R / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dpitch_rollaccel_down_active(:,end+1) = Dpitch_wb_L / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;

                Ddstroke_rollaccel_active(:,end+1) = (Dstroke_wb_R-Dstroke_wb_L) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Ddpitch_rollaccel_active(:,end+1) = (Dpitch_wb_R-Dpitch_wb_L) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
                Dddev_rollaccel_active(:,end+1) = (Ddev_wb_R-Ddev_wb_L) / abs(rollaccel_mean_wb_now_active) * rollaccel_3std;
            end
        end
        
        if abs(yawaccel_mean_wb_now) > abs(yawaccel_limit)
            if yawaccel_mean_wb_now > 0
                t_yaw(:,end+1) = t;

                Dstroke_yawaccel_right(:,end+1) = Dstroke_wb_L / yawaccel_mean_wb_now * yawaccel_3std;
                Dstroke_yawaccel_left(:,end+1) = Dstroke_wb_R / yawaccel_mean_wb_now * yawaccel_3std;
                Ddev_yawaccel_right(:,end+1) = Ddev_wb_L / yawaccel_mean_wb_now * yawaccel_3std;
                Ddev_yawaccel_left(:,end+1) = Ddev_wb_R / yawaccel_mean_wb_now * yawaccel_3std;
                Dpitch_yawaccel_right(:,end+1) = Dpitch_wb_L / yawaccel_mean_wb_now * yawaccel_3std;
                Dpitch_yawaccel_left(:,end+1) = Dpitch_wb_R / yawaccel_mean_wb_now * yawaccel_3std;

                Ddstroke_yawaccel(:,end+1) = (Dstroke_wb_L-Dstroke_wb_R) / yawaccel_mean_wb_now * yawaccel_3std;
                Ddpitch_yawaccel(:,end+1) = (Dpitch_wb_L-Dpitch_wb_R) / yawaccel_mean_wb_now * yawaccel_3std;
                Dddev_yawaccel(:,end+1) = (Ddev_wb_L-Ddev_wb_R) / yawaccel_mean_wb_now * yawaccel_3std;
            else
                t_yaw(:,end+1) = t;

                Dstroke_yawaccel_right(:,end+1) = Dstroke_wb_R / yawaccel_mean_wb_now * yawaccel_3std;
                Dstroke_yawaccel_left(:,end+1) = Dstroke_wb_L / yawaccel_mean_wb_now * yawaccel_3std;
                Ddev_yawaccel_right(:,end+1) = Ddev_wb_R / yawaccel_mean_wb_now * yawaccel_3std;
                Ddev_yawaccel_left(:,end+1) = Ddev_wb_L / yawaccel_mean_wb_now * yawaccel_3std;
                Dpitch_yawaccel_right(:,end+1) = Dpitch_wb_R / yawaccel_mean_wb_now * yawaccel_3std;
                Dpitch_yawaccel_left(:,end+1) = Dpitch_wb_L / yawaccel_mean_wb_now * yawaccel_3std;

                Ddstroke_yawaccel(:,end+1) = (Dstroke_wb_R-Dstroke_wb_L) / yawaccel_mean_wb_now * yawaccel_3std;
                Ddpitch_yawaccel(:,end+1) = (Dpitch_wb_R-Dpitch_wb_L) / yawaccel_mean_wb_now * yawaccel_3std;
                Dddev_yawaccel(:,end+1) = (Ddev_wb_R-Ddev_wb_L) / yawaccel_mean_wb_now * yawaccel_3std;
            end                
        end
        
        if abs(pitchaccel_mean_wb_now) > abs(pitchaccel_limit)
            t_pitch(:,end+1) = t;
            t_pitch(:,end+1) = t;

            Dstroke_pitchaccel(:,end+1) = Dstroke_wb_L / pitchaccel_mean_wb_now * pitchaccel_3std;
            Dstroke_pitchaccel(:,end+1) = Dstroke_wb_R / pitchaccel_mean_wb_now * pitchaccel_3std;
            Ddev_pitchaccel(:,end+1) = Ddev_wb_L / pitchaccel_mean_wb_now * pitchaccel_3std;
            Ddev_pitchaccel(:,end+1) = Ddev_wb_R / pitchaccel_mean_wb_now * pitchaccel_3std;
            Dpitch_pitchaccel(:,end+1) = Dpitch_wb_L / pitchaccel_mean_wb_now * pitchaccel_3std;
            Dpitch_pitchaccel(:,end+1) = Dpitch_wb_R / pitchaccel_mean_wb_now * pitchaccel_3std;
        end
        
        if abs(Fenhance_mean_wb_now) > (1+Fenhance_limit)
            t_F(:,end+1) = t;
            t_F(:,end+1) = t;

            Dstroke_F(:,end+1) = Dstroke_wb_L / Fenhance_mean_wb_now * (1+Fenhance_3std);
            Dstroke_F(:,end+1) = Dstroke_wb_R / Fenhance_mean_wb_now * (1+Fenhance_3std);
            Ddev_F(:,end+1) = Ddev_wb_L / Fenhance_mean_wb_now * (1+Fenhance_3std);
            Ddev_F(:,end+1) = Ddev_wb_R / Fenhance_mean_wb_now * (1+Fenhance_3std);
            Dpitch_F(:,end+1) = Dpitch_wb_L / Fenhance_mean_wb_now * (1+Fenhance_3std);
            Dpitch_F(:,end+1) = Dpitch_wb_R / Fenhance_mean_wb_now * (1+Fenhance_3std);
        end
        
        
        end
    end
% end

% save('WBmod_data.mat','stroke_wb_steady','pitch_wb_steady','dev_wb_steady',...
%     'rollaccel_3std','yawaccel_3std','pitchaccel_3std','Fenhance_3std',...
%     't_roll','t_roll_LR','t_yaw','t_yaw_LR','t_pitch','t_pitch_LR','t_F','t_F_LR',...
% 	'Dstroke_rollaccel_L','Dstroke_yawaccel_L','Dstroke_rollaccel_R','Dstroke_yawaccel_R','Dstroke_pitchaccel','Dstroke_F',...
%     'Dpitch_rollaccel_L','Dpitch_yawaccel_L','Dpitch_rollaccel_R','Dpitch_yawaccel_R','Dpitch_pitchaccel','Dpitch_F',...
%     'Ddev_rollaccel_L','Ddev_yawaccel_L','Ddev_rollaccel_R','Ddev_yawaccel_R','Ddev_pitchaccel','Ddev_F');

%% plot heatmaps

% bins
nx=n;
ny=100;
binx = 0:1/(nx-1):1;

% calc Fenhanced data
calc_WBmod_Fenhance

%% roll accel
figure
plot_WBmod_rollaccel_updown
saveas(gca,'WBmod_rollaccel_updown.fig')
saveas(gca,'WBmod_rollaccel_updown.png')
plot2svg('WBmod_rollaccel_updown.svg')

%% yaw accel
figure
plot_WBmod_yawaccel_lr
saveas(gca,'WBmod_yawaccel_lr.fig')
saveas(gca,'WBmod_yawaccel_lr.png')
plot2svg('WBmod_yawaccel_lr.svg')

%% pitch accel & force
figure
plot_WBmod_pitchaccelNforce
saveas(gca,'WBmod_pitchaccelNforce.fig')
saveas(gca,'WBmod_pitchaccelNforce.png')
plot2svg('WBmod_pitchaccelNforce.svg')


