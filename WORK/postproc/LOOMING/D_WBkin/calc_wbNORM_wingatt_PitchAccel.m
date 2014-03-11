% calc wb MOD Pitch Accel

n_now=0;

%% variables
% wbMOD
freqMOD_wb_PitchAccel = [];

strokeMOD_wb_L_PitchAccel_bins = [];
strokeMOD_wb_R_PitchAccel_bins = [];
strokeMOD_ds_L_PitchAccel_bins = [];
strokeMOD_ds_R_PitchAccel_bins = [];
strokeMOD_us_L_PitchAccel_bins = [];
strokeMOD_us_R_PitchAccel_bins = [];

pitchMOD_wb_L_PitchAccel_bins = [];
pitchMOD_wb_R_PitchAccel_bins = [];
pitchMOD_ds_L_PitchAccel_bins = [];
pitchMOD_ds_R_PitchAccel_bins = [];
pitchMOD_us_L_PitchAccel_bins = [];
pitchMOD_us_R_PitchAccel_bins = [];

devMOD_wb_L_PitchAccel_bins = [];
devMOD_wb_R_PitchAccel_bins = [];
devMOD_ds_L_PitchAccel_bins = [];
devMOD_ds_R_PitchAccel_bins = [];
devMOD_us_L_PitchAccel_bins = [];
devMOD_us_R_PitchAccel_bins = [];

% seq&wb
seq_nr_PitchAccel = [];
wb_nr_PitchAccel = [];

% body kin
V_PitchAccel = [];
pitch_global_PitchAccel = [];
F_PitchAccel = [];

rollaccel_PitchAccel = [];
pitchaccel_PitchAccel = [];
yawaccel_PitchAccel = [];

rollvel_PitchAccel = [];
pitchvel_PitchAccel = [];
yawvel_PitchAccel = [];

% wingbeat kin
dt_wb_PitchAccel = [];
dt_ds_PitchAccel = [];
dt_us_PitchAccel = [];
f_wb_PitchAccel = [];
Rds_PitchAccel = [];

t_wb_L_PitchAccel = [];
t_wb_R_PitchAccel = [];
stroke_wb_L_PitchAccel = [];
stroke_wb_R_PitchAccel = [];
pitch_wb_L_PitchAccel = [];
pitch_wb_R_PitchAccel = [];
dev_wb_L_PitchAccel = [];
dev_wb_R_PitchAccel = [];
aoa_wb_L_PitchAccel = [];
aoa_wb_R_PitchAccel = [];
U_wb_L_PitchAccel = [];
U_wb_R_PitchAccel = [];
Dstroke_wb_PitchAccel = [];
Dpitch_wb_PitchAccel = [];
Ddev_wb_PitchAccel = [];
Daoa_wb_PitchAccel = [];
DU_wb_PitchAccel = [];

t_ds_L_PitchAccel = [];
t_ds_R_PitchAccel = [];
stroke_ds_L_PitchAccel = [];
stroke_ds_R_PitchAccel = [];
pitch_ds_L_PitchAccel = [];
pitch_ds_R_PitchAccel = [];
dev_ds_L_PitchAccel = [];
dev_ds_R_PitchAccel = [];
aoa_ds_L_PitchAccel = [];
aoa_ds_R_PitchAccel = [];
U_ds_L_PitchAccel = [];
U_ds_R_PitchAccel = [];
Dstroke_ds_PitchAccel = [];
Dpitch_ds_PitchAccel = [];
Ddev_ds_PitchAccel = [];
Daoa_ds_PitchAccel = [];
DU_ds_PitchAccel = [];

t_us_L_PitchAccel = [];
t_us_R_PitchAccel = [];
stroke_us_L_PitchAccel = [];
stroke_us_R_PitchAccel = [];
pitch_us_L_PitchAccel = [];
pitch_us_R_PitchAccel = [];
dev_us_L_PitchAccel = [];
dev_us_R_PitchAccel = [];
aoa_us_L_PitchAccel = [];
aoa_us_R_PitchAccel = [];
U_us_L_PitchAccel = [];
U_us_R_PitchAccel = [];
Dstroke_us_PitchAccel = [];
Dpitch_us_PitchAccel = [];
Ddev_us_PitchAccel = [];
Daoa_us_PitchAccel = [];
DU_us_PitchAccel = [];

t_wb_PitchAccel_bins = [];
stroke_wb_L_PitchAccel_bins = [];
stroke_wb_R_PitchAccel_bins = [];
pitch_wb_L_PitchAccel_bins = [];
pitch_wb_R_PitchAccel_bins = [];
dev_wb_L_PitchAccel_bins = [];
dev_wb_R_PitchAccel_bins = [];
aoa_wb_L_PitchAccel_bins = [];
aoa_wb_R_PitchAccel_bins = [];
U_wb_L_PitchAccel_bins = [];
U_wb_R_PitchAccel_bins = [];
Dstroke_wb_PitchAccel_bins = [];
Dpitch_wb_PitchAccel_bins = [];
Ddev_wb_PitchAccel_bins = [];
Daoa_wb_PitchAccel_bins = [];
DU_wb_PitchAccel_bins = [];

t_ds_PitchAccel_bins = [];
stroke_ds_L_PitchAccel_bins = [];
stroke_ds_R_PitchAccel_bins = [];
pitch_ds_L_PitchAccel_bins = [];
pitch_ds_R_PitchAccel_bins = [];
dev_ds_L_PitchAccel_bins = [];
dev_ds_R_PitchAccel_bins = [];
aoa_ds_L_PitchAccel_bins = [];
aoa_ds_R_PitchAccel_bins = [];
U_ds_L_PitchAccel_bins = [];
U_ds_R_PitchAccel_bins = [];
Dstroke_ds_PitchAccel_bins = [];
Dpitch_ds_PitchAccel_bins = [];
Ddev_ds_PitchAccel_bins = [];
Daoa_ds_PitchAccel_bins = [];
DU_ds_PitchAccel_bins = [];

t_us_PitchAccel_bins = [];
stroke_us_L_PitchAccel_bins = [];
stroke_us_R_PitchAccel_bins = [];
pitch_us_L_PitchAccel_bins = [];
pitch_us_R_PitchAccel_bins = [];
dev_us_L_PitchAccel_bins = [];
dev_us_R_PitchAccel_bins = [];
aoa_us_L_PitchAccel_bins = [];
aoa_us_R_PitchAccel_bins = [];
U_us_L_PitchAccel_bins = [];
U_us_R_PitchAccel_bins = [];
Dstroke_us_PitchAccel_bins = [];
Dpitch_us_PitchAccel_bins = [];
Ddev_us_PitchAccel_bins = [];
Daoa_us_PitchAccel_bins = [];
DU_us_PitchAccel_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if abs(pitch_dot_dot_mean_wb(wb)) > pitchaccel_limit_mod
%         if pitch_dot_dot_mean_wb(wb) > pitchaccel_limit_mod
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            F_mean_wb_now = F_mean_wb(wb);

            rollaccel_mean_wb_now = roll_dot_dot_mean_wb(wb);
            pitchaccel_mean_wb_now = pitch_dot_dot_mean_wb(wb);
            yawaccel_mean_wb_now = yaw_dot_dot_mean_wb(wb);

            rollvel_mean_wb_now = roll_dot_mean_wb(wb);
            pitchvel_mean_wb_now = pitch_dot_mean_wb(wb);
            yawvel_mean_wb_now = yaw_dot_mean_wb(wb);

            % wingbeat kin
            dt_ds_now = ([dt_ds_L(wb),dt_ds_R(wb)]);
            dt_us_now = ([dt_us_L(wb),dt_us_R(wb)]);
            f_wb_now = ([f_wb_L(wb),f_wb_R(wb)]);
            Rds_now = dt_ds_now ./ (dt_ds_now + dt_us_now);
            
            stroke_wb_L_now = stroke_wb_L(:,wb);
            stroke_wb_R_now = stroke_wb_R(:,wb);
            stroke_ds_L_now = stroke_ds_L(:,wb);
            stroke_ds_R_now = stroke_ds_R(:,wb);
            stroke_us_L_now = stroke_us_L(:,wb);
            stroke_us_R_now = stroke_us_R(:,wb);
            
            pitch_wb_L_now = pitch_wb_L(:,wb);
            pitch_wb_R_now = pitch_wb_R(:,wb);
            pitch_ds_L_now = pitch_ds_L(:,wb);
            pitch_ds_R_now = pitch_ds_R(:,wb);
            pitch_us_L_now = pitch_us_L(:,wb);
            pitch_us_R_now = pitch_us_R(:,wb);
            
            dev_wb_L_now = dev_wb_L(:,wb);
            dev_wb_R_now = dev_wb_R(:,wb);
            dev_ds_L_now = dev_ds_L(:,wb);
            dev_ds_R_now = dev_ds_R(:,wb);
            dev_us_L_now = dev_us_L(:,wb);
            dev_us_R_now = dev_us_R(:,wb);
            
            U_wb_L_now = U_wb_L(:,wb);
            U_wb_R_now = U_wb_R(:,wb);
            U_ds_L_now = U_ds_L(:,wb);
            U_ds_R_now = U_ds_R(:,wb);
            U_us_L_now = U_us_L(:,wb);
            U_us_R_now = U_us_R(:,wb);
            
            aoa_wb_L_now = aoa_wb_L(:,wb);
            aoa_wb_R_now = aoa_wb_R(:,wb);
            aoa_ds_L_now = aoa_ds_L(:,wb);
            aoa_ds_R_now = aoa_ds_R(:,wb);
            aoa_us_L_now = aoa_us_L(:,wb);
            aoa_us_R_now = aoa_us_R(:,wb);
            
            t_wb_L_now = t_wb_L(:,wb);
            t_wb_R_now = t_wb_R(:,wb);
            t_ds_L_now = t_ds_L(:,wb);
            t_ds_R_now = t_ds_R(:,wb);
            t_us_L_now = t_us_L(:,wb);
            t_us_R_now = t_us_R(:,wb);

            % interpolated data
            stroke_wb_L_interp = stroke_wb_L_bins(:,wb);
            stroke_wb_R_interp = stroke_wb_R_bins(:,wb);
            stroke_ds_L_interp = stroke_ds_L_bins(:,wb);
            stroke_ds_R_interp = stroke_ds_R_bins(:,wb);
            stroke_us_L_interp = stroke_us_L_bins(:,wb);
            stroke_us_R_interp = stroke_us_R_bins(:,wb);
            
            pitch_wb_L_interp = pitch_wb_L_bins(:,wb);
            pitch_wb_R_interp = pitch_wb_R_bins(:,wb);
            pitch_ds_L_interp = pitch_ds_L_bins(:,wb);
            pitch_ds_R_interp = pitch_ds_R_bins(:,wb);
            pitch_us_L_interp = pitch_us_L_bins(:,wb);
            pitch_us_R_interp = pitch_us_R_bins(:,wb);
            
            dev_wb_L_interp = dev_wb_L_bins(:,wb);
            dev_wb_R_interp = dev_wb_R_bins(:,wb);
            dev_ds_L_interp = dev_ds_L_bins(:,wb);
            dev_ds_R_interp = dev_ds_R_bins(:,wb);
            dev_us_L_interp = dev_us_L_bins(:,wb);
            dev_us_R_interp = dev_us_R_bins(:,wb);
            
            U_wb_L_interp = U_wb_L_bins(:,wb);
            U_wb_R_interp = U_wb_R_bins(:,wb);
            U_ds_L_interp = U_ds_L_bins(:,wb);
            U_ds_R_interp = U_ds_R_bins(:,wb);
            U_us_L_interp = U_us_L_bins(:,wb);
            U_us_R_interp = U_us_R_bins(:,wb);
            
            aoa_wb_L_interp = aoa_wb_L_bins(:,wb);
            aoa_wb_R_interp = aoa_wb_R_bins(:,wb);
            aoa_ds_L_interp = aoa_ds_L_bins(:,wb);
            aoa_ds_R_interp = aoa_ds_R_bins(:,wb);
            aoa_us_L_interp = aoa_us_L_bins(:,wb);
            aoa_us_R_interp = aoa_us_R_bins(:,wb);

            Dstroke_wb_interp = Dstroke_wb_bins(:,wb);
            Dstroke_ds_interp = Dstroke_ds_bins(:,wb);
            Dstroke_us_interp = Dstroke_us_bins(:,wb);
            
            Dpitch_wb_interp = Dpitch_wb_bins(:,wb);
            Dpitch_ds_interp = Dpitch_ds_bins(:,wb);
            Dpitch_us_interp = Dpitch_us_bins(:,wb);
            
            Ddev_wb_interp = Ddev_wb_bins(:,wb);
            Ddev_ds_interp = Ddev_ds_bins(:,wb);
            Ddev_us_interp = Ddev_us_bins(:,wb);
            
            DU_wb_interp = DU_wb_bins(:,wb);
            DU_ds_interp = DU_ds_bins(:,wb);
            DU_us_interp = DU_us_bins(:,wb);
            
            Daoa_wb_interp = Daoa_wb_bins(:,wb);
            Daoa_ds_interp = Daoa_ds_bins(:,wb);
            Daoa_us_interp = Daoa_us_bins(:,wb);
            
            % calc wingbeat modifications
            freqMOD_wb_now = (nanmean(f_wb_now) - f_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            
            strokeMOD_wb_L_now = (stroke_wb_L_interp - stroke_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            strokeMOD_wb_R_now = (stroke_wb_R_interp - stroke_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            strokeMOD_ds_L_now = (stroke_ds_L_interp - stroke_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            strokeMOD_ds_R_now = (stroke_ds_R_interp - stroke_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            strokeMOD_us_L_now = (stroke_us_L_interp - stroke_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            strokeMOD_us_R_now = (stroke_us_R_interp - stroke_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            
            pitchMOD_wb_L_now = (pitch_wb_L_interp - pitch_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            pitchMOD_wb_R_now = (pitch_wb_R_interp - pitch_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            pitchMOD_ds_L_now = (pitch_ds_L_interp - pitch_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            pitchMOD_ds_R_now = (pitch_ds_R_interp - pitch_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            pitchMOD_us_L_now = (pitch_us_L_interp - pitch_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            pitchMOD_us_R_now = (pitch_us_R_interp - pitch_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            
            devMOD_wb_L_now = (dev_wb_L_interp - dev_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            devMOD_wb_R_now = (dev_wb_R_interp - dev_wb_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            devMOD_ds_L_now = (dev_ds_L_interp - dev_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            devMOD_ds_R_now = (dev_ds_R_interp - dev_ds_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            devMOD_us_L_now = (dev_us_L_interp - dev_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            devMOD_us_R_now = (dev_us_R_interp - dev_us_steady) / pitchaccel_mean_wb_now * pitchaccel_norm;
            

            %% store data
            n_now=n_now+1;

            seq_nr_PitchAccel(n_now,1) = seq_nr_now;
            wb_nr_PitchAccel(n_now,1) = wb_nr_now;
            
            % body kin
            V_PitchAccel(n_now,1) = V_now;
            pitch_global_PitchAccel(n_now,1) = pitch_global_now;
            F_PitchAccel(n_now,1) = F_mean_wb_now;

            rollaccel_PitchAccel(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_PitchAccel(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_PitchAccel(n_now,1) = yawaccel_mean_wb_now;

            rollvel_PitchAccel(n_now,1) = rollvel_mean_wb_now;
            pitchvel_PitchAccel(n_now,1) = pitchvel_mean_wb_now;
            yawvel_PitchAccel(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_PitchAccel(n_now,1) = nanmean(dt_ds_now);
            dt_us_PitchAccel(n_now,1) = nanmean(dt_us_now);
            f_wb_PitchAccel(n_now,1) = nanmean(f_wb_now);
            Rds_PitchAccel(n_now,1) = nanmean(Rds_now);

            t_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = stroke_wb_L_now;
            stroke_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = stroke_wb_R_now;
            pitch_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = pitch_wb_L_now;
            pitch_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = pitch_wb_R_now;
            dev_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = dev_wb_L_now;
            dev_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = dev_wb_R_now;
            aoa_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = aoa_wb_L_now;
            aoa_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = aoa_wb_R_now;
            U_wb_L_PitchAccel(1:length(t_wb_L_now),n_now) = U_wb_L_now;
            U_wb_R_PitchAccel(1:length(t_wb_R_now),n_now) = U_wb_R_now;
            
            t_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = stroke_ds_L_now;
            stroke_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = stroke_ds_R_now;
            pitch_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = pitch_ds_L_now;
            pitch_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = pitch_ds_R_now;
            dev_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = dev_ds_L_now;
            dev_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = dev_ds_R_now;
            aoa_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = aoa_ds_L_now;
            aoa_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = aoa_ds_R_now;
            U_ds_L_PitchAccel(1:length(t_ds_L_now),n_now) = U_ds_L_now;
            U_ds_R_PitchAccel(1:length(t_ds_R_now),n_now) = U_ds_R_now;
            
            t_us_L_PitchAccel(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R_PitchAccel(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L_PitchAccel(1:length(t_us_L_now),n_now) = stroke_us_L_now;
            stroke_us_R_PitchAccel(1:length(t_us_R_now),n_now) = stroke_us_R_now;
            pitch_us_L_PitchAccel(1:length(t_us_L_now),n_now) = pitch_us_L_now;
            pitch_us_R_PitchAccel(1:length(t_us_R_now),n_now) = pitch_us_R_now;
            dev_us_L_PitchAccel(1:length(t_us_L_now),n_now) = dev_us_L_now;
            dev_us_R_PitchAccel(1:length(t_us_R_now),n_now) = dev_us_R_now;
            aoa_us_L_PitchAccel(1:length(t_us_L_now),n_now) = aoa_us_L_now;
            aoa_us_R_PitchAccel(1:length(t_us_R_now),n_now) = aoa_us_R_now;
            U_us_L_PitchAccel(1:length(t_us_L_now),n_now) = U_us_L_now;
            U_us_R_PitchAccel(1:length(t_us_R_now),n_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_PitchAccel_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_PitchAccel_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_PitchAccel_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_PitchAccel_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_PitchAccel_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_PitchAccel_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_PitchAccel_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_PitchAccel_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_PitchAccel_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_PitchAccel_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_PitchAccel_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_PitchAccel_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_PitchAccel_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_PitchAccel_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_PitchAccel_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_PitchAccel_bins(:,n_now) = DU_wb_interp;
            
            t_ds_PitchAccel_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_PitchAccel_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_PitchAccel_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_PitchAccel_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_PitchAccel_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_PitchAccel_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_PitchAccel_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_PitchAccel_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_PitchAccel_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_PitchAccel_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_PitchAccel_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_PitchAccel_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_PitchAccel_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_PitchAccel_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_PitchAccel_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_PitchAccel_bins(:,n_now) = DU_ds_interp;
            
            t_us_PitchAccel_bins(:,n_now) = t_us_bin;
            stroke_us_L_PitchAccel_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_PitchAccel_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_PitchAccel_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_PitchAccel_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_PitchAccel_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_PitchAccel_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_PitchAccel_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_PitchAccel_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_PitchAccel_bins(:,n_now) = U_us_L_interp;
            U_us_R_PitchAccel_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_PitchAccel_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_PitchAccel_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_PitchAccel_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_PitchAccel_bins(:,n_now) = Daoa_us_interp;
            DU_us_PitchAccel_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_PitchAccel(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_L_PitchAccel_bins(:,n_now) = strokeMOD_wb_L_now;
            strokeMOD_wb_R_PitchAccel_bins(:,n_now) = strokeMOD_wb_R_now;
            strokeMOD_ds_L_PitchAccel_bins(:,n_now) = strokeMOD_ds_L_now;
            strokeMOD_ds_R_PitchAccel_bins(:,n_now) = strokeMOD_ds_R_now;
            strokeMOD_us_L_PitchAccel_bins(:,n_now) = strokeMOD_us_L_now;
            strokeMOD_us_R_PitchAccel_bins(:,n_now) = strokeMOD_us_R_now;

            pitchMOD_wb_L_PitchAccel_bins(:,n_now) = pitchMOD_wb_L_now;
            pitchMOD_wb_R_PitchAccel_bins(:,n_now) = pitchMOD_wb_R_now;
            pitchMOD_ds_L_PitchAccel_bins(:,n_now) = pitchMOD_ds_L_now;
            pitchMOD_ds_R_PitchAccel_bins(:,n_now) = pitchMOD_ds_R_now;
            pitchMOD_us_L_PitchAccel_bins(:,n_now) = pitchMOD_us_L_now;
            pitchMOD_us_R_PitchAccel_bins(:,n_now) = pitchMOD_us_R_now;

            devMOD_wb_L_PitchAccel_bins(:,n_now) = devMOD_wb_L_now;
            devMOD_wb_R_PitchAccel_bins(:,n_now) = devMOD_wb_R_now;
            devMOD_ds_L_PitchAccel_bins(:,n_now) = devMOD_ds_L_now;
            devMOD_ds_R_PitchAccel_bins(:,n_now) = devMOD_ds_R_now;
            devMOD_us_L_PitchAccel_bins(:,n_now) = devMOD_us_L_now;
            devMOD_us_R_PitchAccel_bins(:,n_now) = devMOD_us_R_now;
        end
end

% mean & 95%CI
V_PitchAccel_meanCIstd = [nanmean(V_PitchAccel) 1.96*nanstd(V_PitchAccel)/sqrt(length(V_PitchAccel)) nanstd(V_PitchAccel)];
pitch_global_PitchAccel_meanCIstd = [nanmean(pitch_global_PitchAccel) 1.96*nanstd(pitch_global_PitchAccel)/sqrt(length(pitch_global_PitchAccel)) nanstd(pitch_global_PitchAccel)];

dt_ds_PitchAccel_meanCIstd = [nanmean(dt_ds_PitchAccel) 1.96*nanstd(dt_ds_PitchAccel)/sqrt(length(dt_ds_PitchAccel)) nanstd(dt_ds_PitchAccel)];
dt_us_PitchAccel_meanCIstd = [nanmean(dt_us_PitchAccel) 1.96*nanstd(dt_us_PitchAccel)/sqrt(length(dt_us_PitchAccel)) nanstd(dt_us_PitchAccel)];
f_wb_PitchAccel_meanCIstd = [nanmean(f_wb_PitchAccel) 1.96*nanstd(f_wb_PitchAccel)/sqrt(length(f_wb_PitchAccel)) nanstd(f_wb_PitchAccel)];
Rds_PitchAccel_meanCIstd = [nanmean(Rds_PitchAccel) 1.96*nanstd(Rds_PitchAccel)/sqrt(length(Rds_PitchAccel)) nanstd(Rds_PitchAccel)];

calc_WBfunc_PitchAccel_circmeanCIstd

% wbMOD means & 95%CI
freqMOD_wb_PitchAccel_meanCIstd = [nanmean(freqMOD_wb_PitchAccel) 1.96*nanstd(freqMOD_wb_PitchAccel)/sqrt(length(freqMOD_wb_PitchAccel)) nanstd(freqMOD_wb_PitchAccel)];
calc_WBmod_PitchAccel_circmeanCIstd

%% fit for average wb (bins)
    t_loc = t_wb_PitchAccel_bins(:,1);
    
    stroke_loc = stroke_wb_PitchAccel_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_PitchAccel_bins_meanCIstd(:,1);
    dev_loc = dev_wb_PitchAccel_bins_meanCIstd(:,1);
    Rds_loc = Rds_PitchAccel_meanCIstd(1);

    strokeMOD_loc = strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_PitchAccel_bins_meanCIstd(:,1);
    Rds_loc = Rds_steady_meanCIstd(1);


%% legendre polynomials
% wingbeats
    [stroke_PitchAccel_fit_binmean, stroke_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_PitchAccel_fit_binmean, pitch_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_PitchAccel_fit_binmean, dev_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeats
    [strokeMOD_PitchAccel_fit_binmean, strokeMOD_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_PitchAccel_fit_binmean, pitchMOD_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_PitchAccel_fit_binmean, devMOD_PitchAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

    % fourier series
        plotting = 0;
%         plotting = 1;

        % wingbeats
        [stroke_PitchAccel_fourier_fit_binmean, stroke_PitchAccel_fourier_gof_binmean, stroke_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_PitchAccel_fourier_fit_binmean, pitch_PitchAccel_fourier_gof_binmean, pitch_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_PitchAccel_fourier_fit_binmean, dev_PitchAccel_fourier_gof_binmean, dev_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_PitchAccel_fourier_fit_binmean, strokeMOD_PitchAccel_fourier_gof_binmean, strokeMOD_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_PitchAccel_fourier_fit_binmean, pitchMOD_PitchAccel_fourier_gof_binmean, pitchMOD_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_PitchAccel_fourier_fit_binmean, devMOD_PitchAccel_fourier_gof_binmean, devMOD_PitchAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);

    %% save PitchAccel wb data
n_PitchAccel = n_now   
    save(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs.mat'],...
'n_PitchAccel',...
'seq_nr_PitchAccel',...
'wb_nr_PitchAccel',...
...
'limit_mod',...
'rollaccel_limit_mod',...
'pitchaccel_limit_mod',...
'yawaccel_limit_mod',...
'Fenhance_limit_mod',...
...
'norm',...
'rollaccel_norm',...
'pitchaccel_norm',...
'yawaccel_norm',...
'Fenhance_norm',...
...
'n_pol_stroke',...
'n_pol_pitch',...
'n_pol_dev',...
'n_pol_MOD',...
...
'stroke_PitchAccel_fit_binmean_periodic',...
'pitch_PitchAccel_fit_binmean_periodic',...
'dev_PitchAccel_fit_binmean_periodic',...
...
'strokeMOD_PitchAccel_fit_binmean_periodic',...
'pitchMOD_PitchAccel_fit_binmean_periodic',...
'devMOD_PitchAccel_fit_binmean_periodic',...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_PitchAccel_fourier_coeffs_binmean',...
'pitch_PitchAccel_fourier_coeffs_binmean',...
'dev_PitchAccel_fourier_coeffs_binmean',...
...
'stroke_PitchAccel_fourier_gof_binmean',...
'pitch_PitchAccel_fourier_gof_binmean',...
'dev_PitchAccel_fourier_gof_binmean',...
...
'strokeMOD_PitchAccel_fourier_coeffs_binmean',...
'pitchMOD_PitchAccel_fourier_coeffs_binmean',...
'devMOD_PitchAccel_fourier_coeffs_binmean',...
...
'strokeMOD_PitchAccel_fourier_gof_binmean',...
'pitchMOD_PitchAccel_fourier_gof_binmean',...
'devMOD_PitchAccel_fourier_gof_binmean',...
...
'stroke_wb_PitchAccel_bins_meanCIstd',...
'stroke_ds_PitchAccel_bins_meanCIstd',...
'stroke_us_PitchAccel_bins_meanCIstd',...
...
'pitch_wb_PitchAccel_bins_meanCIstd',...
'pitch_ds_PitchAccel_bins_meanCIstd',...
'pitch_us_PitchAccel_bins_meanCIstd',...
...
'dev_wb_PitchAccel_bins_meanCIstd',...
'dev_ds_PitchAccel_bins_meanCIstd',...
'dev_us_PitchAccel_bins_meanCIstd',...
...
'strokeMOD_wb_PitchAccel_bins_meanCIstd',...
'strokeMOD_ds_PitchAccel_bins_meanCIstd',...
'strokeMOD_us_PitchAccel_bins_meanCIstd',...
...
'pitchMOD_wb_PitchAccel_bins_meanCIstd',...
'pitchMOD_ds_PitchAccel_bins_meanCIstd',...
'pitchMOD_us_PitchAccel_bins_meanCIstd',...
...
'devMOD_wb_PitchAccel_bins_meanCIstd',...
'devMOD_ds_PitchAccel_bins_meanCIstd',...
'devMOD_us_PitchAccel_bins_meanCIstd',...
...
'V_PitchAccel',...
'pitch_global_PitchAccel',...
'dt_ds_PitchAccel',...
'dt_us_PitchAccel',...
'f_wb_PitchAccel',...
'Rds_PitchAccel',...
...
'V_PitchAccel_meanCIstd',...
'pitch_global_PitchAccel_meanCIstd',...
'F_PitchAccel',...
...
'rollaccel_PitchAccel',...
'pitchaccel_PitchAccel',...
'yawaccel_PitchAccel',...
...
'rollvel_PitchAccel',...
'pitchvel_PitchAccel',...
'yawvel_PitchAccel',...
...
'dt_ds_PitchAccel_meanCIstd',...
'dt_us_PitchAccel_meanCIstd',...
'f_wb_PitchAccel_meanCIstd',...
'Rds_PitchAccel_meanCIstd',...
...
'freqMOD_wb_PitchAccel_meanCIstd',...
'freqMOD_wb_PitchAccel',...
...
'strokeMOD_wb_L_PitchAccel_bins',...
'strokeMOD_wb_R_PitchAccel_bins',...
'strokeMOD_ds_L_PitchAccel_bins',...
'strokeMOD_ds_R_PitchAccel_bins',...
'strokeMOD_us_L_PitchAccel_bins',...
'strokeMOD_us_R_PitchAccel_bins',...
...
'pitchMOD_wb_L_PitchAccel_bins',...
'pitchMOD_wb_R_PitchAccel_bins',...
'pitchMOD_ds_L_PitchAccel_bins',...
'pitchMOD_ds_R_PitchAccel_bins',...
'pitchMOD_us_L_PitchAccel_bins',...
'pitchMOD_us_R_PitchAccel_bins',...
...
'devMOD_wb_L_PitchAccel_bins',...
'devMOD_wb_R_PitchAccel_bins',...
'devMOD_ds_L_PitchAccel_bins',...
'devMOD_ds_R_PitchAccel_bins',...
'devMOD_us_L_PitchAccel_bins',...
'devMOD_us_R_PitchAccel_bins',...
...
't_wb_L_PitchAccel',...
't_wb_R_PitchAccel',...
'stroke_wb_L_PitchAccel',...
'stroke_wb_R_PitchAccel',...
'pitch_wb_L_PitchAccel',...
'pitch_wb_R_PitchAccel',...
'dev_wb_L_PitchAccel',...
'dev_wb_R_PitchAccel',...
'aoa_wb_L_PitchAccel',...
'aoa_wb_R_PitchAccel',...
'U_wb_L_PitchAccel',...
'U_wb_R_PitchAccel',...
...
't_ds_L_PitchAccel',...
't_ds_R_PitchAccel',...
'stroke_ds_L_PitchAccel',...
'stroke_ds_R_PitchAccel',...
'pitch_ds_L_PitchAccel',...
'pitch_ds_R_PitchAccel',...
'dev_ds_L_PitchAccel',...
'dev_ds_R_PitchAccel',...
'aoa_ds_L_PitchAccel',...
'aoa_ds_R_PitchAccel',...
'U_ds_L_PitchAccel',...
'U_ds_R_PitchAccel',...
...
't_us_L_PitchAccel',...
't_us_R_PitchAccel',...
'stroke_us_L_PitchAccel',...
'stroke_us_R_PitchAccel',...
'pitch_us_L_PitchAccel',...
'pitch_us_R_PitchAccel',...
'dev_us_L_PitchAccel',...
'dev_us_R_PitchAccel',...
'aoa_us_L_PitchAccel',...
'aoa_us_R_PitchAccel',...
'U_us_L_PitchAccel',...
'U_us_R_PitchAccel',...
...
't_wb_PitchAccel_bins',...
'stroke_wb_L_PitchAccel_bins',...
'stroke_wb_R_PitchAccel_bins',...
'pitch_wb_L_PitchAccel_bins',...
'pitch_wb_R_PitchAccel_bins',...
'dev_wb_L_PitchAccel_bins',...
'dev_wb_R_PitchAccel_bins',...
'aoa_wb_L_PitchAccel_bins',...
'aoa_wb_R_PitchAccel_bins',...
'U_wb_L_PitchAccel_bins',...
'U_wb_R_PitchAccel_bins',...
'Dstroke_wb_PitchAccel_bins',...
'Dpitch_wb_PitchAccel_bins',...
'Ddev_wb_PitchAccel_bins',...
'Daoa_wb_PitchAccel_bins',...
'DU_wb_PitchAccel_bins',...
...
't_ds_PitchAccel_bins',...
'stroke_ds_L_PitchAccel_bins',...
'stroke_ds_R_PitchAccel_bins',...
'pitch_ds_L_PitchAccel_bins',...
'pitch_ds_R_PitchAccel_bins',...
'dev_ds_L_PitchAccel_bins',...
'dev_ds_R_PitchAccel_bins',...
'aoa_ds_L_PitchAccel_bins',...
'aoa_ds_R_PitchAccel_bins',...
'U_ds_L_PitchAccel_bins',...
'U_ds_R_PitchAccel_bins',...
'Dstroke_ds_PitchAccel_bins',...
'Dpitch_ds_PitchAccel_bins',...
'Ddev_ds_PitchAccel_bins',...
'Daoa_ds_PitchAccel_bins',...
'DU_ds_PitchAccel_bins',...
...
't_us_PitchAccel_bins',...
'stroke_us_L_PitchAccel_bins',...
'stroke_us_R_PitchAccel_bins',...
'pitch_us_L_PitchAccel_bins',...
'pitch_us_R_PitchAccel_bins',...
'dev_us_L_PitchAccel_bins',...
'dev_us_R_PitchAccel_bins',...
'aoa_us_L_PitchAccel_bins',...
'aoa_us_R_PitchAccel_bins',...
'U_us_L_PitchAccel_bins',...
'U_us_R_PitchAccel_bins',...
'Dstroke_us_PitchAccel_bins',...
'Dpitch_us_PitchAccel_bins',...
'Ddev_us_PitchAccel_bins',...
'Daoa_us_PitchAccel_bins',...
'DU_us_PitchAccel_bins');



