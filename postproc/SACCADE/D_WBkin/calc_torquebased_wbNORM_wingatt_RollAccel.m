% calc wb MOD Roll Accel

n_now=0;

%% variables
% wbMOD
freqMOD_wb_RollAccel = [];

strokeMOD_wb_up_RollAccel_bins = [];
strokeMOD_wb_down_RollAccel_bins = [];
strokeMOD_ds_up_RollAccel_bins = [];
strokeMOD_ds_down_RollAccel_bins = [];
strokeMOD_us_up_RollAccel_bins = [];
strokeMOD_us_down_RollAccel_bins = [];

pitchMOD_wb_up_RollAccel_bins = [];
pitchMOD_wb_down_RollAccel_bins = [];
pitchMOD_ds_up_RollAccel_bins = [];
pitchMOD_ds_down_RollAccel_bins = [];
pitchMOD_us_up_RollAccel_bins = [];
pitchMOD_us_down_RollAccel_bins = [];

devMOD_wb_up_RollAccel_bins = [];
devMOD_wb_down_RollAccel_bins = [];
devMOD_ds_up_RollAccel_bins = [];
devMOD_ds_down_RollAccel_bins = [];
devMOD_us_up_RollAccel_bins = [];
devMOD_us_down_RollAccel_bins = [];

DstrokeMOD_wb_RollAccel_bins = [];
DstrokeMOD_ds_RollAccel_bins = [];
DstrokeMOD_us_RollAccel_bins = [];

DpitchMOD_wb_RollAccel_bins = [];
DpitchMOD_ds_RollAccel_bins = [];
DpitchMOD_us_RollAccel_bins = [];

DdevMOD_wb_RollAccel_bins = [];
DdevMOD_ds_RollAccel_bins = [];
DdevMOD_us_RollAccel_bins = [];


% seq&wb
seq_nr_RollAccel = [];
wb_nr_RollAccel = [];

% body kin
V_RollAccel = [];
pitch_global_RollAccel = [];
F_RollAccel = [];

Mroll_RollAccel = [];
Mpitch_RollAccel = [];
Myaw_RollAccel = [];

rollaccel_RollAccel = [];
pitchaccel_RollAccel = [];
yawaccel_RollAccel = [];

rollvel_RollAccel = [];
pitchvel_RollAccel = [];
yawvel_RollAccel = [];

rolldamp_RollAccel = [];
pitchdamp_RollAccel = [];
yawdamp_RollAccel = [];

% wingbeat kin
dt_wb_RollAccel = [];
dt_ds_RollAccel = [];
dt_us_RollAccel = [];
f_wb_RollAccel = [];
Rds_RollAccel = [];

t_wb_up_RollAccel = [];
t_wb_down_RollAccel = [];
stroke_wb_up_RollAccel = [];
stroke_wb_down_RollAccel = [];
pitch_wb_up_RollAccel = [];
pitch_wb_down_RollAccel = [];
dev_wb_up_RollAccel = [];
dev_wb_down_RollAccel = [];
aoa_wb_up_RollAccel = [];
aoa_wb_down_RollAccel = [];
U_wb_up_RollAccel = [];
U_wb_down_RollAccel = [];
Dstroke_wb_RollAccel = [];
Dpitch_wb_RollAccel = [];
Ddev_wb_RollAccel = [];
Daoa_wb_RollAccel = [];
DU_wb_RollAccel = [];

t_ds_up_RollAccel = [];
t_ds_down_RollAccel = [];
stroke_ds_up_RollAccel = [];
stroke_ds_down_RollAccel = [];
pitch_ds_up_RollAccel = [];
pitch_ds_down_RollAccel = [];
dev_ds_up_RollAccel = [];
dev_ds_down_RollAccel = [];
aoa_ds_up_RollAccel = [];
aoa_ds_down_RollAccel = [];
U_ds_up_RollAccel = [];
U_ds_down_RollAccel = [];
Dstroke_ds_RollAccel = [];
Dpitch_ds_RollAccel = [];
Ddev_ds_RollAccel = [];
Daoa_ds_RollAccel = [];
DU_ds_RollAccel = [];

t_us_up_RollAccel = [];
t_us_down_RollAccel = [];
stroke_us_up_RollAccel = [];
stroke_us_down_RollAccel = [];
pitch_us_up_RollAccel = [];
pitch_us_down_RollAccel = [];
dev_us_up_RollAccel = [];
dev_us_down_RollAccel = [];
aoa_us_up_RollAccel = [];
aoa_us_down_RollAccel = [];
U_us_up_RollAccel = [];
U_us_down_RollAccel = [];
Dstroke_us_RollAccel = [];
Dpitch_us_RollAccel = [];
Ddev_us_RollAccel = [];
Daoa_us_RollAccel = [];
DU_us_RollAccel = [];

t_wb_RollAccel_bins = [];
stroke_wb_up_RollAccel_bins = [];
stroke_wb_down_RollAccel_bins = [];
pitch_wb_up_RollAccel_bins = [];
pitch_wb_down_RollAccel_bins = [];
dev_wb_up_RollAccel_bins = [];
dev_wb_down_RollAccel_bins = [];
aoa_wb_up_RollAccel_bins = [];
aoa_wb_down_RollAccel_bins = [];
U_wb_up_RollAccel_bins = [];
U_wb_down_RollAccel_bins = [];
Dstroke_wb_RollAccel_bins = [];
Dpitch_wb_RollAccel_bins = [];
Ddev_wb_RollAccel_bins = [];
Daoa_wb_RollAccel_bins = [];
DU_wb_RollAccel_bins = [];

t_ds_RollAccel_bins = [];
stroke_ds_up_RollAccel_bins = [];
stroke_ds_down_RollAccel_bins = [];
pitch_ds_up_RollAccel_bins = [];
pitch_ds_down_RollAccel_bins = [];
dev_ds_up_RollAccel_bins = [];
dev_ds_down_RollAccel_bins = [];
aoa_ds_up_RollAccel_bins = [];
aoa_ds_down_RollAccel_bins = [];
U_ds_up_RollAccel_bins = [];
U_ds_down_RollAccel_bins = [];
Dstroke_ds_RollAccel_bins = [];
Dpitch_ds_RollAccel_bins = [];
Ddev_ds_RollAccel_bins = [];
Daoa_ds_RollAccel_bins = [];
DU_ds_RollAccel_bins = [];

t_us_RollAccel_bins = [];
stroke_us_up_RollAccel_bins = [];
stroke_us_down_RollAccel_bins = [];
pitch_us_up_RollAccel_bins = [];
pitch_us_down_RollAccel_bins = [];
dev_us_up_RollAccel_bins = [];
dev_us_down_RollAccel_bins = [];
aoa_us_up_RollAccel_bins = [];
aoa_us_down_RollAccel_bins = [];
U_us_up_RollAccel_bins = [];
U_us_down_RollAccel_bins = [];
Dstroke_us_RollAccel_bins = [];
Dpitch_us_RollAccel_bins = [];
Ddev_us_RollAccel_bins = [];
Daoa_us_RollAccel_bins = [];
DU_us_RollAccel_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if Mroll_mean_wb(wb) > Mroll_limit_mod
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            F_mean_wb_now = F_mean_wb(wb);

            Mroll_mean_wb_now = Mroll_mean_wb(wb);
            Mpitch_mean_wb_now = Mpitch_mean_wb(wb);
            Myaw_mean_wb_now = Myaw_mean_wb(wb);

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
            
            stroke_wb_up_now = stroke_wb_L(:,wb);
            stroke_wb_down_now = stroke_wb_R(:,wb);
            stroke_ds_up_now = stroke_ds_L(:,wb);
            stroke_ds_down_now = stroke_ds_R(:,wb);
            stroke_us_up_now = stroke_us_L(:,wb);
            stroke_us_down_now = stroke_us_R(:,wb);
            
            pitch_wb_up_now = pitch_wb_L(:,wb);
            pitch_wb_down_now = pitch_wb_R(:,wb);
            pitch_ds_up_now = pitch_ds_L(:,wb);
            pitch_ds_down_now = pitch_ds_R(:,wb);
            pitch_us_up_now = pitch_us_L(:,wb);
            pitch_us_down_now = pitch_us_R(:,wb);
            
            dev_wb_up_now = dev_wb_L(:,wb);
            dev_wb_down_now = dev_wb_R(:,wb);
            dev_ds_up_now = dev_ds_L(:,wb);
            dev_ds_down_now = dev_ds_R(:,wb);
            dev_us_up_now = dev_us_L(:,wb);
            dev_us_down_now = dev_us_R(:,wb);
            
            U_wb_up_now = U_wb_L(:,wb);
            U_wb_down_now = U_wb_R(:,wb);
            U_ds_up_now = U_ds_L(:,wb);
            U_ds_down_now = U_ds_R(:,wb);
            U_us_up_now = U_us_L(:,wb);
            U_us_down_now = U_us_R(:,wb);
            
            aoa_wb_up_now = aoa_wb_L(:,wb);
            aoa_wb_down_now = aoa_wb_R(:,wb);
            aoa_ds_up_now = aoa_ds_L(:,wb);
            aoa_ds_down_now = aoa_ds_R(:,wb);
            aoa_us_up_now = aoa_us_L(:,wb);
            aoa_us_down_now = aoa_us_R(:,wb);
            
            t_wb_up_now = t_wb_L(:,wb);
            t_wb_down_now = t_wb_R(:,wb);
            t_ds_up_now = t_ds_L(:,wb);
            t_ds_down_now = t_ds_R(:,wb);
            t_us_up_now = t_us_L(:,wb);
            t_us_down_now = t_us_R(:,wb);

            % interpolated data
            stroke_wb_up_interp = stroke_wb_L_bins(:,wb);
            stroke_wb_down_interp = stroke_wb_R_bins(:,wb);
            stroke_ds_up_interp = stroke_ds_L_bins(:,wb);
            stroke_ds_down_interp = stroke_ds_R_bins(:,wb);
            stroke_us_up_interp = stroke_us_L_bins(:,wb);
            stroke_us_down_interp = stroke_us_R_bins(:,wb);
            
            pitch_wb_up_interp = pitch_wb_L_bins(:,wb);
            pitch_wb_down_interp = pitch_wb_R_bins(:,wb);
            pitch_ds_up_interp = pitch_ds_L_bins(:,wb);
            pitch_ds_down_interp = pitch_ds_R_bins(:,wb);
            pitch_us_up_interp = pitch_us_L_bins(:,wb);
            pitch_us_down_interp = pitch_us_R_bins(:,wb);
            
            dev_wb_up_interp = dev_wb_L_bins(:,wb);
            dev_wb_down_interp = dev_wb_R_bins(:,wb);
            dev_ds_up_interp = dev_ds_L_bins(:,wb);
            dev_ds_down_interp = dev_ds_R_bins(:,wb);
            dev_us_up_interp = dev_us_L_bins(:,wb);
            dev_us_down_interp = dev_us_R_bins(:,wb);
            
            U_wb_up_interp = U_wb_L_bins(:,wb);
            U_wb_down_interp = U_wb_R_bins(:,wb);
            U_ds_up_interp = U_ds_L_bins(:,wb);
            U_ds_down_interp = U_ds_R_bins(:,wb);
            U_us_up_interp = U_us_L_bins(:,wb);
            U_us_down_interp = U_us_R_bins(:,wb);
            
            aoa_wb_up_interp = aoa_wb_L_bins(:,wb);
            aoa_wb_down_interp = aoa_wb_R_bins(:,wb);
            aoa_ds_up_interp = aoa_ds_L_bins(:,wb);
            aoa_ds_down_interp = aoa_ds_R_bins(:,wb);
            aoa_us_up_interp = aoa_us_L_bins(:,wb);
            aoa_us_down_interp = aoa_us_R_bins(:,wb);

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
            freqMOD_wb_now = (nanmean(f_wb_now) - f_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            strokeMOD_wb_up_now = (stroke_wb_up_interp - stroke_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_wb_down_now = (stroke_wb_down_interp - stroke_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_ds_up_now = (stroke_ds_up_interp - stroke_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_ds_down_now = (stroke_ds_down_interp - stroke_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_us_up_now = (stroke_us_up_interp - stroke_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_us_down_now = (stroke_us_down_interp - stroke_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            pitchMOD_wb_up_now = (pitch_wb_up_interp - pitch_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_wb_down_now = (pitch_wb_down_interp - pitch_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_ds_up_now = (pitch_ds_up_interp - pitch_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_ds_down_now = (pitch_ds_down_interp - pitch_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_us_up_now = (pitch_us_up_interp - pitch_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_us_down_now = (pitch_us_down_interp - pitch_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            devMOD_wb_up_now = (dev_wb_up_interp - dev_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_wb_down_now = (dev_wb_down_interp - dev_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_ds_up_now = (dev_ds_up_interp - dev_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_ds_down_now = (dev_ds_down_interp - dev_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_us_up_now = (dev_us_up_interp - dev_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_us_down_now = (dev_us_down_interp - dev_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            DstrokeMOD_wb_now = Dstroke_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DstrokeMOD_ds_now = Dstroke_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DstrokeMOD_us_now = Dstroke_us_interp / Mroll_mean_wb_now * Mroll_norm;

            DpitchMOD_wb_now = Dpitch_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DpitchMOD_ds_now = Dpitch_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DpitchMOD_us_now = Dpitch_us_interp / Mroll_mean_wb_now * Mroll_norm;

            DdevMOD_wb_now = Ddev_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DdevMOD_ds_now = Ddev_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DdevMOD_us_now = Ddev_us_interp / Mroll_mean_wb_now * Mroll_norm;

            %% store data
            n_now=n_now+1;

            seq_nr_RollAccel(n_now,1) = seq_nr_now;
            wb_nr_RollAccel(n_now,1) = wb_nr_now;
            
            % body kin
            V_RollAccel(n_now,1) = V_now;
            pitch_global_RollAccel(n_now,1) = pitch_global_now;
            F_RollAccel(n_now,1) = F_mean_wb_now;

            Mroll_RollAccel(n_now,1) = Mroll_mean_wb_now;
            Mpitch_RollAccel(n_now,1) = Mpitch_mean_wb_now;
            Myaw_RollAccel(n_now,1) = Myaw_mean_wb_now;

            rollaccel_RollAccel(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_RollAccel(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_RollAccel(n_now,1) = yawaccel_mean_wb_now;

            rollvel_RollAccel(n_now,1) = rollvel_mean_wb_now;
            pitchvel_RollAccel(n_now,1) = pitchvel_mean_wb_now;
            yawvel_RollAccel(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_RollAccel(n_now,1) = nanmean(dt_ds_now);
            dt_us_RollAccel(n_now,1) = nanmean(dt_us_now);
            f_wb_RollAccel(n_now,1) = nanmean(f_wb_now);
            Rds_RollAccel(n_now,1) = nanmean(Rds_now);

            t_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = t_wb_up_now;
            t_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = t_wb_down_now;
            stroke_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = stroke_wb_up_now;
            stroke_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = stroke_wb_down_now;
            pitch_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = pitch_wb_up_now;
            pitch_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = pitch_wb_down_now;
            dev_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = dev_wb_up_now;
            dev_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = dev_wb_down_now;
            aoa_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = aoa_wb_up_now;
            aoa_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = aoa_wb_down_now;
            U_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = U_wb_up_now;
            U_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = U_wb_down_now;
            
            t_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = t_ds_up_now;
            t_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = t_ds_down_now;
            stroke_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = stroke_ds_up_now;
            stroke_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = stroke_ds_down_now;
            pitch_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = pitch_ds_up_now;
            pitch_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = pitch_ds_down_now;
            dev_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = dev_ds_up_now;
            dev_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = dev_ds_down_now;
            aoa_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = aoa_ds_up_now;
            aoa_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = aoa_ds_down_now;
            U_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = U_ds_up_now;
            U_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = U_ds_down_now;
            
            t_us_up_RollAccel(1:length(t_us_up_now),n_now) = t_us_up_now;
            t_us_down_RollAccel(1:length(t_us_down_now),n_now) = t_us_down_now;
            stroke_us_up_RollAccel(1:length(t_us_up_now),n_now) = stroke_us_up_now;
            stroke_us_down_RollAccel(1:length(t_us_down_now),n_now) = stroke_us_down_now;
            pitch_us_up_RollAccel(1:length(t_us_up_now),n_now) = pitch_us_up_now;
            pitch_us_down_RollAccel(1:length(t_us_down_now),n_now) = pitch_us_down_now;
            dev_us_up_RollAccel(1:length(t_us_up_now),n_now) = dev_us_up_now;
            dev_us_down_RollAccel(1:length(t_us_down_now),n_now) = dev_us_down_now;
            aoa_us_up_RollAccel(1:length(t_us_up_now),n_now) = aoa_us_up_now;
            aoa_us_down_RollAccel(1:length(t_us_down_now),n_now) = aoa_us_down_now;
            U_us_up_RollAccel(1:length(t_us_up_now),n_now) = U_us_up_now;
            U_us_down_RollAccel(1:length(t_us_down_now),n_now) = U_us_down_now;
            
            % store interp binned data separate rows
            t_wb_RollAccel_bins(:,n_now) = t_wb_bin;
            stroke_wb_up_RollAccel_bins(:,n_now) = stroke_wb_up_interp;
            stroke_wb_down_RollAccel_bins(:,n_now) = stroke_wb_down_interp;
            pitch_wb_up_RollAccel_bins(:,n_now) = pitch_wb_up_interp;
            pitch_wb_down_RollAccel_bins(:,n_now) = pitch_wb_down_interp;
            dev_wb_up_RollAccel_bins(:,n_now) = dev_wb_up_interp;
            dev_wb_down_RollAccel_bins(:,n_now) = dev_wb_down_interp;
            aoa_wb_up_RollAccel_bins(:,n_now) = aoa_wb_up_interp;
            aoa_wb_down_RollAccel_bins(:,n_now) = aoa_wb_down_interp;
            U_wb_up_RollAccel_bins(:,n_now) = U_wb_up_interp;
            U_wb_down_RollAccel_bins(:,n_now) = U_wb_down_interp;
            Dstroke_wb_RollAccel_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_RollAccel_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_RollAccel_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_RollAccel_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_RollAccel_bins(:,n_now) = DU_wb_interp;
            
            t_ds_RollAccel_bins(:,n_now) = t_ds_bin;
            stroke_ds_up_RollAccel_bins(:,n_now) = stroke_ds_up_interp;
            stroke_ds_down_RollAccel_bins(:,n_now) = stroke_ds_down_interp;
            pitch_ds_up_RollAccel_bins(:,n_now) = pitch_ds_up_interp;
            pitch_ds_down_RollAccel_bins(:,n_now) = pitch_ds_down_interp;
            dev_ds_up_RollAccel_bins(:,n_now) = dev_ds_up_interp;
            dev_ds_down_RollAccel_bins(:,n_now) = dev_ds_down_interp;
            aoa_ds_up_RollAccel_bins(:,n_now) = aoa_ds_up_interp;
            aoa_ds_down_RollAccel_bins(:,n_now) = aoa_ds_down_interp;
            U_ds_up_RollAccel_bins(:,n_now) = U_ds_up_interp;
            U_ds_down_RollAccel_bins(:,n_now) = U_ds_down_interp;
            Dstroke_ds_RollAccel_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_RollAccel_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_RollAccel_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_RollAccel_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_RollAccel_bins(:,n_now) = DU_ds_interp;
            
            t_us_RollAccel_bins(:,n_now) = t_us_bin;
            stroke_us_up_RollAccel_bins(:,n_now) = stroke_us_up_interp;
            stroke_us_down_RollAccel_bins(:,n_now) = stroke_us_down_interp;
            pitch_us_up_RollAccel_bins(:,n_now) = pitch_us_up_interp;
            pitch_us_down_RollAccel_bins(:,n_now) = pitch_us_down_interp;
            dev_us_up_RollAccel_bins(:,n_now) = dev_us_up_interp;
            dev_us_down_RollAccel_bins(:,n_now) = dev_us_down_interp;
            aoa_us_up_RollAccel_bins(:,n_now) = aoa_us_up_interp;
            aoa_us_down_RollAccel_bins(:,n_now) = aoa_us_down_interp;
            U_us_up_RollAccel_bins(:,n_now) = U_us_up_interp;
            U_us_down_RollAccel_bins(:,n_now) = U_us_down_interp;
            Dstroke_us_RollAccel_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_RollAccel_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_RollAccel_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_RollAccel_bins(:,n_now) = Daoa_us_interp;
            DU_us_RollAccel_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_RollAccel(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_up_RollAccel_bins(:,n_now) = strokeMOD_wb_up_now;
            strokeMOD_wb_down_RollAccel_bins(:,n_now) = strokeMOD_wb_down_now;
            strokeMOD_ds_up_RollAccel_bins(:,n_now) = strokeMOD_ds_up_now;
            strokeMOD_ds_down_RollAccel_bins(:,n_now) = strokeMOD_ds_down_now;
            strokeMOD_us_up_RollAccel_bins(:,n_now) = strokeMOD_us_up_now;
            strokeMOD_us_down_RollAccel_bins(:,n_now) = strokeMOD_us_down_now;

            pitchMOD_wb_up_RollAccel_bins(:,n_now) = pitchMOD_wb_up_now;
            pitchMOD_wb_down_RollAccel_bins(:,n_now) = pitchMOD_wb_down_now;
            pitchMOD_ds_up_RollAccel_bins(:,n_now) = pitchMOD_ds_up_now;
            pitchMOD_ds_down_RollAccel_bins(:,n_now) = pitchMOD_ds_down_now;
            pitchMOD_us_up_RollAccel_bins(:,n_now) = pitchMOD_us_up_now;
            pitchMOD_us_down_RollAccel_bins(:,n_now) = pitchMOD_us_down_now;

            devMOD_wb_up_RollAccel_bins(:,n_now) = devMOD_wb_up_now;
            devMOD_wb_down_RollAccel_bins(:,n_now) = devMOD_wb_down_now;
            devMOD_ds_up_RollAccel_bins(:,n_now) = devMOD_ds_up_now;
            devMOD_ds_down_RollAccel_bins(:,n_now) = devMOD_ds_down_now;
            devMOD_us_up_RollAccel_bins(:,n_now) = devMOD_us_up_now;
            devMOD_us_down_RollAccel_bins(:,n_now) = devMOD_us_down_now;
            
            DstrokeMOD_wb_RollAccel_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_RollAccel_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_RollAccel_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_RollAccel_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_RollAccel_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_RollAccel_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_RollAccel_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_RollAccel_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_RollAccel_bins(:,n_now) = DdevMOD_us_now;
            
        elseif Mroll_mean_wb(wb) < -Mroll_limit_mod
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            F_mean_wb_now = F_mean_wb(wb);

            Mroll_mean_wb_now = Mroll_mean_wb(wb);
            Mpitch_mean_wb_now = Mpitch_mean_wb(wb);
            Myaw_mean_wb_now = Myaw_mean_wb(wb);

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
            
            stroke_wb_up_now = stroke_wb_R(:,wb);
            stroke_wb_down_now = stroke_wb_L(:,wb);
            stroke_ds_up_now = stroke_ds_R(:,wb);
            stroke_ds_down_now = stroke_ds_L(:,wb);
            stroke_us_up_now = stroke_us_R(:,wb);
            stroke_us_down_now = stroke_us_L(:,wb);
            
            pitch_wb_up_now = pitch_wb_R(:,wb);
            pitch_wb_down_now = pitch_wb_L(:,wb);
            pitch_ds_up_now = pitch_ds_R(:,wb);
            pitch_ds_down_now = pitch_ds_L(:,wb);
            pitch_us_up_now = pitch_us_R(:,wb);
            pitch_us_down_now = pitch_us_L(:,wb);
            
            dev_wb_up_now = dev_wb_R(:,wb);
            dev_wb_down_now = dev_wb_L(:,wb);
            dev_ds_up_now = dev_ds_R(:,wb);
            dev_ds_down_now = dev_ds_L(:,wb);
            dev_us_up_now = dev_us_R(:,wb);
            dev_us_down_now = dev_us_L(:,wb);
            
            U_wb_up_now = U_wb_R(:,wb);
            U_wb_down_now = U_wb_L(:,wb);
            U_ds_up_now = U_ds_R(:,wb);
            U_ds_down_now = U_ds_L(:,wb);
            U_us_up_now = U_us_R(:,wb);
            U_us_down_now = U_us_L(:,wb);
            
            aoa_wb_up_now = aoa_wb_R(:,wb);
            aoa_wb_down_now = aoa_wb_L(:,wb);
            aoa_ds_up_now = aoa_ds_R(:,wb);
            aoa_ds_down_now = aoa_ds_L(:,wb);
            aoa_us_up_now = aoa_us_R(:,wb);
            aoa_us_down_now = aoa_us_L(:,wb);
            
            t_wb_up_now = t_wb_R(:,wb);
            t_wb_down_now = t_wb_L(:,wb);
            t_ds_up_now = t_ds_R(:,wb);
            t_ds_down_now = t_ds_L(:,wb);
            t_us_up_now = t_us_R(:,wb);
            t_us_down_now = t_us_L(:,wb);

            % interpolated data
            stroke_wb_up_interp = stroke_wb_R_bins(:,wb);
            stroke_wb_down_interp = stroke_wb_L_bins(:,wb);
            stroke_ds_up_interp = stroke_ds_R_bins(:,wb);
            stroke_ds_down_interp = stroke_ds_L_bins(:,wb);
            stroke_us_up_interp = stroke_us_R_bins(:,wb);
            stroke_us_down_interp = stroke_us_L_bins(:,wb);
            
            pitch_wb_up_interp = pitch_wb_R_bins(:,wb);
            pitch_wb_down_interp = pitch_wb_L_bins(:,wb);
            pitch_ds_up_interp = pitch_ds_R_bins(:,wb);
            pitch_ds_down_interp = pitch_ds_L_bins(:,wb);
            pitch_us_up_interp = pitch_us_R_bins(:,wb);
            pitch_us_down_interp = pitch_us_L_bins(:,wb);
            
            dev_wb_up_interp = dev_wb_R_bins(:,wb);
            dev_wb_down_interp = dev_wb_L_bins(:,wb);
            dev_ds_up_interp = dev_ds_R_bins(:,wb);
            dev_ds_down_interp = dev_ds_L_bins(:,wb);
            dev_us_up_interp = dev_us_R_bins(:,wb);
            dev_us_down_interp = dev_us_L_bins(:,wb);
            
            U_wb_up_interp = U_wb_R_bins(:,wb);
            U_wb_down_interp = U_wb_L_bins(:,wb);
            U_ds_up_interp = U_ds_R_bins(:,wb);
            U_ds_down_interp = U_ds_L_bins(:,wb);
            U_us_up_interp = U_us_R_bins(:,wb);
            U_us_down_interp = U_us_L_bins(:,wb);
            
            aoa_wb_up_interp = aoa_wb_R_bins(:,wb);
            aoa_wb_down_interp = aoa_wb_L_bins(:,wb);
            aoa_ds_up_interp = aoa_ds_R_bins(:,wb);
            aoa_ds_down_interp = aoa_ds_L_bins(:,wb);
            aoa_us_up_interp = aoa_us_R_bins(:,wb);
            aoa_us_down_interp = aoa_us_L_bins(:,wb);

            Dstroke_wb_interp = -Dstroke_wb_bins(:,wb);
            Dstroke_ds_interp = -Dstroke_ds_bins(:,wb);
            Dstroke_us_interp = -Dstroke_us_bins(:,wb);
            
            Dpitch_wb_interp = -Dpitch_wb_bins(:,wb);
            Dpitch_ds_interp = -Dpitch_ds_bins(:,wb);
            Dpitch_us_interp = -Dpitch_us_bins(:,wb);
            
            Ddev_wb_interp = -Ddev_wb_bins(:,wb);
            Ddev_ds_interp = -Ddev_ds_bins(:,wb);
            Ddev_us_interp = -Ddev_us_bins(:,wb);
            
            DU_wb_interp = -DU_wb_bins(:,wb);
            DU_ds_interp = -DU_ds_bins(:,wb);
            DU_us_interp = -DU_us_bins(:,wb);
            
            Daoa_wb_interp = -Daoa_wb_bins(:,wb);
            Daoa_ds_interp = -Daoa_ds_bins(:,wb);
            Daoa_us_interp = -Daoa_us_bins(:,wb);
            
            % calc wingbeat modifications
            freqMOD_wb_now = -(nanmean(f_wb_now) - f_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            strokeMOD_wb_up_now = -(stroke_wb_up_interp - stroke_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_wb_down_now = -(stroke_wb_down_interp - stroke_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_ds_up_now = -(stroke_ds_up_interp - stroke_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_ds_down_now = -(stroke_ds_down_interp - stroke_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_us_up_now = -(stroke_us_up_interp - stroke_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            strokeMOD_us_down_now = -(stroke_us_down_interp - stroke_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            pitchMOD_wb_up_now = -(pitch_wb_up_interp - pitch_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_wb_down_now = -(pitch_wb_down_interp - pitch_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_ds_up_now = -(pitch_ds_up_interp - pitch_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_ds_down_now = -(pitch_ds_down_interp - pitch_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_us_up_now = -(pitch_us_up_interp - pitch_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            pitchMOD_us_down_now = -(pitch_us_down_interp - pitch_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            devMOD_wb_up_now = -(dev_wb_up_interp - dev_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_wb_down_now = -(dev_wb_down_interp - dev_wb_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_ds_up_now = -(dev_ds_up_interp - dev_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_ds_down_now = -(dev_ds_down_interp - dev_ds_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_us_up_now = -(dev_us_up_interp - dev_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            devMOD_us_down_now = -(dev_us_down_interp - dev_us_steady) / Mroll_mean_wb_now * Mroll_norm;
            
            DstrokeMOD_wb_now = -Dstroke_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DstrokeMOD_ds_now = -Dstroke_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DstrokeMOD_us_now = -Dstroke_us_interp / Mroll_mean_wb_now * Mroll_norm;

            DpitchMOD_wb_now = -Dpitch_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DpitchMOD_ds_now = -Dpitch_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DpitchMOD_us_now = -Dpitch_us_interp / Mroll_mean_wb_now * Mroll_norm;

            DdevMOD_wb_now = -Ddev_wb_interp / Mroll_mean_wb_now * Mroll_norm;
            DdevMOD_ds_now = -Ddev_ds_interp / Mroll_mean_wb_now * Mroll_norm;
            DdevMOD_us_now = -Ddev_us_interp / Mroll_mean_wb_now * Mroll_norm;

            %% store data
            n_now=n_now+1;

            seq_nr_RollAccel(n_now,1) = seq_nr_now;
            wb_nr_RollAccel(n_now,1) = wb_nr_now;
            
            % body kin
            V_RollAccel(n_now,1) = V_now;
            pitch_global_RollAccel(n_now,1) = pitch_global_now;
            F_RollAccel(n_now,1) = F_mean_wb_now;

            Mroll_RollAccel(n_now,1) = Mroll_mean_wb_now;
            Mpitch_RollAccel(n_now,1) = Mpitch_mean_wb_now;
            Myaw_RollAccel(n_now,1) = Myaw_mean_wb_now;

            rollaccel_RollAccel(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_RollAccel(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_RollAccel(n_now,1) = yawaccel_mean_wb_now;

            rollvel_RollAccel(n_now,1) = rollvel_mean_wb_now;
            pitchvel_RollAccel(n_now,1) = pitchvel_mean_wb_now;
            yawvel_RollAccel(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_RollAccel(n_now,1) = nanmean(dt_ds_now);
            dt_us_RollAccel(n_now,1) = nanmean(dt_us_now);
            f_wb_RollAccel(n_now,1) = nanmean(f_wb_now);
            Rds_RollAccel(n_now,1) = nanmean(Rds_now);

            t_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = t_wb_up_now;
            t_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = t_wb_down_now;
            stroke_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = stroke_wb_up_now;
            stroke_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = stroke_wb_down_now;
            pitch_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = pitch_wb_up_now;
            pitch_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = pitch_wb_down_now;
            dev_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = dev_wb_up_now;
            dev_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = dev_wb_down_now;
            aoa_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = aoa_wb_up_now;
            aoa_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = aoa_wb_down_now;
            U_wb_up_RollAccel(1:length(t_wb_up_now),n_now) = U_wb_up_now;
            U_wb_down_RollAccel(1:length(t_wb_down_now),n_now) = U_wb_down_now;
            
            t_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = t_ds_up_now;
            t_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = t_ds_down_now;
            stroke_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = stroke_ds_up_now;
            stroke_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = stroke_ds_down_now;
            pitch_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = pitch_ds_up_now;
            pitch_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = pitch_ds_down_now;
            dev_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = dev_ds_up_now;
            dev_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = dev_ds_down_now;
            aoa_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = aoa_ds_up_now;
            aoa_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = aoa_ds_down_now;
            U_ds_up_RollAccel(1:length(t_ds_up_now),n_now) = U_ds_up_now;
            U_ds_down_RollAccel(1:length(t_ds_down_now),n_now) = U_ds_down_now;
            
            t_us_up_RollAccel(1:length(t_us_up_now),n_now) = t_us_up_now;
            t_us_down_RollAccel(1:length(t_us_down_now),n_now) = t_us_down_now;
            stroke_us_up_RollAccel(1:length(t_us_up_now),n_now) = stroke_us_up_now;
            stroke_us_down_RollAccel(1:length(t_us_down_now),n_now) = stroke_us_down_now;
            pitch_us_up_RollAccel(1:length(t_us_up_now),n_now) = pitch_us_up_now;
            pitch_us_down_RollAccel(1:length(t_us_down_now),n_now) = pitch_us_down_now;
            dev_us_up_RollAccel(1:length(t_us_up_now),n_now) = dev_us_up_now;
            dev_us_down_RollAccel(1:length(t_us_down_now),n_now) = dev_us_down_now;
            aoa_us_up_RollAccel(1:length(t_us_up_now),n_now) = aoa_us_up_now;
            aoa_us_down_RollAccel(1:length(t_us_down_now),n_now) = aoa_us_down_now;
            U_us_up_RollAccel(1:length(t_us_up_now),n_now) = U_us_up_now;
            U_us_down_RollAccel(1:length(t_us_down_now),n_now) = U_us_down_now;
            
            % store interp binned data separate rows
            t_wb_RollAccel_bins(:,n_now) = t_wb_bin;
            stroke_wb_up_RollAccel_bins(:,n_now) = stroke_wb_up_interp;
            stroke_wb_down_RollAccel_bins(:,n_now) = stroke_wb_down_interp;
            pitch_wb_up_RollAccel_bins(:,n_now) = pitch_wb_up_interp;
            pitch_wb_down_RollAccel_bins(:,n_now) = pitch_wb_down_interp;
            dev_wb_up_RollAccel_bins(:,n_now) = dev_wb_up_interp;
            dev_wb_down_RollAccel_bins(:,n_now) = dev_wb_down_interp;
            aoa_wb_up_RollAccel_bins(:,n_now) = aoa_wb_up_interp;
            aoa_wb_down_RollAccel_bins(:,n_now) = aoa_wb_down_interp;
            U_wb_up_RollAccel_bins(:,n_now) = U_wb_up_interp;
            U_wb_down_RollAccel_bins(:,n_now) = U_wb_down_interp;
            Dstroke_wb_RollAccel_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_RollAccel_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_RollAccel_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_RollAccel_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_RollAccel_bins(:,n_now) = DU_wb_interp;
            
            t_ds_RollAccel_bins(:,n_now) = t_ds_bin;
            stroke_ds_up_RollAccel_bins(:,n_now) = stroke_ds_up_interp;
            stroke_ds_down_RollAccel_bins(:,n_now) = stroke_ds_down_interp;
            pitch_ds_up_RollAccel_bins(:,n_now) = pitch_ds_up_interp;
            pitch_ds_down_RollAccel_bins(:,n_now) = pitch_ds_down_interp;
            dev_ds_up_RollAccel_bins(:,n_now) = dev_ds_up_interp;
            dev_ds_down_RollAccel_bins(:,n_now) = dev_ds_down_interp;
            aoa_ds_up_RollAccel_bins(:,n_now) = aoa_ds_up_interp;
            aoa_ds_down_RollAccel_bins(:,n_now) = aoa_ds_down_interp;
            U_ds_up_RollAccel_bins(:,n_now) = U_ds_up_interp;
            U_ds_down_RollAccel_bins(:,n_now) = U_ds_down_interp;
            Dstroke_ds_RollAccel_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_RollAccel_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_RollAccel_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_RollAccel_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_RollAccel_bins(:,n_now) = DU_ds_interp;
            
            t_us_RollAccel_bins(:,n_now) = t_us_bin;
            stroke_us_up_RollAccel_bins(:,n_now) = stroke_us_up_interp;
            stroke_us_down_RollAccel_bins(:,n_now) = stroke_us_down_interp;
            pitch_us_up_RollAccel_bins(:,n_now) = pitch_us_up_interp;
            pitch_us_down_RollAccel_bins(:,n_now) = pitch_us_down_interp;
            dev_us_up_RollAccel_bins(:,n_now) = dev_us_up_interp;
            dev_us_down_RollAccel_bins(:,n_now) = dev_us_down_interp;
            aoa_us_up_RollAccel_bins(:,n_now) = aoa_us_up_interp;
            aoa_us_down_RollAccel_bins(:,n_now) = aoa_us_down_interp;
            U_us_up_RollAccel_bins(:,n_now) = U_us_up_interp;
            U_us_down_RollAccel_bins(:,n_now) = U_us_down_interp;
            Dstroke_us_RollAccel_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_RollAccel_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_RollAccel_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_RollAccel_bins(:,n_now) = Daoa_us_interp;
            DU_us_RollAccel_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_RollAccel(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_up_RollAccel_bins(:,n_now) = strokeMOD_wb_up_now;
            strokeMOD_wb_down_RollAccel_bins(:,n_now) = strokeMOD_wb_down_now;
            strokeMOD_ds_up_RollAccel_bins(:,n_now) = strokeMOD_ds_up_now;
            strokeMOD_ds_down_RollAccel_bins(:,n_now) = strokeMOD_ds_down_now;
            strokeMOD_us_up_RollAccel_bins(:,n_now) = strokeMOD_us_up_now;
            strokeMOD_us_down_RollAccel_bins(:,n_now) = strokeMOD_us_down_now;

            pitchMOD_wb_up_RollAccel_bins(:,n_now) = pitchMOD_wb_up_now;
            pitchMOD_wb_down_RollAccel_bins(:,n_now) = pitchMOD_wb_down_now;
            pitchMOD_ds_up_RollAccel_bins(:,n_now) = pitchMOD_ds_up_now;
            pitchMOD_ds_down_RollAccel_bins(:,n_now) = pitchMOD_ds_down_now;
            pitchMOD_us_up_RollAccel_bins(:,n_now) = pitchMOD_us_up_now;
            pitchMOD_us_down_RollAccel_bins(:,n_now) = pitchMOD_us_down_now;

            devMOD_wb_up_RollAccel_bins(:,n_now) = devMOD_wb_up_now;
            devMOD_wb_down_RollAccel_bins(:,n_now) = devMOD_wb_down_now;
            devMOD_ds_up_RollAccel_bins(:,n_now) = devMOD_ds_up_now;
            devMOD_ds_down_RollAccel_bins(:,n_now) = devMOD_ds_down_now;
            devMOD_us_up_RollAccel_bins(:,n_now) = devMOD_us_up_now;
            devMOD_us_down_RollAccel_bins(:,n_now) = devMOD_us_down_now;
            
            DstrokeMOD_wb_RollAccel_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_RollAccel_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_RollAccel_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_RollAccel_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_RollAccel_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_RollAccel_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_RollAccel_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_RollAccel_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_RollAccel_bins(:,n_now) = DdevMOD_us_now;
        end
end

% mean & 95%CI
V_RollAccel_meanCIstd = [nanmean(V_RollAccel) 1.96*nanstd(V_RollAccel)/sqrt(length(V_RollAccel)) nanstd(V_RollAccel)];
pitch_global_RollAccel_meanCIstd = [nanmean(pitch_global_RollAccel) 1.96*nanstd(pitch_global_RollAccel)/sqrt(length(pitch_global_RollAccel)) nanstd(pitch_global_RollAccel)];

dt_ds_RollAccel_meanCIstd = [nanmean(dt_ds_RollAccel) 1.96*nanstd(dt_ds_RollAccel)/sqrt(length(dt_ds_RollAccel)) nanstd(dt_ds_RollAccel)];
dt_us_RollAccel_meanCIstd = [nanmean(dt_us_RollAccel) 1.96*nanstd(dt_us_RollAccel)/sqrt(length(dt_us_RollAccel)) nanstd(dt_us_RollAccel)];
f_wb_RollAccel_meanCIstd = [nanmean(f_wb_RollAccel) 1.96*nanstd(f_wb_RollAccel)/sqrt(length(f_wb_RollAccel)) nanstd(f_wb_RollAccel)];
Rds_RollAccel_meanCIstd = [nanmean(Rds_RollAccel) 1.96*nanstd(Rds_RollAccel)/sqrt(length(Rds_RollAccel)) nanstd(Rds_RollAccel)];

calc_WBfunc_RollAccel_circmeanCIstd

% wbMOD means & 95%CI
freqMOD_wb_RollAccel_meanCIstd = [nanmean(freqMOD_wb_RollAccel) 1.96*nanstd(freqMOD_wb_RollAccel)/sqrt(length(freqMOD_wb_RollAccel)) nanstd(freqMOD_wb_RollAccel)];
calc_WBmod_RollAccel_circmeanCIstd

%% WBfits
    t_loc = t_wb_RollAccel_bins(:,1);
    Rds_loc = Rds_RollAccel_meanCIstd(1);

    plotting = 0;
%         plotting = 1;

    %% fit for upwards moving wing (bins)
    stroke_loc = stroke_wb_up_RollAccel_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_up_RollAccel_bins_meanCIstd(:,1);
    dev_loc = dev_wb_up_RollAccel_bins_meanCIstd(:,1);
    
    strokeMOD_loc = strokeMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_up_RollAccel_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_up_RollAccel_fit_binmean, stroke_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_up_RollAccel_fit_binmean, pitch_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_up_RollAccel_fit_binmean, dev_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_up_RollAccel_fit_binmean, strokeMOD_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_up_RollAccel_fit_binmean, pitchMOD_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_up_RollAccel_fit_binmean, devMOD_up_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_up_RollAccel_fourier_fit_binmean, stroke_up_RollAccel_fourier_gof_binmean,stroke_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_up_RollAccel_fourier_fit_binmean, pitch_up_RollAccel_fourier_gof_binmean,pitch_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_up_RollAccel_fourier_fit_binmean, dev_up_RollAccel_fourier_gof_binmean,dev_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_up_RollAccel_fourier_fit_binmean, strokeMOD_up_RollAccel_fourier_gof_binmean,strokeMOD_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_up_RollAccel_fourier_fit_binmean, pitchMOD_up_RollAccel_fourier_gof_binmean,pitchMOD_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_up_RollAccel_fourier_fit_binmean, devMOD_up_RollAccel_fourier_gof_binmean,devMOD_up_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% downwards
    stroke_loc = stroke_wb_down_RollAccel_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_down_RollAccel_bins_meanCIstd(:,1);
    dev_loc = dev_wb_down_RollAccel_bins_meanCIstd(:,1);

    strokeMOD_loc = strokeMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_down_RollAccel_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_down_RollAccel_fit_binmean, stroke_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_down_RollAccel_fit_binmean, pitch_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_down_RollAccel_fit_binmean, dev_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_down_RollAccel_fit_binmean, strokeMOD_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_down_RollAccel_fit_binmean, pitchMOD_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_down_RollAccel_fit_binmean, devMOD_down_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_down_RollAccel_fourier_fit_binmean, stroke_down_RollAccel_fourier_gof_binmean,stroke_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_down_RollAccel_fourier_fit_binmean, pitch_down_RollAccel_fourier_gof_binmean,pitch_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_down_RollAccel_fourier_fit_binmean, dev_down_RollAccel_fourier_gof_binmean,dev_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_down_RollAccel_fourier_fit_binmean, strokeMOD_down_RollAccel_fourier_gof_binmean,strokeMOD_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_down_RollAccel_fourier_fit_binmean, pitchMOD_down_RollAccel_fourier_gof_binmean,pitchMOD_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_down_RollAccel_fourier_fit_binmean, devMOD_down_RollAccel_fourier_gof_binmean,devMOD_down_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% upwards - downwards
    stroke_loc = Dstroke_wb_RollAccel_bins_meanCIstd(:,1);
    pitch_loc = Dpitch_wb_RollAccel_bins_meanCIstd(:,1);
    dev_loc = Ddev_wb_RollAccel_bins_meanCIstd(:,1);

    strokeMOD_loc = DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1);
    pitchMOD_loc = DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1);
    devMOD_loc = DdevMOD_wb_RollAccel_bins_meanCIstd(:,1);
    
%% legendre polynomials
    % wingbeats
    [Dstroke_RollAccel_fit_binmean, Dstroke_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [Dpitch_RollAccel_fit_binmean, Dpitch_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [Ddev_RollAccel_fit_binmean, Ddev_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [DstrokeMOD_RollAccel_fit_binmean, DstrokeMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [DpitchMOD_RollAccel_fit_binmean, DpitchMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [DdevMOD_RollAccel_fit_binmean, DdevMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [Dstroke_RollAccel_fourier_fit_binmean, Dstroke_RollAccel_fourier_gof_binmean,Dstroke_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [Dpitch_RollAccel_fourier_fit_binmean, Dpitch_RollAccel_fourier_gof_binmean,Dpitch_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [Ddev_RollAccel_fourier_fit_binmean, Ddev_RollAccel_fourier_gof_binmean,Ddev_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [DstrokeMOD_RollAccel_fourier_fit_binmean, DstrokeMOD_RollAccel_fourier_gof_binmean,DstrokeMOD_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [DpitchMOD_RollAccel_fourier_fit_binmean, DpitchMOD_RollAccel_fourier_gof_binmean,DpitchMOD_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [DdevMOD_RollAccel_fourier_fit_binmean, DdevMOD_RollAccel_fourier_gof_binmean,DdevMOD_RollAccel_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    

    
    
    
    
    %% save RollAccel wb data
n_RollAccel = n_now   

if save_on == 1

save(['WBmod_torquebased_RollAccel_',num2str(n_RollAccel),'WBs.mat'],...
'n_RollAccel',...
'seq_nr_RollAccel',...
'wb_nr_RollAccel',...
...
'limit_mod',...
'Mroll_limit_mod',...
'Mpitch_limit_mod',...
'Myaw_limit_mod',...
'Fenhance_limit_mod',...
...
'norm',...
'Mroll_norm',...
'Mpitch_norm',...
'Myaw_norm',...
'Fenhance_norm',...
...
'n_pol_stroke',...
'n_pol_pitch',...
'n_pol_dev',...
'n_pol_MOD',...
...
'stroke_up_RollAccel_fit_binmean_periodic',...
'pitch_up_RollAccel_fit_binmean_periodic',...
'dev_up_RollAccel_fit_binmean_periodic',...
...
'stroke_down_RollAccel_fit_binmean_periodic',...
'pitch_down_RollAccel_fit_binmean_periodic',...
'dev_down_RollAccel_fit_binmean_periodic',...
...
'strokeMOD_up_RollAccel_fit_binmean_periodic',...
'pitchMOD_up_RollAccel_fit_binmean_periodic',...
'devMOD_up_RollAccel_fit_binmean_periodic',...
...
'strokeMOD_down_RollAccel_fit_binmean_periodic',...
'pitchMOD_down_RollAccel_fit_binmean_periodic',...
'devMOD_down_RollAccel_fit_binmean_periodic',...
...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_up_RollAccel_fourier_coeffs_binmean',...
'pitch_up_RollAccel_fourier_coeffs_binmean',...
'dev_up_RollAccel_fourier_coeffs_binmean',...
...
'stroke_up_RollAccel_fourier_gof_binmean',...
'pitch_up_RollAccel_fourier_gof_binmean',...
'dev_up_RollAccel_fourier_gof_binmean',...
...
'stroke_down_RollAccel_fourier_coeffs_binmean',...
'pitch_down_RollAccel_fourier_coeffs_binmean',...
'dev_down_RollAccel_fourier_coeffs_binmean',...
...
'stroke_down_RollAccel_fourier_gof_binmean',...
'pitch_down_RollAccel_fourier_gof_binmean',...
'dev_down_RollAccel_fourier_gof_binmean',...
...
'Dstroke_RollAccel_fourier_coeffs_binmean',...
'Dpitch_RollAccel_fourier_coeffs_binmean',...
'Ddev_RollAccel_fourier_coeffs_binmean',...
...
'Dstroke_RollAccel_fourier_gof_binmean',...
'Dpitch_RollAccel_fourier_gof_binmean',...
'Ddev_RollAccel_fourier_gof_binmean',...
...
'strokeMOD_up_RollAccel_fourier_coeffs_binmean',...
'pitchMOD_up_RollAccel_fourier_coeffs_binmean',...
'devMOD_up_RollAccel_fourier_coeffs_binmean',...
...
'strokeMOD_up_RollAccel_fourier_gof_binmean',...
'pitchMOD_up_RollAccel_fourier_gof_binmean',...
'devMOD_up_RollAccel_fourier_gof_binmean',...
...
'strokeMOD_down_RollAccel_fourier_coeffs_binmean',...
'pitchMOD_down_RollAccel_fourier_coeffs_binmean',...
'devMOD_down_RollAccel_fourier_coeffs_binmean',...
...
'strokeMOD_down_RollAccel_fourier_gof_binmean',...
'pitchMOD_down_RollAccel_fourier_gof_binmean',...
'devMOD_down_RollAccel_fourier_gof_binmean',...
...
'DstrokeMOD_RollAccel_fourier_coeffs_binmean',...
'DpitchMOD_RollAccel_fourier_coeffs_binmean',...
'DdevMOD_RollAccel_fourier_coeffs_binmean',...
...
'DstrokeMOD_RollAccel_fourier_gof_binmean',...
'DpitchMOD_RollAccel_fourier_gof_binmean',...
'DdevMOD_RollAccel_fourier_gof_binmean',...
...
'stroke_wb_up_RollAccel_bins_meanCIstd',...
'stroke_ds_up_RollAccel_bins_meanCIstd',...
'stroke_us_up_RollAccel_bins_meanCIstd',...
...
'pitch_wb_up_RollAccel_bins_meanCIstd',...
'pitch_ds_up_RollAccel_bins_meanCIstd',...
'pitch_us_up_RollAccel_bins_meanCIstd',...
...
'dev_wb_up_RollAccel_bins_meanCIstd',...
'dev_ds_up_RollAccel_bins_meanCIstd',...
'dev_us_up_RollAccel_bins_meanCIstd',...
...
'strokeMOD_wb_up_RollAccel_bins_meanCIstd',...
'strokeMOD_ds_up_RollAccel_bins_meanCIstd',...
'strokeMOD_us_up_RollAccel_bins_meanCIstd',...
...
'pitchMOD_wb_up_RollAccel_bins_meanCIstd',...
'pitchMOD_ds_up_RollAccel_bins_meanCIstd',...
'pitchMOD_us_up_RollAccel_bins_meanCIstd',...
...
'devMOD_wb_up_RollAccel_bins_meanCIstd',...
'devMOD_ds_up_RollAccel_bins_meanCIstd',...
'devMOD_us_up_RollAccel_bins_meanCIstd',...
...
'stroke_wb_down_RollAccel_bins_meanCIstd',...
'stroke_ds_down_RollAccel_bins_meanCIstd',...
'stroke_us_down_RollAccel_bins_meanCIstd',...
...
'pitch_wb_down_RollAccel_bins_meanCIstd',...
'pitch_ds_down_RollAccel_bins_meanCIstd',...
'pitch_us_down_RollAccel_bins_meanCIstd',...
...
'dev_wb_down_RollAccel_bins_meanCIstd',...
'dev_ds_down_RollAccel_bins_meanCIstd',...
'dev_us_down_RollAccel_bins_meanCIstd',...
...
'strokeMOD_wb_down_RollAccel_bins_meanCIstd',...
'strokeMOD_ds_down_RollAccel_bins_meanCIstd',...
'strokeMOD_us_down_RollAccel_bins_meanCIstd',...
...
'pitchMOD_wb_down_RollAccel_bins_meanCIstd',...
'pitchMOD_ds_down_RollAccel_bins_meanCIstd',...
'pitchMOD_us_down_RollAccel_bins_meanCIstd',...
...
'devMOD_wb_down_RollAccel_bins_meanCIstd',...
'devMOD_ds_down_RollAccel_bins_meanCIstd',...
'devMOD_us_down_RollAccel_bins_meanCIstd',...
...
'DstrokeMOD_wb_RollAccel_bins_meanCIstd',...
'DstrokeMOD_ds_RollAccel_bins_meanCIstd',...
'DstrokeMOD_us_RollAccel_bins_meanCIstd',...
...
'DpitchMOD_wb_RollAccel_bins_meanCIstd',...
'DpitchMOD_ds_RollAccel_bins_meanCIstd',...
'DpitchMOD_us_RollAccel_bins_meanCIstd',...
...
'DdevMOD_wb_RollAccel_bins_meanCIstd',...
'DdevMOD_ds_RollAccel_bins_meanCIstd',...
'DdevMOD_us_RollAccel_bins_meanCIstd',...
...
'V_RollAccel',...
'pitch_global_RollAccel',...
'dt_ds_RollAccel',...
'dt_us_RollAccel',...
'f_wb_RollAccel',...
'Rds_RollAccel',...
...
'V_RollAccel_meanCIstd',...
'pitch_global_RollAccel_meanCIstd',...
'F_RollAccel',...
...
'Mroll_RollAccel',...
'Mpitch_RollAccel',...
'Myaw_RollAccel',...
...
'rollaccel_RollAccel',...
'pitchaccel_RollAccel',...
'yawaccel_RollAccel',...
...
'rollvel_RollAccel',...
'pitchvel_RollAccel',...
'yawvel_RollAccel',...
...
'dt_ds_RollAccel_meanCIstd',...
'dt_us_RollAccel_meanCIstd',...
'f_wb_RollAccel_meanCIstd',...
'Rds_RollAccel_meanCIstd',...
...
'freqMOD_wb_RollAccel_meanCIstd',...
'freqMOD_wb_RollAccel',...
...
'strokeMOD_wb_up_RollAccel_bins',...
'strokeMOD_wb_down_RollAccel_bins',...
'strokeMOD_ds_up_RollAccel_bins',...
'strokeMOD_ds_down_RollAccel_bins',...
'strokeMOD_us_up_RollAccel_bins',...
'strokeMOD_us_down_RollAccel_bins',...
...
'pitchMOD_wb_up_RollAccel_bins',...
'pitchMOD_wb_down_RollAccel_bins',...
'pitchMOD_ds_up_RollAccel_bins',...
'pitchMOD_ds_down_RollAccel_bins',...
'pitchMOD_us_up_RollAccel_bins',...
'pitchMOD_us_down_RollAccel_bins',...
...
'devMOD_wb_up_RollAccel_bins',...
'devMOD_wb_down_RollAccel_bins',...
'devMOD_ds_up_RollAccel_bins',...
'devMOD_ds_down_RollAccel_bins',...
'devMOD_us_up_RollAccel_bins',...
'devMOD_us_down_RollAccel_bins',...
...
'DstrokeMOD_wb_RollAccel_bins',...
'DstrokeMOD_ds_RollAccel_bins',...
'DstrokeMOD_us_RollAccel_bins',...
...
'DpitchMOD_wb_RollAccel_bins',...
'DpitchMOD_ds_RollAccel_bins',...
'DpitchMOD_us_RollAccel_bins',...
...
'DdevMOD_wb_RollAccel_bins',...
'DdevMOD_ds_RollAccel_bins',...
'DdevMOD_us_RollAccel_bins',...
...
't_wb_up_RollAccel',...
't_wb_down_RollAccel',...
'stroke_wb_up_RollAccel',...
'stroke_wb_down_RollAccel',...
'pitch_wb_up_RollAccel',...
'pitch_wb_down_RollAccel',...
'dev_wb_up_RollAccel',...
'dev_wb_down_RollAccel',...
'aoa_wb_up_RollAccel',...
'aoa_wb_down_RollAccel',...
'U_wb_up_RollAccel',...
'U_wb_down_RollAccel',...
...
't_ds_up_RollAccel',...
't_ds_down_RollAccel',...
'stroke_ds_up_RollAccel',...
'stroke_ds_down_RollAccel',...
'pitch_ds_up_RollAccel',...
'pitch_ds_down_RollAccel',...
'dev_ds_up_RollAccel',...
'dev_ds_down_RollAccel',...
'aoa_ds_up_RollAccel',...
'aoa_ds_down_RollAccel',...
'U_ds_up_RollAccel',...
'U_ds_down_RollAccel',...
...
't_us_up_RollAccel',...
't_us_down_RollAccel',...
'stroke_us_up_RollAccel',...
'stroke_us_down_RollAccel',...
'pitch_us_up_RollAccel',...
'pitch_us_down_RollAccel',...
'dev_us_up_RollAccel',...
'dev_us_down_RollAccel',...
'aoa_us_up_RollAccel',...
'aoa_us_down_RollAccel',...
'U_us_up_RollAccel',...
'U_us_down_RollAccel',...
...
't_wb_RollAccel_bins',...
'stroke_wb_up_RollAccel_bins',...
'stroke_wb_down_RollAccel_bins',...
'pitch_wb_up_RollAccel_bins',...
'pitch_wb_down_RollAccel_bins',...
'dev_wb_up_RollAccel_bins',...
'dev_wb_down_RollAccel_bins',...
'aoa_wb_up_RollAccel_bins',...
'aoa_wb_down_RollAccel_bins',...
'U_wb_up_RollAccel_bins',...
'U_wb_down_RollAccel_bins',...
'Dstroke_wb_RollAccel_bins',...
'Dpitch_wb_RollAccel_bins',...
'Ddev_wb_RollAccel_bins',...
'Daoa_wb_RollAccel_bins',...
'DU_wb_RollAccel_bins',...
...
't_ds_RollAccel_bins',...
'stroke_ds_up_RollAccel_bins',...
'stroke_ds_down_RollAccel_bins',...
'pitch_ds_up_RollAccel_bins',...
'pitch_ds_down_RollAccel_bins',...
'dev_ds_up_RollAccel_bins',...
'dev_ds_down_RollAccel_bins',...
'aoa_ds_up_RollAccel_bins',...
'aoa_ds_down_RollAccel_bins',...
'U_ds_up_RollAccel_bins',...
'U_ds_down_RollAccel_bins',...
'Dstroke_ds_RollAccel_bins',...
'Dpitch_ds_RollAccel_bins',...
'Ddev_ds_RollAccel_bins',...
'Daoa_ds_RollAccel_bins',...
'DU_ds_RollAccel_bins',...
...
't_us_RollAccel_bins',...
'stroke_us_up_RollAccel_bins',...
'stroke_us_down_RollAccel_bins',...
'pitch_us_up_RollAccel_bins',...
'pitch_us_down_RollAccel_bins',...
'dev_us_up_RollAccel_bins',...
'dev_us_down_RollAccel_bins',...
'aoa_us_up_RollAccel_bins',...
'aoa_us_down_RollAccel_bins',...
'U_us_up_RollAccel_bins',...
'U_us_down_RollAccel_bins',...
'Dstroke_us_RollAccel_bins',...
'Dpitch_us_RollAccel_bins',...
'Ddev_us_RollAccel_bins',...
'Daoa_us_RollAccel_bins',...
'DU_us_RollAccel_bins');

end
