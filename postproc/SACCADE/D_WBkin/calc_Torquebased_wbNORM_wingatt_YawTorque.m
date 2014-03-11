% calc wb MOD Yaw Torque

n_now=0;

%% variables
% wbMOD
freqMOD_wb_YawTorque = [];

strokeMOD_wb_fwd_YawTorque_bins = [];
strokeMOD_wb_rwd_YawTorque_bins = [];
strokeMOD_ds_fwd_YawTorque_bins = [];
strokeMOD_ds_rwd_YawTorque_bins = [];
strokeMOD_us_fwd_YawTorque_bins = [];
strokeMOD_us_rwd_YawTorque_bins = [];

pitchMOD_wb_fwd_YawTorque_bins = [];
pitchMOD_wb_rwd_YawTorque_bins = [];
pitchMOD_ds_fwd_YawTorque_bins = [];
pitchMOD_ds_rwd_YawTorque_bins = [];
pitchMOD_us_fwd_YawTorque_bins = [];
pitchMOD_us_rwd_YawTorque_bins = [];

devMOD_wb_fwd_YawTorque_bins = [];
devMOD_wb_rwd_YawTorque_bins = [];
devMOD_ds_fwd_YawTorque_bins = [];
devMOD_ds_rwd_YawTorque_bins = [];
devMOD_us_fwd_YawTorque_bins = [];
devMOD_us_rwd_YawTorque_bins = [];

DstrokeMOD_wb_YawTorque_bins = [];
DstrokeMOD_ds_YawTorque_bins = [];
DstrokeMOD_us_YawTorque_bins = [];

DpitchMOD_wb_YawTorque_bins = [];
DpitchMOD_ds_YawTorque_bins = [];
DpitchMOD_us_YawTorque_bins = [];

DdevMOD_wb_YawTorque_bins = [];
DdevMOD_ds_YawTorque_bins = [];
DdevMOD_us_YawTorque_bins = [];


% seq&wb
seq_nr_YawTorque = [];
wb_nr_YawTorque = [];

% body kin
V_YawTorque = [];
pitch_global_YawTorque = [];
F_YawTorque = [];

Mroll_accel_YawTorque = [];
Mpitch_accel_YawTorque = [];
Myaw_accel_YawTorque = [];
M_L_accel_YawTorque = [];
M_R_accel_YawTorque = [];

Mroll_damp_YawTorque = [];
Mpitch_damp_YawTorque = [];
Myaw_damp_YawTorque = [];
M_L_damp_YawTorque = [];
M_R_damp_YawTorque = [];

Mroll_YawTorque = [];
Mpitch_YawTorque = [];
Myaw_YawTorque = [];
M_L_YawTorque = [];
M_R_YawTorque = [];

rollaccel_YawTorque = [];
pitchaccel_YawTorque = [];
yawaccel_YawTorque = [];

rollvel_YawTorque = [];
pitchvel_YawTorque = [];
yawvel_YawTorque = [];

% wingbeat kin
dt_wb_YawTorque = [];
dt_ds_YawTorque = [];
dt_us_YawTorque = [];
f_wb_YawTorque = [];
Rds_YawTorque = [];

t_wb_fwd_YawTorque = [];
t_wb_rwd_YawTorque = [];
stroke_wb_fwd_YawTorque = [];
stroke_wb_rwd_YawTorque = [];
pitch_wb_fwd_YawTorque = [];
pitch_wb_rwd_YawTorque = [];
dev_wb_fwd_YawTorque = [];
dev_wb_rwd_YawTorque = [];
aoa_wb_fwd_YawTorque = [];
aoa_wb_rwd_YawTorque = [];
U_wb_fwd_YawTorque = [];
U_wb_rwd_YawTorque = [];
Dstroke_wb_YawTorque = [];
Dpitch_wb_YawTorque = [];
Ddev_wb_YawTorque = [];
Daoa_wb_YawTorque = [];
DU_wb_YawTorque = [];

t_ds_fwd_YawTorque = [];
t_ds_rwd_YawTorque = [];
stroke_ds_fwd_YawTorque = [];
stroke_ds_rwd_YawTorque = [];
pitch_ds_fwd_YawTorque = [];
pitch_ds_rwd_YawTorque = [];
dev_ds_fwd_YawTorque = [];
dev_ds_rwd_YawTorque = [];
aoa_ds_fwd_YawTorque = [];
aoa_ds_rwd_YawTorque = [];
U_ds_fwd_YawTorque = [];
U_ds_rwd_YawTorque = [];
Dstroke_ds_YawTorque = [];
Dpitch_ds_YawTorque = [];
Ddev_ds_YawTorque = [];
Daoa_ds_YawTorque = [];
DU_ds_YawTorque = [];

t_us_fwd_YawTorque = [];
t_us_rwd_YawTorque = [];
stroke_us_fwd_YawTorque = [];
stroke_us_rwd_YawTorque = [];
pitch_us_fwd_YawTorque = [];
pitch_us_rwd_YawTorque = [];
dev_us_fwd_YawTorque = [];
dev_us_rwd_YawTorque = [];
aoa_us_fwd_YawTorque = [];
aoa_us_rwd_YawTorque = [];
U_us_fwd_YawTorque = [];
U_us_rwd_YawTorque = [];
Dstroke_us_YawTorque = [];
Dpitch_us_YawTorque = [];
Ddev_us_YawTorque = [];
Daoa_us_YawTorque = [];
DU_us_YawTorque = [];

t_wb_YawTorque_bins = [];
stroke_wb_fwd_YawTorque_bins = [];
stroke_wb_rwd_YawTorque_bins = [];
pitch_wb_fwd_YawTorque_bins = [];
pitch_wb_rwd_YawTorque_bins = [];
dev_wb_fwd_YawTorque_bins = [];
dev_wb_rwd_YawTorque_bins = [];
aoa_wb_fwd_YawTorque_bins = [];
aoa_wb_rwd_YawTorque_bins = [];
U_wb_fwd_YawTorque_bins = [];
U_wb_rwd_YawTorque_bins = [];
Dstroke_wb_YawTorque_bins = [];
Dpitch_wb_YawTorque_bins = [];
Ddev_wb_YawTorque_bins = [];
Daoa_wb_YawTorque_bins = [];
DU_wb_YawTorque_bins = [];

t_ds_YawTorque_bins = [];
stroke_ds_fwd_YawTorque_bins = [];
stroke_ds_rwd_YawTorque_bins = [];
pitch_ds_fwd_YawTorque_bins = [];
pitch_ds_rwd_YawTorque_bins = [];
dev_ds_fwd_YawTorque_bins = [];
dev_ds_rwd_YawTorque_bins = [];
aoa_ds_fwd_YawTorque_bins = [];
aoa_ds_rwd_YawTorque_bins = [];
U_ds_fwd_YawTorque_bins = [];
U_ds_rwd_YawTorque_bins = [];
Dstroke_ds_YawTorque_bins = [];
Dpitch_ds_YawTorque_bins = [];
Ddev_ds_YawTorque_bins = [];
Daoa_ds_YawTorque_bins = [];
DU_ds_YawTorque_bins = [];

t_us_YawTorque_bins = [];
stroke_us_fwd_YawTorque_bins = [];
stroke_us_rwd_YawTorque_bins = [];
pitch_us_fwd_YawTorque_bins = [];
pitch_us_rwd_YawTorque_bins = [];
dev_us_fwd_YawTorque_bins = [];
dev_us_rwd_YawTorque_bins = [];
aoa_us_fwd_YawTorque_bins = [];
aoa_us_rwd_YawTorque_bins = [];
U_us_fwd_YawTorque_bins = [];
U_us_rwd_YawTorque_bins = [];
Dstroke_us_YawTorque_bins = [];
Dpitch_us_YawTorque_bins = [];
Ddev_us_YawTorque_bins = [];
Daoa_us_YawTorque_bins = [];
DU_us_YawTorque_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if Myaw_mean_wb(wb) > Myaw_limit_mod
%         if yaw_dot_dot_mean_wb(wb) > yawaccel_limit_mod
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            F_mean_wb_now = F_mean_wb(wb);

            Mroll_mean_wb_accel_now = Mroll_mean_wb_accel(wb);
            Mpitch_mean_wb_accel_now = Mpitch_mean_wb_accel(wb);
            Myaw_mean_wb_accel_now = Myaw_mean_wb_accel(wb);
            M_L_mean_wb_accel_now = M_L_mean_wb_accel(wb);
            M_R_mean_wb_accel_now = M_R_mean_wb_accel(wb);

            Mroll_mean_wb_damp_now = Mroll_mean_wb_damp(wb);
            Mpitch_mean_wb_damp_now = Mpitch_mean_wb_damp(wb);
            Myaw_mean_wb_damp_now = Myaw_mean_wb_damp(wb);
            M_L_mean_wb_damp_now = M_L_mean_wb_damp(wb);
            M_R_mean_wb_damp_now = M_R_mean_wb_damp(wb);

            Mroll_mean_wb_now = Mroll_mean_wb(wb);
            Mpitch_mean_wb_now = Mpitch_mean_wb(wb);
            Myaw_mean_wb_now = Myaw_mean_wb(wb);
            M_L_mean_wb_now = M_L_mean_wb(wb);
            M_R_mean_wb_now = M_R_mean_wb(wb);

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
            
            stroke_wb_fwd_now = stroke_wb_L(:,wb);
            stroke_wb_rwd_now = stroke_wb_R(:,wb);
            stroke_ds_fwd_now = stroke_ds_L(:,wb);
            stroke_ds_rwd_now = stroke_ds_R(:,wb);
            stroke_us_fwd_now = stroke_us_L(:,wb);
            stroke_us_rwd_now = stroke_us_R(:,wb);
            
            pitch_wb_fwd_now = pitch_wb_L(:,wb);
            pitch_wb_rwd_now = pitch_wb_R(:,wb);
            pitch_ds_fwd_now = pitch_ds_L(:,wb);
            pitch_ds_rwd_now = pitch_ds_R(:,wb);
            pitch_us_fwd_now = pitch_us_L(:,wb);
            pitch_us_rwd_now = pitch_us_R(:,wb);
            
            dev_wb_fwd_now = dev_wb_L(:,wb);
            dev_wb_rwd_now = dev_wb_R(:,wb);
            dev_ds_fwd_now = dev_ds_L(:,wb);
            dev_ds_rwd_now = dev_ds_R(:,wb);
            dev_us_fwd_now = dev_us_L(:,wb);
            dev_us_rwd_now = dev_us_R(:,wb);
            
            U_wb_fwd_now = U_wb_L(:,wb);
            U_wb_rwd_now = U_wb_R(:,wb);
            U_ds_fwd_now = U_ds_L(:,wb);
            U_ds_rwd_now = U_ds_R(:,wb);
            U_us_fwd_now = U_us_L(:,wb);
            U_us_rwd_now = U_us_R(:,wb);
            
            aoa_wb_fwd_now = aoa_wb_L(:,wb);
            aoa_wb_rwd_now = aoa_wb_R(:,wb);
            aoa_ds_fwd_now = aoa_ds_L(:,wb);
            aoa_ds_rwd_now = aoa_ds_R(:,wb);
            aoa_us_fwd_now = aoa_us_L(:,wb);
            aoa_us_rwd_now = aoa_us_R(:,wb);
            
            t_wb_fwd_now = t_wb_L(:,wb);
            t_wb_rwd_now = t_wb_R(:,wb);
            t_ds_fwd_now = t_ds_L(:,wb);
            t_ds_rwd_now = t_ds_R(:,wb);
            t_us_fwd_now = t_us_L(:,wb);
            t_us_rwd_now = t_us_R(:,wb);

            % interpolated data
            stroke_wb_fwd_interp = stroke_wb_L_bins(:,wb);
            stroke_wb_rwd_interp = stroke_wb_R_bins(:,wb);
            stroke_ds_fwd_interp = stroke_ds_L_bins(:,wb);
            stroke_ds_rwd_interp = stroke_ds_R_bins(:,wb);
            stroke_us_fwd_interp = stroke_us_L_bins(:,wb);
            stroke_us_rwd_interp = stroke_us_R_bins(:,wb);
            
            pitch_wb_fwd_interp = pitch_wb_L_bins(:,wb);
            pitch_wb_rwd_interp = pitch_wb_R_bins(:,wb);
            pitch_ds_fwd_interp = pitch_ds_L_bins(:,wb);
            pitch_ds_rwd_interp = pitch_ds_R_bins(:,wb);
            pitch_us_fwd_interp = pitch_us_L_bins(:,wb);
            pitch_us_rwd_interp = pitch_us_R_bins(:,wb);
            
            dev_wb_fwd_interp = dev_wb_L_bins(:,wb);
            dev_wb_rwd_interp = dev_wb_R_bins(:,wb);
            dev_ds_fwd_interp = dev_ds_L_bins(:,wb);
            dev_ds_rwd_interp = dev_ds_R_bins(:,wb);
            dev_us_fwd_interp = dev_us_L_bins(:,wb);
            dev_us_rwd_interp = dev_us_R_bins(:,wb);
            
            U_wb_fwd_interp = U_wb_L_bins(:,wb);
            U_wb_rwd_interp = U_wb_R_bins(:,wb);
            U_ds_fwd_interp = U_ds_L_bins(:,wb);
            U_ds_rwd_interp = U_ds_R_bins(:,wb);
            U_us_fwd_interp = U_us_L_bins(:,wb);
            U_us_rwd_interp = U_us_R_bins(:,wb);
            
            aoa_wb_fwd_interp = aoa_wb_L_bins(:,wb);
            aoa_wb_rwd_interp = aoa_wb_R_bins(:,wb);
            aoa_ds_fwd_interp = aoa_ds_L_bins(:,wb);
            aoa_ds_rwd_interp = aoa_ds_R_bins(:,wb);
            aoa_us_fwd_interp = aoa_us_L_bins(:,wb);
            aoa_us_rwd_interp = aoa_us_R_bins(:,wb);

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
            freqMOD_wb_now = (nanmean(f_wb_now) - f_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            strokeMOD_wb_fwd_now = (stroke_wb_fwd_interp - stroke_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_wb_rwd_now = (stroke_wb_rwd_interp - stroke_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_ds_fwd_now = (stroke_ds_fwd_interp - stroke_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_ds_rwd_now = (stroke_ds_rwd_interp - stroke_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_us_fwd_now = (stroke_us_fwd_interp - stroke_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_us_rwd_now = (stroke_us_rwd_interp - stroke_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            pitchMOD_wb_fwd_now = (pitch_wb_fwd_interp - pitch_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_wb_rwd_now = (pitch_wb_rwd_interp - pitch_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_ds_fwd_now = (pitch_ds_fwd_interp - pitch_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_ds_rwd_now = (pitch_ds_rwd_interp - pitch_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_us_fwd_now = (pitch_us_fwd_interp - pitch_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_us_rwd_now = (pitch_us_rwd_interp - pitch_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            devMOD_wb_fwd_now = (dev_wb_fwd_interp - dev_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_wb_rwd_now = (dev_wb_rwd_interp - dev_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_ds_fwd_now = (dev_ds_fwd_interp - dev_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_ds_rwd_now = (dev_ds_rwd_interp - dev_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_us_fwd_now = (dev_us_fwd_interp - dev_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_us_rwd_now = (dev_us_rwd_interp - dev_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            DstrokeMOD_wb_now = Dstroke_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DstrokeMOD_ds_now = Dstroke_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DstrokeMOD_us_now = Dstroke_us_interp / Myaw_mean_wb_now * Myaw_norm;

            DpitchMOD_wb_now = Dpitch_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DpitchMOD_ds_now = Dpitch_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DpitchMOD_us_now = Dpitch_us_interp / Myaw_mean_wb_now * Myaw_norm;

            DdevMOD_wb_now = Ddev_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DdevMOD_ds_now = Ddev_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DdevMOD_us_now = Ddev_us_interp / Myaw_mean_wb_now * Myaw_norm;

            %% store data
            n_now=n_now+1;

            seq_nr_YawTorque(n_now,1) = seq_nr_now;
            wb_nr_YawTorque(n_now,1) = wb_nr_now;
            
            % body kin
            V_YawTorque(n_now,1) = V_now;
            pitch_global_YawTorque(n_now,1) = pitch_global_now;
            F_YawTorque(n_now,1) = F_mean_wb_now;

            Mroll_accel_YawTorque(n_now,1) = Mroll_mean_wb_accel_now;
            Mpitch_accel_YawTorque(n_now,1) = Mpitch_mean_wb_accel_now;
            Myaw_accel_YawTorque(n_now,1) = Myaw_mean_wb_accel_now;
            M_L_accel_YawTorque(n_now,1) = M_L_mean_wb_accel_now;
            M_R_accel_YawTorque(n_now,1) = M_R_mean_wb_accel_now;

            Mroll_damp_YawTorque(n_now,1) = Mroll_mean_wb_damp_now;
            Mpitch_damp_YawTorque(n_now,1) = Mpitch_mean_wb_damp_now;
            Myaw_damp_YawTorque(n_now,1) = Myaw_mean_wb_damp_now;
            M_L_damp_YawTorque(n_now,1) = M_L_mean_wb_damp_now;
            M_R_damp_YawTorque(n_now,1) = M_R_mean_wb_damp_now;

            Mroll_YawTorque(n_now,1) = Mroll_mean_wb_now;
            Mpitch_YawTorque(n_now,1) = Mpitch_mean_wb_now;
            Myaw_YawTorque(n_now,1) = Myaw_mean_wb_now;
            M_L_YawTorque(n_now,1) = M_L_mean_wb_now;
            M_R_YawTorque(n_now,1) = M_R_mean_wb_now;

            rollaccel_YawTorque(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_YawTorque(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_YawTorque(n_now,1) = yawaccel_mean_wb_now;

            rollvel_YawTorque(n_now,1) = rollvel_mean_wb_now;
            pitchvel_YawTorque(n_now,1) = pitchvel_mean_wb_now;
            yawvel_YawTorque(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_YawTorque(n_now,1) = nanmean(dt_ds_now);
            dt_us_YawTorque(n_now,1) = nanmean(dt_us_now);
            f_wb_YawTorque(n_now,1) = nanmean(f_wb_now);
            Rds_YawTorque(n_now,1) = nanmean(Rds_now);

            t_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = t_wb_fwd_now;
            t_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = t_wb_rwd_now;
            stroke_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = stroke_wb_fwd_now;
            stroke_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = stroke_wb_rwd_now;
            pitch_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = pitch_wb_fwd_now;
            pitch_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = pitch_wb_rwd_now;
            dev_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = dev_wb_fwd_now;
            dev_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = dev_wb_rwd_now;
            aoa_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = aoa_wb_fwd_now;
            aoa_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = aoa_wb_rwd_now;
            U_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = U_wb_fwd_now;
            U_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = U_wb_rwd_now;
            
            t_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = t_ds_fwd_now;
            t_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = t_ds_rwd_now;
            stroke_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = stroke_ds_fwd_now;
            stroke_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = stroke_ds_rwd_now;
            pitch_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = pitch_ds_fwd_now;
            pitch_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = pitch_ds_rwd_now;
            dev_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = dev_ds_fwd_now;
            dev_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = dev_ds_rwd_now;
            aoa_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = aoa_ds_fwd_now;
            aoa_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = aoa_ds_rwd_now;
            U_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = U_ds_fwd_now;
            U_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = U_ds_rwd_now;
            
            t_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = t_us_fwd_now;
            t_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = t_us_rwd_now;
            stroke_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = stroke_us_fwd_now;
            stroke_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = stroke_us_rwd_now;
            pitch_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = pitch_us_fwd_now;
            pitch_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = pitch_us_rwd_now;
            dev_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = dev_us_fwd_now;
            dev_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = dev_us_rwd_now;
            aoa_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = aoa_us_fwd_now;
            aoa_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = aoa_us_rwd_now;
            U_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = U_us_fwd_now;
            U_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = U_us_rwd_now;
            
            % store interp binned data separate rows
            t_wb_YawTorque_bins(:,n_now) = t_wb_bin;
            stroke_wb_fwd_YawTorque_bins(:,n_now) = stroke_wb_fwd_interp;
            stroke_wb_rwd_YawTorque_bins(:,n_now) = stroke_wb_rwd_interp;
            pitch_wb_fwd_YawTorque_bins(:,n_now) = pitch_wb_fwd_interp;
            pitch_wb_rwd_YawTorque_bins(:,n_now) = pitch_wb_rwd_interp;
            dev_wb_fwd_YawTorque_bins(:,n_now) = dev_wb_fwd_interp;
            dev_wb_rwd_YawTorque_bins(:,n_now) = dev_wb_rwd_interp;
            aoa_wb_fwd_YawTorque_bins(:,n_now) = aoa_wb_fwd_interp;
            aoa_wb_rwd_YawTorque_bins(:,n_now) = aoa_wb_rwd_interp;
            U_wb_fwd_YawTorque_bins(:,n_now) = U_wb_fwd_interp;
            U_wb_rwd_YawTorque_bins(:,n_now) = U_wb_rwd_interp;
            Dstroke_wb_YawTorque_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_YawTorque_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_YawTorque_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_YawTorque_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_YawTorque_bins(:,n_now) = DU_wb_interp;
            
            t_ds_YawTorque_bins(:,n_now) = t_ds_bin;
            stroke_ds_fwd_YawTorque_bins(:,n_now) = stroke_ds_fwd_interp;
            stroke_ds_rwd_YawTorque_bins(:,n_now) = stroke_ds_rwd_interp;
            pitch_ds_fwd_YawTorque_bins(:,n_now) = pitch_ds_fwd_interp;
            pitch_ds_rwd_YawTorque_bins(:,n_now) = pitch_ds_rwd_interp;
            dev_ds_fwd_YawTorque_bins(:,n_now) = dev_ds_fwd_interp;
            dev_ds_rwd_YawTorque_bins(:,n_now) = dev_ds_rwd_interp;
            aoa_ds_fwd_YawTorque_bins(:,n_now) = aoa_ds_fwd_interp;
            aoa_ds_rwd_YawTorque_bins(:,n_now) = aoa_ds_rwd_interp;
            U_ds_fwd_YawTorque_bins(:,n_now) = U_ds_fwd_interp;
            U_ds_rwd_YawTorque_bins(:,n_now) = U_ds_rwd_interp;
            Dstroke_ds_YawTorque_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_YawTorque_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_YawTorque_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_YawTorque_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_YawTorque_bins(:,n_now) = DU_ds_interp;
            
            t_us_YawTorque_bins(:,n_now) = t_us_bin;
            stroke_us_fwd_YawTorque_bins(:,n_now) = stroke_us_fwd_interp;
            stroke_us_rwd_YawTorque_bins(:,n_now) = stroke_us_rwd_interp;
            pitch_us_fwd_YawTorque_bins(:,n_now) = pitch_us_fwd_interp;
            pitch_us_rwd_YawTorque_bins(:,n_now) = pitch_us_rwd_interp;
            dev_us_fwd_YawTorque_bins(:,n_now) = dev_us_fwd_interp;
            dev_us_rwd_YawTorque_bins(:,n_now) = dev_us_rwd_interp;
            aoa_us_fwd_YawTorque_bins(:,n_now) = aoa_us_fwd_interp;
            aoa_us_rwd_YawTorque_bins(:,n_now) = aoa_us_rwd_interp;
            U_us_fwd_YawTorque_bins(:,n_now) = U_us_fwd_interp;
            U_us_rwd_YawTorque_bins(:,n_now) = U_us_rwd_interp;
            Dstroke_us_YawTorque_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_YawTorque_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_YawTorque_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_YawTorque_bins(:,n_now) = Daoa_us_interp;
            DU_us_YawTorque_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_YawTorque(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_fwd_YawTorque_bins(:,n_now) = strokeMOD_wb_fwd_now;
            strokeMOD_wb_rwd_YawTorque_bins(:,n_now) = strokeMOD_wb_rwd_now;
            strokeMOD_ds_fwd_YawTorque_bins(:,n_now) = strokeMOD_ds_fwd_now;
            strokeMOD_ds_rwd_YawTorque_bins(:,n_now) = strokeMOD_ds_rwd_now;
            strokeMOD_us_fwd_YawTorque_bins(:,n_now) = strokeMOD_us_fwd_now;
            strokeMOD_us_rwd_YawTorque_bins(:,n_now) = strokeMOD_us_rwd_now;

            pitchMOD_wb_fwd_YawTorque_bins(:,n_now) = pitchMOD_wb_fwd_now;
            pitchMOD_wb_rwd_YawTorque_bins(:,n_now) = pitchMOD_wb_rwd_now;
            pitchMOD_ds_fwd_YawTorque_bins(:,n_now) = pitchMOD_ds_fwd_now;
            pitchMOD_ds_rwd_YawTorque_bins(:,n_now) = pitchMOD_ds_rwd_now;
            pitchMOD_us_fwd_YawTorque_bins(:,n_now) = pitchMOD_us_fwd_now;
            pitchMOD_us_rwd_YawTorque_bins(:,n_now) = pitchMOD_us_rwd_now;

            devMOD_wb_fwd_YawTorque_bins(:,n_now) = devMOD_wb_fwd_now;
            devMOD_wb_rwd_YawTorque_bins(:,n_now) = devMOD_wb_rwd_now;
            devMOD_ds_fwd_YawTorque_bins(:,n_now) = devMOD_ds_fwd_now;
            devMOD_ds_rwd_YawTorque_bins(:,n_now) = devMOD_ds_rwd_now;
            devMOD_us_fwd_YawTorque_bins(:,n_now) = devMOD_us_fwd_now;
            devMOD_us_rwd_YawTorque_bins(:,n_now) = devMOD_us_rwd_now;
            
            DstrokeMOD_wb_YawTorque_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_YawTorque_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_YawTorque_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_YawTorque_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_YawTorque_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_YawTorque_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_YawTorque_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_YawTorque_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_YawTorque_bins(:,n_now) = DdevMOD_us_now;
            
        elseif Myaw_mean_wb(wb) < -Myaw_limit_mod
%         elseif yaw_dot_dot_mean_wb(wb) < -yawaccel_limit_mod
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            F_mean_wb_now = F_mean_wb(wb);

            Mroll_mean_wb_accel_now = Mroll_mean_wb_accel(wb);
            Mpitch_mean_wb_accel_now = Mpitch_mean_wb_accel(wb);
            Myaw_mean_wb_accel_now = Myaw_mean_wb_accel(wb);
            M_L_mean_wb_accel_now = M_L_mean_wb_accel(wb);
            M_R_mean_wb_accel_now = M_R_mean_wb_accel(wb);

            Mroll_mean_wb_damp_now = Mroll_mean_wb_damp(wb);
            Mpitch_mean_wb_damp_now = Mpitch_mean_wb_damp(wb);
            Myaw_mean_wb_damp_now = Myaw_mean_wb_damp(wb);
            M_L_mean_wb_damp_now = M_L_mean_wb_damp(wb);
            M_R_mean_wb_damp_now = M_R_mean_wb_damp(wb);

            Mroll_mean_wb_now = Mroll_mean_wb(wb);
            Mpitch_mean_wb_now = Mpitch_mean_wb(wb);
            Myaw_mean_wb_now = Myaw_mean_wb(wb);
            M_L_mean_wb_now = M_L_mean_wb(wb);
            M_R_mean_wb_now = M_R_mean_wb(wb);

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
            
            stroke_wb_fwd_now = stroke_wb_R(:,wb);
            stroke_wb_rwd_now = stroke_wb_L(:,wb);
            stroke_ds_fwd_now = stroke_ds_R(:,wb);
            stroke_ds_rwd_now = stroke_ds_L(:,wb);
            stroke_us_fwd_now = stroke_us_R(:,wb);
            stroke_us_rwd_now = stroke_us_L(:,wb);
            
            pitch_wb_fwd_now = pitch_wb_R(:,wb);
            pitch_wb_rwd_now = pitch_wb_L(:,wb);
            pitch_ds_fwd_now = pitch_ds_R(:,wb);
            pitch_ds_rwd_now = pitch_ds_L(:,wb);
            pitch_us_fwd_now = pitch_us_R(:,wb);
            pitch_us_rwd_now = pitch_us_L(:,wb);
            
            dev_wb_fwd_now = dev_wb_R(:,wb);
            dev_wb_rwd_now = dev_wb_L(:,wb);
            dev_ds_fwd_now = dev_ds_R(:,wb);
            dev_ds_rwd_now = dev_ds_L(:,wb);
            dev_us_fwd_now = dev_us_R(:,wb);
            dev_us_rwd_now = dev_us_L(:,wb);
            
            U_wb_fwd_now = U_wb_R(:,wb);
            U_wb_rwd_now = U_wb_L(:,wb);
            U_ds_fwd_now = U_ds_R(:,wb);
            U_ds_rwd_now = U_ds_L(:,wb);
            U_us_fwd_now = U_us_R(:,wb);
            U_us_rwd_now = U_us_L(:,wb);
            
            aoa_wb_fwd_now = aoa_wb_R(:,wb);
            aoa_wb_rwd_now = aoa_wb_L(:,wb);
            aoa_ds_fwd_now = aoa_ds_R(:,wb);
            aoa_ds_rwd_now = aoa_ds_L(:,wb);
            aoa_us_fwd_now = aoa_us_R(:,wb);
            aoa_us_rwd_now = aoa_us_L(:,wb);
            
            t_wb_fwd_now = t_wb_R(:,wb);
            t_wb_rwd_now = t_wb_L(:,wb);
            t_ds_fwd_now = t_ds_R(:,wb);
            t_ds_rwd_now = t_ds_L(:,wb);
            t_us_fwd_now = t_us_R(:,wb);
            t_us_rwd_now = t_us_L(:,wb);

            % interpolated data
            stroke_wb_fwd_interp = stroke_wb_R_bins(:,wb);
            stroke_wb_rwd_interp = stroke_wb_L_bins(:,wb);
            stroke_ds_fwd_interp = stroke_ds_R_bins(:,wb);
            stroke_ds_rwd_interp = stroke_ds_L_bins(:,wb);
            stroke_us_fwd_interp = stroke_us_R_bins(:,wb);
            stroke_us_rwd_interp = stroke_us_L_bins(:,wb);
            
            pitch_wb_fwd_interp = pitch_wb_R_bins(:,wb);
            pitch_wb_rwd_interp = pitch_wb_L_bins(:,wb);
            pitch_ds_fwd_interp = pitch_ds_R_bins(:,wb);
            pitch_ds_rwd_interp = pitch_ds_L_bins(:,wb);
            pitch_us_fwd_interp = pitch_us_R_bins(:,wb);
            pitch_us_rwd_interp = pitch_us_L_bins(:,wb);
            
            dev_wb_fwd_interp = dev_wb_R_bins(:,wb);
            dev_wb_rwd_interp = dev_wb_L_bins(:,wb);
            dev_ds_fwd_interp = dev_ds_R_bins(:,wb);
            dev_ds_rwd_interp = dev_ds_L_bins(:,wb);
            dev_us_fwd_interp = dev_us_R_bins(:,wb);
            dev_us_rwd_interp = dev_us_L_bins(:,wb);
            
            U_wb_fwd_interp = U_wb_R_bins(:,wb);
            U_wb_rwd_interp = U_wb_L_bins(:,wb);
            U_ds_fwd_interp = U_ds_R_bins(:,wb);
            U_ds_rwd_interp = U_ds_L_bins(:,wb);
            U_us_fwd_interp = U_us_R_bins(:,wb);
            U_us_rwd_interp = U_us_L_bins(:,wb);
            
            aoa_wb_fwd_interp = aoa_wb_R_bins(:,wb);
            aoa_wb_rwd_interp = aoa_wb_L_bins(:,wb);
            aoa_ds_fwd_interp = aoa_ds_R_bins(:,wb);
            aoa_ds_rwd_interp = aoa_ds_L_bins(:,wb);
            aoa_us_fwd_interp = aoa_us_R_bins(:,wb);
            aoa_us_rwd_interp = aoa_us_L_bins(:,wb);

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
            freqMOD_wb_now = -(nanmean(f_wb_now) - f_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            strokeMOD_wb_fwd_now = -(stroke_wb_fwd_interp - stroke_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_wb_rwd_now = -(stroke_wb_rwd_interp - stroke_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_ds_fwd_now = -(stroke_ds_fwd_interp - stroke_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_ds_rwd_now = -(stroke_ds_rwd_interp - stroke_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_us_fwd_now = -(stroke_us_fwd_interp - stroke_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            strokeMOD_us_rwd_now = -(stroke_us_rwd_interp - stroke_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            pitchMOD_wb_fwd_now = -(pitch_wb_fwd_interp - pitch_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_wb_rwd_now = -(pitch_wb_rwd_interp - pitch_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_ds_fwd_now = -(pitch_ds_fwd_interp - pitch_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_ds_rwd_now = -(pitch_ds_rwd_interp - pitch_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_us_fwd_now = -(pitch_us_fwd_interp - pitch_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            pitchMOD_us_rwd_now = -(pitch_us_rwd_interp - pitch_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            devMOD_wb_fwd_now = -(dev_wb_fwd_interp - dev_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_wb_rwd_now = -(dev_wb_rwd_interp - dev_wb_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_ds_fwd_now = -(dev_ds_fwd_interp - dev_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_ds_rwd_now = -(dev_ds_rwd_interp - dev_ds_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_us_fwd_now = -(dev_us_fwd_interp - dev_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            devMOD_us_rwd_now = -(dev_us_rwd_interp - dev_us_steady) / Myaw_mean_wb_now * Myaw_norm;
            
            DstrokeMOD_wb_now = -Dstroke_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DstrokeMOD_ds_now = -Dstroke_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DstrokeMOD_us_now = -Dstroke_us_interp / Myaw_mean_wb_now * Myaw_norm;

            DpitchMOD_wb_now = -Dpitch_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DpitchMOD_ds_now = -Dpitch_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DpitchMOD_us_now = -Dpitch_us_interp / Myaw_mean_wb_now * Myaw_norm;

            DdevMOD_wb_now = -Ddev_wb_interp / Myaw_mean_wb_now * Myaw_norm;
            DdevMOD_ds_now = -Ddev_ds_interp / Myaw_mean_wb_now * Myaw_norm;
            DdevMOD_us_now = -Ddev_us_interp / Myaw_mean_wb_now * Myaw_norm;

            %% store data
            n_now=n_now+1;

            seq_nr_YawTorque(n_now,1) = seq_nr_now;
            wb_nr_YawTorque(n_now,1) = wb_nr_now;
            
            % body kin
            V_YawTorque(n_now,1) = V_now;
            pitch_global_YawTorque(n_now,1) = pitch_global_now;
            F_YawTorque(n_now,1) = F_mean_wb_now;

            Mroll_accel_YawTorque(n_now,1) = Mroll_mean_wb_accel_now;
            Mpitch_accel_YawTorque(n_now,1) = Mpitch_mean_wb_accel_now;
            Myaw_accel_YawTorque(n_now,1) = Myaw_mean_wb_accel_now;
            M_L_accel_YawTorque(n_now,1) = M_L_mean_wb_accel_now;
            M_R_accel_YawTorque(n_now,1) = M_R_mean_wb_accel_now;

            Mroll_damp_YawTorque(n_now,1) = Mroll_mean_wb_damp_now;
            Mpitch_damp_YawTorque(n_now,1) = Mpitch_mean_wb_damp_now;
            Myaw_damp_YawTorque(n_now,1) = Myaw_mean_wb_damp_now;
            M_L_damp_YawTorque(n_now,1) = M_L_mean_wb_damp_now;
            M_R_damp_YawTorque(n_now,1) = M_R_mean_wb_damp_now;

            Mroll_YawTorque(n_now,1) = Mroll_mean_wb_now;
            Mpitch_YawTorque(n_now,1) = Mpitch_mean_wb_now;
            Myaw_YawTorque(n_now,1) = Myaw_mean_wb_now;
            M_L_YawTorque(n_now,1) = M_L_mean_wb_now;
            M_R_YawTorque(n_now,1) = M_R_mean_wb_now;

            rollaccel_YawTorque(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_YawTorque(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_YawTorque(n_now,1) = yawaccel_mean_wb_now;

            rollvel_YawTorque(n_now,1) = rollvel_mean_wb_now;
            pitchvel_YawTorque(n_now,1) = pitchvel_mean_wb_now;
            yawvel_YawTorque(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_YawTorque(n_now,1) = nanmean(dt_ds_now);
            dt_us_YawTorque(n_now,1) = nanmean(dt_us_now);
            f_wb_YawTorque(n_now,1) = nanmean(f_wb_now);
            Rds_YawTorque(n_now,1) = nanmean(Rds_now);

            t_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = t_wb_fwd_now;
            t_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = t_wb_rwd_now;
            stroke_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = stroke_wb_fwd_now;
            stroke_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = stroke_wb_rwd_now;
            pitch_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = pitch_wb_fwd_now;
            pitch_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = pitch_wb_rwd_now;
            dev_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = dev_wb_fwd_now;
            dev_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = dev_wb_rwd_now;
            aoa_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = aoa_wb_fwd_now;
            aoa_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = aoa_wb_rwd_now;
            U_wb_fwd_YawTorque(1:length(t_wb_fwd_now),n_now) = U_wb_fwd_now;
            U_wb_rwd_YawTorque(1:length(t_wb_rwd_now),n_now) = U_wb_rwd_now;
            
            t_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = t_ds_fwd_now;
            t_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = t_ds_rwd_now;
            stroke_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = stroke_ds_fwd_now;
            stroke_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = stroke_ds_rwd_now;
            pitch_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = pitch_ds_fwd_now;
            pitch_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = pitch_ds_rwd_now;
            dev_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = dev_ds_fwd_now;
            dev_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = dev_ds_rwd_now;
            aoa_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = aoa_ds_fwd_now;
            aoa_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = aoa_ds_rwd_now;
            U_ds_fwd_YawTorque(1:length(t_ds_fwd_now),n_now) = U_ds_fwd_now;
            U_ds_rwd_YawTorque(1:length(t_ds_rwd_now),n_now) = U_ds_rwd_now;
            
            t_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = t_us_fwd_now;
            t_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = t_us_rwd_now;
            stroke_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = stroke_us_fwd_now;
            stroke_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = stroke_us_rwd_now;
            pitch_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = pitch_us_fwd_now;
            pitch_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = pitch_us_rwd_now;
            dev_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = dev_us_fwd_now;
            dev_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = dev_us_rwd_now;
            aoa_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = aoa_us_fwd_now;
            aoa_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = aoa_us_rwd_now;
            U_us_fwd_YawTorque(1:length(t_us_fwd_now),n_now) = U_us_fwd_now;
            U_us_rwd_YawTorque(1:length(t_us_rwd_now),n_now) = U_us_rwd_now;
            
            % store interp binned data separate rows
            t_wb_YawTorque_bins(:,n_now) = t_wb_bin;
            stroke_wb_fwd_YawTorque_bins(:,n_now) = stroke_wb_fwd_interp;
            stroke_wb_rwd_YawTorque_bins(:,n_now) = stroke_wb_rwd_interp;
            pitch_wb_fwd_YawTorque_bins(:,n_now) = pitch_wb_fwd_interp;
            pitch_wb_rwd_YawTorque_bins(:,n_now) = pitch_wb_rwd_interp;
            dev_wb_fwd_YawTorque_bins(:,n_now) = dev_wb_fwd_interp;
            dev_wb_rwd_YawTorque_bins(:,n_now) = dev_wb_rwd_interp;
            aoa_wb_fwd_YawTorque_bins(:,n_now) = aoa_wb_fwd_interp;
            aoa_wb_rwd_YawTorque_bins(:,n_now) = aoa_wb_rwd_interp;
            U_wb_fwd_YawTorque_bins(:,n_now) = U_wb_fwd_interp;
            U_wb_rwd_YawTorque_bins(:,n_now) = U_wb_rwd_interp;
            Dstroke_wb_YawTorque_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_YawTorque_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_YawTorque_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_YawTorque_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_YawTorque_bins(:,n_now) = DU_wb_interp;
            
            t_ds_YawTorque_bins(:,n_now) = t_ds_bin;
            stroke_ds_fwd_YawTorque_bins(:,n_now) = stroke_ds_fwd_interp;
            stroke_ds_rwd_YawTorque_bins(:,n_now) = stroke_ds_rwd_interp;
            pitch_ds_fwd_YawTorque_bins(:,n_now) = pitch_ds_fwd_interp;
            pitch_ds_rwd_YawTorque_bins(:,n_now) = pitch_ds_rwd_interp;
            dev_ds_fwd_YawTorque_bins(:,n_now) = dev_ds_fwd_interp;
            dev_ds_rwd_YawTorque_bins(:,n_now) = dev_ds_rwd_interp;
            aoa_ds_fwd_YawTorque_bins(:,n_now) = aoa_ds_fwd_interp;
            aoa_ds_rwd_YawTorque_bins(:,n_now) = aoa_ds_rwd_interp;
            U_ds_fwd_YawTorque_bins(:,n_now) = U_ds_fwd_interp;
            U_ds_rwd_YawTorque_bins(:,n_now) = U_ds_rwd_interp;
            Dstroke_ds_YawTorque_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_YawTorque_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_YawTorque_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_YawTorque_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_YawTorque_bins(:,n_now) = DU_ds_interp;
            
            t_us_YawTorque_bins(:,n_now) = t_us_bin;
            stroke_us_fwd_YawTorque_bins(:,n_now) = stroke_us_fwd_interp;
            stroke_us_rwd_YawTorque_bins(:,n_now) = stroke_us_rwd_interp;
            pitch_us_fwd_YawTorque_bins(:,n_now) = pitch_us_fwd_interp;
            pitch_us_rwd_YawTorque_bins(:,n_now) = pitch_us_rwd_interp;
            dev_us_fwd_YawTorque_bins(:,n_now) = dev_us_fwd_interp;
            dev_us_rwd_YawTorque_bins(:,n_now) = dev_us_rwd_interp;
            aoa_us_fwd_YawTorque_bins(:,n_now) = aoa_us_fwd_interp;
            aoa_us_rwd_YawTorque_bins(:,n_now) = aoa_us_rwd_interp;
            U_us_fwd_YawTorque_bins(:,n_now) = U_us_fwd_interp;
            U_us_rwd_YawTorque_bins(:,n_now) = U_us_rwd_interp;
            Dstroke_us_YawTorque_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_YawTorque_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_YawTorque_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_YawTorque_bins(:,n_now) = Daoa_us_interp;
            DU_us_YawTorque_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_YawTorque(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_fwd_YawTorque_bins(:,n_now) = strokeMOD_wb_fwd_now;
            strokeMOD_wb_rwd_YawTorque_bins(:,n_now) = strokeMOD_wb_rwd_now;
            strokeMOD_ds_fwd_YawTorque_bins(:,n_now) = strokeMOD_ds_fwd_now;
            strokeMOD_ds_rwd_YawTorque_bins(:,n_now) = strokeMOD_ds_rwd_now;
            strokeMOD_us_fwd_YawTorque_bins(:,n_now) = strokeMOD_us_fwd_now;
            strokeMOD_us_rwd_YawTorque_bins(:,n_now) = strokeMOD_us_rwd_now;

            pitchMOD_wb_fwd_YawTorque_bins(:,n_now) = pitchMOD_wb_fwd_now;
            pitchMOD_wb_rwd_YawTorque_bins(:,n_now) = pitchMOD_wb_rwd_now;
            pitchMOD_ds_fwd_YawTorque_bins(:,n_now) = pitchMOD_ds_fwd_now;
            pitchMOD_ds_rwd_YawTorque_bins(:,n_now) = pitchMOD_ds_rwd_now;
            pitchMOD_us_fwd_YawTorque_bins(:,n_now) = pitchMOD_us_fwd_now;
            pitchMOD_us_rwd_YawTorque_bins(:,n_now) = pitchMOD_us_rwd_now;

            devMOD_wb_fwd_YawTorque_bins(:,n_now) = devMOD_wb_fwd_now;
            devMOD_wb_rwd_YawTorque_bins(:,n_now) = devMOD_wb_rwd_now;
            devMOD_ds_fwd_YawTorque_bins(:,n_now) = devMOD_ds_fwd_now;
            devMOD_ds_rwd_YawTorque_bins(:,n_now) = devMOD_ds_rwd_now;
            devMOD_us_fwd_YawTorque_bins(:,n_now) = devMOD_us_fwd_now;
            devMOD_us_rwd_YawTorque_bins(:,n_now) = devMOD_us_rwd_now;
            
            DstrokeMOD_wb_YawTorque_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_YawTorque_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_YawTorque_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_YawTorque_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_YawTorque_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_YawTorque_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_YawTorque_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_YawTorque_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_YawTorque_bins(:,n_now) = DdevMOD_us_now;
        end
end

% mean & 95%CI
V_YawTorque_meanCIstd = [nanmean(V_YawTorque) 1.96*nanstd(V_YawTorque)/sqrt(length(V_YawTorque)) nanstd(V_YawTorque)];
pitch_global_YawTorque_meanCIstd = [nanmean(pitch_global_YawTorque) 1.96*nanstd(pitch_global_YawTorque)/sqrt(length(pitch_global_YawTorque)) nanstd(pitch_global_YawTorque)];

dt_ds_YawTorque_meanCIstd = [nanmean(dt_ds_YawTorque) 1.96*nanstd(dt_ds_YawTorque)/sqrt(length(dt_ds_YawTorque)) nanstd(dt_ds_YawTorque)];
dt_us_YawTorque_meanCIstd = [nanmean(dt_us_YawTorque) 1.96*nanstd(dt_us_YawTorque)/sqrt(length(dt_us_YawTorque)) nanstd(dt_us_YawTorque)];
f_wb_YawTorque_meanCIstd = [nanmean(f_wb_YawTorque) 1.96*nanstd(f_wb_YawTorque)/sqrt(length(f_wb_YawTorque)) nanstd(f_wb_YawTorque)];
Rds_YawTorque_meanCIstd = [nanmean(Rds_YawTorque) 1.96*nanstd(Rds_YawTorque)/sqrt(length(Rds_YawTorque)) nanstd(Rds_YawTorque)];

calc_WBfunc_YawTorque_circmeanCIstd

% wbMOD means & 95%CI
freqMOD_wb_YawTorque_meanCIstd = [nanmean(freqMOD_wb_YawTorque) 1.96*nanstd(freqMOD_wb_YawTorque)/sqrt(length(freqMOD_wb_YawTorque)) nanstd(freqMOD_wb_YawTorque)];
calc_WBmod_YawTorque_circmeanCIstd

%% WBfits
    t_loc = t_wb_YawTorque_bins(:,1);
    Rds_loc = Rds_YawTorque_meanCIstd(1);

    plotting = 0;
%         plotting = 1;

    %% fit for forwards moving wing (bins)
    stroke_loc = stroke_wb_fwd_YawTorque_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_fwd_YawTorque_bins_meanCIstd(:,1);
    dev_loc = dev_wb_fwd_YawTorque_bins_meanCIstd(:,1);
    
    strokeMOD_loc = strokeMOD_wb_fwd_YawTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_fwd_YawTorque_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_fwd_YawTorque_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_fwd_YawTorque_fit_binmean, stroke_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_fwd_YawTorque_fit_binmean, pitch_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_fwd_YawTorque_fit_binmean, dev_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_fwd_YawTorque_fit_binmean, strokeMOD_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_fwd_YawTorque_fit_binmean, pitchMOD_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_fwd_YawTorque_fit_binmean, devMOD_fwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_fwd_YawTorque_fourier_fit_binmean, stroke_fwd_YawTorque_fourier_gof_binmean,stroke_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_fwd_YawTorque_fourier_fit_binmean, pitch_fwd_YawTorque_fourier_gof_binmean,pitch_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_fwd_YawTorque_fourier_fit_binmean, dev_fwd_YawTorque_fourier_gof_binmean,dev_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_fwd_YawTorque_fourier_fit_binmean, strokeMOD_fwd_YawTorque_fourier_gof_binmean,strokeMOD_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_fwd_YawTorque_fourier_fit_binmean, pitchMOD_fwd_YawTorque_fourier_gof_binmean,pitchMOD_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_fwd_YawTorque_fourier_fit_binmean, devMOD_fwd_YawTorque_fourier_gof_binmean,devMOD_fwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% rearwards
    stroke_loc = stroke_wb_rwd_YawTorque_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_rwd_YawTorque_bins_meanCIstd(:,1);
    dev_loc = dev_wb_rwd_YawTorque_bins_meanCIstd(:,1);

    strokeMOD_loc = strokeMOD_wb_rwd_YawTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_rwd_YawTorque_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_rwd_YawTorque_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_rwd_YawTorque_fit_binmean, stroke_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_rwd_YawTorque_fit_binmean, pitch_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_rwd_YawTorque_fit_binmean, dev_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_rwd_YawTorque_fit_binmean, strokeMOD_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_rwd_YawTorque_fit_binmean, pitchMOD_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_rwd_YawTorque_fit_binmean, devMOD_rwd_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_rwd_YawTorque_fourier_fit_binmean, stroke_rwd_YawTorque_fourier_gof_binmean,stroke_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_rwd_YawTorque_fourier_fit_binmean, pitch_rwd_YawTorque_fourier_gof_binmean,pitch_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_rwd_YawTorque_fourier_fit_binmean, dev_rwd_YawTorque_fourier_gof_binmean,dev_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_rwd_YawTorque_fourier_fit_binmean, strokeMOD_rwd_YawTorque_fourier_gof_binmean,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_rwd_YawTorque_fourier_fit_binmean, pitchMOD_rwd_YawTorque_fourier_gof_binmean,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_rwd_YawTorque_fourier_fit_binmean, devMOD_rwd_YawTorque_fourier_gof_binmean,devMOD_rwd_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% forwards - rearwards
    stroke_loc = Dstroke_wb_YawTorque_bins_meanCIstd(:,1);
    pitch_loc = Dpitch_wb_YawTorque_bins_meanCIstd(:,1);
    dev_loc = Ddev_wb_YawTorque_bins_meanCIstd(:,1);

    strokeMOD_loc = DstrokeMOD_wb_YawTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = DpitchMOD_wb_YawTorque_bins_meanCIstd(:,1);
    devMOD_loc = DdevMOD_wb_YawTorque_bins_meanCIstd(:,1);
    
%% legendre polynomials
    % wingbeats
    [Dstroke_YawTorque_fit_binmean, Dstroke_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [Dpitch_YawTorque_fit_binmean, Dpitch_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [Ddev_YawTorque_fit_binmean, Ddev_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [DstrokeMOD_YawTorque_fit_binmean, DstrokeMOD_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [DpitchMOD_YawTorque_fit_binmean, DpitchMOD_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [DdevMOD_YawTorque_fit_binmean, DdevMOD_YawTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [Dstroke_YawTorque_fourier_fit_binmean, Dstroke_YawTorque_fourier_gof_binmean,Dstroke_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [Dpitch_YawTorque_fourier_fit_binmean, Dpitch_YawTorque_fourier_gof_binmean,Dpitch_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [Ddev_YawTorque_fourier_fit_binmean, Ddev_YawTorque_fourier_gof_binmean,Ddev_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [DstrokeMOD_YawTorque_fourier_fit_binmean, DstrokeMOD_YawTorque_fourier_gof_binmean,DstrokeMOD_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [DpitchMOD_YawTorque_fourier_fit_binmean, DpitchMOD_YawTorque_fourier_gof_binmean,DpitchMOD_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [DdevMOD_YawTorque_fourier_fit_binmean, DdevMOD_YawTorque_fourier_gof_binmean,DdevMOD_YawTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    

    
    
    
    
    %% save YawTorque wb data
n_YawTorque = n_now   

if save_on == 1

save(['WBmod_torquebased_YawTorque_',num2str(n_YawTorque),'WBs.mat'],...
'n_YawTorque',...
'seq_nr_YawTorque',...
'wb_nr_YawTorque',...
...
'limit_mod',...
'Mroll_limit_mod',...
'Mpitch_limit_mod',...
'Myaw_limit_mod',...
'M_L_limit_mod',...
'M_R_limit_mod',...
'Fenhance_limit_mod',...
...
'norm',...
'Mroll_norm',...
'Mpitch_norm',...
'Myaw_norm',...
'M_L_norm',...
'M_R_norm',...
'Fenhance_norm',...
...
'n_pol_stroke',...
'n_pol_pitch',...
'n_pol_dev',...
'n_pol_MOD',...
...
'stroke_fwd_YawTorque_fit_binmean_periodic',...
'pitch_fwd_YawTorque_fit_binmean_periodic',...
'dev_fwd_YawTorque_fit_binmean_periodic',...
...
'stroke_rwd_YawTorque_fit_binmean_periodic',...
'pitch_rwd_YawTorque_fit_binmean_periodic',...
'dev_rwd_YawTorque_fit_binmean_periodic',...
...
'strokeMOD_fwd_YawTorque_fit_binmean_periodic',...
'pitchMOD_fwd_YawTorque_fit_binmean_periodic',...
'devMOD_fwd_YawTorque_fit_binmean_periodic',...
...
'strokeMOD_rwd_YawTorque_fit_binmean_periodic',...
'pitchMOD_rwd_YawTorque_fit_binmean_periodic',...
'devMOD_rwd_YawTorque_fit_binmean_periodic',...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_fwd_YawTorque_fourier_coeffs_binmean',...
'pitch_fwd_YawTorque_fourier_coeffs_binmean',...
'dev_fwd_YawTorque_fourier_coeffs_binmean',...
...
'stroke_fwd_YawTorque_fourier_gof_binmean',...
'pitch_fwd_YawTorque_fourier_gof_binmean',...
'dev_fwd_YawTorque_fourier_gof_binmean',...
...
'stroke_rwd_YawTorque_fourier_coeffs_binmean',...
'pitch_rwd_YawTorque_fourier_coeffs_binmean',...
'dev_rwd_YawTorque_fourier_coeffs_binmean',...
...
'stroke_rwd_YawTorque_fourier_gof_binmean',...
'pitch_rwd_YawTorque_fourier_gof_binmean',...
'dev_rwd_YawTorque_fourier_gof_binmean',...
...
'Dstroke_YawTorque_fourier_coeffs_binmean',...
'Dpitch_YawTorque_fourier_coeffs_binmean',...
'Ddev_YawTorque_fourier_coeffs_binmean',...
...
'Dstroke_YawTorque_fourier_gof_binmean',...
'Dpitch_YawTorque_fourier_gof_binmean',...
'Ddev_YawTorque_fourier_gof_binmean',...
...
'strokeMOD_fwd_YawTorque_fourier_coeffs_binmean',...
'pitchMOD_fwd_YawTorque_fourier_coeffs_binmean',...
'devMOD_fwd_YawTorque_fourier_coeffs_binmean',...
...
'strokeMOD_fwd_YawTorque_fourier_gof_binmean',...
'pitchMOD_fwd_YawTorque_fourier_gof_binmean',...
'devMOD_fwd_YawTorque_fourier_gof_binmean',...
...
'strokeMOD_rwd_YawTorque_fourier_coeffs_binmean',...
'pitchMOD_rwd_YawTorque_fourier_coeffs_binmean',...
'devMOD_rwd_YawTorque_fourier_coeffs_binmean',...
...
'strokeMOD_rwd_YawTorque_fourier_gof_binmean',...
'pitchMOD_rwd_YawTorque_fourier_gof_binmean',...
'devMOD_rwd_YawTorque_fourier_gof_binmean',...
...
'DstrokeMOD_YawTorque_fourier_coeffs_binmean',...
'DpitchMOD_YawTorque_fourier_coeffs_binmean',...
'DdevMOD_YawTorque_fourier_coeffs_binmean',...
...
'DstrokeMOD_YawTorque_fourier_gof_binmean',...
'DpitchMOD_YawTorque_fourier_gof_binmean',...
'DdevMOD_YawTorque_fourier_gof_binmean',...
...
'stroke_wb_fwd_YawTorque_bins_meanCIstd',...
'stroke_ds_fwd_YawTorque_bins_meanCIstd',...
'stroke_us_fwd_YawTorque_bins_meanCIstd',...
...
'pitch_wb_fwd_YawTorque_bins_meanCIstd',...
'pitch_ds_fwd_YawTorque_bins_meanCIstd',...
'pitch_us_fwd_YawTorque_bins_meanCIstd',...
...
'dev_wb_fwd_YawTorque_bins_meanCIstd',...
'dev_ds_fwd_YawTorque_bins_meanCIstd',...
'dev_us_fwd_YawTorque_bins_meanCIstd',...
...
'strokeMOD_wb_fwd_YawTorque_bins_meanCIstd',...
'strokeMOD_ds_fwd_YawTorque_bins_meanCIstd',...
'strokeMOD_us_fwd_YawTorque_bins_meanCIstd',...
...
'pitchMOD_wb_fwd_YawTorque_bins_meanCIstd',...
'pitchMOD_ds_fwd_YawTorque_bins_meanCIstd',...
'pitchMOD_us_fwd_YawTorque_bins_meanCIstd',...
...
'devMOD_wb_fwd_YawTorque_bins_meanCIstd',...
'devMOD_ds_fwd_YawTorque_bins_meanCIstd',...
'devMOD_us_fwd_YawTorque_bins_meanCIstd',...
...
'stroke_wb_rwd_YawTorque_bins_meanCIstd',...
'stroke_ds_rwd_YawTorque_bins_meanCIstd',...
'stroke_us_rwd_YawTorque_bins_meanCIstd',...
...
'pitch_wb_rwd_YawTorque_bins_meanCIstd',...
'pitch_ds_rwd_YawTorque_bins_meanCIstd',...
'pitch_us_rwd_YawTorque_bins_meanCIstd',...
...
'dev_wb_rwd_YawTorque_bins_meanCIstd',...
'dev_ds_rwd_YawTorque_bins_meanCIstd',...
'dev_us_rwd_YawTorque_bins_meanCIstd',...
...
'strokeMOD_wb_rwd_YawTorque_bins_meanCIstd',...
'strokeMOD_ds_rwd_YawTorque_bins_meanCIstd',...
'strokeMOD_us_rwd_YawTorque_bins_meanCIstd',...
...
'pitchMOD_wb_rwd_YawTorque_bins_meanCIstd',...
'pitchMOD_ds_rwd_YawTorque_bins_meanCIstd',...
'pitchMOD_us_rwd_YawTorque_bins_meanCIstd',...
...
'devMOD_wb_rwd_YawTorque_bins_meanCIstd',...
'devMOD_ds_rwd_YawTorque_bins_meanCIstd',...
'devMOD_us_rwd_YawTorque_bins_meanCIstd',...
...
'DstrokeMOD_wb_YawTorque_bins_meanCIstd',...
'DstrokeMOD_ds_YawTorque_bins_meanCIstd',...
'DstrokeMOD_us_YawTorque_bins_meanCIstd',...
...
'DpitchMOD_wb_YawTorque_bins_meanCIstd',...
'DpitchMOD_ds_YawTorque_bins_meanCIstd',...
'DpitchMOD_us_YawTorque_bins_meanCIstd',...
...
'DdevMOD_wb_YawTorque_bins_meanCIstd',...
'DdevMOD_ds_YawTorque_bins_meanCIstd',...
'DdevMOD_us_YawTorque_bins_meanCIstd',...
...
'V_YawTorque',...
'pitch_global_YawTorque',...
'dt_ds_YawTorque',...
'dt_us_YawTorque',...
'f_wb_YawTorque',...
'Rds_YawTorque',...
...
'V_YawTorque_meanCIstd',...
'pitch_global_YawTorque_meanCIstd',...
'F_YawTorque',...
...
'Mroll_accel_YawTorque',...
'Mpitch_accel_YawTorque',...
'Myaw_accel_YawTorque',...
'M_L_accel_YawTorque',...
'M_R_accel_YawTorque',...
...
'Mroll_damp_YawTorque',...
'Mpitch_damp_YawTorque',...
'Myaw_damp_YawTorque',...
'M_L_damp_YawTorque',...
'M_R_damp_YawTorque',...
...
'Mroll_YawTorque',...
'Mpitch_YawTorque',...
'Myaw_YawTorque',...
'M_L_YawTorque',...
'M_R_YawTorque',...
...
'rollaccel_YawTorque',...
'pitchaccel_YawTorque',...
'yawaccel_YawTorque',...
...
'rollvel_YawTorque',...
'pitchvel_YawTorque',...
'yawvel_YawTorque',...
...
'dt_ds_YawTorque_meanCIstd',...
'dt_us_YawTorque_meanCIstd',...
'f_wb_YawTorque_meanCIstd',...
'Rds_YawTorque_meanCIstd',...
...
'freqMOD_wb_YawTorque_meanCIstd',...
'freqMOD_wb_YawTorque',...
...
'strokeMOD_wb_fwd_YawTorque_bins',...
'strokeMOD_wb_rwd_YawTorque_bins',...
'strokeMOD_ds_fwd_YawTorque_bins',...
'strokeMOD_ds_rwd_YawTorque_bins',...
'strokeMOD_us_fwd_YawTorque_bins',...
'strokeMOD_us_rwd_YawTorque_bins',...
...
'pitchMOD_wb_fwd_YawTorque_bins',...
'pitchMOD_wb_rwd_YawTorque_bins',...
'pitchMOD_ds_fwd_YawTorque_bins',...
'pitchMOD_ds_rwd_YawTorque_bins',...
'pitchMOD_us_fwd_YawTorque_bins',...
'pitchMOD_us_rwd_YawTorque_bins',...
...
'devMOD_wb_fwd_YawTorque_bins',...
'devMOD_wb_rwd_YawTorque_bins',...
'devMOD_ds_fwd_YawTorque_bins',...
'devMOD_ds_rwd_YawTorque_bins',...
'devMOD_us_fwd_YawTorque_bins',...
'devMOD_us_rwd_YawTorque_bins',...
...
'DstrokeMOD_wb_YawTorque_bins',...
'DstrokeMOD_ds_YawTorque_bins',...
'DstrokeMOD_us_YawTorque_bins',...
...
'DpitchMOD_wb_YawTorque_bins',...
'DpitchMOD_ds_YawTorque_bins',...
'DpitchMOD_us_YawTorque_bins',...
...
'DdevMOD_wb_YawTorque_bins',...
'DdevMOD_ds_YawTorque_bins',...
'DdevMOD_us_YawTorque_bins',...
...
't_wb_fwd_YawTorque',...
't_wb_rwd_YawTorque',...
'stroke_wb_fwd_YawTorque',...
'stroke_wb_rwd_YawTorque',...
'pitch_wb_fwd_YawTorque',...
'pitch_wb_rwd_YawTorque',...
'dev_wb_fwd_YawTorque',...
'dev_wb_rwd_YawTorque',...
'aoa_wb_fwd_YawTorque',...
'aoa_wb_rwd_YawTorque',...
'U_wb_fwd_YawTorque',...
'U_wb_rwd_YawTorque',...
...
't_ds_fwd_YawTorque',...
't_ds_rwd_YawTorque',...
'stroke_ds_fwd_YawTorque',...
'stroke_ds_rwd_YawTorque',...
'pitch_ds_fwd_YawTorque',...
'pitch_ds_rwd_YawTorque',...
'dev_ds_fwd_YawTorque',...
'dev_ds_rwd_YawTorque',...
'aoa_ds_fwd_YawTorque',...
'aoa_ds_rwd_YawTorque',...
'U_ds_fwd_YawTorque',...
'U_ds_rwd_YawTorque',...
...
't_us_fwd_YawTorque',...
't_us_rwd_YawTorque',...
'stroke_us_fwd_YawTorque',...
'stroke_us_rwd_YawTorque',...
'pitch_us_fwd_YawTorque',...
'pitch_us_rwd_YawTorque',...
'dev_us_fwd_YawTorque',...
'dev_us_rwd_YawTorque',...
'aoa_us_fwd_YawTorque',...
'aoa_us_rwd_YawTorque',...
'U_us_fwd_YawTorque',...
'U_us_rwd_YawTorque',...
...
't_wb_YawTorque_bins',...
'stroke_wb_fwd_YawTorque_bins',...
'stroke_wb_rwd_YawTorque_bins',...
'pitch_wb_fwd_YawTorque_bins',...
'pitch_wb_rwd_YawTorque_bins',...
'dev_wb_fwd_YawTorque_bins',...
'dev_wb_rwd_YawTorque_bins',...
'aoa_wb_fwd_YawTorque_bins',...
'aoa_wb_rwd_YawTorque_bins',...
'U_wb_fwd_YawTorque_bins',...
'U_wb_rwd_YawTorque_bins',...
'Dstroke_wb_YawTorque_bins',...
'Dpitch_wb_YawTorque_bins',...
'Ddev_wb_YawTorque_bins',...
'Daoa_wb_YawTorque_bins',...
'DU_wb_YawTorque_bins',...
...
't_ds_YawTorque_bins',...
'stroke_ds_fwd_YawTorque_bins',...
'stroke_ds_rwd_YawTorque_bins',...
'pitch_ds_fwd_YawTorque_bins',...
'pitch_ds_rwd_YawTorque_bins',...
'dev_ds_fwd_YawTorque_bins',...
'dev_ds_rwd_YawTorque_bins',...
'aoa_ds_fwd_YawTorque_bins',...
'aoa_ds_rwd_YawTorque_bins',...
'U_ds_fwd_YawTorque_bins',...
'U_ds_rwd_YawTorque_bins',...
'Dstroke_ds_YawTorque_bins',...
'Dpitch_ds_YawTorque_bins',...
'Ddev_ds_YawTorque_bins',...
'Daoa_ds_YawTorque_bins',...
'DU_ds_YawTorque_bins',...
...
't_us_YawTorque_bins',...
'stroke_us_fwd_YawTorque_bins',...
'stroke_us_rwd_YawTorque_bins',...
'pitch_us_fwd_YawTorque_bins',...
'pitch_us_rwd_YawTorque_bins',...
'dev_us_fwd_YawTorque_bins',...
'dev_us_rwd_YawTorque_bins',...
'aoa_us_fwd_YawTorque_bins',...
'aoa_us_rwd_YawTorque_bins',...
'U_us_fwd_YawTorque_bins',...
'U_us_rwd_YawTorque_bins',...
'Dstroke_us_YawTorque_bins',...
'Dpitch_us_YawTorque_bins',...
'Ddev_us_YawTorque_bins',...
'Daoa_us_YawTorque_bins',...
'DU_us_YawTorque_bins');

end
