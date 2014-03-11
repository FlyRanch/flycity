% calc wb MOD Roll Torque

n_now=0;

%% variables
% wbMOD
freqMOD_wb_RollTorque = [];

strokeMOD_wb_up_RollTorque_bins = [];
strokeMOD_wb_down_RollTorque_bins = [];
strokeMOD_ds_up_RollTorque_bins = [];
strokeMOD_ds_down_RollTorque_bins = [];
strokeMOD_us_up_RollTorque_bins = [];
strokeMOD_us_down_RollTorque_bins = [];

pitchMOD_wb_up_RollTorque_bins = [];
pitchMOD_wb_down_RollTorque_bins = [];
pitchMOD_ds_up_RollTorque_bins = [];
pitchMOD_ds_down_RollTorque_bins = [];
pitchMOD_us_up_RollTorque_bins = [];
pitchMOD_us_down_RollTorque_bins = [];

devMOD_wb_up_RollTorque_bins = [];
devMOD_wb_down_RollTorque_bins = [];
devMOD_ds_up_RollTorque_bins = [];
devMOD_ds_down_RollTorque_bins = [];
devMOD_us_up_RollTorque_bins = [];
devMOD_us_down_RollTorque_bins = [];

DstrokeMOD_wb_RollTorque_bins = [];
DstrokeMOD_ds_RollTorque_bins = [];
DstrokeMOD_us_RollTorque_bins = [];

DpitchMOD_wb_RollTorque_bins = [];
DpitchMOD_ds_RollTorque_bins = [];
DpitchMOD_us_RollTorque_bins = [];

DdevMOD_wb_RollTorque_bins = [];
DdevMOD_ds_RollTorque_bins = [];
DdevMOD_us_RollTorque_bins = [];

% seq&wb
seq_nr_RollTorque = [];
wb_nr_RollTorque = [];

% body kin
V_RollTorque = [];
pitch_global_RollTorque = [];
F_RollTorque = [];

Mroll_accel_RollTorque = [];
Mpitch_accel_RollTorque = [];
Myaw_accel_RollTorque = [];
M_L_accel_RollTorque = [];
M_R_accel_RollTorque = [];

Mroll_damp_RollTorque = [];
Mpitch_damp_RollTorque = [];
Myaw_damp_RollTorque = [];
M_L_damp_RollTorque = [];
M_R_damp_RollTorque = [];

Mroll_RollTorque = [];
Mpitch_RollTorque = [];
Myaw_RollTorque = [];
M_L_RollTorque = [];
M_R_RollTorque = [];

rollaccel_RollTorque = [];
pitchaccel_RollTorque = [];
yawaccel_RollTorque = [];
rotRaccel_RollTorque = [];

rollvel_RollTorque = [];
pitchvel_RollTorque = [];
yawvel_RollTorque = [];
rotRvel_RollTorque = [];

% wingbeat kin
dt_wb_RollTorque = [];
dt_ds_RollTorque = [];
dt_us_RollTorque = [];
f_wb_RollTorque = [];
Rds_RollTorque = [];

t_wb_up_RollTorque = [];
t_wb_down_RollTorque = [];
stroke_wb_up_RollTorque = [];
stroke_wb_down_RollTorque = [];
pitch_wb_up_RollTorque = [];
pitch_wb_down_RollTorque = [];
dev_wb_up_RollTorque = [];
dev_wb_down_RollTorque = [];
aoa_wb_up_RollTorque = [];
aoa_wb_down_RollTorque = [];
U_wb_up_RollTorque = [];
U_wb_down_RollTorque = [];
Dstroke_wb_RollTorque = [];
Dpitch_wb_RollTorque = [];
Ddev_wb_RollTorque = [];
Daoa_wb_RollTorque = [];
DU_wb_RollTorque = [];

t_ds_up_RollTorque = [];
t_ds_down_RollTorque = [];
stroke_ds_up_RollTorque = [];
stroke_ds_down_RollTorque = [];
pitch_ds_up_RollTorque = [];
pitch_ds_down_RollTorque = [];
dev_ds_up_RollTorque = [];
dev_ds_down_RollTorque = [];
aoa_ds_up_RollTorque = [];
aoa_ds_down_RollTorque = [];
U_ds_up_RollTorque = [];
U_ds_down_RollTorque = [];
Dstroke_ds_RollTorque = [];
Dpitch_ds_RollTorque = [];
Ddev_ds_RollTorque = [];
Daoa_ds_RollTorque = [];
DU_ds_RollTorque = [];

t_us_up_RollTorque = [];
t_us_down_RollTorque = [];
stroke_us_up_RollTorque = [];
stroke_us_down_RollTorque = [];
pitch_us_up_RollTorque = [];
pitch_us_down_RollTorque = [];
dev_us_up_RollTorque = [];
dev_us_down_RollTorque = [];
aoa_us_up_RollTorque = [];
aoa_us_down_RollTorque = [];
U_us_up_RollTorque = [];
U_us_down_RollTorque = [];
Dstroke_us_RollTorque = [];
Dpitch_us_RollTorque = [];
Ddev_us_RollTorque = [];
Daoa_us_RollTorque = [];
DU_us_RollTorque = [];

t_wb_RollTorque_bins = [];
stroke_wb_up_RollTorque_bins = [];
stroke_wb_down_RollTorque_bins = [];
pitch_wb_up_RollTorque_bins = [];
pitch_wb_down_RollTorque_bins = [];
dev_wb_up_RollTorque_bins = [];
dev_wb_down_RollTorque_bins = [];
aoa_wb_up_RollTorque_bins = [];
aoa_wb_down_RollTorque_bins = [];
U_wb_up_RollTorque_bins = [];
U_wb_down_RollTorque_bins = [];
Dstroke_wb_RollTorque_bins = [];
Dpitch_wb_RollTorque_bins = [];
Ddev_wb_RollTorque_bins = [];
Daoa_wb_RollTorque_bins = [];
DU_wb_RollTorque_bins = [];

t_ds_RollTorque_bins = [];
stroke_ds_up_RollTorque_bins = [];
stroke_ds_down_RollTorque_bins = [];
pitch_ds_up_RollTorque_bins = [];
pitch_ds_down_RollTorque_bins = [];
dev_ds_up_RollTorque_bins = [];
dev_ds_down_RollTorque_bins = [];
aoa_ds_up_RollTorque_bins = [];
aoa_ds_down_RollTorque_bins = [];
U_ds_up_RollTorque_bins = [];
U_ds_down_RollTorque_bins = [];
Dstroke_ds_RollTorque_bins = [];
Dpitch_ds_RollTorque_bins = [];
Ddev_ds_RollTorque_bins = [];
Daoa_ds_RollTorque_bins = [];
DU_ds_RollTorque_bins = [];

t_us_RollTorque_bins = [];
stroke_us_up_RollTorque_bins = [];
stroke_us_down_RollTorque_bins = [];
pitch_us_up_RollTorque_bins = [];
pitch_us_down_RollTorque_bins = [];
dev_us_up_RollTorque_bins = [];
dev_us_down_RollTorque_bins = [];
aoa_us_up_RollTorque_bins = [];
aoa_us_down_RollTorque_bins = [];
U_us_up_RollTorque_bins = [];
U_us_down_RollTorque_bins = [];
Dstroke_us_RollTorque_bins = [];
Dpitch_us_RollTorque_bins = [];
Ddev_us_RollTorque_bins = [];
Daoa_us_RollTorque_bins = [];
DU_us_RollTorque_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if Mroll_mean_wb(wb) > Mroll_limit_mod
%         if roll_dot_dot_mean_wb(wb) > rollaccel_limit_mod

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
            rotRaccel_mean_wb_now = rot_dot_dot_R_mean_wb(wb);

            rollvel_mean_wb_now = roll_dot_mean_wb(wb);
            pitchvel_mean_wb_now = pitch_dot_mean_wb(wb);
            yawvel_mean_wb_now = yaw_dot_mean_wb(wb);
            rotRvel_mean_wb_now = rot_dot_R_mean_wb(wb);

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

            seq_nr_RollTorque(n_now,1) = seq_nr_now;
            wb_nr_RollTorque(n_now,1) = wb_nr_now;
            
            % body kin
            V_RollTorque(n_now,1) = V_now;
            pitch_global_RollTorque(n_now,1) = pitch_global_now;
            F_RollTorque(n_now,1) = F_mean_wb_now;

            Mroll_accel_RollTorque(n_now,1) = Mroll_mean_wb_accel_now;
            Mpitch_accel_RollTorque(n_now,1) = Mpitch_mean_wb_accel_now;
            Myaw_accel_RollTorque(n_now,1) = Myaw_mean_wb_accel_now;
            M_L_accel_RollTorque(n_now,1) = M_L_mean_wb_accel_now;
            M_R_accel_RollTorque(n_now,1) = M_R_mean_wb_accel_now;

            Mroll_damp_RollTorque(n_now,1) = Mroll_mean_wb_damp_now;
            Mpitch_damp_RollTorque(n_now,1) = Mpitch_mean_wb_damp_now;
            Myaw_damp_RollTorque(n_now,1) = Myaw_mean_wb_damp_now;
            M_L_damp_RollTorque(n_now,1) = M_L_mean_wb_damp_now;
            M_R_damp_RollTorque(n_now,1) = M_R_mean_wb_damp_now;

            Mroll_RollTorque(n_now,1) = Mroll_mean_wb_now;
            Mpitch_RollTorque(n_now,1) = Mpitch_mean_wb_now;
            Myaw_RollTorque(n_now,1) = Myaw_mean_wb_now;
            M_L_RollTorque(n_now,1) = M_L_mean_wb_now;
            M_R_RollTorque(n_now,1) = M_R_mean_wb_now;

            rollaccel_RollTorque(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_RollTorque(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_RollTorque(n_now,1) = yawaccel_mean_wb_now;
            rotRaccel_RollTorque(n_now,1) = rotRaccel_mean_wb_now;

            rollvel_RollTorque(n_now,1) = rollvel_mean_wb_now;
            pitchvel_RollTorque(n_now,1) = pitchvel_mean_wb_now;
            yawvel_RollTorque(n_now,1) = yawvel_mean_wb_now;
            rotRvel_RollTorque(n_now,1) = rotRvel_mean_wb_now;

            % wingbeat kin
            dt_ds_RollTorque(n_now,1) = nanmean(dt_ds_now);
            dt_us_RollTorque(n_now,1) = nanmean(dt_us_now);
            f_wb_RollTorque(n_now,1) = nanmean(f_wb_now);
            Rds_RollTorque(n_now,1) = nanmean(Rds_now);

            t_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = t_wb_up_now;
            t_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = t_wb_down_now;
            stroke_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = stroke_wb_up_now;
            stroke_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = stroke_wb_down_now;
            pitch_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = pitch_wb_up_now;
            pitch_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = pitch_wb_down_now;
            dev_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = dev_wb_up_now;
            dev_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = dev_wb_down_now;
            aoa_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = aoa_wb_up_now;
            aoa_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = aoa_wb_down_now;
            U_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = U_wb_up_now;
            U_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = U_wb_down_now;
            
            t_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = t_ds_up_now;
            t_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = t_ds_down_now;
            stroke_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = stroke_ds_up_now;
            stroke_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = stroke_ds_down_now;
            pitch_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = pitch_ds_up_now;
            pitch_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = pitch_ds_down_now;
            dev_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = dev_ds_up_now;
            dev_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = dev_ds_down_now;
            aoa_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = aoa_ds_up_now;
            aoa_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = aoa_ds_down_now;
            U_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = U_ds_up_now;
            U_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = U_ds_down_now;
            
            t_us_up_RollTorque(1:length(t_us_up_now),n_now) = t_us_up_now;
            t_us_down_RollTorque(1:length(t_us_down_now),n_now) = t_us_down_now;
            stroke_us_up_RollTorque(1:length(t_us_up_now),n_now) = stroke_us_up_now;
            stroke_us_down_RollTorque(1:length(t_us_down_now),n_now) = stroke_us_down_now;
            pitch_us_up_RollTorque(1:length(t_us_up_now),n_now) = pitch_us_up_now;
            pitch_us_down_RollTorque(1:length(t_us_down_now),n_now) = pitch_us_down_now;
            dev_us_up_RollTorque(1:length(t_us_up_now),n_now) = dev_us_up_now;
            dev_us_down_RollTorque(1:length(t_us_down_now),n_now) = dev_us_down_now;
            aoa_us_up_RollTorque(1:length(t_us_up_now),n_now) = aoa_us_up_now;
            aoa_us_down_RollTorque(1:length(t_us_down_now),n_now) = aoa_us_down_now;
            U_us_up_RollTorque(1:length(t_us_up_now),n_now) = U_us_up_now;
            U_us_down_RollTorque(1:length(t_us_down_now),n_now) = U_us_down_now;
            
            % store interp binned data separate rows
            t_wb_RollTorque_bins(:,n_now) = t_wb_bin;
            stroke_wb_up_RollTorque_bins(:,n_now) = stroke_wb_up_interp;
            stroke_wb_down_RollTorque_bins(:,n_now) = stroke_wb_down_interp;
            pitch_wb_up_RollTorque_bins(:,n_now) = pitch_wb_up_interp;
            pitch_wb_down_RollTorque_bins(:,n_now) = pitch_wb_down_interp;
            dev_wb_up_RollTorque_bins(:,n_now) = dev_wb_up_interp;
            dev_wb_down_RollTorque_bins(:,n_now) = dev_wb_down_interp;
            aoa_wb_up_RollTorque_bins(:,n_now) = aoa_wb_up_interp;
            aoa_wb_down_RollTorque_bins(:,n_now) = aoa_wb_down_interp;
            U_wb_up_RollTorque_bins(:,n_now) = U_wb_up_interp;
            U_wb_down_RollTorque_bins(:,n_now) = U_wb_down_interp;
            Dstroke_wb_RollTorque_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_RollTorque_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_RollTorque_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_RollTorque_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_RollTorque_bins(:,n_now) = DU_wb_interp;
            
            t_ds_RollTorque_bins(:,n_now) = t_ds_bin;
            stroke_ds_up_RollTorque_bins(:,n_now) = stroke_ds_up_interp;
            stroke_ds_down_RollTorque_bins(:,n_now) = stroke_ds_down_interp;
            pitch_ds_up_RollTorque_bins(:,n_now) = pitch_ds_up_interp;
            pitch_ds_down_RollTorque_bins(:,n_now) = pitch_ds_down_interp;
            dev_ds_up_RollTorque_bins(:,n_now) = dev_ds_up_interp;
            dev_ds_down_RollTorque_bins(:,n_now) = dev_ds_down_interp;
            aoa_ds_up_RollTorque_bins(:,n_now) = aoa_ds_up_interp;
            aoa_ds_down_RollTorque_bins(:,n_now) = aoa_ds_down_interp;
            U_ds_up_RollTorque_bins(:,n_now) = U_ds_up_interp;
            U_ds_down_RollTorque_bins(:,n_now) = U_ds_down_interp;
            Dstroke_ds_RollTorque_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_RollTorque_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_RollTorque_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_RollTorque_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_RollTorque_bins(:,n_now) = DU_ds_interp;
            
            t_us_RollTorque_bins(:,n_now) = t_us_bin;
            stroke_us_up_RollTorque_bins(:,n_now) = stroke_us_up_interp;
            stroke_us_down_RollTorque_bins(:,n_now) = stroke_us_down_interp;
            pitch_us_up_RollTorque_bins(:,n_now) = pitch_us_up_interp;
            pitch_us_down_RollTorque_bins(:,n_now) = pitch_us_down_interp;
            dev_us_up_RollTorque_bins(:,n_now) = dev_us_up_interp;
            dev_us_down_RollTorque_bins(:,n_now) = dev_us_down_interp;
            aoa_us_up_RollTorque_bins(:,n_now) = aoa_us_up_interp;
            aoa_us_down_RollTorque_bins(:,n_now) = aoa_us_down_interp;
            U_us_up_RollTorque_bins(:,n_now) = U_us_up_interp;
            U_us_down_RollTorque_bins(:,n_now) = U_us_down_interp;
            Dstroke_us_RollTorque_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_RollTorque_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_RollTorque_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_RollTorque_bins(:,n_now) = Daoa_us_interp;
            DU_us_RollTorque_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_RollTorque(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_up_RollTorque_bins(:,n_now) = strokeMOD_wb_up_now;
            strokeMOD_wb_down_RollTorque_bins(:,n_now) = strokeMOD_wb_down_now;
            strokeMOD_ds_up_RollTorque_bins(:,n_now) = strokeMOD_ds_up_now;
            strokeMOD_ds_down_RollTorque_bins(:,n_now) = strokeMOD_ds_down_now;
            strokeMOD_us_up_RollTorque_bins(:,n_now) = strokeMOD_us_up_now;
            strokeMOD_us_down_RollTorque_bins(:,n_now) = strokeMOD_us_down_now;

            pitchMOD_wb_up_RollTorque_bins(:,n_now) = pitchMOD_wb_up_now;
            pitchMOD_wb_down_RollTorque_bins(:,n_now) = pitchMOD_wb_down_now;
            pitchMOD_ds_up_RollTorque_bins(:,n_now) = pitchMOD_ds_up_now;
            pitchMOD_ds_down_RollTorque_bins(:,n_now) = pitchMOD_ds_down_now;
            pitchMOD_us_up_RollTorque_bins(:,n_now) = pitchMOD_us_up_now;
            pitchMOD_us_down_RollTorque_bins(:,n_now) = pitchMOD_us_down_now;

            devMOD_wb_up_RollTorque_bins(:,n_now) = devMOD_wb_up_now;
            devMOD_wb_down_RollTorque_bins(:,n_now) = devMOD_wb_down_now;
            devMOD_ds_up_RollTorque_bins(:,n_now) = devMOD_ds_up_now;
            devMOD_ds_down_RollTorque_bins(:,n_now) = devMOD_ds_down_now;
            devMOD_us_up_RollTorque_bins(:,n_now) = devMOD_us_up_now;
            devMOD_us_down_RollTorque_bins(:,n_now) = devMOD_us_down_now;
            
            DstrokeMOD_wb_RollTorque_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_RollTorque_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_RollTorque_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_RollTorque_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_RollTorque_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_RollTorque_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_RollTorque_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_RollTorque_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_RollTorque_bins(:,n_now) = DdevMOD_us_now;
            
        elseif Mroll_mean_wb(wb) < -Mroll_limit_mod
%         elseif roll_dot_dot_mean_wb(wb) < -rollaccel_limit_mod
            
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

            seq_nr_RollTorque(n_now,1) = seq_nr_now;
            wb_nr_RollTorque(n_now,1) = wb_nr_now;
            
            % body kin
            V_RollTorque(n_now,1) = V_now;
            pitch_global_RollTorque(n_now,1) = pitch_global_now;
            F_RollTorque(n_now,1) = F_mean_wb_now;

            Mroll_accel_RollTorque(n_now,1) = Mroll_mean_wb_accel_now;
            Mpitch_accel_RollTorque(n_now,1) = Mpitch_mean_wb_accel_now;
            Myaw_accel_RollTorque(n_now,1) = Myaw_mean_wb_accel_now;
            M_L_accel_RollTorque(n_now,1) = M_L_mean_wb_accel_now;
            M_R_accel_RollTorque(n_now,1) = M_R_mean_wb_accel_now;

            Mroll_damp_RollTorque(n_now,1) = Mroll_mean_wb_damp_now;
            Mpitch_damp_RollTorque(n_now,1) = Mpitch_mean_wb_damp_now;
            Myaw_damp_RollTorque(n_now,1) = Myaw_mean_wb_damp_now;
            M_L_damp_RollTorque(n_now,1) = M_L_mean_wb_damp_now;
            M_R_damp_RollTorque(n_now,1) = M_R_mean_wb_damp_now;

            Mroll_RollTorque(n_now,1) = Mroll_mean_wb_now;
            Mpitch_RollTorque(n_now,1) = Mpitch_mean_wb_now;
            Myaw_RollTorque(n_now,1) = Myaw_mean_wb_now;
            M_L_RollTorque(n_now,1) = M_L_mean_wb_now;
            M_R_RollTorque(n_now,1) = M_R_mean_wb_now;

            rollaccel_RollTorque(n_now,1) = rollaccel_mean_wb_now;
            pitchaccel_RollTorque(n_now,1) = pitchaccel_mean_wb_now;
            yawaccel_RollTorque(n_now,1) = yawaccel_mean_wb_now;

            rollvel_RollTorque(n_now,1) = rollvel_mean_wb_now;
            pitchvel_RollTorque(n_now,1) = pitchvel_mean_wb_now;
            yawvel_RollTorque(n_now,1) = yawvel_mean_wb_now;
            
            % wingbeat kin
            dt_ds_RollTorque(n_now,1) = nanmean(dt_ds_now);
            dt_us_RollTorque(n_now,1) = nanmean(dt_us_now);
            f_wb_RollTorque(n_now,1) = nanmean(f_wb_now);
            Rds_RollTorque(n_now,1) = nanmean(Rds_now);

            t_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = t_wb_up_now;
            t_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = t_wb_down_now;
            stroke_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = stroke_wb_up_now;
            stroke_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = stroke_wb_down_now;
            pitch_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = pitch_wb_up_now;
            pitch_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = pitch_wb_down_now;
            dev_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = dev_wb_up_now;
            dev_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = dev_wb_down_now;
            aoa_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = aoa_wb_up_now;
            aoa_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = aoa_wb_down_now;
            U_wb_up_RollTorque(1:length(t_wb_up_now),n_now) = U_wb_up_now;
            U_wb_down_RollTorque(1:length(t_wb_down_now),n_now) = U_wb_down_now;
            
            t_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = t_ds_up_now;
            t_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = t_ds_down_now;
            stroke_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = stroke_ds_up_now;
            stroke_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = stroke_ds_down_now;
            pitch_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = pitch_ds_up_now;
            pitch_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = pitch_ds_down_now;
            dev_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = dev_ds_up_now;
            dev_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = dev_ds_down_now;
            aoa_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = aoa_ds_up_now;
            aoa_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = aoa_ds_down_now;
            U_ds_up_RollTorque(1:length(t_ds_up_now),n_now) = U_ds_up_now;
            U_ds_down_RollTorque(1:length(t_ds_down_now),n_now) = U_ds_down_now;
            
            t_us_up_RollTorque(1:length(t_us_up_now),n_now) = t_us_up_now;
            t_us_down_RollTorque(1:length(t_us_down_now),n_now) = t_us_down_now;
            stroke_us_up_RollTorque(1:length(t_us_up_now),n_now) = stroke_us_up_now;
            stroke_us_down_RollTorque(1:length(t_us_down_now),n_now) = stroke_us_down_now;
            pitch_us_up_RollTorque(1:length(t_us_up_now),n_now) = pitch_us_up_now;
            pitch_us_down_RollTorque(1:length(t_us_down_now),n_now) = pitch_us_down_now;
            dev_us_up_RollTorque(1:length(t_us_up_now),n_now) = dev_us_up_now;
            dev_us_down_RollTorque(1:length(t_us_down_now),n_now) = dev_us_down_now;
            aoa_us_up_RollTorque(1:length(t_us_up_now),n_now) = aoa_us_up_now;
            aoa_us_down_RollTorque(1:length(t_us_down_now),n_now) = aoa_us_down_now;
            U_us_up_RollTorque(1:length(t_us_up_now),n_now) = U_us_up_now;
            U_us_down_RollTorque(1:length(t_us_down_now),n_now) = U_us_down_now;
            
            % store interp binned data separate rows
            t_wb_RollTorque_bins(:,n_now) = t_wb_bin;
            stroke_wb_up_RollTorque_bins(:,n_now) = stroke_wb_up_interp;
            stroke_wb_down_RollTorque_bins(:,n_now) = stroke_wb_down_interp;
            pitch_wb_up_RollTorque_bins(:,n_now) = pitch_wb_up_interp;
            pitch_wb_down_RollTorque_bins(:,n_now) = pitch_wb_down_interp;
            dev_wb_up_RollTorque_bins(:,n_now) = dev_wb_up_interp;
            dev_wb_down_RollTorque_bins(:,n_now) = dev_wb_down_interp;
            aoa_wb_up_RollTorque_bins(:,n_now) = aoa_wb_up_interp;
            aoa_wb_down_RollTorque_bins(:,n_now) = aoa_wb_down_interp;
            U_wb_up_RollTorque_bins(:,n_now) = U_wb_up_interp;
            U_wb_down_RollTorque_bins(:,n_now) = U_wb_down_interp;
            Dstroke_wb_RollTorque_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_RollTorque_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_RollTorque_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_RollTorque_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_RollTorque_bins(:,n_now) = DU_wb_interp;
            
            t_ds_RollTorque_bins(:,n_now) = t_ds_bin;
            stroke_ds_up_RollTorque_bins(:,n_now) = stroke_ds_up_interp;
            stroke_ds_down_RollTorque_bins(:,n_now) = stroke_ds_down_interp;
            pitch_ds_up_RollTorque_bins(:,n_now) = pitch_ds_up_interp;
            pitch_ds_down_RollTorque_bins(:,n_now) = pitch_ds_down_interp;
            dev_ds_up_RollTorque_bins(:,n_now) = dev_ds_up_interp;
            dev_ds_down_RollTorque_bins(:,n_now) = dev_ds_down_interp;
            aoa_ds_up_RollTorque_bins(:,n_now) = aoa_ds_up_interp;
            aoa_ds_down_RollTorque_bins(:,n_now) = aoa_ds_down_interp;
            U_ds_up_RollTorque_bins(:,n_now) = U_ds_up_interp;
            U_ds_down_RollTorque_bins(:,n_now) = U_ds_down_interp;
            Dstroke_ds_RollTorque_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_RollTorque_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_RollTorque_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_RollTorque_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_RollTorque_bins(:,n_now) = DU_ds_interp;
            
            t_us_RollTorque_bins(:,n_now) = t_us_bin;
            stroke_us_up_RollTorque_bins(:,n_now) = stroke_us_up_interp;
            stroke_us_down_RollTorque_bins(:,n_now) = stroke_us_down_interp;
            pitch_us_up_RollTorque_bins(:,n_now) = pitch_us_up_interp;
            pitch_us_down_RollTorque_bins(:,n_now) = pitch_us_down_interp;
            dev_us_up_RollTorque_bins(:,n_now) = dev_us_up_interp;
            dev_us_down_RollTorque_bins(:,n_now) = dev_us_down_interp;
            aoa_us_up_RollTorque_bins(:,n_now) = aoa_us_up_interp;
            aoa_us_down_RollTorque_bins(:,n_now) = aoa_us_down_interp;
            U_us_up_RollTorque_bins(:,n_now) = U_us_up_interp;
            U_us_down_RollTorque_bins(:,n_now) = U_us_down_interp;
            Dstroke_us_RollTorque_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_RollTorque_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_RollTorque_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_RollTorque_bins(:,n_now) = Daoa_us_interp;
            DU_us_RollTorque_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_RollTorque(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_up_RollTorque_bins(:,n_now) = strokeMOD_wb_up_now;
            strokeMOD_wb_down_RollTorque_bins(:,n_now) = strokeMOD_wb_down_now;
            strokeMOD_ds_up_RollTorque_bins(:,n_now) = strokeMOD_ds_up_now;
            strokeMOD_ds_down_RollTorque_bins(:,n_now) = strokeMOD_ds_down_now;
            strokeMOD_us_up_RollTorque_bins(:,n_now) = strokeMOD_us_up_now;
            strokeMOD_us_down_RollTorque_bins(:,n_now) = strokeMOD_us_down_now;

            pitchMOD_wb_up_RollTorque_bins(:,n_now) = pitchMOD_wb_up_now;
            pitchMOD_wb_down_RollTorque_bins(:,n_now) = pitchMOD_wb_down_now;
            pitchMOD_ds_up_RollTorque_bins(:,n_now) = pitchMOD_ds_up_now;
            pitchMOD_ds_down_RollTorque_bins(:,n_now) = pitchMOD_ds_down_now;
            pitchMOD_us_up_RollTorque_bins(:,n_now) = pitchMOD_us_up_now;
            pitchMOD_us_down_RollTorque_bins(:,n_now) = pitchMOD_us_down_now;

            devMOD_wb_up_RollTorque_bins(:,n_now) = devMOD_wb_up_now;
            devMOD_wb_down_RollTorque_bins(:,n_now) = devMOD_wb_down_now;
            devMOD_ds_up_RollTorque_bins(:,n_now) = devMOD_ds_up_now;
            devMOD_ds_down_RollTorque_bins(:,n_now) = devMOD_ds_down_now;
            devMOD_us_up_RollTorque_bins(:,n_now) = devMOD_us_up_now;
            devMOD_us_down_RollTorque_bins(:,n_now) = devMOD_us_down_now;
            
            DstrokeMOD_wb_RollTorque_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_RollTorque_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_RollTorque_bins(:,n_now) = DstrokeMOD_us_now;
            
            DpitchMOD_wb_RollTorque_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_RollTorque_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_RollTorque_bins(:,n_now) = DpitchMOD_us_now;
            
            DdevMOD_wb_RollTorque_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_RollTorque_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_RollTorque_bins(:,n_now) = DdevMOD_us_now;
        end
end

% mean & 95%CI
V_RollTorque_meanCIstd = [nanmean(V_RollTorque) 1.96*nanstd(V_RollTorque)/sqrt(length(V_RollTorque)) nanstd(V_RollTorque)];
pitch_global_RollTorque_meanCIstd = [nanmean(pitch_global_RollTorque) 1.96*nanstd(pitch_global_RollTorque)/sqrt(length(pitch_global_RollTorque)) nanstd(pitch_global_RollTorque)];

dt_ds_RollTorque_meanCIstd = [nanmean(dt_ds_RollTorque) 1.96*nanstd(dt_ds_RollTorque)/sqrt(length(dt_ds_RollTorque)) nanstd(dt_ds_RollTorque)];
dt_us_RollTorque_meanCIstd = [nanmean(dt_us_RollTorque) 1.96*nanstd(dt_us_RollTorque)/sqrt(length(dt_us_RollTorque)) nanstd(dt_us_RollTorque)];
f_wb_RollTorque_meanCIstd = [nanmean(f_wb_RollTorque) 1.96*nanstd(f_wb_RollTorque)/sqrt(length(f_wb_RollTorque)) nanstd(f_wb_RollTorque)];
Rds_RollTorque_meanCIstd = [nanmean(Rds_RollTorque) 1.96*nanstd(Rds_RollTorque)/sqrt(length(Rds_RollTorque)) nanstd(Rds_RollTorque)];

calc_WBfunc_RollTorque_circmeanCIstd

% wbMOD means & 95%CI
freqMOD_wb_RollTorque_meanCIstd = [nanmean(freqMOD_wb_RollTorque) 1.96*nanstd(freqMOD_wb_RollTorque)/sqrt(length(freqMOD_wb_RollTorque)) nanstd(freqMOD_wb_RollTorque)];
calc_WBmod_RollTorque_circmeanCIstd

%% WBfits
    t_loc = t_wb_RollTorque_bins(:,1);
    Rds_loc = Rds_RollTorque_meanCIstd(1);

    plotting = 0;
%         plotting = 1;

    %% fit for upwards moving wing (bins)
    stroke_loc = stroke_wb_up_RollTorque_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_up_RollTorque_bins_meanCIstd(:,1);
    dev_loc = dev_wb_up_RollTorque_bins_meanCIstd(:,1);
    
    strokeMOD_loc = strokeMOD_wb_up_RollTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_up_RollTorque_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_up_RollTorque_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_up_RollTorque_fit_binmean, stroke_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_up_RollTorque_fit_binmean, pitch_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_up_RollTorque_fit_binmean, dev_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_up_RollTorque_fit_binmean, strokeMOD_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_up_RollTorque_fit_binmean, pitchMOD_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_up_RollTorque_fit_binmean, devMOD_up_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_up_RollTorque_fourier_fit_binmean, stroke_up_RollTorque_fourier_gof_binmean,stroke_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_up_RollTorque_fourier_fit_binmean, pitch_up_RollTorque_fourier_gof_binmean,pitch_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_up_RollTorque_fourier_fit_binmean, dev_up_RollTorque_fourier_gof_binmean,dev_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_up_RollTorque_fourier_fit_binmean, strokeMOD_up_RollTorque_fourier_gof_binmean,strokeMOD_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_up_RollTorque_fourier_fit_binmean, pitchMOD_up_RollTorque_fourier_gof_binmean,pitchMOD_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_up_RollTorque_fourier_fit_binmean, devMOD_up_RollTorque_fourier_gof_binmean,devMOD_up_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% downwards
    stroke_loc = stroke_wb_down_RollTorque_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_down_RollTorque_bins_meanCIstd(:,1);
    dev_loc = dev_wb_down_RollTorque_bins_meanCIstd(:,1);

    strokeMOD_loc = strokeMOD_wb_down_RollTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_down_RollTorque_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_down_RollTorque_bins_meanCIstd(:,1);

%% legendre polynomials
% wingbeats
    [stroke_down_RollTorque_fit_binmean, stroke_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_down_RollTorque_fit_binmean, pitch_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_down_RollTorque_fit_binmean, dev_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [strokeMOD_down_RollTorque_fit_binmean, strokeMOD_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_down_RollTorque_fit_binmean, pitchMOD_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_down_RollTorque_fit_binmean, devMOD_down_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [stroke_down_RollTorque_fourier_fit_binmean, stroke_down_RollTorque_fourier_gof_binmean,stroke_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [pitch_down_RollTorque_fourier_fit_binmean, pitch_down_RollTorque_fourier_gof_binmean,pitch_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [dev_down_RollTorque_fourier_fit_binmean, dev_down_RollTorque_fourier_gof_binmean,dev_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [strokeMOD_down_RollTorque_fourier_fit_binmean, strokeMOD_down_RollTorque_fourier_gof_binmean,strokeMOD_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [pitchMOD_down_RollTorque_fourier_fit_binmean, pitchMOD_down_RollTorque_fourier_gof_binmean,pitchMOD_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [devMOD_down_RollTorque_fourier_fit_binmean, devMOD_down_RollTorque_fourier_gof_binmean,devMOD_down_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    
%% upwards - downwards
    stroke_loc = Dstroke_wb_RollTorque_bins_meanCIstd(:,1);
    pitch_loc = Dpitch_wb_RollTorque_bins_meanCIstd(:,1);
    dev_loc = Ddev_wb_RollTorque_bins_meanCIstd(:,1);

    strokeMOD_loc = DstrokeMOD_wb_RollTorque_bins_meanCIstd(:,1);
    pitchMOD_loc = DpitchMOD_wb_RollTorque_bins_meanCIstd(:,1);
    devMOD_loc = DdevMOD_wb_RollTorque_bins_meanCIstd(:,1);
    
%% legendre polynomials
    % wingbeats
    [Dstroke_RollTorque_fit_binmean, Dstroke_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [Dpitch_RollTorque_fit_binmean, Dpitch_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [Ddev_RollTorque_fit_binmean, Ddev_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeatMODs
    [DstrokeMOD_RollTorque_fit_binmean, DstrokeMOD_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
    [DpitchMOD_RollTorque_fit_binmean, DpitchMOD_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
    [DdevMOD_RollTorque_fit_binmean, DdevMOD_RollTorque_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

%% fourier series
        % wingbeats
        [Dstroke_RollTorque_fourier_fit_binmean, Dstroke_RollTorque_fourier_gof_binmean,Dstroke_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
        [Dpitch_RollTorque_fourier_fit_binmean, Dpitch_RollTorque_fourier_gof_binmean,Dpitch_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
        [Ddev_RollTorque_fourier_fit_binmean, Ddev_RollTorque_fourier_gof_binmean,Ddev_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);

        % wingbeatMODs
        [DstrokeMOD_RollTorque_fourier_fit_binmean, DstrokeMOD_RollTorque_fourier_gof_binmean,DstrokeMOD_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plotting);
        [DpitchMOD_RollTorque_fourier_fit_binmean, DpitchMOD_RollTorque_fourier_gof_binmean,DpitchMOD_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plotting);
        [DdevMOD_RollTorque_fourier_fit_binmean, DdevMOD_RollTorque_fourier_gof_binmean,DdevMOD_RollTorque_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plotting);
    

    
    
    
    
    %% save RollTorque wb data
n_RollTorque = n_now   

if save_on == 1

save(['WBmod_torquebased_RollTorque_',num2str(n_RollTorque),'WBs.mat'],...
'n_RollTorque',...
'seq_nr_RollTorque',...
'wb_nr_RollTorque',...
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
'stroke_up_RollTorque_fit_binmean_periodic',...
'pitch_up_RollTorque_fit_binmean_periodic',...
'dev_up_RollTorque_fit_binmean_periodic',...
...
'stroke_down_RollTorque_fit_binmean_periodic',...
'pitch_down_RollTorque_fit_binmean_periodic',...
'dev_down_RollTorque_fit_binmean_periodic',...
...
'strokeMOD_up_RollTorque_fit_binmean_periodic',...
'pitchMOD_up_RollTorque_fit_binmean_periodic',...
'devMOD_up_RollTorque_fit_binmean_periodic',...
...
'strokeMOD_down_RollTorque_fit_binmean_periodic',...
'pitchMOD_down_RollTorque_fit_binmean_periodic',...
'devMOD_down_RollTorque_fit_binmean_periodic',...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_up_RollTorque_fourier_coeffs_binmean',...
'pitch_up_RollTorque_fourier_coeffs_binmean',...
'dev_up_RollTorque_fourier_coeffs_binmean',...
...
'stroke_up_RollTorque_fourier_gof_binmean',...
'pitch_up_RollTorque_fourier_gof_binmean',...
'dev_up_RollTorque_fourier_gof_binmean',...
...
'stroke_down_RollTorque_fourier_coeffs_binmean',...
'pitch_down_RollTorque_fourier_coeffs_binmean',...
'dev_down_RollTorque_fourier_coeffs_binmean',...
...
'stroke_down_RollTorque_fourier_gof_binmean',...
'pitch_down_RollTorque_fourier_gof_binmean',...
'dev_down_RollTorque_fourier_gof_binmean',...
...
'Dstroke_RollTorque_fourier_coeffs_binmean',...
'Dpitch_RollTorque_fourier_coeffs_binmean',...
'Ddev_RollTorque_fourier_coeffs_binmean',...
...
'Dstroke_RollTorque_fourier_gof_binmean',...
'Dpitch_RollTorque_fourier_gof_binmean',...
'Ddev_RollTorque_fourier_gof_binmean',...
...
'strokeMOD_up_RollTorque_fourier_coeffs_binmean',...
'pitchMOD_up_RollTorque_fourier_coeffs_binmean',...
'devMOD_up_RollTorque_fourier_coeffs_binmean',...
...
'strokeMOD_up_RollTorque_fourier_gof_binmean',...
'pitchMOD_up_RollTorque_fourier_gof_binmean',...
'devMOD_up_RollTorque_fourier_gof_binmean',...
...
'strokeMOD_down_RollTorque_fourier_coeffs_binmean',...
'pitchMOD_down_RollTorque_fourier_coeffs_binmean',...
'devMOD_down_RollTorque_fourier_coeffs_binmean',...
...
'strokeMOD_down_RollTorque_fourier_gof_binmean',...
'pitchMOD_down_RollTorque_fourier_gof_binmean',...
'devMOD_down_RollTorque_fourier_gof_binmean',...
...
'DstrokeMOD_RollTorque_fourier_coeffs_binmean',...
'DpitchMOD_RollTorque_fourier_coeffs_binmean',...
'DdevMOD_RollTorque_fourier_coeffs_binmean',...
...
'DstrokeMOD_RollTorque_fourier_gof_binmean',...
'DpitchMOD_RollTorque_fourier_gof_binmean',...
'DdevMOD_RollTorque_fourier_gof_binmean',...
...
'stroke_wb_up_RollTorque_bins_meanCIstd',...
'stroke_ds_up_RollTorque_bins_meanCIstd',...
'stroke_us_up_RollTorque_bins_meanCIstd',...
...
'pitch_wb_up_RollTorque_bins_meanCIstd',...
'pitch_ds_up_RollTorque_bins_meanCIstd',...
'pitch_us_up_RollTorque_bins_meanCIstd',...
...
'dev_wb_up_RollTorque_bins_meanCIstd',...
'dev_ds_up_RollTorque_bins_meanCIstd',...
'dev_us_up_RollTorque_bins_meanCIstd',...
...
'strokeMOD_wb_up_RollTorque_bins_meanCIstd',...
'strokeMOD_ds_up_RollTorque_bins_meanCIstd',...
'strokeMOD_us_up_RollTorque_bins_meanCIstd',...
...
'pitchMOD_wb_up_RollTorque_bins_meanCIstd',...
'pitchMOD_ds_up_RollTorque_bins_meanCIstd',...
'pitchMOD_us_up_RollTorque_bins_meanCIstd',...
...
'devMOD_wb_up_RollTorque_bins_meanCIstd',...
'devMOD_ds_up_RollTorque_bins_meanCIstd',...
'devMOD_us_up_RollTorque_bins_meanCIstd',...
...
'stroke_wb_down_RollTorque_bins_meanCIstd',...
'stroke_ds_down_RollTorque_bins_meanCIstd',...
'stroke_us_down_RollTorque_bins_meanCIstd',...
...
'pitch_wb_down_RollTorque_bins_meanCIstd',...
'pitch_ds_down_RollTorque_bins_meanCIstd',...
'pitch_us_down_RollTorque_bins_meanCIstd',...
...
'dev_wb_down_RollTorque_bins_meanCIstd',...
'dev_ds_down_RollTorque_bins_meanCIstd',...
'dev_us_down_RollTorque_bins_meanCIstd',...
...
'strokeMOD_wb_down_RollTorque_bins_meanCIstd',...
'strokeMOD_ds_down_RollTorque_bins_meanCIstd',...
'strokeMOD_us_down_RollTorque_bins_meanCIstd',...
...
'pitchMOD_wb_down_RollTorque_bins_meanCIstd',...
'pitchMOD_ds_down_RollTorque_bins_meanCIstd',...
'pitchMOD_us_down_RollTorque_bins_meanCIstd',...
...
'devMOD_wb_down_RollTorque_bins_meanCIstd',...
'devMOD_ds_down_RollTorque_bins_meanCIstd',...
'devMOD_us_down_RollTorque_bins_meanCIstd',...
...
'DstrokeMOD_wb_RollTorque_bins_meanCIstd',...
'DstrokeMOD_ds_RollTorque_bins_meanCIstd',...
'DstrokeMOD_us_RollTorque_bins_meanCIstd',...
...
'DpitchMOD_wb_RollTorque_bins_meanCIstd',...
'DpitchMOD_ds_RollTorque_bins_meanCIstd',...
'DpitchMOD_us_RollTorque_bins_meanCIstd',...
...
'DdevMOD_wb_RollTorque_bins_meanCIstd',...
'DdevMOD_ds_RollTorque_bins_meanCIstd',...
'DdevMOD_us_RollTorque_bins_meanCIstd',...
...
'V_RollTorque',...
'pitch_global_RollTorque',...
'dt_ds_RollTorque',...
'dt_us_RollTorque',...
'f_wb_RollTorque',...
'Rds_RollTorque',...
...
'V_RollTorque_meanCIstd',...
'pitch_global_RollTorque_meanCIstd',...
'F_RollTorque',...
...
'Mroll_accel_RollTorque',...
'Mpitch_accel_RollTorque',...
'Myaw_accel_RollTorque',...
'M_L_accel_RollTorque',...
'M_R_accel_RollTorque',...
...
'Mroll_damp_RollTorque',...
'Mpitch_damp_RollTorque',...
'Myaw_damp_RollTorque',...
'M_L_damp_RollTorque',...
'M_R_damp_RollTorque',...
...
'Mroll_RollTorque',...
'Mpitch_RollTorque',...
'Myaw_RollTorque',...
'M_L_RollTorque',...
'M_R_RollTorque',...
...
'rollaccel_RollTorque',...
'pitchaccel_RollTorque',...
'yawaccel_RollTorque',...
'rotRaccel_RollTorque',...
...
'rollvel_RollTorque',...
'pitchvel_RollTorque',...
'yawvel_RollTorque',...
'rotRvel_RollTorque',...
...
'dt_ds_RollTorque_meanCIstd',...
'dt_us_RollTorque_meanCIstd',...
'f_wb_RollTorque_meanCIstd',...
'Rds_RollTorque_meanCIstd',...
...
'freqMOD_wb_RollTorque_meanCIstd',...
'freqMOD_wb_RollTorque',...
...
'strokeMOD_wb_up_RollTorque_bins',...
'strokeMOD_wb_down_RollTorque_bins',...
'strokeMOD_ds_up_RollTorque_bins',...
'strokeMOD_ds_down_RollTorque_bins',...
'strokeMOD_us_up_RollTorque_bins',...
'strokeMOD_us_down_RollTorque_bins',...
...
'pitchMOD_wb_up_RollTorque_bins',...
'pitchMOD_wb_down_RollTorque_bins',...
'pitchMOD_ds_up_RollTorque_bins',...
'pitchMOD_ds_down_RollTorque_bins',...
'pitchMOD_us_up_RollTorque_bins',...
'pitchMOD_us_down_RollTorque_bins',...
...
'devMOD_wb_up_RollTorque_bins',...
'devMOD_wb_down_RollTorque_bins',...
'devMOD_ds_up_RollTorque_bins',...
'devMOD_ds_down_RollTorque_bins',...
'devMOD_us_up_RollTorque_bins',...
'devMOD_us_down_RollTorque_bins',...
...
'DstrokeMOD_wb_RollTorque_bins',...
'DstrokeMOD_ds_RollTorque_bins',...
'DstrokeMOD_us_RollTorque_bins',...
...
'DpitchMOD_wb_RollTorque_bins',...
'DpitchMOD_ds_RollTorque_bins',...
'DpitchMOD_us_RollTorque_bins',...
...
'DdevMOD_wb_RollTorque_bins',...
'DdevMOD_ds_RollTorque_bins',...
'DdevMOD_us_RollTorque_bins',...
...
't_wb_up_RollTorque',...
't_wb_down_RollTorque',...
'stroke_wb_up_RollTorque',...
'stroke_wb_down_RollTorque',...
'pitch_wb_up_RollTorque',...
'pitch_wb_down_RollTorque',...
'dev_wb_up_RollTorque',...
'dev_wb_down_RollTorque',...
'aoa_wb_up_RollTorque',...
'aoa_wb_down_RollTorque',...
'U_wb_up_RollTorque',...
'U_wb_down_RollTorque',...
...
't_ds_up_RollTorque',...
't_ds_down_RollTorque',...
'stroke_ds_up_RollTorque',...
'stroke_ds_down_RollTorque',...
'pitch_ds_up_RollTorque',...
'pitch_ds_down_RollTorque',...
'dev_ds_up_RollTorque',...
'dev_ds_down_RollTorque',...
'aoa_ds_up_RollTorque',...
'aoa_ds_down_RollTorque',...
'U_ds_up_RollTorque',...
'U_ds_down_RollTorque',...
...
't_us_up_RollTorque',...
't_us_down_RollTorque',...
'stroke_us_up_RollTorque',...
'stroke_us_down_RollTorque',...
'pitch_us_up_RollTorque',...
'pitch_us_down_RollTorque',...
'dev_us_up_RollTorque',...
'dev_us_down_RollTorque',...
'aoa_us_up_RollTorque',...
'aoa_us_down_RollTorque',...
'U_us_up_RollTorque',...
'U_us_down_RollTorque',...
...
't_wb_RollTorque_bins',...
'stroke_wb_up_RollTorque_bins',...
'stroke_wb_down_RollTorque_bins',...
'pitch_wb_up_RollTorque_bins',...
'pitch_wb_down_RollTorque_bins',...
'dev_wb_up_RollTorque_bins',...
'dev_wb_down_RollTorque_bins',...
'aoa_wb_up_RollTorque_bins',...
'aoa_wb_down_RollTorque_bins',...
'U_wb_up_RollTorque_bins',...
'U_wb_down_RollTorque_bins',...
'Dstroke_wb_RollTorque_bins',...
'Dpitch_wb_RollTorque_bins',...
'Ddev_wb_RollTorque_bins',...
'Daoa_wb_RollTorque_bins',...
'DU_wb_RollTorque_bins',...
...
't_ds_RollTorque_bins',...
'stroke_ds_up_RollTorque_bins',...
'stroke_ds_down_RollTorque_bins',...
'pitch_ds_up_RollTorque_bins',...
'pitch_ds_down_RollTorque_bins',...
'dev_ds_up_RollTorque_bins',...
'dev_ds_down_RollTorque_bins',...
'aoa_ds_up_RollTorque_bins',...
'aoa_ds_down_RollTorque_bins',...
'U_ds_up_RollTorque_bins',...
'U_ds_down_RollTorque_bins',...
'Dstroke_ds_RollTorque_bins',...
'Dpitch_ds_RollTorque_bins',...
'Ddev_ds_RollTorque_bins',...
'Daoa_ds_RollTorque_bins',...
'DU_ds_RollTorque_bins',...
...
't_us_RollTorque_bins',...
'stroke_us_up_RollTorque_bins',...
'stroke_us_down_RollTorque_bins',...
'pitch_us_up_RollTorque_bins',...
'pitch_us_down_RollTorque_bins',...
'dev_us_up_RollTorque_bins',...
'dev_us_down_RollTorque_bins',...
'aoa_us_up_RollTorque_bins',...
'aoa_us_down_RollTorque_bins',...
'U_us_up_RollTorque_bins',...
'U_us_down_RollTorque_bins',...
'Dstroke_us_RollTorque_bins',...
'Dpitch_us_RollTorque_bins',...
'Ddev_us_RollTorque_bins',...
'Daoa_us_RollTorque_bins',...
'DU_us_RollTorque_bins');

end
