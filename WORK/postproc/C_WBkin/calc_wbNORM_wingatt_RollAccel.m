
n_now=0;

% steady wb data
f_wb_steady = f_wb_steady_meanCIstd(:,1);

stroke_wb_steady = stroke_wb_steady_bins_meanCIstd(:,1);
stroke_ds_steady = stroke_ds_steady_bins_meanCIstd(:,1);
stroke_us_steady = stroke_us_steady_bins_meanCIstd(:,1);

pitch_wb_steady = pitch_wb_steady_bins_meanCIstd(:,1);
pitch_ds_steady = pitch_ds_steady_bins_meanCIstd(:,1);
pitch_us_steady = pitch_us_steady_bins_meanCIstd(:,1);

dev_wb_steady = dev_wb_steady_bins_meanCIstd(:,1);
dev_ds_steady = dev_ds_steady_bins_meanCIstd(:,1);
dev_us_steady = dev_us_steady_bins_meanCIstd(:,1);

% limits
limit = .5;
rollaccel_limit = limit*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_limit = limit*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_limit = limit*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_limit = limit*nanstd(F_mean_wb(:));

% normalization value
norm = 3;
rollaccel_norm = norm*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_norm = norm*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_norm = norm*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_norm = norm*nanstd(F_mean_wb(:));

%% variables
% wbMOD
freqMOD_wb_RollAccel = [];

strokeMOD_wb_L_RollAccel_bins = [];
strokeMOD_wb_R_RollAccel_bins = [];
strokeMOD_ds_L_RollAccel_bins = [];
strokeMOD_ds_R_RollAccel_bins = [];
strokeMOD_us_L_RollAccel_bins = [];
strokeMOD_us_R_RollAccel_bins = [];

pitchMOD_wb_L_RollAccel_bins = [];
pitchMOD_wb_R_RollAccel_bins = [];
pitchMOD_ds_L_RollAccel_bins = [];
pitchMOD_ds_R_RollAccel_bins = [];
pitchMOD_us_L_RollAccel_bins = [];
pitchMOD_us_R_RollAccel_bins = [];

devMOD_wb_L_RollAccel_bins = [];
devMOD_wb_R_RollAccel_bins = [];
devMOD_ds_L_RollAccel_bins = [];
devMOD_ds_R_RollAccel_bins = [];
devMOD_us_L_RollAccel_bins = [];
devMOD_us_R_RollAccel_bins = [];

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
Fenhance_RollAccel = [];

rollaccel_RollAccel = [];
pitchaccel_RollAccel = [];
yawaccel_RollAccel = [];

rollvel_RollAccel = [];
pitchvel_RollAccel = [];
yawvel_RollAccel = [];

% wingbeat kin
dt_wb_RollAccel = [];
dt_ds_RollAccel = [];
dt_us_RollAccel = [];
f_wb_RollAccel = [];
Rds_RollAccel = [];

t_wb_L_RollAccel = [];
t_wb_R_RollAccel = [];
stroke_wb_L_RollAccel = [];
stroke_wb_R_RollAccel = [];
pitch_wb_L_RollAccel = [];
pitch_wb_R_RollAccel = [];
dev_wb_L_RollAccel = [];
dev_wb_R_RollAccel = [];
aoa_wb_L_RollAccel = [];
aoa_wb_R_RollAccel = [];
U_wb_L_RollAccel = [];
U_wb_R_RollAccel = [];
Dstroke_wb_RollAccel = [];
Dpitch_wb_RollAccel = [];
Ddev_wb_RollAccel = [];
Daoa_wb_RollAccel = [];
DU_wb_RollAccel = [];

t_ds_L_RollAccel = [];
t_ds_R_RollAccel = [];
stroke_ds_L_RollAccel = [];
stroke_ds_R_RollAccel = [];
pitch_ds_L_RollAccel = [];
pitch_ds_R_RollAccel = [];
dev_ds_L_RollAccel = [];
dev_ds_R_RollAccel = [];
aoa_ds_L_RollAccel = [];
aoa_ds_R_RollAccel = [];
U_ds_L_RollAccel = [];
U_ds_R_RollAccel = [];
Dstroke_ds_RollAccel = [];
Dpitch_ds_RollAccel = [];
Ddev_ds_RollAccel = [];
Daoa_ds_RollAccel = [];
DU_ds_RollAccel = [];

t_us_L_RollAccel = [];
t_us_R_RollAccel = [];
stroke_us_L_RollAccel = [];
stroke_us_R_RollAccel = [];
pitch_us_L_RollAccel = [];
pitch_us_R_RollAccel = [];
dev_us_L_RollAccel = [];
dev_us_R_RollAccel = [];
aoa_us_L_RollAccel = [];
aoa_us_R_RollAccel = [];
U_us_L_RollAccel = [];
U_us_R_RollAccel = [];
Dstroke_us_RollAccel = [];
Dpitch_us_RollAccel = [];
Ddev_us_RollAccel = [];
Daoa_us_RollAccel = [];
DU_us_RollAccel = [];

t_wb_RollAccel_bins = [];
stroke_wb_L_RollAccel_bins = [];
stroke_wb_R_RollAccel_bins = [];
pitch_wb_L_RollAccel_bins = [];
pitch_wb_R_RollAccel_bins = [];
dev_wb_L_RollAccel_bins = [];
dev_wb_R_RollAccel_bins = [];
aoa_wb_L_RollAccel_bins = [];
aoa_wb_R_RollAccel_bins = [];
U_wb_L_RollAccel_bins = [];
U_wb_R_RollAccel_bins = [];
Dstroke_wb_RollAccel_bins = [];
Dpitch_wb_RollAccel_bins = [];
Ddev_wb_RollAccel_bins = [];
Daoa_wb_RollAccel_bins = [];
DU_wb_RollAccel_bins = [];

t_ds_RollAccel_bins = [];
stroke_ds_L_RollAccel_bins = [];
stroke_ds_R_RollAccel_bins = [];
pitch_ds_L_RollAccel_bins = [];
pitch_ds_R_RollAccel_bins = [];
dev_ds_L_RollAccel_bins = [];
dev_ds_R_RollAccel_bins = [];
aoa_ds_L_RollAccel_bins = [];
aoa_ds_R_RollAccel_bins = [];
U_ds_L_RollAccel_bins = [];
U_ds_R_RollAccel_bins = [];
Dstroke_ds_RollAccel_bins = [];
Dpitch_ds_RollAccel_bins = [];
Ddev_ds_RollAccel_bins = [];
Daoa_ds_RollAccel_bins = [];
DU_ds_RollAccel_bins = [];

t_us_RollAccel_bins = [];
stroke_us_L_RollAccel_bins = [];
stroke_us_R_RollAccel_bins = [];
pitch_us_L_RollAccel_bins = [];
pitch_us_R_RollAccel_bins = [];
dev_us_L_RollAccel_bins = [];
dev_us_R_RollAccel_bins = [];
aoa_us_L_RollAccel_bins = [];
aoa_us_R_RollAccel_bins = [];
U_us_L_RollAccel_bins = [];
U_us_R_RollAccel_bins = [];
Dstroke_us_RollAccel_bins = [];
Dpitch_us_RollAccel_bins = [];
Ddev_us_RollAccel_bins = [];
Daoa_us_RollAccel_bins = [];
DU_us_RollAccel_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if abs(roll_dot_dot_mean_wb(wb)) > rollaccel_limit
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            % body kin
            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            RollAccel_mean_wb_now = F_mean_wb(wb);

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
            freqMOD_wb_now = (nanmean(f_wb_now) - f_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            
            strokeMOD_wb_L_now = (stroke_wb_L_interp - stroke_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            strokeMOD_wb_R_now = (stroke_wb_L_interp - stroke_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            strokeMOD_ds_L_now = (stroke_ds_L_interp - stroke_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            strokeMOD_ds_R_now = (stroke_ds_L_interp - stroke_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            strokeMOD_us_L_now = (stroke_us_L_interp - stroke_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            strokeMOD_us_R_now = (stroke_us_L_interp - stroke_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            
            pitchMOD_wb_L_now = (pitch_wb_L_interp - pitch_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            pitchMOD_wb_R_now = (pitch_wb_L_interp - pitch_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            pitchMOD_ds_L_now = (pitch_ds_L_interp - pitch_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            pitchMOD_ds_R_now = (pitch_ds_L_interp - pitch_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            pitchMOD_us_L_now = (pitch_us_L_interp - pitch_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            pitchMOD_us_R_now = (pitch_us_L_interp - pitch_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            
            devMOD_wb_L_now = (dev_wb_L_interp - dev_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            devMOD_wb_R_now = (dev_wb_L_interp - dev_wb_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            devMOD_ds_L_now = (dev_ds_L_interp - dev_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            devMOD_ds_R_now = (dev_ds_L_interp - dev_ds_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            devMOD_us_L_now = (dev_us_L_interp - dev_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            devMOD_us_R_now = (dev_us_L_interp - dev_us_steady) / RollAccel_mean_wb_now * rollaccel_norm;
            
            DstrokeMOD_wb_now = Dstroke_wb_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DstrokeMOD_ds_now = Dstroke_ds_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DstrokeMOD_us_now = Dstroke_us_interp / RollAccel_mean_wb_now * rollaccel_norm;

            DpitchMOD_wb_now = Dpitch_wb_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DpitchMOD_ds_now = Dpitch_ds_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DpitchMOD_us_now = Dpitch_us_interp / RollAccel_mean_wb_now * rollaccel_norm;

            DdevMOD_wb_now = Ddev_wb_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DdevMOD_ds_now = Ddev_ds_interp / RollAccel_mean_wb_now * rollaccel_norm;
            DdevMOD_us_now = Ddev_us_interp / RollAccel_mean_wb_now * rollaccel_norm;

            %% store data
            n_now=n_now+1;

            seq_nr_RollAccel(n_now,1) = seq_nr_now;
            wb_nr_RollAccel(n_now,1) = wb_nr_now;
            
            % body kin
            V_RollAccel(n_now,1) = V_now;
            pitch_global_RollAccel(n_now,1) = pitch_global_now;
            Fenhance_RollAccel(n_now,1) = RollAccel_mean_wb_now;

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

            t_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = stroke_wb_L_now;
            stroke_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = stroke_wb_R_now;
            pitch_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = pitch_wb_L_now;
            pitch_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = pitch_wb_R_now;
            dev_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = dev_wb_L_now;
            dev_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = dev_wb_R_now;
            aoa_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = aoa_wb_L_now;
            aoa_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = aoa_wb_R_now;
            U_wb_L_RollAccel(1:length(t_wb_L_now),n_now) = U_wb_L_now;
            U_wb_R_RollAccel(1:length(t_wb_R_now),n_now) = U_wb_R_now;
            
            t_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = stroke_ds_L_now;
            stroke_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = stroke_ds_R_now;
            pitch_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = pitch_ds_L_now;
            pitch_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = pitch_ds_R_now;
            dev_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = dev_ds_L_now;
            dev_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = dev_ds_R_now;
            aoa_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = aoa_ds_L_now;
            aoa_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = aoa_ds_R_now;
            U_ds_L_RollAccel(1:length(t_ds_L_now),n_now) = U_ds_L_now;
            U_ds_R_RollAccel(1:length(t_ds_R_now),n_now) = U_ds_R_now;
            
            t_us_L_RollAccel(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R_RollAccel(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L_RollAccel(1:length(t_us_L_now),n_now) = stroke_us_L_now;
            stroke_us_R_RollAccel(1:length(t_us_R_now),n_now) = stroke_us_R_now;
            pitch_us_L_RollAccel(1:length(t_us_L_now),n_now) = pitch_us_L_now;
            pitch_us_R_RollAccel(1:length(t_us_R_now),n_now) = pitch_us_R_now;
            dev_us_L_RollAccel(1:length(t_us_L_now),n_now) = dev_us_L_now;
            dev_us_R_RollAccel(1:length(t_us_R_now),n_now) = dev_us_R_now;
            aoa_us_L_RollAccel(1:length(t_us_L_now),n_now) = aoa_us_L_now;
            aoa_us_R_RollAccel(1:length(t_us_R_now),n_now) = aoa_us_R_now;
            U_us_L_RollAccel(1:length(t_us_L_now),n_now) = U_us_L_now;
            U_us_R_RollAccel(1:length(t_us_R_now),n_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_RollAccel_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_RollAccel_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_RollAccel_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_RollAccel_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_RollAccel_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_RollAccel_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_RollAccel_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_RollAccel_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_RollAccel_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_RollAccel_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_RollAccel_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_RollAccel_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_RollAccel_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_RollAccel_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_RollAccel_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_RollAccel_bins(:,n_now) = DU_wb_interp;
            
            t_ds_RollAccel_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_RollAccel_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_RollAccel_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_RollAccel_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_RollAccel_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_RollAccel_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_RollAccel_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_RollAccel_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_RollAccel_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_RollAccel_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_RollAccel_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_RollAccel_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_RollAccel_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_RollAccel_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_RollAccel_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_RollAccel_bins(:,n_now) = DU_ds_interp;
            
            t_us_RollAccel_bins(:,n_now) = t_us_bin;
            stroke_us_L_RollAccel_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_RollAccel_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_RollAccel_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_RollAccel_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_RollAccel_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_RollAccel_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_RollAccel_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_RollAccel_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_RollAccel_bins(:,n_now) = U_us_L_interp;
            U_us_R_RollAccel_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_RollAccel_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_RollAccel_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_RollAccel_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_RollAccel_bins(:,n_now) = Daoa_us_interp;
            DU_us_RollAccel_bins(:,n_now) = DU_us_interp;
            
            freqMOD_wb_RollAccel(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_L_RollAccel_bins(:,n_now) = strokeMOD_wb_L_now;
            strokeMOD_wb_R_RollAccel_bins(:,n_now) = strokeMOD_wb_R_now;
            strokeMOD_ds_L_RollAccel_bins(:,n_now) = strokeMOD_ds_L_now;
            strokeMOD_ds_R_RollAccel_bins(:,n_now) = strokeMOD_ds_R_now;
            strokeMOD_us_L_RollAccel_bins(:,n_now) = strokeMOD_us_L_now;
            strokeMOD_us_R_RollAccel_bins(:,n_now) = strokeMOD_us_R_now;

            pitchMOD_wb_L_RollAccel_bins(:,n_now) = pitchMOD_wb_L_now;
            pitchMOD_wb_R_RollAccel_bins(:,n_now) = pitchMOD_wb_R_now;
            pitchMOD_ds_L_RollAccel_bins(:,n_now) = pitchMOD_ds_L_now;
            pitchMOD_ds_R_RollAccel_bins(:,n_now) = pitchMOD_ds_R_now;
            pitchMOD_us_L_RollAccel_bins(:,n_now) = pitchMOD_us_L_now;
            pitchMOD_us_R_RollAccel_bins(:,n_now) = pitchMOD_us_R_now;

            devMOD_wb_L_RollAccel_bins(:,n_now) = devMOD_wb_L_now;
            devMOD_wb_R_RollAccel_bins(:,n_now) = devMOD_wb_R_now;
            devMOD_ds_L_RollAccel_bins(:,n_now) = devMOD_ds_L_now;
            devMOD_ds_R_RollAccel_bins(:,n_now) = devMOD_ds_R_now;
            devMOD_us_L_RollAccel_bins(:,n_now) = devMOD_us_L_now;
            devMOD_us_R_RollAccel_bins(:,n_now) = devMOD_us_R_now;
            
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

%% legendre polynomials
% wingbeats
    n_pol_stroke = 12; % Order of used polynomials
    n_pol_pitch = 14; % Order of used polynomials
    n_pol_dev = 10; % Order of used polynomials 

    % fit for average wb (bins)
    t_loc = t_wb_RollAccel_bins(:,1);
    stroke_loc = stroke_wb_RollAccel_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_RollAccel_bins_meanCIstd(:,1);
    dev_loc = dev_wb_RollAccel_bins_meanCIstd(:,1);
    Rds_loc = Rds_RollAccel_meanCIstd(1);

    [stroke_RollAccel_fit_binmean, stroke_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
    [pitch_RollAccel_fit_binmean, pitch_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
    [dev_RollAccel_fit_binmean, dev_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

% wingbeats
    n_pol_strokeMOD = 10; % Order of used polynomials
    n_pol_pitchMOD = 10; % Order of used polynomials
    n_pol_devMOD = 10; % Order of used polynomials 

    % fit for average wb (bins)
    t_loc = t_wb_RollAccel_bins(:,1);
    strokeMOD_loc = strokeMOD_wb_RollAccel_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_RollAccel_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_RollAccel_bins_meanCIstd(:,1);
    Rds_loc = Rds_steady_meanCIstd(1);

    [strokeMOD_RollAccel_fit_binmean, strokeMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_strokeMOD,t_loc,strokeMOD_loc,Rds_loc);
    [pitchMOD_RollAccel_fit_binmean, pitchMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitchMOD,t_loc,pitchMOD_loc,Rds_loc);
    [devMOD_RollAccel_fit_binmean, devMOD_RollAccel_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_devMOD,t_loc,devMOD_loc,Rds_loc);

n_RollAccel = n_now   

    % save RollAccel wb data
    cd ..
    save(['WBmod_RollAccel_',num2str(n_RollAccel),'WBs.mat'],...
'n_RollAccel',...
'seq_nr_RollAccel',...
'wb_nr_RollAccel',...
...
'limit',...
'rollaccel_limit',...
'pitchaccel_limit',...
'yawaccel_limit',...
'Fenhance_limit',...
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
...
'stroke_RollAccel_fit_binmean_periodic',...
'pitch_RollAccel_fit_binmean_periodic',...
'dev_RollAccel_fit_binmean_periodic',...
...
'n_pol_strokeMOD',...
'n_pol_pitchMOD',...
'n_pol_devMOD',...
...
'strokeMOD_RollAccel_fit_binmean_periodic',...
'pitchMOD_RollAccel_fit_binmean_periodic',...
'devMOD_RollAccel_fit_binmean_periodic',...
...
'stroke_wb_RollAccel_bins_meanCIstd',...
'stroke_ds_RollAccel_bins_meanCIstd',...
'stroke_us_RollAccel_bins_meanCIstd',...
...
'pitch_wb_RollAccel_bins_meanCIstd',...
'pitch_ds_RollAccel_bins_meanCIstd',...
'pitch_us_RollAccel_bins_meanCIstd',...
...
'dev_wb_RollAccel_bins_meanCIstd',...
'dev_ds_RollAccel_bins_meanCIstd',...
'dev_us_RollAccel_bins_meanCIstd',...
...
'strokeMOD_wb_RollAccel_bins_meanCIstd',...
'strokeMOD_ds_RollAccel_bins_meanCIstd',...
'strokeMOD_us_RollAccel_bins_meanCIstd',...
...
'pitchMOD_wb_RollAccel_bins_meanCIstd',...
'pitchMOD_ds_RollAccel_bins_meanCIstd',...
'pitchMOD_us_RollAccel_bins_meanCIstd',...
...
'devMOD_wb_RollAccel_bins_meanCIstd',...
'devMOD_ds_RollAccel_bins_meanCIstd',...
'devMOD_us_RollAccel_bins_meanCIstd',...
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
'V_RollAccel',...
'pitch_global_RollAccel',...
'dt_ds_RollAccel',...
'dt_us_RollAccel',...
'f_wb_RollAccel',...
'Rds_RollAccel',...
...
'V_RollAccel_meanCIstd',...
'pitch_global_RollAccel_meanCIstd',...
'Fenhance_RollAccel',...
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
'strokeMOD_wb_L_RollAccel_bins',...
'strokeMOD_wb_R_RollAccel_bins',...
'strokeMOD_ds_L_RollAccel_bins',...
'strokeMOD_ds_R_RollAccel_bins',...
'strokeMOD_us_L_RollAccel_bins',...
'strokeMOD_us_R_RollAccel_bins',...
...
'pitchMOD_wb_L_RollAccel_bins',...
'pitchMOD_wb_R_RollAccel_bins',...
'pitchMOD_ds_L_RollAccel_bins',...
'pitchMOD_ds_R_RollAccel_bins',...
'pitchMOD_us_L_RollAccel_bins',...
'pitchMOD_us_R_RollAccel_bins',...
...
'devMOD_wb_L_RollAccel_bins',...
'devMOD_wb_R_RollAccel_bins',...
'devMOD_ds_L_RollAccel_bins',...
'devMOD_ds_R_RollAccel_bins',...
'devMOD_us_L_RollAccel_bins',...
'devMOD_us_R_RollAccel_bins',...
...
't_wb_L_RollAccel',...
't_wb_R_RollAccel',...
'stroke_wb_L_RollAccel',...
'stroke_wb_R_RollAccel',...
'pitch_wb_L_RollAccel',...
'pitch_wb_R_RollAccel',...
'dev_wb_L_RollAccel',...
'dev_wb_R_RollAccel',...
'aoa_wb_L_RollAccel',...
'aoa_wb_R_RollAccel',...
'U_wb_L_RollAccel',...
'U_wb_R_RollAccel',...
...
't_ds_L_RollAccel',...
't_ds_R_RollAccel',...
'stroke_ds_L_RollAccel',...
'stroke_ds_R_RollAccel',...
'pitch_ds_L_RollAccel',...
'pitch_ds_R_RollAccel',...
'dev_ds_L_RollAccel',...
'dev_ds_R_RollAccel',...
'aoa_ds_L_RollAccel',...
'aoa_ds_R_RollAccel',...
'U_ds_L_RollAccel',...
'U_ds_R_RollAccel',...
...
't_us_L_RollAccel',...
't_us_R_RollAccel',...
'stroke_us_L_RollAccel',...
'stroke_us_R_RollAccel',...
'pitch_us_L_RollAccel',...
'pitch_us_R_RollAccel',...
'dev_us_L_RollAccel',...
'dev_us_R_RollAccel',...
'aoa_us_L_RollAccel',...
'aoa_us_R_RollAccel',...
'U_us_L_RollAccel',...
'U_us_R_RollAccel',...
...
't_wb_RollAccel_bins',...
'stroke_wb_L_RollAccel_bins',...
'stroke_wb_R_RollAccel_bins',...
'pitch_wb_L_RollAccel_bins',...
'pitch_wb_R_RollAccel_bins',...
'dev_wb_L_RollAccel_bins',...
'dev_wb_R_RollAccel_bins',...
'aoa_wb_L_RollAccel_bins',...
'aoa_wb_R_RollAccel_bins',...
'U_wb_L_RollAccel_bins',...
'U_wb_R_RollAccel_bins',...
'Dstroke_wb_RollAccel_bins',...
'Dpitch_wb_RollAccel_bins',...
'Ddev_wb_RollAccel_bins',...
'Daoa_wb_RollAccel_bins',...
'DU_wb_RollAccel_bins',...
...
't_ds_RollAccel_bins',...
'stroke_ds_L_RollAccel_bins',...
'stroke_ds_R_RollAccel_bins',...
'pitch_ds_L_RollAccel_bins',...
'pitch_ds_R_RollAccel_bins',...
'dev_ds_L_RollAccel_bins',...
'dev_ds_R_RollAccel_bins',...
'aoa_ds_L_RollAccel_bins',...
'aoa_ds_R_RollAccel_bins',...
'U_ds_L_RollAccel_bins',...
'U_ds_R_RollAccel_bins',...
'Dstroke_ds_RollAccel_bins',...
'Dpitch_ds_RollAccel_bins',...
'Ddev_ds_RollAccel_bins',...
'Daoa_ds_RollAccel_bins',...
'DU_ds_RollAccel_bins',...
...
't_us_RollAccel_bins',...
'stroke_us_L_RollAccel_bins',...
'stroke_us_R_RollAccel_bins',...
'pitch_us_L_RollAccel_bins',...
'pitch_us_R_RollAccel_bins',...
'dev_us_L_RollAccel_bins',...
'dev_us_R_RollAccel_bins',...
'aoa_us_L_RollAccel_bins',...
'aoa_us_R_RollAccel_bins',...
'U_us_L_RollAccel_bins',...
'U_us_R_RollAccel_bins',...
'Dstroke_us_RollAccel_bins',...
'Dpitch_us_RollAccel_bins',...
'Ddev_us_RollAccel_bins',...
'Daoa_us_RollAccel_bins',...
'DU_us_RollAccel_bins');

cd('WBmean_figs')


