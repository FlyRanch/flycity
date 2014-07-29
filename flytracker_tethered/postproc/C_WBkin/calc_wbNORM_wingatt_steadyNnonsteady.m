% steady wingbeat

rollaccel_limit_steady = limit_steady*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_limit_steady = limit_steady*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_limit_steady = limit_steady*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_limit_steady = limit_steady*nanstd(F_mean_wb(:));

n_now=0;
m_now=0;

%% steady
seq_nr_steady = [];
wb_nr_steady = [];

t_wb_L_steady = [];
t_wb_R_steady = [];
stroke_wb_L_steady = [];
stroke_wb_R_steady = [];
pitch_wb_L_steady = [];
pitch_wb_R_steady = [];
dev_wb_L_steady = [];
dev_wb_R_steady = [];
aoa_wb_L_steady = [];
aoa_wb_R_steady = [];
U_wb_L_steady = [];
U_wb_R_steady = [];
Dstroke_wb_steady = [];
Dpitch_wb_steady = [];
Ddev_wb_steady = [];
Daoa_wb_steady = [];
DU_wb_steady = [];

t_ds_L_steady = [];
t_ds_R_steady = [];
stroke_ds_L_steady = [];
stroke_ds_R_steady = [];
pitch_ds_L_steady = [];
pitch_ds_R_steady = [];
dev_ds_L_steady = [];
dev_ds_R_steady = [];
aoa_ds_L_steady = [];
aoa_ds_R_steady = [];
U_ds_L_steady = [];
U_ds_R_steady = [];
Dstroke_ds_steady = [];
Dpitch_ds_steady = [];
Ddev_ds_steady = [];
Daoa_ds_steady = [];
DU_ds_steady = [];

t_us_L_steady = [];
t_us_R_steady = [];
stroke_us_L_steady = [];
stroke_us_R_steady = [];
pitch_us_L_steady = [];
pitch_us_R_steady = [];
dev_us_L_steady = [];
dev_us_R_steady = [];
aoa_us_L_steady = [];
aoa_us_R_steady = [];
U_us_L_steady = [];
U_us_R_steady = [];
Dstroke_us_steady = [];
Dpitch_us_steady = [];
Ddev_us_steady = [];
Daoa_us_steady = [];
DU_us_steady = [];

t_wb_steady_bins = [];
stroke_wb_L_steady_bins = [];
stroke_wb_R_steady_bins = [];
pitch_wb_L_steady_bins = [];
pitch_wb_R_steady_bins = [];
dev_wb_L_steady_bins = [];
dev_wb_R_steady_bins = [];
aoa_wb_L_steady_bins = [];
aoa_wb_R_steady_bins = [];
U_wb_L_steady_bins = [];
U_wb_R_steady_bins = [];
Dstroke_wb_steady_bins = [];
Dpitch_wb_steady_bins = [];
Ddev_wb_steady_bins = [];
Daoa_wb_steady_bins = [];
DU_wb_steady_bins = [];

t_ds_steady_bins = [];
stroke_ds_L_steady_bins = [];
stroke_ds_R_steady_bins = [];
pitch_ds_L_steady_bins = [];
pitch_ds_R_steady_bins = [];
dev_ds_L_steady_bins = [];
dev_ds_R_steady_bins = [];
aoa_ds_L_steady_bins = [];
aoa_ds_R_steady_bins = [];
U_ds_L_steady_bins = [];
U_ds_R_steady_bins = [];
Dstroke_ds_steady_bins = [];
Dpitch_ds_steady_bins = [];
Ddev_ds_steady_bins = [];
Daoa_ds_steady_bins = [];
DU_ds_steady_bins = [];

t_us_steady_bins = [];
stroke_us_L_steady_bins = [];
stroke_us_R_steady_bins = [];
pitch_us_L_steady_bins = [];
pitch_us_R_steady_bins = [];
dev_us_L_steady_bins = [];
dev_us_R_steady_bins = [];
aoa_us_L_steady_bins = [];
aoa_us_R_steady_bins = [];
U_us_L_steady_bins = [];
U_us_R_steady_bins = [];
Dstroke_us_steady_bins = [];
Dpitch_us_steady_bins = [];
Ddev_us_steady_bins = [];
Daoa_us_steady_bins = [];
DU_us_steady_bins = [];

dt_wb_steady = [];
dt_ds_steady = [];
dt_us_steady = [];
f_wb_steady = [];
Rds_steady = [];

V_steady = [];
pitch_global_steady = [];

%% NONsteady
seq_nr_NONsteady = [];
wb_nr_NONsteady = [];

t_wb_L_NONsteady = [];
t_wb_R_NONsteady = [];
stroke_wb_L_NONsteady = [];
stroke_wb_R_NONsteady = [];
pitch_wb_L_NONsteady = [];
pitch_wb_R_NONsteady = [];
dev_wb_L_NONsteady = [];
dev_wb_R_NONsteady = [];
aoa_wb_L_NONsteady = [];
aoa_wb_R_NONsteady = [];
U_wb_L_NONsteady = [];
U_wb_R_NONsteady = [];
Dstroke_wb_NONsteady = [];
Dpitch_wb_NONsteady = [];
Ddev_wb_NONsteady = [];
Daoa_wb_NONsteady = [];
DU_wb_NONsteady = [];

t_ds_L_NONsteady = [];
t_ds_R_NONsteady = [];
stroke_ds_L_NONsteady = [];
stroke_ds_R_NONsteady = [];
pitch_ds_L_NONsteady = [];
pitch_ds_R_NONsteady = [];
dev_ds_L_NONsteady = [];
dev_ds_R_NONsteady = [];
aoa_ds_L_NONsteady = [];
aoa_ds_R_NONsteady = [];
U_ds_L_NONsteady = [];
U_ds_R_NONsteady = [];
Dstroke_ds_NONsteady = [];
Dpitch_ds_NONsteady = [];
Ddev_ds_NONsteady = [];
Daoa_ds_NONsteady = [];
DU_ds_NONsteady = [];

t_us_L_NONsteady = [];
t_us_R_NONsteady = [];
stroke_us_L_NONsteady = [];
stroke_us_R_NONsteady = [];
pitch_us_L_NONsteady = [];
pitch_us_R_NONsteady = [];
dev_us_L_NONsteady = [];
dev_us_R_NONsteady = [];
aoa_us_L_NONsteady = [];
aoa_us_R_NONsteady = [];
U_us_L_NONsteady = [];
U_us_R_NONsteady = [];
Dstroke_us_NONsteady = [];
Dpitch_us_NONsteady = [];
Ddev_us_NONsteady = [];
Daoa_us_NONsteady = [];
DU_us_NONsteady = [];

t_wb_NONsteady_bins = [];
stroke_wb_L_NONsteady_bins = [];
stroke_wb_R_NONsteady_bins = [];
pitch_wb_L_NONsteady_bins = [];
pitch_wb_R_NONsteady_bins = [];
dev_wb_L_NONsteady_bins = [];
dev_wb_R_NONsteady_bins = [];
aoa_wb_L_NONsteady_bins = [];
aoa_wb_R_NONsteady_bins = [];
U_wb_L_NONsteady_bins = [];
U_wb_R_NONsteady_bins = [];
Dstroke_wb_NONsteady_bins = [];
Dpitch_wb_NONsteady_bins = [];
Ddev_wb_NONsteady_bins = [];
Daoa_wb_NONsteady_bins = [];
DU_wb_NONsteady_bins = [];

t_ds_NONsteady_bins = [];
stroke_ds_L_NONsteady_bins = [];
stroke_ds_R_NONsteady_bins = [];
pitch_ds_L_NONsteady_bins = [];
pitch_ds_R_NONsteady_bins = [];
dev_ds_L_NONsteady_bins = [];
dev_ds_R_NONsteady_bins = [];
aoa_ds_L_NONsteady_bins = [];
aoa_ds_R_NONsteady_bins = [];
U_ds_L_NONsteady_bins = [];
U_ds_R_NONsteady_bins = [];
Dstroke_ds_NONsteady_bins = [];
Dpitch_ds_NONsteady_bins = [];
Ddev_ds_NONsteady_bins = [];
Daoa_ds_NONsteady_bins = [];
DU_ds_NONsteady_bins = [];

t_us_NONsteady_bins = [];
stroke_us_L_NONsteady_bins = [];
stroke_us_R_NONsteady_bins = [];
pitch_us_L_NONsteady_bins = [];
pitch_us_R_NONsteady_bins = [];
dev_us_L_NONsteady_bins = [];
dev_us_R_NONsteady_bins = [];
aoa_us_L_NONsteady_bins = [];
aoa_us_R_NONsteady_bins = [];
U_us_L_NONsteady_bins = [];
U_us_R_NONsteady_bins = [];
Dstroke_us_NONsteady_bins = [];
Dpitch_us_NONsteady_bins = [];
Ddev_us_NONsteady_bins = [];
Daoa_us_NONsteady_bins = [];
DU_us_NONsteady_bins = [];

dt_wb_NONsteady = [];
dt_ds_NONsteady = [];
dt_us_NONsteady = [];
f_wb_NONsteady = [];
Rds_NONsteady = [];

V_NONsteady = [];
pitch_global_NONsteady = [];

%% loop
for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        %% steady wb
        if abs(roll_dot_dot_mean_wb(wb)) < rollaccel_limit_steady &&...
                abs(pitch_dot_dot_mean_wb(wb)) < pitchaccel_limit_steady &&...
                abs(yaw_dot_dot_mean_wb(wb)) < yawaccel_limit_steady &&...
                abs(F_mean_wb(wb)-1) < Fenhance_limit_steady
            
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            
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
            
            
            %% store steady data
            n_now=n_now+1;

            seq_nr_steady(n_now,1) = seq_nr_now;
            wb_nr_steady(n_now,1) = wb_nr_now;
            
            V_steady(n_now,1) = V_now;
            pitch_global_steady(n_now,1) = pitch_global_now;
            
            dt_ds_steady(n_now,1) = nanmean(dt_ds_now);
            dt_us_steady(n_now,1) = nanmean(dt_us_now);
            f_wb_steady(n_now,1) = nanmean(f_wb_now);
            Rds_steady(n_now,1) = nanmean(Rds_now);

            
            t_wb_L_steady(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R_steady(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L_steady(1:length(t_wb_L_now),n_now) = stroke_wb_L_now;
            stroke_wb_R_steady(1:length(t_wb_R_now),n_now) = stroke_wb_R_now;
            pitch_wb_L_steady(1:length(t_wb_L_now),n_now) = pitch_wb_L_now;
            pitch_wb_R_steady(1:length(t_wb_R_now),n_now) = pitch_wb_R_now;
            dev_wb_L_steady(1:length(t_wb_L_now),n_now) = dev_wb_L_now;
            dev_wb_R_steady(1:length(t_wb_R_now),n_now) = dev_wb_R_now;
            aoa_wb_L_steady(1:length(t_wb_L_now),n_now) = aoa_wb_L_now;
            aoa_wb_R_steady(1:length(t_wb_R_now),n_now) = aoa_wb_R_now;
            U_wb_L_steady(1:length(t_wb_L_now),n_now) = U_wb_L_now;
            U_wb_R_steady(1:length(t_wb_R_now),n_now) = U_wb_R_now;
            
            t_ds_L_steady(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R_steady(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L_steady(1:length(t_ds_L_now),n_now) = stroke_ds_L_now;
            stroke_ds_R_steady(1:length(t_ds_R_now),n_now) = stroke_ds_R_now;
            pitch_ds_L_steady(1:length(t_ds_L_now),n_now) = pitch_ds_L_now;
            pitch_ds_R_steady(1:length(t_ds_R_now),n_now) = pitch_ds_R_now;
            dev_ds_L_steady(1:length(t_ds_L_now),n_now) = dev_ds_L_now;
            dev_ds_R_steady(1:length(t_ds_R_now),n_now) = dev_ds_R_now;
            aoa_ds_L_steady(1:length(t_ds_L_now),n_now) = aoa_ds_L_now;
            aoa_ds_R_steady(1:length(t_ds_R_now),n_now) = aoa_ds_R_now;
            U_ds_L_steady(1:length(t_ds_L_now),n_now) = U_ds_L_now;
            U_ds_R_steady(1:length(t_ds_R_now),n_now) = U_ds_R_now;
            
            t_us_L_steady(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R_steady(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L_steady(1:length(t_us_L_now),n_now) = stroke_us_L_now;
            stroke_us_R_steady(1:length(t_us_R_now),n_now) = stroke_us_R_now;
            pitch_us_L_steady(1:length(t_us_L_now),n_now) = pitch_us_L_now;
            pitch_us_R_steady(1:length(t_us_R_now),n_now) = pitch_us_R_now;
            dev_us_L_steady(1:length(t_us_L_now),n_now) = dev_us_L_now;
            dev_us_R_steady(1:length(t_us_R_now),n_now) = dev_us_R_now;
            aoa_us_L_steady(1:length(t_us_L_now),n_now) = aoa_us_L_now;
            aoa_us_R_steady(1:length(t_us_R_now),n_now) = aoa_us_R_now;
            U_us_L_steady(1:length(t_us_L_now),n_now) = U_us_L_now;
            U_us_R_steady(1:length(t_us_R_now),n_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_steady_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_steady_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_steady_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_steady_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_steady_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_steady_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_steady_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_steady_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_steady_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_steady_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_steady_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_steady_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_steady_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_steady_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_steady_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_steady_bins(:,n_now) = DU_wb_interp;
            
            t_ds_steady_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_steady_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_steady_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_steady_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_steady_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_steady_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_steady_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_steady_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_steady_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_steady_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_steady_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_steady_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_steady_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_steady_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_steady_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_steady_bins(:,n_now) = DU_ds_interp;
            
            t_us_steady_bins(:,n_now) = t_us_bin;
            stroke_us_L_steady_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_steady_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_steady_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_steady_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_steady_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_steady_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_steady_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_steady_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_steady_bins(:,n_now) = U_us_L_interp;
            U_us_R_steady_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_steady_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_steady_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_steady_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_steady_bins(:,n_now) = Daoa_us_interp;
            DU_us_steady_bins(:,n_now) = DU_us_interp;
            
        %% NONsteady wb
        else
            
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);

            V_now = V_mean_wb(wb);
            pitch_global_now = pitch_global_mean_wb(wb);
            
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
            
            
            %% store NONsteady data
            m_now=m_now+1;

            seq_nr_NONsteady(m_now,1) = seq_nr_now;
            wb_nr_NONsteady(m_now,1) = wb_nr_now;
            
            V_NONsteady(m_now,1) = V_now;
            pitch_global_NONsteady(m_now,1) = pitch_global_now;
            
            dt_ds_NONsteady(m_now,1) = nanmean(dt_ds_now);
            dt_us_NONsteady(m_now,1) = nanmean(dt_us_now);
            f_wb_NONsteady(m_now,1) = nanmean(f_wb_now);
            Rds_NONsteady(m_now,1) = nanmean(Rds_now);

            
            t_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = t_wb_L_now;
            t_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = t_wb_R_now;
            stroke_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = stroke_wb_L_now;
            stroke_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = stroke_wb_R_now;
            pitch_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = pitch_wb_L_now;
            pitch_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = pitch_wb_R_now;
            dev_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = dev_wb_L_now;
            dev_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = dev_wb_R_now;
            aoa_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = aoa_wb_L_now;
            aoa_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = aoa_wb_R_now;
            U_wb_L_NONsteady(1:length(t_wb_L_now),m_now) = U_wb_L_now;
            U_wb_R_NONsteady(1:length(t_wb_R_now),m_now) = U_wb_R_now;
            
            t_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = t_ds_L_now;
            t_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = t_ds_R_now;
            stroke_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = stroke_ds_L_now;
            stroke_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = stroke_ds_R_now;
            pitch_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = pitch_ds_L_now;
            pitch_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = pitch_ds_R_now;
            dev_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = dev_ds_L_now;
            dev_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = dev_ds_R_now;
            aoa_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = aoa_ds_L_now;
            aoa_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = aoa_ds_R_now;
            U_ds_L_NONsteady(1:length(t_ds_L_now),m_now) = U_ds_L_now;
            U_ds_R_NONsteady(1:length(t_ds_R_now),m_now) = U_ds_R_now;
            
            t_us_L_NONsteady(1:length(t_us_L_now),m_now) = t_us_L_now;
            t_us_R_NONsteady(1:length(t_us_R_now),m_now) = t_us_R_now;
            stroke_us_L_NONsteady(1:length(t_us_L_now),m_now) = stroke_us_L_now;
            stroke_us_R_NONsteady(1:length(t_us_R_now),m_now) = stroke_us_R_now;
            pitch_us_L_NONsteady(1:length(t_us_L_now),m_now) = pitch_us_L_now;
            pitch_us_R_NONsteady(1:length(t_us_R_now),m_now) = pitch_us_R_now;
            dev_us_L_NONsteady(1:length(t_us_L_now),m_now) = dev_us_L_now;
            dev_us_R_NONsteady(1:length(t_us_R_now),m_now) = dev_us_R_now;
            aoa_us_L_NONsteady(1:length(t_us_L_now),m_now) = aoa_us_L_now;
            aoa_us_R_NONsteady(1:length(t_us_R_now),m_now) = aoa_us_R_now;
            U_us_L_NONsteady(1:length(t_us_L_now),m_now) = U_us_L_now;
            U_us_R_NONsteady(1:length(t_us_R_now),m_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_NONsteady_bins(:,m_now) = t_wb_bin;
            stroke_wb_L_NONsteady_bins(:,m_now) = stroke_wb_L_interp;
            stroke_wb_R_NONsteady_bins(:,m_now) = stroke_wb_R_interp;
            pitch_wb_L_NONsteady_bins(:,m_now) = pitch_wb_L_interp;
            pitch_wb_R_NONsteady_bins(:,m_now) = pitch_wb_R_interp;
            dev_wb_L_NONsteady_bins(:,m_now) = dev_wb_L_interp;
            dev_wb_R_NONsteady_bins(:,m_now) = dev_wb_R_interp;
            aoa_wb_L_NONsteady_bins(:,m_now) = aoa_wb_L_interp;
            aoa_wb_R_NONsteady_bins(:,m_now) = aoa_wb_R_interp;
            U_wb_L_NONsteady_bins(:,m_now) = U_wb_L_interp;
            U_wb_R_NONsteady_bins(:,m_now) = U_wb_R_interp;
            Dstroke_wb_NONsteady_bins(:,m_now) = Dstroke_wb_interp;
            Dpitch_wb_NONsteady_bins(:,m_now) = Dpitch_wb_interp;
            Ddev_wb_NONsteady_bins(:,m_now) = Ddev_wb_interp;
            Daoa_wb_NONsteady_bins(:,m_now) = Daoa_wb_interp;
            DU_wb_NONsteady_bins(:,m_now) = DU_wb_interp;
            
            t_ds_NONsteady_bins(:,m_now) = t_ds_bin;
            stroke_ds_L_NONsteady_bins(:,m_now) = stroke_ds_L_interp;
            stroke_ds_R_NONsteady_bins(:,m_now) = stroke_ds_R_interp;
            pitch_ds_L_NONsteady_bins(:,m_now) = pitch_ds_L_interp;
            pitch_ds_R_NONsteady_bins(:,m_now) = pitch_ds_R_interp;
            dev_ds_L_NONsteady_bins(:,m_now) = dev_ds_L_interp;
            dev_ds_R_NONsteady_bins(:,m_now) = dev_ds_R_interp;
            aoa_ds_L_NONsteady_bins(:,m_now) = aoa_ds_L_interp;
            aoa_ds_R_NONsteady_bins(:,m_now) = aoa_ds_R_interp;
            U_ds_L_NONsteady_bins(:,m_now) = U_ds_L_interp;
            U_ds_R_NONsteady_bins(:,m_now) = U_ds_R_interp;
            Dstroke_ds_NONsteady_bins(:,m_now) = Dstroke_ds_interp;
            Dpitch_ds_NONsteady_bins(:,m_now) = Dpitch_ds_interp;
            Ddev_ds_NONsteady_bins(:,m_now) = Ddev_ds_interp;
            Daoa_ds_NONsteady_bins(:,m_now) = Daoa_ds_interp;
            DU_ds_NONsteady_bins(:,m_now) = DU_ds_interp;
            
            t_us_NONsteady_bins(:,m_now) = t_us_bin;
            stroke_us_L_NONsteady_bins(:,m_now) = stroke_us_L_interp;
            stroke_us_R_NONsteady_bins(:,m_now) = stroke_us_R_interp;
            pitch_us_L_NONsteady_bins(:,m_now) = pitch_us_L_interp;
            pitch_us_R_NONsteady_bins(:,m_now) = pitch_us_R_interp;
            dev_us_L_NONsteady_bins(:,m_now) = dev_us_L_interp;
            dev_us_R_NONsteady_bins(:,m_now) = dev_us_R_interp;
            aoa_us_L_NONsteady_bins(:,m_now) = aoa_us_L_interp;
            aoa_us_R_NONsteady_bins(:,m_now) = aoa_us_R_interp;
            U_us_L_NONsteady_bins(:,m_now) = U_us_L_interp;
            U_us_R_NONsteady_bins(:,m_now) = U_us_R_interp;
            Dstroke_us_NONsteady_bins(:,m_now) = Dstroke_us_interp;
            Dpitch_us_NONsteady_bins(:,m_now) = Dpitch_us_interp;
            Ddev_us_NONsteady_bins(:,m_now) = Ddev_us_interp;
            Daoa_us_NONsteady_bins(:,m_now) = Daoa_us_interp;
            DU_us_NONsteady_bins(:,m_now) = DU_us_interp;
            
    end
end

% % mean steady data
% V_steady_meanCIstd = [nanmean(V_steady) 1.96*nanstd(V_steady)/sqrt(length(V_steady)) nanstd(V_steady)];
% pitch_global_steady_meanCIstd = [nanmean(pitch_global_steady) 1.96*nanstd(pitch_global_steady)/sqrt(length(pitch_global_steady)) nanstd(pitch_global_steady)];
% 
% dt_ds_steady_meanCIstd = [nanmean(dt_ds_steady) 1.96*nanstd(dt_ds_steady)/sqrt(length(dt_ds_steady)) nanstd(dt_ds_steady)];
% dt_us_steady_meanCIstd = [nanmean(dt_us_steady) 1.96*nanstd(dt_us_steady)/sqrt(length(dt_us_steady)) nanstd(dt_us_steady)];
% f_wb_steady_meanCIstd = [nanmean(f_wb_steady) 1.96*nanstd(f_wb_steady)/sqrt(length(f_wb_steady)) nanstd(f_wb_steady)];
% Rds_steady_meanCIstd = [nanmean(Rds_steady) 1.96*nanstd(Rds_steady)/sqrt(length(Rds_steady)) nanstd(Rds_steady)];
% 
% 
% calc_WBfunc_steady_circmeanCIstd
% 
% 
% %% fit for average wb (bins)
% t_loc = t_wb_steady_bins(:,1);
% stroke_loc = stroke_wb_steady_bins_meanCIstd(:,1);
% pitch_loc = pitch_wb_steady_bins_meanCIstd(:,1);
% dev_loc = dev_wb_steady_bins_meanCIstd(:,1);
% Rds_loc = Rds_steady_meanCIstd(1);
% 
% % calc mean steady wb polies
% [stroke_steady_fit_binmean, stroke_steady_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
% [pitch_steady_fit_binmean, pitch_steady_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
% [dev_steady_fit_binmean, dev_steady_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);
% 
% % calc fourier series
% plotting = 0;
% % plotting = 1;
% 
% [stroke_steady_fourier_fit_binmean, stroke_steady_fourier_gof_binmean, stroke_steady_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plotting);
% [pitch_steady_fourier_fit_binmean, pitch_steady_fourier_gof_binmean, pitch_steady_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plotting);
% [dev_steady_fourier_fit_binmean, dev_steady_fourier_gof_binmean, dev_steady_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plotting);
% 
%     
% % save steady wb data
% n_steady = n_now   

if save_on == 1
    save(['WBdataset_steadyNnonsteady.mat']);
end
