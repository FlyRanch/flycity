% calc wb MOD CLIP d2ndMom

n_now=0;

%% variables

% bodyMODs
velMOD_dS3 = [];

slipMOD_dS3 = [];
pitchMOD_dS3 = [];
rollMOD_dS3 = [];

Fsp_pitchMOD_dS3 = [];
Fsp_rollMOD_dS3 = [];

% wbMOD
freqMOD_wb_dS3 = [];

strokeMOD_wb_L_dS3_bins = [];
strokeMOD_wb_R_dS3_bins = [];
strokeMOD_ds_L_dS3_bins = [];
strokeMOD_ds_R_dS3_bins = [];
strokeMOD_us_L_dS3_bins = [];
strokeMOD_us_R_dS3_bins = [];

pitchMOD_wb_L_dS3_bins = [];
pitchMOD_wb_R_dS3_bins = [];
pitchMOD_ds_L_dS3_bins = [];
pitchMOD_ds_R_dS3_bins = [];
pitchMOD_us_L_dS3_bins = [];
pitchMOD_us_R_dS3_bins = [];

devMOD_wb_L_dS3_bins = [];
devMOD_wb_R_dS3_bins = [];
devMOD_ds_L_dS3_bins = [];
devMOD_ds_R_dS3_bins = [];
devMOD_us_L_dS3_bins = [];
devMOD_us_R_dS3_bins = [];

DstrokeMOD_wb_dS3_bins = [];
DstrokeMOD_ds_dS3_bins = [];
DstrokeMOD_us_dS3_bins = [];

DpitchMOD_wb_dS3_bins = [];
DpitchMOD_ds_dS3_bins = [];
DpitchMOD_us_dS3_bins = [];

DdevMOD_wb_dS3_bins = [];
DdevMOD_ds_dS3_bins = [];
DdevMOD_us_dS3_bins = [];

% seq&wb
seq_nr_dS3 = [];
wb_nr_dS3 = [];

% clip data
clip_side_dS3 = [];
clip_type_dS3 = [];

S2Clipped_dS3 = [];
S2Intact_dS3 = [];
S2Ratio_dS3 = [];

S3Clipped_dS3 = [];
S3Intact_dS3 = [];
S3Ratio_dS3 = [];

AreaClipped_dS3 = [];
AreaIntact_dS3 = [];
AreaRatio_dS3 = [];

LengthClipped_dS3 = [];
LengthIntact_dS3 = [];
LengthRatio_dS3 = [];

% body kin
vel_dS3 = [];

slip_dS3 = [];
pitch_dS3 = [];
roll_dS3 = [];

Fsp_pitch_dS3 = [];
Fsp_roll_dS3 = [];

% wingbeat kin
dt_wb_dS3 = [];
dt_ds_dS3 = [];
dt_us_dS3 = [];
f_wb_dS3 = [];
Rds_dS3 = [];

t_wb_L_dS3 = [];
t_wb_R_dS3 = [];
stroke_wb_L_dS3 = [];
stroke_wb_R_dS3 = [];
pitch_wb_L_dS3 = [];
pitch_wb_R_dS3 = [];
dev_wb_L_dS3 = [];
dev_wb_R_dS3 = [];
aoa_wb_L_dS3 = [];
aoa_wb_R_dS3 = [];
U_wb_L_dS3 = [];
U_wb_R_dS3 = [];
Dstroke_wb_dS3 = [];
Dpitch_wb_dS3 = [];
Ddev_wb_dS3 = [];
Daoa_wb_dS3 = [];
DU_wb_dS3 = [];

t_ds_L_dS3 = [];
t_ds_R_dS3 = [];
stroke_ds_L_dS3 = [];
stroke_ds_R_dS3 = [];
pitch_ds_L_dS3 = [];
pitch_ds_R_dS3 = [];
dev_ds_L_dS3 = [];
dev_ds_R_dS3 = [];
aoa_ds_L_dS3 = [];
aoa_ds_R_dS3 = [];
U_ds_L_dS3 = [];
U_ds_R_dS3 = [];
Dstroke_ds_dS3 = [];
Dpitch_ds_dS3 = [];
Ddev_ds_dS3 = [];
Daoa_ds_dS3 = [];
DU_ds_dS3 = [];

t_us_L_dS3 = [];
t_us_R_dS3 = [];
stroke_us_L_dS3 = [];
stroke_us_R_dS3 = [];
pitch_us_L_dS3 = [];
pitch_us_R_dS3 = [];
dev_us_L_dS3 = [];
dev_us_R_dS3 = [];
aoa_us_L_dS3 = [];
aoa_us_R_dS3 = [];
U_us_L_dS3 = [];
U_us_R_dS3 = [];
Dstroke_us_dS3 = [];
Dpitch_us_dS3 = [];
Ddev_us_dS3 = [];
Daoa_us_dS3 = [];
DU_us_dS3 = [];

t_wb_dS3_bins = [];
stroke_wb_L_dS3_bins = [];
stroke_wb_R_dS3_bins = [];
pitch_wb_L_dS3_bins = [];
pitch_wb_R_dS3_bins = [];
dev_wb_L_dS3_bins = [];
dev_wb_R_dS3_bins = [];
aoa_wb_L_dS3_bins = [];
aoa_wb_R_dS3_bins = [];
U_wb_L_dS3_bins = [];
U_wb_R_dS3_bins = [];
Dstroke_wb_dS3_bins = [];
Dpitch_wb_dS3_bins = [];
Ddev_wb_dS3_bins = [];
Daoa_wb_dS3_bins = [];
DU_wb_dS3_bins = [];

t_ds_dS3_bins = [];
stroke_ds_L_dS3_bins = [];
stroke_ds_R_dS3_bins = [];
pitch_ds_L_dS3_bins = [];
pitch_ds_R_dS3_bins = [];
dev_ds_L_dS3_bins = [];
dev_ds_R_dS3_bins = [];
aoa_ds_L_dS3_bins = [];
aoa_ds_R_dS3_bins = [];
U_ds_L_dS3_bins = [];
U_ds_R_dS3_bins = [];
Dstroke_ds_dS3_bins = [];
Dpitch_ds_dS3_bins = [];
Ddev_ds_dS3_bins = [];
Daoa_ds_dS3_bins = [];
DU_ds_dS3_bins = [];

t_us_dS3_bins = [];
stroke_us_L_dS3_bins = [];
stroke_us_R_dS3_bins = [];
pitch_us_L_dS3_bins = [];
pitch_us_R_dS3_bins = [];
dev_us_L_dS3_bins = [];
dev_us_R_dS3_bins = [];
aoa_us_L_dS3_bins = [];
aoa_us_R_dS3_bins = [];
U_us_L_dS3_bins = [];
U_us_R_dS3_bins = [];
Dstroke_us_dS3_bins = [];
Dpitch_us_dS3_bins = [];
Ddev_us_dS3_bins = [];
Daoa_us_dS3_bins = [];
DU_us_dS3_bins = [];

for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if steady_nr_mean_wb(wb) == 1 && S3Ratio(wb) < S3Ratio_cut
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);
            
            % clip data
            clip_side_now = clip_side(wb);
            clip_type_now = clip_type(wb);

            S2Clipped_now = S2Clipped(wb);
            S2Intact_now = S2Intact(wb);
            S2Ratio_now = S2Ratio(wb);

            S3Clipped_now = S3Clipped(wb);
            S3Intact_now = S3Intact(wb);
            S3Ratio_now = S3Ratio(wb);

            AreaClipped_now = AreaClipped(wb);
            AreaIntact_now = AreaIntact(wb);
            AreaRatio_now = AreaRatio(wb);

            LengthClipped_now = LengthClipped(wb);
            LengthIntact_now = LengthIntact(wb);
            LengthRatio_now = LengthRatio(wb);

            % body kin
            vel_now = V_mean_wb(wb);
            
            slip_now = slip_mean_wb(wb);
            pitch_now = pitch_mean_wb(wb);
            roll_now = roll_mean_wb(wb);
            
            Fsp_pitch_now = Fsp_pitch_mean_wb(wb);
            Fsp_roll_now = Fsp_roll_mean_wb(wb);
            
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
            
            
            %% calc bodykin modifications
            velMOD_now   = (nanmean(vel_now) - vel_steady)          / (S3Ratio_now-1);

            slipMOD_now  =  nanmean(slip_now)                       / (S3Ratio_now-1);
            pitchMOD_now = (nanmean(pitch_now) - pitch_body_steady) / (S3Ratio_now-1);
            rollMOD_now  =  nanmean(roll_now)                       / (S3Ratio_now-1);

            Fsp_pitchMOD_now = nanmean(Fsp_pitch_now)               / (S3Ratio_now-1);
            Fsp_rollMOD_now  = nanmean(Fsp_roll_now)                / (S3Ratio_now-1);

            %% calc wingbeat modifications
            freqMOD_wb_now =     (nanmean(f_wb_now)  - f_wb_steady)      / (S3Ratio_now-1);
            
            strokeMOD_wb_L_now = (stroke_wb_L_interp - stroke_wb_steady) / (S3Ratio_now-1);
            strokeMOD_wb_R_now = (stroke_wb_R_interp - stroke_wb_steady) / (S3Ratio_now-1);
            strokeMOD_ds_L_now = (stroke_ds_L_interp - stroke_ds_steady) / (S3Ratio_now-1);
            strokeMOD_ds_R_now = (stroke_ds_R_interp - stroke_ds_steady) / (S3Ratio_now-1);
            strokeMOD_us_L_now = (stroke_us_L_interp - stroke_us_steady) / (S3Ratio_now-1);
            strokeMOD_us_R_now = (stroke_us_R_interp - stroke_us_steady) / (S3Ratio_now-1);
            
            pitchMOD_wb_L_now = (pitch_wb_L_interp - pitch_wb_steady) / (S3Ratio_now-1);
            pitchMOD_wb_R_now = (pitch_wb_R_interp - pitch_wb_steady) / (S3Ratio_now-1);
            pitchMOD_ds_L_now = (pitch_ds_L_interp - pitch_ds_steady) / (S3Ratio_now-1);
            pitchMOD_ds_R_now = (pitch_ds_R_interp - pitch_ds_steady) / (S3Ratio_now-1);
            pitchMOD_us_L_now = (pitch_us_L_interp - pitch_us_steady) / (S3Ratio_now-1);
            pitchMOD_us_R_now = (pitch_us_R_interp - pitch_us_steady) / (S3Ratio_now-1);
            
            devMOD_wb_L_now = (dev_wb_L_interp - dev_wb_steady) / (S3Ratio_now-1);
            devMOD_wb_R_now = (dev_wb_R_interp - dev_wb_steady) / (S3Ratio_now-1);
            devMOD_ds_L_now = (dev_ds_L_interp - dev_ds_steady) / (S3Ratio_now-1);
            devMOD_ds_R_now = (dev_ds_R_interp - dev_ds_steady) / (S3Ratio_now-1);
            devMOD_us_L_now = (dev_us_L_interp - dev_us_steady) / (S3Ratio_now-1);
            devMOD_us_R_now = (dev_us_R_interp - dev_us_steady) / (S3Ratio_now-1);
            
            DstrokeMOD_wb_now = Dstroke_wb_interp / (S3Ratio_now-1);
            DstrokeMOD_ds_now = Dstroke_ds_interp / (S3Ratio_now-1);
            DstrokeMOD_us_now = Dstroke_us_interp / (S3Ratio_now-1);

            DpitchMOD_wb_now = Dpitch_wb_interp / (S3Ratio_now-1);
            DpitchMOD_ds_now = Dpitch_ds_interp / (S3Ratio_now-1);
            DpitchMOD_us_now = Dpitch_us_interp / (S3Ratio_now-1);

            DdevMOD_wb_now = Ddev_wb_interp / (S3Ratio_now-1);
            DdevMOD_ds_now = Ddev_ds_interp / (S3Ratio_now-1);
            DdevMOD_us_now = Ddev_us_interp / (S3Ratio_now-1);
            

            %% store data
            n_now=n_now+1;

            seq_nr_dS3(n_now,1) = seq_nr_now;
            wb_nr_dS3(n_now,1) = wb_nr_now;

            % clip data
            clip_side_dS3(n_now,1) = clip_side_now;
            clip_type_dS3(n_now,1) = clip_type_now;

            S2Clipped_dS3(n_now,1) = S2Clipped_now;
            S2Intact_dS3(n_now,1) = S2Intact_now;
            S2Ratio_dS3(n_now,1) = S2Ratio_now;

            S3Clipped_dS3(n_now,1) = S3Clipped_now;
            S3Intact_dS3(n_now,1) = S3Intact_now;
            S3Ratio_dS3(n_now,1) = S3Ratio_now;

            AreaClipped_dS3(n_now,1) = AreaClipped_now;
            AreaIntact_dS3(n_now,1) = AreaIntact_now;
            AreaRatio_dS3(n_now,1) = AreaRatio_now;

            LengthClipped_dS3(n_now,1) = LengthClipped_now;
            LengthIntact_dS3(n_now,1) = LengthIntact_now;
            LengthRatio_dS3(n_now,1) = LengthRatio_now;

            % body kin
            vel_dS3(n_now,1) = vel_now;

            slip_dS3(n_now,1) = slip_now;
            pitch_dS3(n_now,1) = pitch_now;
            roll_dS3(n_now,1) = roll_now;

            Fsp_pitch_dS3(n_now,1) = Fsp_pitch_now;
            Fsp_roll_dS3(n_now,1) = Fsp_roll_now;

            % wingbeat kin
            dt_ds_dS3(n_now,1) = nanmean(dt_ds_now);
            dt_us_dS3(n_now,1) = nanmean(dt_us_now);
            f_wb_dS3(n_now,1) = nanmean(f_wb_now);
            Rds_dS3(n_now,1) = nanmean(Rds_now);

            t_wb_L_dS3(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R_dS3(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L_dS3(1:length(t_wb_L_now),n_now) = stroke_wb_L_now;
            stroke_wb_R_dS3(1:length(t_wb_R_now),n_now) = stroke_wb_R_now;
            pitch_wb_L_dS3(1:length(t_wb_L_now),n_now) = pitch_wb_L_now;
            pitch_wb_R_dS3(1:length(t_wb_R_now),n_now) = pitch_wb_R_now;
            dev_wb_L_dS3(1:length(t_wb_L_now),n_now) = dev_wb_L_now;
            dev_wb_R_dS3(1:length(t_wb_R_now),n_now) = dev_wb_R_now;
            aoa_wb_L_dS3(1:length(t_wb_L_now),n_now) = aoa_wb_L_now;
            aoa_wb_R_dS3(1:length(t_wb_R_now),n_now) = aoa_wb_R_now;
            U_wb_L_dS3(1:length(t_wb_L_now),n_now) = U_wb_L_now;
            U_wb_R_dS3(1:length(t_wb_R_now),n_now) = U_wb_R_now;
            
            t_ds_L_dS3(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R_dS3(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L_dS3(1:length(t_ds_L_now),n_now) = stroke_ds_L_now;
            stroke_ds_R_dS3(1:length(t_ds_R_now),n_now) = stroke_ds_R_now;
            pitch_ds_L_dS3(1:length(t_ds_L_now),n_now) = pitch_ds_L_now;
            pitch_ds_R_dS3(1:length(t_ds_R_now),n_now) = pitch_ds_R_now;
            dev_ds_L_dS3(1:length(t_ds_L_now),n_now) = dev_ds_L_now;
            dev_ds_R_dS3(1:length(t_ds_R_now),n_now) = dev_ds_R_now;
            aoa_ds_L_dS3(1:length(t_ds_L_now),n_now) = aoa_ds_L_now;
            aoa_ds_R_dS3(1:length(t_ds_R_now),n_now) = aoa_ds_R_now;
            U_ds_L_dS3(1:length(t_ds_L_now),n_now) = U_ds_L_now;
            U_ds_R_dS3(1:length(t_ds_R_now),n_now) = U_ds_R_now;
            
            t_us_L_dS3(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R_dS3(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L_dS3(1:length(t_us_L_now),n_now) = stroke_us_L_now;
            stroke_us_R_dS3(1:length(t_us_R_now),n_now) = stroke_us_R_now;
            pitch_us_L_dS3(1:length(t_us_L_now),n_now) = pitch_us_L_now;
            pitch_us_R_dS3(1:length(t_us_R_now),n_now) = pitch_us_R_now;
            dev_us_L_dS3(1:length(t_us_L_now),n_now) = dev_us_L_now;
            dev_us_R_dS3(1:length(t_us_R_now),n_now) = dev_us_R_now;
            aoa_us_L_dS3(1:length(t_us_L_now),n_now) = aoa_us_L_now;
            aoa_us_R_dS3(1:length(t_us_R_now),n_now) = aoa_us_R_now;
            U_us_L_dS3(1:length(t_us_L_now),n_now) = U_us_L_now;
            U_us_R_dS3(1:length(t_us_R_now),n_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_dS3_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_dS3_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_dS3_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_dS3_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_dS3_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_dS3_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_dS3_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_dS3_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_dS3_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_dS3_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_dS3_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_dS3_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_dS3_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_dS3_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_dS3_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_dS3_bins(:,n_now) = DU_wb_interp;
            
            t_ds_dS3_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_dS3_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_dS3_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_dS3_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_dS3_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_dS3_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_dS3_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_dS3_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_dS3_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_dS3_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_dS3_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_dS3_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_dS3_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_dS3_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_dS3_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_dS3_bins(:,n_now) = DU_ds_interp;
            
            t_us_dS3_bins(:,n_now) = t_us_bin;
            stroke_us_L_dS3_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_dS3_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_dS3_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_dS3_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_dS3_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_dS3_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_dS3_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_dS3_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_dS3_bins(:,n_now) = U_us_L_interp;
            U_us_R_dS3_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_dS3_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_dS3_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_dS3_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_dS3_bins(:,n_now) = Daoa_us_interp;
            DU_us_dS3_bins(:,n_now) = DU_us_interp;
            
            %% store MOD data
            
            %  bodykin modifications
            velMOD_dS3(n_now,1) = velMOD_now;
            
            slipMOD_dS3(n_now,1) = slipMOD_now;
            pitchMOD_dS3(n_now,1) = pitchMOD_now;
            rollMOD_dS3(n_now,1) = rollMOD_now;

            Fsp_pitchMOD_dS3(n_now,1) = Fsp_pitchMOD_now;
            Fsp_rollMOD_dS3(n_now,1) = Fsp_rollMOD_now;

            %  WB modifications
            freqMOD_wb_dS3(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_L_dS3_bins(:,n_now) = strokeMOD_wb_L_now;
            strokeMOD_wb_R_dS3_bins(:,n_now) = strokeMOD_wb_R_now;
            strokeMOD_ds_L_dS3_bins(:,n_now) = strokeMOD_ds_L_now;
            strokeMOD_ds_R_dS3_bins(:,n_now) = strokeMOD_ds_R_now;
            strokeMOD_us_L_dS3_bins(:,n_now) = strokeMOD_us_L_now;
            strokeMOD_us_R_dS3_bins(:,n_now) = strokeMOD_us_R_now;

            pitchMOD_wb_L_dS3_bins(:,n_now) = pitchMOD_wb_L_now;
            pitchMOD_wb_R_dS3_bins(:,n_now) = pitchMOD_wb_R_now;
            pitchMOD_ds_L_dS3_bins(:,n_now) = pitchMOD_ds_L_now;
            pitchMOD_ds_R_dS3_bins(:,n_now) = pitchMOD_ds_R_now;
            pitchMOD_us_L_dS3_bins(:,n_now) = pitchMOD_us_L_now;
            pitchMOD_us_R_dS3_bins(:,n_now) = pitchMOD_us_R_now;

            devMOD_wb_L_dS3_bins(:,n_now) = devMOD_wb_L_now;
            devMOD_wb_R_dS3_bins(:,n_now) = devMOD_wb_R_now;
            devMOD_ds_L_dS3_bins(:,n_now) = devMOD_ds_L_now;
            devMOD_ds_R_dS3_bins(:,n_now) = devMOD_ds_R_now;
            devMOD_us_L_dS3_bins(:,n_now) = devMOD_us_L_now;
            devMOD_us_R_dS3_bins(:,n_now) = devMOD_us_R_now;

            DstrokeMOD_wb_dS3_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_dS3_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_dS3_bins(:,n_now) = DstrokeMOD_us_now;

            DpitchMOD_wb_dS3_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_dS3_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_dS3_bins(:,n_now) = DpitchMOD_us_now;

            DdevMOD_wb_dS3_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_dS3_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_dS3_bins(:,n_now) = DdevMOD_us_now;
        end
end

%% mean & 95%CI
% body kin
vel_dS3_meanCIstd = [nanmean(vel_dS3) 1.96*nanstd(vel_dS3)/sqrt(length(vel_dS3)) nanstd(vel_dS3)];

slip_dS3_meanCIstd = [nanmean(slip_dS3) 1.96*nanstd(slip_dS3)/sqrt(length(slip_dS3)) nanstd(slip_dS3)];
pitch_dS3_meanCIstd = [nanmean(pitch_dS3) 1.96*nanstd(pitch_dS3)/sqrt(length(pitch_dS3)) nanstd(pitch_dS3)];
roll_dS3_meanCIstd = [nanmean(roll_dS3) 1.96*nanstd(roll_dS3)/sqrt(length(roll_dS3)) nanstd(roll_dS3)];

Fsp_pitch_dS3_meanCIstd = [nanmean(Fsp_pitch_dS3) 1.96*nanstd(Fsp_pitch_dS3)/sqrt(length(Fsp_pitch_dS3)) nanstd(Fsp_pitch_dS3)];
Fsp_roll_dS3_meanCIstd = [nanmean(Fsp_roll_dS3) 1.96*nanstd(Fsp_roll_dS3)/sqrt(length(Fsp_roll_dS3)) nanstd(Fsp_roll_dS3)];

% WB kin
dt_ds_dS3_meanCIstd = [nanmean(dt_ds_dS3) 1.96*nanstd(dt_ds_dS3)/sqrt(length(dt_ds_dS3)) nanstd(dt_ds_dS3)];
dt_us_dS3_meanCIstd = [nanmean(dt_us_dS3) 1.96*nanstd(dt_us_dS3)/sqrt(length(dt_us_dS3)) nanstd(dt_us_dS3)];
f_wb_dS3_meanCIstd = [nanmean(f_wb_dS3) 1.96*nanstd(f_wb_dS3)/sqrt(length(f_wb_dS3)) nanstd(f_wb_dS3)];
Rds_dS3_meanCIstd = [nanmean(Rds_dS3) 1.96*nanstd(Rds_dS3)/sqrt(length(Rds_dS3)) nanstd(Rds_dS3)];

calc_WBfunc_dS3_circmeanCIstd

%% wbMOD means & 95%CI
velMOD_dS3_meanCIstd = [nanmean(velMOD_dS3) 1.96*nanstd(velMOD_dS3)/sqrt(length(velMOD_dS3)) nanstd(velMOD_dS3)];

slipMOD_dS3_meanCIstd = [nanmean(slipMOD_dS3) 1.96*nanstd(slipMOD_dS3)/sqrt(length(slipMOD_dS3)) nanstd(slipMOD_dS3)];
pitchMOD_dS3_meanCIstd = [nanmean(pitchMOD_dS3) 1.96*nanstd(pitchMOD_dS3)/sqrt(length(pitchMOD_dS3)) nanstd(pitchMOD_dS3)];
rollMOD_dS3_meanCIstd = [nanmean(rollMOD_dS3) 1.96*nanstd(rollMOD_dS3)/sqrt(length(rollMOD_dS3)) nanstd(rollMOD_dS3)];

Fsp_pitchMOD_dS3_meanCIstd = [nanmean(Fsp_pitchMOD_dS3) 1.96*nanstd(Fsp_pitchMOD_dS3)/sqrt(length(Fsp_pitchMOD_dS3)) nanstd(Fsp_pitchMOD_dS3)];
Fsp_rollMOD_dS3_meanCIstd = [nanmean(Fsp_rollMOD_dS3) 1.96*nanstd(Fsp_rollMOD_dS3)/sqrt(length(Fsp_rollMOD_dS3)) nanstd(Fsp_rollMOD_dS3)];

freqMOD_wb_dS3_meanCIstd = [nanmean(freqMOD_wb_dS3) 1.96*nanstd(freqMOD_wb_dS3)/sqrt(length(freqMOD_wb_dS3)) nanstd(freqMOD_wb_dS3)];

calc_WBmod_dS3_circmeanCIstd

%% WBfits
    t_loc = t_wb_dS3_bins(:,1);
    Rds_loc = Rds_dS3_meanCIstd(1);

%% fit for LEFT wing (bins)
    stroke_loc = stroke_wb_L_dS3_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_L_dS3_bins_meanCIstd(:,1);
    dev_loc = dev_wb_L_dS3_bins_meanCIstd(:,1);
    
    strokeMOD_loc = strokeMOD_wb_L_dS3_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_L_dS3_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_L_dS3_bins_meanCIstd(:,1);

        %% legendre polynomials
        % wingbeats
            [stroke_L_dS3_fit_binmean, stroke_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
            [pitch_L_dS3_fit_binmean, pitch_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
            [dev_L_dS3_fit_binmean, dev_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

        % wingbeatMODs
            [strokeMOD_L_dS3_fit_binmean, strokeMOD_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
            [pitchMOD_L_dS3_fit_binmean, pitchMOD_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
            [devMOD_L_dS3_fit_binmean, devMOD_L_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

        %% fourier series
            % wingbeats
            [stroke_L_dS3_fourier_fit_binmean, stroke_L_dS3_fourier_gof_binmean,stroke_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
            [pitch_L_dS3_fourier_fit_binmean, pitch_L_dS3_fourier_gof_binmean,pitch_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
            [dev_L_dS3_fourier_fit_binmean, dev_L_dS3_fourier_gof_binmean,dev_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);

            % wingbeatMODs
            [strokeMOD_L_dS3_fourier_fit_binmean, strokeMOD_L_dS3_fourier_gof_binmean,strokeMOD_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
            [pitchMOD_L_dS3_fourier_fit_binmean, pitchMOD_L_dS3_fourier_gof_binmean,pitchMOD_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
            [devMOD_L_dS3_fourier_fit_binmean, devMOD_L_dS3_fourier_gof_binmean,devMOD_L_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
%% fit for Right wing (bins)
    stroke_loc = stroke_wb_R_dS3_bins_meanCIstd(:,1);
    pitch_loc = pitch_wb_R_dS3_bins_meanCIstd(:,1);
    dev_loc = dev_wb_R_dS3_bins_meanCIstd(:,1);

    strokeMOD_loc = strokeMOD_wb_R_dS3_bins_meanCIstd(:,1);
    pitchMOD_loc = pitchMOD_wb_R_dS3_bins_meanCIstd(:,1);
    devMOD_loc = devMOD_wb_R_dS3_bins_meanCIstd(:,1);

        %% legendre polynomials
        % wingbeats
            [stroke_R_dS3_fit_binmean, stroke_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
            [pitch_R_dS3_fit_binmean, pitch_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
            [dev_R_dS3_fit_binmean, dev_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

        % wingbeatMODs
            [strokeMOD_R_dS3_fit_binmean, strokeMOD_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
            [pitchMOD_R_dS3_fit_binmean, pitchMOD_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
            [devMOD_R_dS3_fit_binmean, devMOD_R_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

        %% fourier series
            % wingbeats
            [stroke_R_dS3_fourier_fit_binmean, stroke_R_dS3_fourier_gof_binmean,stroke_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
            [pitch_R_dS3_fourier_fit_binmean, pitch_R_dS3_fourier_gof_binmean,pitch_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
            [dev_R_dS3_fourier_fit_binmean, dev_R_dS3_fourier_gof_binmean,dev_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);

            % wingbeatMODs
            [strokeMOD_R_dS3_fourier_fit_binmean, strokeMOD_R_dS3_fourier_gof_binmean,strokeMOD_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
            [pitchMOD_R_dS3_fourier_fit_binmean, pitchMOD_R_dS3_fourier_gof_binmean,pitchMOD_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
            [devMOD_R_dS3_fourier_fit_binmean, devMOD_R_dS3_fourier_gof_binmean,devMOD_R_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
%% L-R
    stroke_loc = Dstroke_wb_dS3_bins_meanCIstd(:,1);
    pitch_loc = Dpitch_wb_dS3_bins_meanCIstd(:,1);
    dev_loc = Ddev_wb_dS3_bins_meanCIstd(:,1);

    strokeMOD_loc = DstrokeMOD_wb_dS3_bins_meanCIstd(:,1);
    pitchMOD_loc = DpitchMOD_wb_dS3_bins_meanCIstd(:,1);
    devMOD_loc = DdevMOD_wb_dS3_bins_meanCIstd(:,1);
    
        %% legendre polynomials
            % wingbeats
            [Dstroke_dS3_fit_binmean, Dstroke_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
            [Dpitch_dS3_fit_binmean, Dpitch_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
            [Ddev_dS3_fit_binmean, Ddev_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

        % wingbeatMODs
            [DstrokeMOD_dS3_fit_binmean, DstrokeMOD_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
            [DpitchMOD_dS3_fit_binmean, DpitchMOD_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
            [DdevMOD_dS3_fit_binmean, DdevMOD_dS3_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

        %% fourier series
            % wingbeats
            [Dstroke_dS3_fourier_fit_binmean, Dstroke_dS3_fourier_gof_binmean,Dstroke_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,MOD_fourier_order,plot_fourier);
            [Dpitch_dS3_fourier_fit_binmean, Dpitch_dS3_fourier_gof_binmean,Dpitch_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,MOD_fourier_order,plot_fourier);
            [Ddev_dS3_fourier_fit_binmean, Ddev_dS3_fourier_gof_binmean,Ddev_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,MOD_fourier_order,plot_fourier);

            % wingbeatMODs
            [DstrokeMOD_dS3_fourier_fit_binmean, DstrokeMOD_dS3_fourier_gof_binmean,DstrokeMOD_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
            [DpitchMOD_dS3_fourier_fit_binmean, DpitchMOD_dS3_fourier_gof_binmean,DpitchMOD_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
            [DdevMOD_dS3_fourier_fit_binmean, DdevMOD_dS3_fourier_gof_binmean,DdevMOD_dS3_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);

%% save CLIPd wb data
n_dS3 = n_now   

if save_on == 1

    save(['WBmod_dS3_',num2str(n_dS3),'WBs.mat'],...
'n_dS3',...
'seq_nr_dS3',...
'wb_nr_dS3',...
'S3Ratio_cut',...
'S3Ratio_plot',...
...
'n_pol_stroke',...
'n_pol_pitch',...
'n_pol_dev',...
'n_pol_MOD',...
...
'clip_side_dS3',...
'clip_type_dS3',...
...
'S2Clipped_dS3',...
'S2Intact_dS3',...
'S2Ratio_dS3',...
...
'S3Clipped_dS3',...
'S3Intact_dS3',...
'S3Ratio_dS3',...
...
'AreaClipped_dS3',...
'AreaIntact_dS3',...
'AreaRatio_dS3',...
...
'LengthClipped_dS3',...
'LengthIntact_dS3',...
'LengthRatio_dS3',...
...
'stroke_R_dS3_fit_binmean_periodic',...
'pitch_R_dS3_fit_binmean_periodic',...
'dev_R_dS3_fit_binmean_periodic',...
...
'stroke_L_dS3_fit_binmean_periodic',...
'pitch_L_dS3_fit_binmean_periodic',...
'dev_L_dS3_fit_binmean_periodic',...
...
'strokeMOD_R_dS3_fit_binmean_periodic',...
'pitchMOD_R_dS3_fit_binmean_periodic',...
'devMOD_R_dS3_fit_binmean_periodic',...
...
'strokeMOD_L_dS3_fit_binmean_periodic',...
'pitchMOD_L_dS3_fit_binmean_periodic',...
'devMOD_L_dS3_fit_binmean_periodic',...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_R_dS3_fourier_coeffs_binmean',...
'pitch_R_dS3_fourier_coeffs_binmean',...
'dev_R_dS3_fourier_coeffs_binmean',...
...
'stroke_R_dS3_fourier_gof_binmean',...
'pitch_R_dS3_fourier_gof_binmean',...
'dev_R_dS3_fourier_gof_binmean',...
...
'stroke_L_dS3_fourier_coeffs_binmean',...
'pitch_L_dS3_fourier_coeffs_binmean',...
'dev_L_dS3_fourier_coeffs_binmean',...
...
'stroke_L_dS3_fourier_gof_binmean',...
'pitch_L_dS3_fourier_gof_binmean',...
'dev_L_dS3_fourier_gof_binmean',...
...
'Dstroke_dS3_fourier_coeffs_binmean',...
'Dpitch_dS3_fourier_coeffs_binmean',...
'Ddev_dS3_fourier_coeffs_binmean',...
...
'Dstroke_dS3_fourier_gof_binmean',...
'Dpitch_dS3_fourier_gof_binmean',...
'Ddev_dS3_fourier_gof_binmean',...
...
'strokeMOD_R_dS3_fourier_coeffs_binmean',...
'pitchMOD_R_dS3_fourier_coeffs_binmean',...
'devMOD_R_dS3_fourier_coeffs_binmean',...
...
'strokeMOD_R_dS3_fourier_gof_binmean',...
'pitchMOD_R_dS3_fourier_gof_binmean',...
'devMOD_R_dS3_fourier_gof_binmean',...
...
'strokeMOD_L_dS3_fourier_coeffs_binmean',...
'pitchMOD_L_dS3_fourier_coeffs_binmean',...
'devMOD_L_dS3_fourier_coeffs_binmean',...
...
'strokeMOD_L_dS3_fourier_gof_binmean',...
'pitchMOD_L_dS3_fourier_gof_binmean',...
'devMOD_L_dS3_fourier_gof_binmean',...
...
'DstrokeMOD_dS3_fourier_coeffs_binmean',...
'DpitchMOD_dS3_fourier_coeffs_binmean',...
'DdevMOD_dS3_fourier_coeffs_binmean',...
...
'DstrokeMOD_dS3_fourier_gof_binmean',...
'DpitchMOD_dS3_fourier_gof_binmean',...
'DdevMOD_dS3_fourier_gof_binmean',...
...
'stroke_wb_R_dS3_bins_meanCIstd',...
'stroke_ds_R_dS3_bins_meanCIstd',...
'stroke_us_R_dS3_bins_meanCIstd',...
...
'pitch_wb_R_dS3_bins_meanCIstd',...
'pitch_ds_R_dS3_bins_meanCIstd',...
'pitch_us_R_dS3_bins_meanCIstd',...
...
'dev_wb_R_dS3_bins_meanCIstd',...
'dev_ds_R_dS3_bins_meanCIstd',...
'dev_us_R_dS3_bins_meanCIstd',...
...
'strokeMOD_wb_R_dS3_bins_meanCIstd',...
'strokeMOD_ds_R_dS3_bins_meanCIstd',...
'strokeMOD_us_R_dS3_bins_meanCIstd',...
...
'pitchMOD_wb_R_dS3_bins_meanCIstd',...
'pitchMOD_ds_R_dS3_bins_meanCIstd',...
'pitchMOD_us_R_dS3_bins_meanCIstd',...
...
'devMOD_wb_R_dS3_bins_meanCIstd',...
'devMOD_ds_R_dS3_bins_meanCIstd',...
'devMOD_us_R_dS3_bins_meanCIstd',...
...
'stroke_wb_L_dS3_bins_meanCIstd',...
'stroke_ds_L_dS3_bins_meanCIstd',...
'stroke_us_L_dS3_bins_meanCIstd',...
...
'pitch_wb_L_dS3_bins_meanCIstd',...
'pitch_ds_L_dS3_bins_meanCIstd',...
'pitch_us_L_dS3_bins_meanCIstd',...
...
'dev_wb_L_dS3_bins_meanCIstd',...
'dev_ds_L_dS3_bins_meanCIstd',...
'dev_us_L_dS3_bins_meanCIstd',...
...
'strokeMOD_wb_L_dS3_bins_meanCIstd',...
'strokeMOD_ds_L_dS3_bins_meanCIstd',...
'strokeMOD_us_L_dS3_bins_meanCIstd',...
...
'pitchMOD_wb_L_dS3_bins_meanCIstd',...
'pitchMOD_ds_L_dS3_bins_meanCIstd',...
'pitchMOD_us_L_dS3_bins_meanCIstd',...
...
'devMOD_wb_L_dS3_bins_meanCIstd',...
'devMOD_ds_L_dS3_bins_meanCIstd',...
'devMOD_us_L_dS3_bins_meanCIstd',...
...
'DstrokeMOD_wb_dS3_bins_meanCIstd',...
'DstrokeMOD_ds_dS3_bins_meanCIstd',...
'DstrokeMOD_us_dS3_bins_meanCIstd',...
...
'DpitchMOD_wb_dS3_bins_meanCIstd',...
'DpitchMOD_ds_dS3_bins_meanCIstd',...
'DpitchMOD_us_dS3_bins_meanCIstd',...
...
'DdevMOD_wb_dS3_bins_meanCIstd',...
'DdevMOD_ds_dS3_bins_meanCIstd',...
'DdevMOD_us_dS3_bins_meanCIstd',...
...
'vel_dS3',...
'slip_dS3',...
'pitch_dS3',...
'roll_dS3',...
'Fsp_pitch_dS3',...
'Fsp_roll_dS3',...
...
'dt_ds_dS3',...
'dt_us_dS3',...
'f_wb_dS3',...
'Rds_dS3',...
...
'vel_dS3_meanCIstd',...
'slip_dS3_meanCIstd',...
'pitch_dS3_meanCIstd',...
'roll_dS3_meanCIstd',...
'Fsp_pitch_dS3_meanCIstd',...
'Fsp_roll_dS3_meanCIstd',...
...
'dt_ds_dS3_meanCIstd',...
'dt_us_dS3_meanCIstd',...
'f_wb_dS3_meanCIstd',...
'Rds_dS3_meanCIstd',...
...
'velMOD_dS3_meanCIstd',...
'velMOD_dS3',...
...
'slipMOD_dS3_meanCIstd',...
'slipMOD_dS3',...
...
'pitchMOD_dS3_meanCIstd',...
'pitchMOD_dS3',...
...
'rollMOD_dS3_meanCIstd',...
'rollMOD_dS3',...
...
'Fsp_pitchMOD_dS3_meanCIstd',...
'Fsp_pitchMOD_dS3',...
...
'Fsp_rollMOD_dS3_meanCIstd',...
'Fsp_rollMOD_dS3',...
...
'freqMOD_wb_dS3_meanCIstd',...
'freqMOD_wb_dS3',...
...
'strokeMOD_wb_L_dS3_bins',...
'strokeMOD_wb_R_dS3_bins',...
'strokeMOD_ds_L_dS3_bins',...
'strokeMOD_ds_R_dS3_bins',...
'strokeMOD_us_L_dS3_bins',...
'strokeMOD_us_R_dS3_bins',...
...
'pitchMOD_wb_L_dS3_bins',...
'pitchMOD_wb_R_dS3_bins',...
'pitchMOD_ds_L_dS3_bins',...
'pitchMOD_ds_R_dS3_bins',...
'pitchMOD_us_L_dS3_bins',...
'pitchMOD_us_R_dS3_bins',...
...
'devMOD_wb_L_dS3_bins',...
'devMOD_wb_R_dS3_bins',...
'devMOD_ds_L_dS3_bins',...
'devMOD_ds_R_dS3_bins',...
'devMOD_us_L_dS3_bins',...
'devMOD_us_R_dS3_bins',...
...
'DstrokeMOD_wb_dS3_bins',...
'DstrokeMOD_ds_dS3_bins',...
'DstrokeMOD_us_dS3_bins',...
...
'DpitchMOD_wb_dS3_bins',...
'DpitchMOD_ds_dS3_bins',...
'DpitchMOD_us_dS3_bins',...
...
'DdevMOD_wb_dS3_bins',...
'DdevMOD_ds_dS3_bins',...
'DdevMOD_us_dS3_bins',...
...
't_wb_L_dS3',...
't_wb_R_dS3',...
'stroke_wb_L_dS3',...
'stroke_wb_R_dS3',...
'pitch_wb_L_dS3',...
'pitch_wb_R_dS3',...
'dev_wb_L_dS3',...
'dev_wb_R_dS3',...
'aoa_wb_L_dS3',...
'aoa_wb_R_dS3',...
'U_wb_L_dS3',...
'U_wb_R_dS3',...
...
't_ds_L_dS3',...
't_ds_R_dS3',...
'stroke_ds_L_dS3',...
'stroke_ds_R_dS3',...
'pitch_ds_L_dS3',...
'pitch_ds_R_dS3',...
'dev_ds_L_dS3',...
'dev_ds_R_dS3',...
'aoa_ds_L_dS3',...
'aoa_ds_R_dS3',...
'U_ds_L_dS3',...
'U_ds_R_dS3',...
...
't_us_L_dS3',...
't_us_R_dS3',...
'stroke_us_L_dS3',...
'stroke_us_R_dS3',...
'pitch_us_L_dS3',...
'pitch_us_R_dS3',...
'dev_us_L_dS3',...
'dev_us_R_dS3',...
'aoa_us_L_dS3',...
'aoa_us_R_dS3',...
'U_us_L_dS3',...
'U_us_R_dS3',...
...
't_wb_dS3_bins',...
'stroke_wb_L_dS3_bins',...
'stroke_wb_R_dS3_bins',...
'pitch_wb_L_dS3_bins',...
'pitch_wb_R_dS3_bins',...
'dev_wb_L_dS3_bins',...
'dev_wb_R_dS3_bins',...
'aoa_wb_L_dS3_bins',...
'aoa_wb_R_dS3_bins',...
'U_wb_L_dS3_bins',...
'U_wb_R_dS3_bins',...
'Dstroke_wb_dS3_bins',...
'Dpitch_wb_dS3_bins',...
'Ddev_wb_dS3_bins',...
'Daoa_wb_dS3_bins',...
'DU_wb_dS3_bins',...
...
't_ds_dS3_bins',...
'stroke_ds_L_dS3_bins',...
'stroke_ds_R_dS3_bins',...
'pitch_ds_L_dS3_bins',...
'pitch_ds_R_dS3_bins',...
'dev_ds_L_dS3_bins',...
'dev_ds_R_dS3_bins',...
'aoa_ds_L_dS3_bins',...
'aoa_ds_R_dS3_bins',...
'U_ds_L_dS3_bins',...
'U_ds_R_dS3_bins',...
'Dstroke_ds_dS3_bins',...
'Dpitch_ds_dS3_bins',...
'Ddev_ds_dS3_bins',...
'Daoa_ds_dS3_bins',...
'DU_ds_dS3_bins',...
...
't_us_dS3_bins',...
'stroke_us_L_dS3_bins',...
'stroke_us_R_dS3_bins',...
'pitch_us_L_dS3_bins',...
'pitch_us_R_dS3_bins',...
'dev_us_L_dS3_bins',...
'dev_us_R_dS3_bins',...
'aoa_us_L_dS3_bins',...
'aoa_us_R_dS3_bins',...
'U_us_L_dS3_bins',...
'U_us_R_dS3_bins',...
'Dstroke_us_dS3_bins',...
'Dpitch_us_dS3_bins',...
'Ddev_us_dS3_bins',...
'Daoa_us_dS3_bins',...
'DU_us_dS3_bins');

end

