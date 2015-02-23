% calc wb MOD CLIP AeroCenterFuncs

n_now=0;

%% variables

% bodyMODs
velMOD_AeroCenterFuncs = [];

slipMOD_AeroCenterFuncs = [];
pitchMOD_AeroCenterFuncs = [];
rollMOD_AeroCenterFuncs = [];

Fsp_pitchMOD_AeroCenterFuncs = [];
Fsp_rollMOD_AeroCenterFuncs = [];

% wbMOD
freqMOD_wb_AeroCenterFuncs = [];

strokeMOD_wb_L_AeroCenterFuncs_bins = [];
strokeMOD_wb_R_AeroCenterFuncs_bins = [];
strokeMOD_ds_L_AeroCenterFuncs_bins = [];
strokeMOD_ds_R_AeroCenterFuncs_bins = [];
strokeMOD_us_L_AeroCenterFuncs_bins = [];
strokeMOD_us_R_AeroCenterFuncs_bins = [];

pitchMOD_wb_L_AeroCenterFuncs_bins = [];
pitchMOD_wb_R_AeroCenterFuncs_bins = [];
pitchMOD_ds_L_AeroCenterFuncs_bins = [];
pitchMOD_ds_R_AeroCenterFuncs_bins = [];
pitchMOD_us_L_AeroCenterFuncs_bins = [];
pitchMOD_us_R_AeroCenterFuncs_bins = [];

devMOD_wb_L_AeroCenterFuncs_bins = [];
devMOD_wb_R_AeroCenterFuncs_bins = [];
devMOD_ds_L_AeroCenterFuncs_bins = [];
devMOD_ds_R_AeroCenterFuncs_bins = [];
devMOD_us_L_AeroCenterFuncs_bins = [];
devMOD_us_R_AeroCenterFuncs_bins = [];

DstrokeMOD_wb_AeroCenterFuncs_bins = [];
DstrokeMOD_ds_AeroCenterFuncs_bins = [];
DstrokeMOD_us_AeroCenterFuncs_bins = [];

DpitchMOD_wb_AeroCenterFuncs_bins = [];
DpitchMOD_ds_AeroCenterFuncs_bins = [];
DpitchMOD_us_AeroCenterFuncs_bins = [];

DdevMOD_wb_AeroCenterFuncs_bins = [];
DdevMOD_ds_AeroCenterFuncs_bins = [];
DdevMOD_us_AeroCenterFuncs_bins = [];

% seq&wb
seq_nr_AeroCenterFuncs = [];
wb_nr_AeroCenterFuncs = [];

% clip data
clip_side_AeroCenterFuncs = [];
clip_type_AeroCenterFuncs = [];

AeroCenterFuncClipped_AeroCenterFuncs = [];
AeroCenterFuncIntact_AeroCenterFuncs = [];
AeroCenterFuncRatio_AeroCenterFuncs = [];

FirstMomentNormClipped_AeroCenterFuncs = [];
FirstMomentNormIntact_AeroCenterFuncs = [];
FirstMomentRatio_AeroCenterFuncs = [];

SecondMomentNormClipped_AeroCenterFuncs = [];
SecondMomentNormIntact_AeroCenterFuncs = [];
SecondMomentRatio_AeroCenterFuncs = [];

ThirdMomentNormClipped_AeroCenterFuncs = [];
ThirdMomentNormIntact_AeroCenterFuncs = [];
ThirdMomentRatio_AeroCenterFuncs = [];

AreaNormClipped_AeroCenterFuncs = [];
AreaNormIntact_AeroCenterFuncs = [];
AreaRatio_AeroCenterFuncs = [];

LengthClipped_pixels_AeroCenterFuncs = [];
LengthIntact_pixels_AeroCenterFuncs = [];
LengthRatio_AeroCenterFuncs = [];

% body kin
vel_AeroCenterFuncs = [];

slip_AeroCenterFuncs = [];
pitch_AeroCenterFuncs = [];
roll_AeroCenterFuncs = [];

Fsp_pitch_AeroCenterFuncs = [];
Fsp_roll_AeroCenterFuncs = [];

% wingbeat kin
dt_wb_AeroCenterFuncs = [];
dt_ds_AeroCenterFuncs = [];
dt_us_AeroCenterFuncs = [];
f_wb_AeroCenterFuncs = [];
Rds_AeroCenterFuncs = [];

t_wb_L_AeroCenterFuncs = [];
t_wb_R_AeroCenterFuncs = [];
stroke_wb_L_AeroCenterFuncs = [];
stroke_wb_R_AeroCenterFuncs = [];
pitch_wb_L_AeroCenterFuncs = [];
pitch_wb_R_AeroCenterFuncs = [];
dev_wb_L_AeroCenterFuncs = [];
dev_wb_R_AeroCenterFuncs = [];
aoa_wb_L_AeroCenterFuncs = [];
aoa_wb_R_AeroCenterFuncs = [];
U_wb_L_AeroCenterFuncs = [];
U_wb_R_AeroCenterFuncs = [];
Dstroke_wb_AeroCenterFuncs = [];
Dpitch_wb_AeroCenterFuncs = [];
Ddev_wb_AeroCenterFuncs = [];
Daoa_wb_AeroCenterFuncs = [];
DU_wb_AeroCenterFuncs = [];

t_ds_L_AeroCenterFuncs = [];
t_ds_R_AeroCenterFuncs = [];
stroke_ds_L_AeroCenterFuncs = [];
stroke_ds_R_AeroCenterFuncs = [];
pitch_ds_L_AeroCenterFuncs = [];
pitch_ds_R_AeroCenterFuncs = [];
dev_ds_L_AeroCenterFuncs = [];
dev_ds_R_AeroCenterFuncs = [];
aoa_ds_L_AeroCenterFuncs = [];
aoa_ds_R_AeroCenterFuncs = [];
U_ds_L_AeroCenterFuncs = [];
U_ds_R_AeroCenterFuncs = [];
Dstroke_ds_AeroCenterFuncs = [];
Dpitch_ds_AeroCenterFuncs = [];
Ddev_ds_AeroCenterFuncs = [];
Daoa_ds_AeroCenterFuncs = [];
DU_ds_AeroCenterFuncs = [];

t_us_L_AeroCenterFuncs = [];
t_us_R_AeroCenterFuncs = [];
stroke_us_L_AeroCenterFuncs = [];
stroke_us_R_AeroCenterFuncs = [];
pitch_us_L_AeroCenterFuncs = [];
pitch_us_R_AeroCenterFuncs = [];
dev_us_L_AeroCenterFuncs = [];
dev_us_R_AeroCenterFuncs = [];
aoa_us_L_AeroCenterFuncs = [];
aoa_us_R_AeroCenterFuncs = [];
U_us_L_AeroCenterFuncs = [];
U_us_R_AeroCenterFuncs = [];
Dstroke_us_AeroCenterFuncs = [];
Dpitch_us_AeroCenterFuncs = [];
Ddev_us_AeroCenterFuncs = [];
Daoa_us_AeroCenterFuncs = [];
DU_us_AeroCenterFuncs = [];

t_wb_AeroCenterFuncs_bins = [];
stroke_wb_L_AeroCenterFuncs_bins = [];
stroke_wb_R_AeroCenterFuncs_bins = [];
pitch_wb_L_AeroCenterFuncs_bins = [];
pitch_wb_R_AeroCenterFuncs_bins = [];
dev_wb_L_AeroCenterFuncs_bins = [];
dev_wb_R_AeroCenterFuncs_bins = [];
aoa_wb_L_AeroCenterFuncs_bins = [];
aoa_wb_R_AeroCenterFuncs_bins = [];
U_wb_L_AeroCenterFuncs_bins = [];
U_wb_R_AeroCenterFuncs_bins = [];
Dstroke_wb_AeroCenterFuncs_bins = [];
Dpitch_wb_AeroCenterFuncs_bins = [];
Ddev_wb_AeroCenterFuncs_bins = [];
Daoa_wb_AeroCenterFuncs_bins = [];
DU_wb_AeroCenterFuncs_bins = [];

t_ds_AeroCenterFuncs_bins = [];
stroke_ds_L_AeroCenterFuncs_bins = [];
stroke_ds_R_AeroCenterFuncs_bins = [];
pitch_ds_L_AeroCenterFuncs_bins = [];
pitch_ds_R_AeroCenterFuncs_bins = [];
dev_ds_L_AeroCenterFuncs_bins = [];
dev_ds_R_AeroCenterFuncs_bins = [];
aoa_ds_L_AeroCenterFuncs_bins = [];
aoa_ds_R_AeroCenterFuncs_bins = [];
U_ds_L_AeroCenterFuncs_bins = [];
U_ds_R_AeroCenterFuncs_bins = [];
Dstroke_ds_AeroCenterFuncs_bins = [];
Dpitch_ds_AeroCenterFuncs_bins = [];
Ddev_ds_AeroCenterFuncs_bins = [];
Daoa_ds_AeroCenterFuncs_bins = [];
DU_ds_AeroCenterFuncs_bins = [];

t_us_AeroCenterFuncs_bins = [];
stroke_us_L_AeroCenterFuncs_bins = [];
stroke_us_R_AeroCenterFuncs_bins = [];
pitch_us_L_AeroCenterFuncs_bins = [];
pitch_us_R_AeroCenterFuncs_bins = [];
dev_us_L_AeroCenterFuncs_bins = [];
dev_us_R_AeroCenterFuncs_bins = [];
aoa_us_L_AeroCenterFuncs_bins = [];
aoa_us_R_AeroCenterFuncs_bins = [];
U_us_L_AeroCenterFuncs_bins = [];
U_us_R_AeroCenterFuncs_bins = [];
Dstroke_us_AeroCenterFuncs_bins = [];
Dpitch_us_AeroCenterFuncs_bins = [];
Ddev_us_AeroCenterFuncs_bins = [];
Daoa_us_AeroCenterFuncs_bins = [];
DU_us_AeroCenterFuncs_bins = [];

%% filter wingbeats based on FirstMomentRatio
for wb = 1:size(wb_nr,1)
    counter = size(wb_nr,1) -wb
        
        if steady_nr_mean_wb(wb) == 1 && FirstMomentRatio(wb) < FirstMomentRatio_cut
            
            % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);
            
            % clip data
            clip_side_now = clip_side(wb);
            clip_type_now = clip_type(wb);

            AeroCenterFuncClipped_now = AeroCenterFuncClipped(wb);
            AeroCenterFuncIntact_now = AeroCenterFuncIntact(wb);
            AeroCenterFuncRatio_now = AeroCenterFuncRatio(wb);

            FirstMomentNormClipped_now = FirstMomentNormClipped(wb);
            FirstMomentNormIntact_now = FirstMomentNormIntact(wb);
            FirstMomentRatio_now = FirstMomentRatio(wb);

            SecondMomentNormClipped_now = SecondMomentNormClipped(wb);
            SecondMomentNormIntact_now = SecondMomentNormIntact(wb);
            SecondMomentRatio_now = SecondMomentRatio(wb);

            ThirdMomentNormClipped_now = ThirdMomentNormClipped(wb);
            ThirdMomentNormIntact_now = ThirdMomentNormIntact(wb);
            ThirdMomentRatio_now = ThirdMomentRatio(wb);

            AreaNormClipped_now = AreaNormClipped(wb);
            AreaNormIntact_now = AreaNormIntact(wb);
            AreaRatio_now = AreaRatio(wb);

            LengthClipped_pixels_now = LengthClipped_pixels(wb);
            LengthIntact_pixels_now = LengthIntact_pixels(wb);
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
            
            
            %% calc bodykin modifications: AeroCenterFuncRatio
            velMOD_now   = (nanmean(vel_now) - vel_steady)          / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);

            slipMOD_now  =  nanmean(slip_now)                       / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            pitchMOD_now = (nanmean(pitch_now) - pitch_body_steady) / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            rollMOD_now  =  nanmean(roll_now)                       / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);

            Fsp_pitchMOD_now = nanmean(Fsp_pitch_now)               / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            Fsp_rollMOD_now  = nanmean(Fsp_roll_now)                / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);

            %% calc wingbeat frequency modifications: AeroCenterFuncRatio
            freqMOD_wb_now =     (nanmean(f_wb_now)  - f_wb_steady)      / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            
            %% Left wing = intact wing: correlate with AeroCenterFuncIntact
            strokeMOD_wb_L_now = (stroke_wb_L_interp - stroke_wb_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            strokeMOD_ds_L_now = (stroke_ds_L_interp - stroke_ds_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            strokeMOD_us_L_now = (stroke_us_L_interp - stroke_us_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            
            pitchMOD_wb_L_now = (pitch_wb_L_interp - pitch_wb_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            pitchMOD_ds_L_now = (pitch_ds_L_interp - pitch_ds_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            pitchMOD_us_L_now = (pitch_us_L_interp - pitch_us_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            
            devMOD_wb_L_now = (dev_wb_L_interp - dev_wb_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            devMOD_ds_L_now = (dev_ds_L_interp - dev_ds_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            devMOD_us_L_now = (dev_us_L_interp - dev_us_steady) / (AeroCenterFuncIntact_now-AeroCenterFuncIntact_steady);
            
            %% Right wing = clipped wing: correlate with AeroCenterFuncClipped
            strokeMOD_wb_R_now = (stroke_wb_R_interp - stroke_wb_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            strokeMOD_ds_R_now = (stroke_ds_R_interp - stroke_ds_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            strokeMOD_us_R_now = (stroke_us_R_interp - stroke_us_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            
            pitchMOD_wb_R_now = (pitch_wb_R_interp - pitch_wb_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            pitchMOD_ds_R_now = (pitch_ds_R_interp - pitch_ds_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            pitchMOD_us_R_now = (pitch_us_R_interp - pitch_us_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            
            devMOD_wb_R_now = (dev_wb_R_interp - dev_wb_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            devMOD_ds_R_now = (dev_ds_R_interp - dev_ds_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            devMOD_us_R_now = (dev_us_R_interp - dev_us_steady) / (AeroCenterFuncClipped_now-AeroCenterFuncClipped_steady);
            
            %% Left - Right wingbeat kin: correlate with AeroCenterFuncRatio
            DstrokeMOD_wb_now = Dstroke_wb_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DstrokeMOD_ds_now = Dstroke_ds_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DstrokeMOD_us_now = Dstroke_us_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);

            DpitchMOD_wb_now = Dpitch_wb_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DpitchMOD_ds_now = Dpitch_ds_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DpitchMOD_us_now = Dpitch_us_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);

            DdevMOD_wb_now = Ddev_wb_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DdevMOD_ds_now = Ddev_ds_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            DdevMOD_us_now = Ddev_us_interp / (AeroCenterFuncRatio_now-AeroCenterFuncRatio_steady);
            
            %% store data
            n_now=n_now+1;

            seq_nr_AeroCenterFuncs(n_now,1) = seq_nr_now;
            wb_nr_AeroCenterFuncs(n_now,1) = wb_nr_now;

            % clip data
            clip_side_AeroCenterFuncs(n_now,1) = clip_side_now;
            clip_type_AeroCenterFuncs(n_now,1) = clip_type_now;

            AeroCenterFuncClipped_AeroCenterFuncs(n_now,1) = AeroCenterFuncClipped_now;
            AeroCenterFuncIntact_AeroCenterFuncs(n_now,1) = AeroCenterFuncIntact_now;
            AeroCenterFuncRatio_AeroCenterFuncs(n_now,1) = AeroCenterFuncRatio_now;

            FirstMomentNormClipped_AeroCenterFuncs(n_now,1) = FirstMomentNormClipped_now;
            FirstMomentNormIntact_AeroCenterFuncs(n_now,1) = FirstMomentNormIntact_now;
            FirstMomentRatio_AeroCenterFuncs(n_now,1) = FirstMomentRatio_now;

            SecondMomentNormClipped_AeroCenterFuncs(n_now,1) = SecondMomentNormClipped_now;
            SecondMomentNormIntact_AeroCenterFuncs(n_now,1) = SecondMomentNormIntact_now;
            SecondMomentRatio_AeroCenterFuncs(n_now,1) = SecondMomentRatio_now;

            ThirdMomentNormClipped_AeroCenterFuncs(n_now,1) = ThirdMomentNormClipped_now;
            ThirdMomentNormIntact_AeroCenterFuncs(n_now,1) = ThirdMomentNormIntact_now;
            ThirdMomentRatio_AeroCenterFuncs(n_now,1) = ThirdMomentRatio_now;

            AreaNormClipped_AeroCenterFuncs(n_now,1) = AreaNormClipped_now;
            AreaNormIntact_AeroCenterFuncs(n_now,1) = AreaNormIntact_now;
            AreaRatio_AeroCenterFuncs(n_now,1) = AreaRatio_now;

            LengthClipped_pixels_AeroCenterFuncs(n_now,1) = LengthClipped_pixels_now;
            LengthIntact_pixels_AeroCenterFuncs(n_now,1) = LengthIntact_pixels_now;
            LengthRatio_AeroCenterFuncs(n_now,1) = LengthRatio_now;

            % body kin
            vel_AeroCenterFuncs(n_now,1) = vel_now;

            slip_AeroCenterFuncs(n_now,1) = slip_now;
            pitch_AeroCenterFuncs(n_now,1) = pitch_now;
            roll_AeroCenterFuncs(n_now,1) = roll_now;

            Fsp_pitch_AeroCenterFuncs(n_now,1) = Fsp_pitch_now;
            Fsp_roll_AeroCenterFuncs(n_now,1) = Fsp_roll_now;

            % wingbeat kin
            dt_ds_AeroCenterFuncs(n_now,1) = nanmean(dt_ds_now);
            dt_us_AeroCenterFuncs(n_now,1) = nanmean(dt_us_now);
            f_wb_AeroCenterFuncs(n_now,1) = nanmean(f_wb_now);
            Rds_AeroCenterFuncs(n_now,1) = nanmean(Rds_now);

            t_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = stroke_wb_L_now;
            stroke_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = stroke_wb_R_now;
            pitch_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = pitch_wb_L_now;
            pitch_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = pitch_wb_R_now;
            dev_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = dev_wb_L_now;
            dev_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = dev_wb_R_now;
            aoa_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = aoa_wb_L_now;
            aoa_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = aoa_wb_R_now;
            U_wb_L_AeroCenterFuncs(1:length(t_wb_L_now),n_now) = U_wb_L_now;
            U_wb_R_AeroCenterFuncs(1:length(t_wb_R_now),n_now) = U_wb_R_now;
            
            t_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = stroke_ds_L_now;
            stroke_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = stroke_ds_R_now;
            pitch_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = pitch_ds_L_now;
            pitch_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = pitch_ds_R_now;
            dev_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = dev_ds_L_now;
            dev_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = dev_ds_R_now;
            aoa_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = aoa_ds_L_now;
            aoa_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = aoa_ds_R_now;
            U_ds_L_AeroCenterFuncs(1:length(t_ds_L_now),n_now) = U_ds_L_now;
            U_ds_R_AeroCenterFuncs(1:length(t_ds_R_now),n_now) = U_ds_R_now;
            
            t_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = stroke_us_L_now;
            stroke_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = stroke_us_R_now;
            pitch_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = pitch_us_L_now;
            pitch_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = pitch_us_R_now;
            dev_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = dev_us_L_now;
            dev_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = dev_us_R_now;
            aoa_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = aoa_us_L_now;
            aoa_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = aoa_us_R_now;
            U_us_L_AeroCenterFuncs(1:length(t_us_L_now),n_now) = U_us_L_now;
            U_us_R_AeroCenterFuncs(1:length(t_us_R_now),n_now) = U_us_R_now;
            
            % store interp binned data separate rows
            t_wb_AeroCenterFuncs_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_AeroCenterFuncs_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_AeroCenterFuncs_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_AeroCenterFuncs_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_AeroCenterFuncs_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_AeroCenterFuncs_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_AeroCenterFuncs_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_AeroCenterFuncs_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_AeroCenterFuncs_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_AeroCenterFuncs_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_AeroCenterFuncs_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_AeroCenterFuncs_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_AeroCenterFuncs_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_AeroCenterFuncs_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_AeroCenterFuncs_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_AeroCenterFuncs_bins(:,n_now) = DU_wb_interp;
            
            t_ds_AeroCenterFuncs_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_AeroCenterFuncs_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_AeroCenterFuncs_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_AeroCenterFuncs_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_AeroCenterFuncs_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_AeroCenterFuncs_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_AeroCenterFuncs_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_AeroCenterFuncs_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_AeroCenterFuncs_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_AeroCenterFuncs_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_AeroCenterFuncs_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_AeroCenterFuncs_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_AeroCenterFuncs_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_AeroCenterFuncs_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_AeroCenterFuncs_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_AeroCenterFuncs_bins(:,n_now) = DU_ds_interp;
            
            t_us_AeroCenterFuncs_bins(:,n_now) = t_us_bin;
            stroke_us_L_AeroCenterFuncs_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_AeroCenterFuncs_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_AeroCenterFuncs_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_AeroCenterFuncs_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_AeroCenterFuncs_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_AeroCenterFuncs_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_AeroCenterFuncs_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_AeroCenterFuncs_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_AeroCenterFuncs_bins(:,n_now) = U_us_L_interp;
            U_us_R_AeroCenterFuncs_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_AeroCenterFuncs_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_AeroCenterFuncs_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_AeroCenterFuncs_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_AeroCenterFuncs_bins(:,n_now) = Daoa_us_interp;
            DU_us_AeroCenterFuncs_bins(:,n_now) = DU_us_interp;
            
            %% store MOD data
            
            %  bodykin modifications
            velMOD_AeroCenterFuncs(n_now,1) = velMOD_now;
            
            slipMOD_AeroCenterFuncs(n_now,1) = slipMOD_now;
            pitchMOD_AeroCenterFuncs(n_now,1) = pitchMOD_now;
            rollMOD_AeroCenterFuncs(n_now,1) = rollMOD_now;

            Fsp_pitchMOD_AeroCenterFuncs(n_now,1) = Fsp_pitchMOD_now;
            Fsp_rollMOD_AeroCenterFuncs(n_now,1) = Fsp_rollMOD_now;

            %  WB modifications
            freqMOD_wb_AeroCenterFuncs(n_now,1) = freqMOD_wb_now;
            
            strokeMOD_wb_L_AeroCenterFuncs_bins(:,n_now) = strokeMOD_wb_L_now;
            strokeMOD_wb_R_AeroCenterFuncs_bins(:,n_now) = strokeMOD_wb_R_now;
            strokeMOD_ds_L_AeroCenterFuncs_bins(:,n_now) = strokeMOD_ds_L_now;
            strokeMOD_ds_R_AeroCenterFuncs_bins(:,n_now) = strokeMOD_ds_R_now;
            strokeMOD_us_L_AeroCenterFuncs_bins(:,n_now) = strokeMOD_us_L_now;
            strokeMOD_us_R_AeroCenterFuncs_bins(:,n_now) = strokeMOD_us_R_now;

            pitchMOD_wb_L_AeroCenterFuncs_bins(:,n_now) = pitchMOD_wb_L_now;
            pitchMOD_wb_R_AeroCenterFuncs_bins(:,n_now) = pitchMOD_wb_R_now;
            pitchMOD_ds_L_AeroCenterFuncs_bins(:,n_now) = pitchMOD_ds_L_now;
            pitchMOD_ds_R_AeroCenterFuncs_bins(:,n_now) = pitchMOD_ds_R_now;
            pitchMOD_us_L_AeroCenterFuncs_bins(:,n_now) = pitchMOD_us_L_now;
            pitchMOD_us_R_AeroCenterFuncs_bins(:,n_now) = pitchMOD_us_R_now;

            devMOD_wb_L_AeroCenterFuncs_bins(:,n_now) = devMOD_wb_L_now;
            devMOD_wb_R_AeroCenterFuncs_bins(:,n_now) = devMOD_wb_R_now;
            devMOD_ds_L_AeroCenterFuncs_bins(:,n_now) = devMOD_ds_L_now;
            devMOD_ds_R_AeroCenterFuncs_bins(:,n_now) = devMOD_ds_R_now;
            devMOD_us_L_AeroCenterFuncs_bins(:,n_now) = devMOD_us_L_now;
            devMOD_us_R_AeroCenterFuncs_bins(:,n_now) = devMOD_us_R_now;

            DstrokeMOD_wb_AeroCenterFuncs_bins(:,n_now) = DstrokeMOD_wb_now;
            DstrokeMOD_ds_AeroCenterFuncs_bins(:,n_now) = DstrokeMOD_ds_now;
            DstrokeMOD_us_AeroCenterFuncs_bins(:,n_now) = DstrokeMOD_us_now;

            DpitchMOD_wb_AeroCenterFuncs_bins(:,n_now) = DpitchMOD_wb_now;
            DpitchMOD_ds_AeroCenterFuncs_bins(:,n_now) = DpitchMOD_ds_now;
            DpitchMOD_us_AeroCenterFuncs_bins(:,n_now) = DpitchMOD_us_now;

            DdevMOD_wb_AeroCenterFuncs_bins(:,n_now) = DdevMOD_wb_now;
            DdevMOD_ds_AeroCenterFuncs_bins(:,n_now) = DdevMOD_ds_now;
            DdevMOD_us_AeroCenterFuncs_bins(:,n_now) = DdevMOD_us_now;
        end
end

%% mean & 95%CI
% body kin
vel_AeroCenterFuncs_meanCIstd = [nanmean(vel_AeroCenterFuncs) 1.96*nanstd(vel_AeroCenterFuncs)/sqrt(length(vel_AeroCenterFuncs)) nanstd(vel_AeroCenterFuncs)];

slip_AeroCenterFuncs_meanCIstd = [nanmean(slip_AeroCenterFuncs) 1.96*nanstd(slip_AeroCenterFuncs)/sqrt(length(slip_AeroCenterFuncs)) nanstd(slip_AeroCenterFuncs)];
pitch_AeroCenterFuncs_meanCIstd = [nanmean(pitch_AeroCenterFuncs) 1.96*nanstd(pitch_AeroCenterFuncs)/sqrt(length(pitch_AeroCenterFuncs)) nanstd(pitch_AeroCenterFuncs)];
roll_AeroCenterFuncs_meanCIstd = [nanmean(roll_AeroCenterFuncs) 1.96*nanstd(roll_AeroCenterFuncs)/sqrt(length(roll_AeroCenterFuncs)) nanstd(roll_AeroCenterFuncs)];

Fsp_pitch_AeroCenterFuncs_meanCIstd = [nanmean(Fsp_pitch_AeroCenterFuncs) 1.96*nanstd(Fsp_pitch_AeroCenterFuncs)/sqrt(length(Fsp_pitch_AeroCenterFuncs)) nanstd(Fsp_pitch_AeroCenterFuncs)];
Fsp_roll_AeroCenterFuncs_meanCIstd = [nanmean(Fsp_roll_AeroCenterFuncs) 1.96*nanstd(Fsp_roll_AeroCenterFuncs)/sqrt(length(Fsp_roll_AeroCenterFuncs)) nanstd(Fsp_roll_AeroCenterFuncs)];

% WB kin
dt_ds_AeroCenterFuncs_meanCIstd = [nanmean(dt_ds_AeroCenterFuncs) 1.96*nanstd(dt_ds_AeroCenterFuncs)/sqrt(length(dt_ds_AeroCenterFuncs)) nanstd(dt_ds_AeroCenterFuncs)];
dt_us_AeroCenterFuncs_meanCIstd = [nanmean(dt_us_AeroCenterFuncs) 1.96*nanstd(dt_us_AeroCenterFuncs)/sqrt(length(dt_us_AeroCenterFuncs)) nanstd(dt_us_AeroCenterFuncs)];
f_wb_AeroCenterFuncs_meanCIstd = [nanmean(f_wb_AeroCenterFuncs) 1.96*nanstd(f_wb_AeroCenterFuncs)/sqrt(length(f_wb_AeroCenterFuncs)) nanstd(f_wb_AeroCenterFuncs)];
Rds_AeroCenterFuncs_meanCIstd = [nanmean(Rds_AeroCenterFuncs) 1.96*nanstd(Rds_AeroCenterFuncs)/sqrt(length(Rds_AeroCenterFuncs)) nanstd(Rds_AeroCenterFuncs)];

% calc_WBfunc_AeroCenterFuncs_circmeanCIstd
calc_WBfunc_AeroCenterFuncs_NOcircmeanCIstd

%% wbMOD means & 95%CI
velMOD_AeroCenterFuncs_meanCIstd = [nanmean(velMOD_AeroCenterFuncs) 1.96*nanstd(velMOD_AeroCenterFuncs)/sqrt(length(velMOD_AeroCenterFuncs)) nanstd(velMOD_AeroCenterFuncs)];

slipMOD_AeroCenterFuncs_meanCIstd = [nanmean(slipMOD_AeroCenterFuncs) 1.96*nanstd(slipMOD_AeroCenterFuncs)/sqrt(length(slipMOD_AeroCenterFuncs)) nanstd(slipMOD_AeroCenterFuncs)];
pitchMOD_AeroCenterFuncs_meanCIstd = [nanmean(pitchMOD_AeroCenterFuncs) 1.96*nanstd(pitchMOD_AeroCenterFuncs)/sqrt(length(pitchMOD_AeroCenterFuncs)) nanstd(pitchMOD_AeroCenterFuncs)];
rollMOD_AeroCenterFuncs_meanCIstd = [nanmean(rollMOD_AeroCenterFuncs) 1.96*nanstd(rollMOD_AeroCenterFuncs)/sqrt(length(rollMOD_AeroCenterFuncs)) nanstd(rollMOD_AeroCenterFuncs)];

Fsp_pitchMOD_AeroCenterFuncs_meanCIstd = [nanmean(Fsp_pitchMOD_AeroCenterFuncs) 1.96*nanstd(Fsp_pitchMOD_AeroCenterFuncs)/sqrt(length(Fsp_pitchMOD_AeroCenterFuncs)) nanstd(Fsp_pitchMOD_AeroCenterFuncs)];
Fsp_rollMOD_AeroCenterFuncs_meanCIstd = [nanmean(Fsp_rollMOD_AeroCenterFuncs) 1.96*nanstd(Fsp_rollMOD_AeroCenterFuncs)/sqrt(length(Fsp_rollMOD_AeroCenterFuncs)) nanstd(Fsp_rollMOD_AeroCenterFuncs)];

freqMOD_wb_AeroCenterFuncs_meanCIstd = [nanmean(freqMOD_wb_AeroCenterFuncs) 1.96*nanstd(freqMOD_wb_AeroCenterFuncs)/sqrt(length(freqMOD_wb_AeroCenterFuncs)) nanstd(freqMOD_wb_AeroCenterFuncs)];

% calc_WBmod_AeroCenterFuncs_circmeanCIstd
calc_WBmod_AeroCenterFuncs_NOcircmeanCIstd

%% WBfits
%% fit for LEFT wing (bins)
        %% legendre polynomials
            t_loc = t_wb_AeroCenterFuncs_bins(:,1);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);

            stroke_loc = stroke_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
            pitch_loc = pitch_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
            dev_loc = dev_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);

            strokeMOD_loc = strokeMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
            pitchMOD_loc = pitchMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
            devMOD_loc = devMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);

            [stroke_L_AeroCenterFuncs_fit_binmean, stroke_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
            [pitch_L_AeroCenterFuncs_fit_binmean, pitch_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
            [dev_L_AeroCenterFuncs_fit_binmean, dev_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

            [strokeMOD_L_AeroCenterFuncs_fit_binmean, strokeMOD_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
            [pitchMOD_L_AeroCenterFuncs_fit_binmean, pitchMOD_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
            [devMOD_L_AeroCenterFuncs_fit_binmean, devMOD_L_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

        %% fourier series ALL data
            t_loc = t_wb_AeroCenterFuncs_bins(:);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
            
            stroke_loc = stroke_wb_L_AeroCenterFuncs_bins(:);
            pitch_loc = pitch_wb_L_AeroCenterFuncs_bins(:);
            dev_loc = dev_wb_L_AeroCenterFuncs_bins(:);

            strokeMOD_loc = strokeMOD_wb_L_AeroCenterFuncs_bins(:);
            pitchMOD_loc = pitchMOD_wb_L_AeroCenterFuncs_bins(:);
            devMOD_loc = devMOD_wb_L_AeroCenterFuncs_bins(:);

            % wingbeats
            [stroke_L_AeroCenterFuncs_fourier_fit_binmean, stroke_L_AeroCenterFuncs_fourier_gof_binmean,stroke_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
            [pitch_L_AeroCenterFuncs_fourier_fit_binmean, pitch_L_AeroCenterFuncs_fourier_gof_binmean,pitch_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
            [dev_L_AeroCenterFuncs_fourier_fit_binmean, dev_L_AeroCenterFuncs_fourier_gof_binmean,dev_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);

            % wingbeatMODs
            [strokeMOD_L_AeroCenterFuncs_fourier_fit_binmean, strokeMOD_L_AeroCenterFuncs_fourier_gof_binmean,strokeMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
            [pitchMOD_L_AeroCenterFuncs_fourier_fit_binmean, pitchMOD_L_AeroCenterFuncs_fourier_gof_binmean,pitchMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
            [devMOD_L_AeroCenterFuncs_fourier_fit_binmean, devMOD_L_AeroCenterFuncs_fourier_gof_binmean,devMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
        %% fourier series circmeans
%             t_loc = t_wb_AeroCenterFuncs_bins(:,1);
%             Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
% 
%             stroke_loc = stroke_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
%             pitch_loc = pitch_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
%             dev_loc = dev_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             strokeMOD_loc = strokeMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
%             pitchMOD_loc = pitchMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
%             devMOD_loc = devMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             % wingbeats
%             [stroke_L_AeroCenterFuncs_fourier_fit_binmean, stroke_L_AeroCenterFuncs_fourier_gof_binmean,stroke_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
%             [pitch_L_AeroCenterFuncs_fourier_fit_binmean, pitch_L_AeroCenterFuncs_fourier_gof_binmean,pitch_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
%             [dev_L_AeroCenterFuncs_fourier_fit_binmean, dev_L_AeroCenterFuncs_fourier_gof_binmean,dev_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);
% 
%             % wingbeatMODs
%             [strokeMOD_L_AeroCenterFuncs_fourier_fit_binmean, strokeMOD_L_AeroCenterFuncs_fourier_gof_binmean,strokeMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
%             [pitchMOD_L_AeroCenterFuncs_fourier_fit_binmean, pitchMOD_L_AeroCenterFuncs_fourier_gof_binmean,pitchMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
%             [devMOD_L_AeroCenterFuncs_fourier_fit_binmean, devMOD_L_AeroCenterFuncs_fourier_gof_binmean,devMOD_L_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
%% fit for Right wing (bins)
        %% legendre polynomials
            t_loc = t_wb_AeroCenterFuncs_bins(:,1);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);

            stroke_loc = stroke_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
            pitch_loc = pitch_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
            dev_loc = dev_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);

            strokeMOD_loc = strokeMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
            pitchMOD_loc = pitchMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
            devMOD_loc = devMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);

            [stroke_R_AeroCenterFuncs_fit_binmean, stroke_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,stroke_loc,Rds_loc);
            [pitch_R_AeroCenterFuncs_fit_binmean, pitch_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,pitch_loc,Rds_loc);
            [dev_R_AeroCenterFuncs_fit_binmean, dev_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,dev_loc,Rds_loc);

            [strokeMOD_R_AeroCenterFuncs_fit_binmean, strokeMOD_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,strokeMOD_loc,Rds_loc);
            [pitchMOD_R_AeroCenterFuncs_fit_binmean, pitchMOD_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,pitchMOD_loc,Rds_loc);
            [devMOD_R_AeroCenterFuncs_fit_binmean, devMOD_R_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,devMOD_loc,Rds_loc);

        %% fourier series ALL data
            t_loc = t_wb_AeroCenterFuncs_bins(:);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
            
            stroke_loc = stroke_wb_R_AeroCenterFuncs_bins(:);
            pitch_loc = pitch_wb_R_AeroCenterFuncs_bins(:);
            dev_loc = dev_wb_R_AeroCenterFuncs_bins(:);

            strokeMOD_loc = strokeMOD_wb_R_AeroCenterFuncs_bins(:);
            pitchMOD_loc = pitchMOD_wb_R_AeroCenterFuncs_bins(:);
            devMOD_loc = devMOD_wb_R_AeroCenterFuncs_bins(:);

            % wingbeats
            [stroke_R_AeroCenterFuncs_fourier_fit_binmean, stroke_R_AeroCenterFuncs_fourier_gof_binmean,stroke_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
            [pitch_R_AeroCenterFuncs_fourier_fit_binmean, pitch_R_AeroCenterFuncs_fourier_gof_binmean,pitch_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
            [dev_R_AeroCenterFuncs_fourier_fit_binmean, dev_R_AeroCenterFuncs_fourier_gof_binmean,dev_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);

            % wingbeatMODs
            [strokeMOD_R_AeroCenterFuncs_fourier_fit_binmean, strokeMOD_R_AeroCenterFuncs_fourier_gof_binmean,strokeMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
            [pitchMOD_R_AeroCenterFuncs_fourier_fit_binmean, pitchMOD_R_AeroCenterFuncs_fourier_gof_binmean,pitchMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
            [devMOD_R_AeroCenterFuncs_fourier_fit_binmean, devMOD_R_AeroCenterFuncs_fourier_gof_binmean,devMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
        %% fourier series circmeans
%             t_loc = t_wb_AeroCenterFuncs_bins(:,1);
%             Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
% 
%             stroke_loc = stroke_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
%             pitch_loc = pitch_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
%             dev_loc = dev_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             strokeMOD_loc = strokeMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
%             pitchMOD_loc = pitchMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
%             devMOD_loc = devMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             % wingbeats
%             [stroke_R_AeroCenterFuncs_fourier_fit_binmean, stroke_R_AeroCenterFuncs_fourier_gof_binmean,stroke_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, stroke_loc,stroke_fourier_order,plot_fourier);
%             [pitch_R_AeroCenterFuncs_fourier_fit_binmean, pitch_R_AeroCenterFuncs_fourier_gof_binmean,pitch_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitch_loc,pitch_fourier_order,plot_fourier);
%             [dev_R_AeroCenterFuncs_fourier_fit_binmean, dev_R_AeroCenterFuncs_fourier_gof_binmean,dev_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, dev_loc,dev_fourier_order,plot_fourier);
% 
%             % wingbeatMODs
%             [strokeMOD_R_AeroCenterFuncs_fourier_fit_binmean, strokeMOD_R_AeroCenterFuncs_fourier_gof_binmean,strokeMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, strokeMOD_loc,MOD_fourier_order,plot_fourier);
%             [pitchMOD_R_AeroCenterFuncs_fourier_fit_binmean, pitchMOD_R_AeroCenterFuncs_fourier_gof_binmean,pitchMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, pitchMOD_loc,MOD_fourier_order,plot_fourier);
%             [devMOD_R_AeroCenterFuncs_fourier_fit_binmean, devMOD_R_AeroCenterFuncs_fourier_gof_binmean,devMOD_R_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, devMOD_loc,MOD_fourier_order,plot_fourier);
    
    
%% L-R

    stroke_loc = Dstroke_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
    pitch_loc = Dpitch_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
    dev_loc = Ddev_wb_AeroCenterFuncs_bins_meanCIstd(:,1);

    strokeMOD_loc = DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
    pitchMOD_loc = DpitchMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
    devMOD_loc = DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);

        %% legendre polynomials
            t_loc = t_wb_AeroCenterFuncs_bins(:,1);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);

            Dstroke_loc = Dstroke_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
            Dpitch_loc = Dpitch_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
            Ddev_loc = Ddev_wb_AeroCenterFuncs_bins_meanCIstd(:,1);

            DstrokeMOD_loc = DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
            DpitchMOD_loc = DpitchMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
            DdevMOD_loc = DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);

            [Dstroke_AeroCenterFuncs_fit_binmean, Dstroke_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_stroke,t_loc,Dstroke_loc,Rds_loc);
            [Dpitch_AeroCenterFuncs_fit_binmean, Dpitch_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_pitch,t_loc,Dpitch_loc,Rds_loc);
            [Ddev_AeroCenterFuncs_fit_binmean, Ddev_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_dev,t_loc,Ddev_loc,Rds_loc);

            [DstrokeMOD_AeroCenterFuncs_fit_binmean, DstrokeMOD_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,DstrokeMOD_loc,Rds_loc);
            [DpitchMOD_AeroCenterFuncs_fit_binmean, DpitchMOD_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,DpitchMOD_loc,Rds_loc);
            [DdevMOD_AeroCenterFuncs_fit_binmean, DdevMOD_AeroCenterFuncs_fit_binmean_periodic] = WBfitting_singlevar_updowncut(n_pol_MOD,t_loc,DdevMOD_loc,Rds_loc);

        %% fourier series ALL data
            t_loc = t_wb_AeroCenterFuncs_bins(:);
            Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
            
            Dstroke_loc = Dstroke_wb_AeroCenterFuncs_bins(:);
            Dpitch_loc = Dpitch_wb_AeroCenterFuncs_bins(:);
            Ddev_loc = Ddev_wb_AeroCenterFuncs_bins(:);

            DstrokeMOD_loc = DstrokeMOD_wb_AeroCenterFuncs_bins(:);
            DpitchMOD_loc = DpitchMOD_wb_AeroCenterFuncs_bins(:);
            DdevMOD_loc = DdevMOD_wb_AeroCenterFuncs_bins(:);

            % wingbeats
            [Dstroke_AeroCenterFuncs_fourier_fit_binmean, Dstroke_AeroCenterFuncs_fourier_gof_binmean,Dstroke_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Dstroke_loc,stroke_fourier_order,plot_fourier);
            [Dpitch_AeroCenterFuncs_fourier_fit_binmean, Dpitch_AeroCenterFuncs_fourier_gof_binmean,Dpitch_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Dpitch_loc,pitch_fourier_order,plot_fourier);
            [Ddev_AeroCenterFuncs_fourier_fit_binmean, Ddev_AeroCenterFuncs_fourier_gof_binmean,Ddev_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Ddev_loc,dev_fourier_order,plot_fourier);

            % wingbeatMODs
            [DstrokeMOD_AeroCenterFuncs_fourier_fit_binmean, DstrokeMOD_AeroCenterFuncs_fourier_gof_binmean,DstrokeMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DstrokeMOD_loc,MOD_fourier_order,plot_fourier);
            [DpitchMOD_AeroCenterFuncs_fourier_fit_binmean, DpitchMOD_AeroCenterFuncs_fourier_gof_binmean,DpitchMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DpitchMOD_loc,MOD_fourier_order,plot_fourier);
            [DdevMOD_AeroCenterFuncs_fourier_fit_binmean, DdevMOD_AeroCenterFuncs_fourier_gof_binmean,DdevMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DdevMOD_loc,MOD_fourier_order,plot_fourier);
    
        %% fourier series circmeans
%             t_loc = t_wb_AeroCenterFuncs_bins(:,1);
%             Rds_loc = Rds_AeroCenterFuncs_meanCIstd(1);
% 
%             Dstroke_loc = Dstroke_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
%             Dpitch_loc = Dpitch_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
%             Ddev_loc = Ddev_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             DstrokeMOD_loc = DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
%             DpitchMOD_loc = DpitchMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
%             DdevMOD_loc = DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd(:,1);
% 
%             % wingbeats
%             [Dstroke_AeroCenterFuncs_fourier_fit_binmean, Dstroke_AeroCenterFuncs_fourier_gof_binmean,Dstroke_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Dstroke_loc,stroke_fourier_order,plot_fourier);
%             [Dpitch_AeroCenterFuncs_fourier_fit_binmean, Dpitch_AeroCenterFuncs_fourier_gof_binmean,Dpitch_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Dpitch_loc,pitch_fourier_order,plot_fourier);
%             [Ddev_AeroCenterFuncs_fourier_fit_binmean, Ddev_AeroCenterFuncs_fourier_gof_binmean,Ddev_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, Ddev_loc,dev_fourier_order,plot_fourier);
% 
%             % wingbeatMODs
%             [DstrokeMOD_AeroCenterFuncs_fourier_fit_binmean, DstrokeMOD_AeroCenterFuncs_fourier_gof_binmean,DstrokeMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DstrokeMOD_loc,MOD_fourier_order,plot_fourier);
%             [DpitchMOD_AeroCenterFuncs_fourier_fit_binmean, DpitchMOD_AeroCenterFuncs_fourier_gof_binmean,DpitchMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DpitchMOD_loc,MOD_fourier_order,plot_fourier);
%             [DdevMOD_AeroCenterFuncs_fourier_fit_binmean, DdevMOD_AeroCenterFuncs_fourier_gof_binmean,DdevMOD_AeroCenterFuncs_fourier_coeffs_binmean] = create_fourier_series_4thN8th_order(t_loc, DdevMOD_loc,MOD_fourier_order,plot_fourier);

%% save CLIPd wb data
n_AeroCenterFuncs = n_now   

if save_on == 1

    save(['WBmod_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.mat'],...
'n_AeroCenterFuncs',...
'seq_nr_AeroCenterFuncs',...
'wb_nr_AeroCenterFuncs',...
'FirstMomentRatio_cut',...
'FirstMomentRatio_plot',...
...
'n_pol_stroke',...
'n_pol_pitch',...
'n_pol_dev',...
'n_pol_MOD',...
...
'clip_side_AeroCenterFuncs',...
'clip_type_AeroCenterFuncs',...
...
'AeroCenterFuncClipped',...
'AeroCenterFuncIntact',...
'AeroCenterFuncRatio',...
...
'FirstMomentNormClipped',...
'FirstMomentNormIntact',...
'FirstMomentRatio',...
...
'AeroCenterFuncClipped_AeroCenterFuncs',...
'AeroCenterFuncIntact_AeroCenterFuncs',...
'AeroCenterFuncRatio_AeroCenterFuncs',...
...
'FirstMomentNormClipped_AeroCenterFuncs',...
'FirstMomentNormIntact_AeroCenterFuncs',...
'FirstMomentRatio_AeroCenterFuncs',...
...
'SecondMomentNormClipped_AeroCenterFuncs',...
'SecondMomentNormIntact_AeroCenterFuncs',...
'SecondMomentRatio_AeroCenterFuncs',...
...
'ThirdMomentNormClipped_AeroCenterFuncs',...
'ThirdMomentNormIntact_AeroCenterFuncs',...
'ThirdMomentRatio_AeroCenterFuncs',...
...
'AreaNormClipped_AeroCenterFuncs',...
'AreaNormIntact_AeroCenterFuncs',...
'AreaRatio_AeroCenterFuncs',...
...
'LengthClipped_pixels_AeroCenterFuncs',...
'LengthIntact_pixels_AeroCenterFuncs',...
'LengthRatio_AeroCenterFuncs',...
...
'stroke_R_AeroCenterFuncs_fit_binmean_periodic',...
'pitch_R_AeroCenterFuncs_fit_binmean_periodic',...
'dev_R_AeroCenterFuncs_fit_binmean_periodic',...
...
'stroke_L_AeroCenterFuncs_fit_binmean_periodic',...
'pitch_L_AeroCenterFuncs_fit_binmean_periodic',...
'dev_L_AeroCenterFuncs_fit_binmean_periodic',...
...
'strokeMOD_R_AeroCenterFuncs_fit_binmean_periodic',...
'pitchMOD_R_AeroCenterFuncs_fit_binmean_periodic',...
'devMOD_R_AeroCenterFuncs_fit_binmean_periodic',...
...
'strokeMOD_L_AeroCenterFuncs_fit_binmean_periodic',...
'pitchMOD_L_AeroCenterFuncs_fit_binmean_periodic',...
'devMOD_L_AeroCenterFuncs_fit_binmean_periodic',...
...
'stroke_fourier_order',...
'pitch_fourier_order',...
'dev_fourier_order',...
'MOD_fourier_order',...
...
'stroke_R_AeroCenterFuncs_fourier_coeffs_binmean',...
'pitch_R_AeroCenterFuncs_fourier_coeffs_binmean',...
'dev_R_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'stroke_R_AeroCenterFuncs_fourier_gof_binmean',...
'pitch_R_AeroCenterFuncs_fourier_gof_binmean',...
'dev_R_AeroCenterFuncs_fourier_gof_binmean',...
...
'stroke_L_AeroCenterFuncs_fourier_coeffs_binmean',...
'pitch_L_AeroCenterFuncs_fourier_coeffs_binmean',...
'dev_L_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'stroke_L_AeroCenterFuncs_fourier_gof_binmean',...
'pitch_L_AeroCenterFuncs_fourier_gof_binmean',...
'dev_L_AeroCenterFuncs_fourier_gof_binmean',...
...
'Dstroke_AeroCenterFuncs_fourier_coeffs_binmean',...
'Dpitch_AeroCenterFuncs_fourier_coeffs_binmean',...
'Ddev_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'Dstroke_AeroCenterFuncs_fourier_gof_binmean',...
'Dpitch_AeroCenterFuncs_fourier_gof_binmean',...
'Ddev_AeroCenterFuncs_fourier_gof_binmean',...
...
'strokeMOD_R_AeroCenterFuncs_fourier_coeffs_binmean',...
'pitchMOD_R_AeroCenterFuncs_fourier_coeffs_binmean',...
'devMOD_R_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'strokeMOD_R_AeroCenterFuncs_fourier_gof_binmean',...
'pitchMOD_R_AeroCenterFuncs_fourier_gof_binmean',...
'devMOD_R_AeroCenterFuncs_fourier_gof_binmean',...
...
'strokeMOD_L_AeroCenterFuncs_fourier_coeffs_binmean',...
'pitchMOD_L_AeroCenterFuncs_fourier_coeffs_binmean',...
'devMOD_L_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'strokeMOD_L_AeroCenterFuncs_fourier_gof_binmean',...
'pitchMOD_L_AeroCenterFuncs_fourier_gof_binmean',...
'devMOD_L_AeroCenterFuncs_fourier_gof_binmean',...
...
'DstrokeMOD_AeroCenterFuncs_fourier_coeffs_binmean',...
'DpitchMOD_AeroCenterFuncs_fourier_coeffs_binmean',...
'DdevMOD_AeroCenterFuncs_fourier_coeffs_binmean',...
...
'DstrokeMOD_AeroCenterFuncs_fourier_gof_binmean',...
'DpitchMOD_AeroCenterFuncs_fourier_gof_binmean',...
'DdevMOD_AeroCenterFuncs_fourier_gof_binmean',...
...
'stroke_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'stroke_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'stroke_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'pitch_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'pitch_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'pitch_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'dev_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'dev_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'dev_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'strokeMOD_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'strokeMOD_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'strokeMOD_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'pitchMOD_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'pitchMOD_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'pitchMOD_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'devMOD_wb_R_AeroCenterFuncs_bins_meanCIstd',...
'devMOD_ds_R_AeroCenterFuncs_bins_meanCIstd',...
'devMOD_us_R_AeroCenterFuncs_bins_meanCIstd',...
...
'stroke_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'stroke_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'stroke_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'pitch_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'pitch_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'pitch_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'dev_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'dev_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'dev_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'strokeMOD_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'strokeMOD_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'strokeMOD_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'pitchMOD_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'pitchMOD_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'pitchMOD_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'devMOD_wb_L_AeroCenterFuncs_bins_meanCIstd',...
'devMOD_ds_L_AeroCenterFuncs_bins_meanCIstd',...
'devMOD_us_L_AeroCenterFuncs_bins_meanCIstd',...
...
'DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd',...
'DstrokeMOD_ds_AeroCenterFuncs_bins_meanCIstd',...
'DstrokeMOD_us_AeroCenterFuncs_bins_meanCIstd',...
...
'DpitchMOD_wb_AeroCenterFuncs_bins_meanCIstd',...
'DpitchMOD_ds_AeroCenterFuncs_bins_meanCIstd',...
'DpitchMOD_us_AeroCenterFuncs_bins_meanCIstd',...
...
'DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd',...
'DdevMOD_ds_AeroCenterFuncs_bins_meanCIstd',...
'DdevMOD_us_AeroCenterFuncs_bins_meanCIstd',...
...
'vel_AeroCenterFuncs',...
'slip_AeroCenterFuncs',...
'pitch_AeroCenterFuncs',...
'roll_AeroCenterFuncs',...
'Fsp_pitch_AeroCenterFuncs',...
'Fsp_roll_AeroCenterFuncs',...
...
'dt_ds_AeroCenterFuncs',...
'dt_us_AeroCenterFuncs',...
'f_wb_AeroCenterFuncs',...
'Rds_AeroCenterFuncs',...
...
'vel_AeroCenterFuncs_meanCIstd',...
'slip_AeroCenterFuncs_meanCIstd',...
'pitch_AeroCenterFuncs_meanCIstd',...
'roll_AeroCenterFuncs_meanCIstd',...
'Fsp_pitch_AeroCenterFuncs_meanCIstd',...
'Fsp_roll_AeroCenterFuncs_meanCIstd',...
...
'dt_ds_AeroCenterFuncs_meanCIstd',...
'dt_us_AeroCenterFuncs_meanCIstd',...
'f_wb_AeroCenterFuncs_meanCIstd',...
'Rds_AeroCenterFuncs_meanCIstd',...
...
'velMOD_AeroCenterFuncs_meanCIstd',...
'velMOD_AeroCenterFuncs',...
...
'slipMOD_AeroCenterFuncs_meanCIstd',...
'slipMOD_AeroCenterFuncs',...
...
'pitchMOD_AeroCenterFuncs_meanCIstd',...
'pitchMOD_AeroCenterFuncs',...
...
'rollMOD_AeroCenterFuncs_meanCIstd',...
'rollMOD_AeroCenterFuncs',...
...
'Fsp_pitchMOD_AeroCenterFuncs_meanCIstd',...
'Fsp_pitchMOD_AeroCenterFuncs',...
...
'Fsp_rollMOD_AeroCenterFuncs_meanCIstd',...
'Fsp_rollMOD_AeroCenterFuncs',...
...
'freqMOD_wb_AeroCenterFuncs_meanCIstd',...
'freqMOD_wb_AeroCenterFuncs',...
...
'strokeMOD_wb_L_AeroCenterFuncs_bins',...
'strokeMOD_wb_R_AeroCenterFuncs_bins',...
'strokeMOD_ds_L_AeroCenterFuncs_bins',...
'strokeMOD_ds_R_AeroCenterFuncs_bins',...
'strokeMOD_us_L_AeroCenterFuncs_bins',...
'strokeMOD_us_R_AeroCenterFuncs_bins',...
...
'pitchMOD_wb_L_AeroCenterFuncs_bins',...
'pitchMOD_wb_R_AeroCenterFuncs_bins',...
'pitchMOD_ds_L_AeroCenterFuncs_bins',...
'pitchMOD_ds_R_AeroCenterFuncs_bins',...
'pitchMOD_us_L_AeroCenterFuncs_bins',...
'pitchMOD_us_R_AeroCenterFuncs_bins',...
...
'devMOD_wb_L_AeroCenterFuncs_bins',...
'devMOD_wb_R_AeroCenterFuncs_bins',...
'devMOD_ds_L_AeroCenterFuncs_bins',...
'devMOD_ds_R_AeroCenterFuncs_bins',...
'devMOD_us_L_AeroCenterFuncs_bins',...
'devMOD_us_R_AeroCenterFuncs_bins',...
...
'DstrokeMOD_wb_AeroCenterFuncs_bins',...
'DstrokeMOD_ds_AeroCenterFuncs_bins',...
'DstrokeMOD_us_AeroCenterFuncs_bins',...
...
'DpitchMOD_wb_AeroCenterFuncs_bins',...
'DpitchMOD_ds_AeroCenterFuncs_bins',...
'DpitchMOD_us_AeroCenterFuncs_bins',...
...
'DdevMOD_wb_AeroCenterFuncs_bins',...
'DdevMOD_ds_AeroCenterFuncs_bins',...
'DdevMOD_us_AeroCenterFuncs_bins',...
...
't_wb_L_AeroCenterFuncs',...
't_wb_R_AeroCenterFuncs',...
'stroke_wb_L_AeroCenterFuncs',...
'stroke_wb_R_AeroCenterFuncs',...
'pitch_wb_L_AeroCenterFuncs',...
'pitch_wb_R_AeroCenterFuncs',...
'dev_wb_L_AeroCenterFuncs',...
'dev_wb_R_AeroCenterFuncs',...
'aoa_wb_L_AeroCenterFuncs',...
'aoa_wb_R_AeroCenterFuncs',...
'U_wb_L_AeroCenterFuncs',...
'U_wb_R_AeroCenterFuncs',...
...
't_ds_L_AeroCenterFuncs',...
't_ds_R_AeroCenterFuncs',...
'stroke_ds_L_AeroCenterFuncs',...
'stroke_ds_R_AeroCenterFuncs',...
'pitch_ds_L_AeroCenterFuncs',...
'pitch_ds_R_AeroCenterFuncs',...
'dev_ds_L_AeroCenterFuncs',...
'dev_ds_R_AeroCenterFuncs',...
'aoa_ds_L_AeroCenterFuncs',...
'aoa_ds_R_AeroCenterFuncs',...
'U_ds_L_AeroCenterFuncs',...
'U_ds_R_AeroCenterFuncs',...
...
't_us_L_AeroCenterFuncs',...
't_us_R_AeroCenterFuncs',...
'stroke_us_L_AeroCenterFuncs',...
'stroke_us_R_AeroCenterFuncs',...
'pitch_us_L_AeroCenterFuncs',...
'pitch_us_R_AeroCenterFuncs',...
'dev_us_L_AeroCenterFuncs',...
'dev_us_R_AeroCenterFuncs',...
'aoa_us_L_AeroCenterFuncs',...
'aoa_us_R_AeroCenterFuncs',...
'U_us_L_AeroCenterFuncs',...
'U_us_R_AeroCenterFuncs',...
...
't_wb_AeroCenterFuncs_bins',...
'stroke_wb_L_AeroCenterFuncs_bins',...
'stroke_wb_R_AeroCenterFuncs_bins',...
'pitch_wb_L_AeroCenterFuncs_bins',...
'pitch_wb_R_AeroCenterFuncs_bins',...
'dev_wb_L_AeroCenterFuncs_bins',...
'dev_wb_R_AeroCenterFuncs_bins',...
'aoa_wb_L_AeroCenterFuncs_bins',...
'aoa_wb_R_AeroCenterFuncs_bins',...
'U_wb_L_AeroCenterFuncs_bins',...
'U_wb_R_AeroCenterFuncs_bins',...
'Dstroke_wb_AeroCenterFuncs_bins',...
'Dpitch_wb_AeroCenterFuncs_bins',...
'Ddev_wb_AeroCenterFuncs_bins',...
'Daoa_wb_AeroCenterFuncs_bins',...
'DU_wb_AeroCenterFuncs_bins',...
...
't_ds_AeroCenterFuncs_bins',...
'stroke_ds_L_AeroCenterFuncs_bins',...
'stroke_ds_R_AeroCenterFuncs_bins',...
'pitch_ds_L_AeroCenterFuncs_bins',...
'pitch_ds_R_AeroCenterFuncs_bins',...
'dev_ds_L_AeroCenterFuncs_bins',...
'dev_ds_R_AeroCenterFuncs_bins',...
'aoa_ds_L_AeroCenterFuncs_bins',...
'aoa_ds_R_AeroCenterFuncs_bins',...
'U_ds_L_AeroCenterFuncs_bins',...
'U_ds_R_AeroCenterFuncs_bins',...
'Dstroke_ds_AeroCenterFuncs_bins',...
'Dpitch_ds_AeroCenterFuncs_bins',...
'Ddev_ds_AeroCenterFuncs_bins',...
'Daoa_ds_AeroCenterFuncs_bins',...
'DU_ds_AeroCenterFuncs_bins',...
...
't_us_AeroCenterFuncs_bins',...
'stroke_us_L_AeroCenterFuncs_bins',...
'stroke_us_R_AeroCenterFuncs_bins',...
'pitch_us_L_AeroCenterFuncs_bins',...
'pitch_us_R_AeroCenterFuncs_bins',...
'dev_us_L_AeroCenterFuncs_bins',...
'dev_us_R_AeroCenterFuncs_bins',...
'aoa_us_L_AeroCenterFuncs_bins',...
'aoa_us_R_AeroCenterFuncs_bins',...
'U_us_L_AeroCenterFuncs_bins',...
'U_us_R_AeroCenterFuncs_bins',...
'Dstroke_us_AeroCenterFuncs_bins',...
'Dpitch_us_AeroCenterFuncs_bins',...
'Ddev_us_AeroCenterFuncs_bins',...
'Daoa_us_AeroCenterFuncs_bins',...
'DU_us_AeroCenterFuncs_bins');

end

