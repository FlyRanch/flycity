%% cal WBmod kin data
% RollAccel
Dstroke_RollAccel_max = DstrokeMOD_RollAccel_max * roll_dot_dot_norm_mirror_max;
Dpitch_RollAccel_max = DpitchMOD_RollAccel_max * roll_dot_dot_norm_mirror_max;
Ddev_RollAccel_max = DdevMOD_RollAccel_max * roll_dot_dot_norm_mirror_max;

Dstroke_RollAccel_min = DstrokeMOD_RollAccel_max * roll_dot_dot_norm_mirror_min;
Dpitch_RollAccel_min = DpitchMOD_RollAccel_max * roll_dot_dot_norm_mirror_min;
Ddev_RollAccel_min = DdevMOD_RollAccel_max * roll_dot_dot_norm_mirror_min;

% YawAccel
Dstroke_YawAccel_max = DstrokeMOD_YawAccel_max * yaw_dot_dot_norm_mirror_max;
Dpitch_YawAccel_max = DpitchMOD_YawAccel_max * yaw_dot_dot_norm_mirror_max;
Ddev_YawAccel_max = DdevMOD_YawAccel_max * yaw_dot_dot_norm_mirror_max;

Dstroke_YawAccel_min = DstrokeMOD_YawAccel_max * yaw_dot_dot_norm_mirror_min;
Dpitch_YawAccel_min = DpitchMOD_YawAccel_max * yaw_dot_dot_norm_mirror_min;
Ddev_YawAccel_min = DdevMOD_YawAccel_max * yaw_dot_dot_norm_mirror_min;

% PitchAccel
stroke_PitchAccel_max = strokeMOD_PitchAccel_max * pitch_dot_dot_norm_max;
pitch_PitchAccel_max = pitchMOD_PitchAccel_max * pitch_dot_dot_norm_max;
dev_PitchAccel_max = devMOD_PitchAccel_max * pitch_dot_dot_norm_max;

stroke_PitchAccel_min = strokeMOD_PitchAccel_max * pitch_dot_dot_norm_min;
pitch_PitchAccel_min = pitchMOD_PitchAccel_max * pitch_dot_dot_norm_min;
dev_PitchAccel_min = devMOD_PitchAccel_max * pitch_dot_dot_norm_min;

% Fenhance
stroke_Fenhance_max = strokeMOD_Fenhance_max * F_norm_max;
pitch_Fenhance_max = pitchMOD_Fenhance_max * F_norm_max;
dev_Fenhance_max = devMOD_Fenhance_max * F_norm_max;
freq_Fenhance_max = freqMOD_Fenhance * F_norm_max;

%% cal WBmod kin data WBmean
% RollAccel_WBmean
Dstroke_RollAccel_WBmean_max = DstrokeMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_max;
Dpitch_RollAccel_WBmean_max = DpitchMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_max;
Ddev_RollAccel_WBmean_max = DdevMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_max;

Dstroke_RollAccel_WBmean_min = DstrokeMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_min;
Dpitch_RollAccel_WBmean_min = DpitchMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_min;
Ddev_RollAccel_WBmean_min = DdevMOD_RollAccel_max * roll_dot_dot_WBmean_norm_mirror_min;

% YawAccel_WBmean
Dstroke_YawAccel_WBmean_max = DstrokeMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_max;
Dpitch_YawAccel_WBmean_max = DpitchMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_max;
Ddev_YawAccel_WBmean_max = DdevMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_max;

Dstroke_YawAccel_WBmean_min = DstrokeMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_min;
Dpitch_YawAccel_WBmean_min = DpitchMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_min;
Ddev_YawAccel_WBmean_min = DdevMOD_YawAccel_max * yaw_dot_dot_WBmean_norm_mirror_min;

% PitchAccel_WBmean
stroke_PitchAccel_WBmean_max = strokeMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_max;
pitch_PitchAccel_WBmean_max = pitchMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_max;
dev_PitchAccel_WBmean_max = devMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_max;

stroke_PitchAccel_WBmean_min = strokeMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_min;
pitch_PitchAccel_WBmean_min = pitchMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_min;
dev_PitchAccel_WBmean_min = devMOD_PitchAccel_max * pitch_dot_dot_WBmean_norm_min;

% Fenhance
stroke_Fenhance_WBmean_max = strokeMOD_Fenhance_max * F_WBmean_norm_max;
pitch_Fenhance_WBmean_max = pitchMOD_Fenhance_max * F_WBmean_norm_max;
dev_Fenhance_WBmean_max = devMOD_Fenhance_max * F_WBmean_norm_max;
freq_Fenhance_WBmean_max = freqMOD_Fenhance * F_WBmean_norm_max;
