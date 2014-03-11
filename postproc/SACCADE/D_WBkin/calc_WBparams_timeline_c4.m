
%% stroke
% L
stroke_wb_L_seq_bins_MAX1_c4 = max(stroke_wb_L_seq_bins_mean_c4(1:100,:))';
stroke_wb_L_seq_bins_MAX2_c4 = max(stroke_wb_L_seq_bins_mean_c4(100:end,:))';
stroke_wb_L_seq_bins_MAX1_c4 = [stroke_wb_L_seq_bins_MAX1_c4;nan];
stroke_wb_L_seq_bins_MAX2_c4 = [nan;stroke_wb_L_seq_bins_MAX2_c4];
stroke_wb_L_seq_bins_MAX_c4 = [stroke_wb_L_seq_bins_MAX1_c4 stroke_wb_L_seq_bins_MAX2_c4];
stroke_wb_L_seq_bins_MAX_c4 = max(stroke_wb_L_seq_bins_MAX_c4(1:end-1,:)')';

stroke_wb_L_seq_bins_MIN_c4 = min(stroke_wb_L_seq_bins_mean_c4)';
Astroke_wb_L_seq_bins_c4 = stroke_wb_L_seq_bins_MAX_c4 - stroke_wb_L_seq_bins_MIN_c4;

% R
stroke_wb_R_seq_bins_MAX1_c4 = max(stroke_wb_R_seq_bins_mean_c4(1:100,:))';
stroke_wb_R_seq_bins_MAX2_c4 = max(stroke_wb_R_seq_bins_mean_c4(100:end,:))';
stroke_wb_R_seq_bins_MAX1_c4 = [stroke_wb_R_seq_bins_MAX1_c4;nan];
stroke_wb_R_seq_bins_MAX2_c4 = [nan;stroke_wb_R_seq_bins_MAX2_c4];
stroke_wb_R_seq_bins_MAX_c4 = [stroke_wb_R_seq_bins_MAX1_c4 stroke_wb_R_seq_bins_MAX2_c4];
stroke_wb_R_seq_bins_MAX_c4 = max(stroke_wb_R_seq_bins_MAX_c4(1:end-1,:)')';

stroke_wb_R_seq_bins_MIN_c4 = min(stroke_wb_R_seq_bins_mean_c4)';
Astroke_wb_R_seq_bins_c4 = stroke_wb_R_seq_bins_MAX_c4 - stroke_wb_R_seq_bins_MIN_c4;

% R+L
Astroke_wb_seq_bins_c4 = nanmean([Astroke_wb_R_seq_bins_c4 Astroke_wb_L_seq_bins_c4]')';
stroke_wb_seq_bins_MAX_c4 = nanmean([stroke_wb_R_seq_bins_MAX_c4 stroke_wb_L_seq_bins_MAX_c4]')';
stroke_wb_seq_bins_MIN_c4 = nanmean([stroke_wb_R_seq_bins_MIN_c4 stroke_wb_L_seq_bins_MIN_c4]')';

% R-L
dAstroke_wb_seq_bins_c4 = Astroke_wb_R_seq_bins_c4 - Astroke_wb_L_seq_bins_c4;
Dstroke_wb_seq_bins_MAX_c4 = stroke_wb_R_seq_bins_MAX_c4 - stroke_wb_L_seq_bins_MAX_c4;
Dstroke_wb_seq_bins_MIN_c4 = stroke_wb_R_seq_bins_MIN_c4 - stroke_wb_L_seq_bins_MIN_c4;

%% wing rotation
% L
pitch_wb_L_seq_bins_MAXmidDS_c4 = max(pitch_wb_L_seq_bins_mean_c4(40:100,:))'-90;
pitch_wb_L_seq_bins_MINmidUS_c4 = min(pitch_wb_L_seq_bins_mean_c4(150:end,:))'-90;
Apitch_wb_L_seq_bins_c4 = pitch_wb_L_seq_bins_MAXmidDS_c4 - pitch_wb_L_seq_bins_MINmidUS_c4;

% R
pitch_wb_R_seq_bins_MAXmidDS_c4 = max(pitch_wb_R_seq_bins_mean_c4(40:100,:))'-90;
pitch_wb_R_seq_bins_MINmidUS_c4 = min(pitch_wb_R_seq_bins_mean_c4(150:end,:))'-90;
Apitch_wb_R_seq_bins_c4 = pitch_wb_R_seq_bins_MAXmidDS_c4 - pitch_wb_R_seq_bins_MINmidUS_c4;

% R+L
Apitch_wb_seq_bins_c4 = nanmean([Apitch_wb_R_seq_bins_c4 Apitch_wb_L_seq_bins_c4]')';
pitch_wb_seq_bins_MAXmidDS_c4 = nanmean([pitch_wb_R_seq_bins_MAXmidDS_c4 pitch_wb_L_seq_bins_MAXmidDS_c4]')';
pitch_wb_seq_bins_MINmidUS_c4 = nanmean([pitch_wb_R_seq_bins_MINmidUS_c4 pitch_wb_L_seq_bins_MINmidUS_c4]')';

% R-L
dApitch_wb_seq_bins_c4 = Apitch_wb_R_seq_bins_c4 - Apitch_wb_L_seq_bins_c4;
Dpitch_wb_seq_bins_MAXmidDS_c4 = pitch_wb_R_seq_bins_MAXmidDS_c4 - pitch_wb_L_seq_bins_MAXmidDS_c4;
Dpitch_wb_seq_bins_MINmidUS_c4 = pitch_wb_R_seq_bins_MINmidUS_c4 - pitch_wb_L_seq_bins_MINmidUS_c4;

%% wing deviation
% L
dev_wb_L_seq_bins_MAXds_c4 = max(dev_wb_L_seq_bins_mean_c4(1:50,:))'-90;
dev_wb_L_seq_bins_MINds_c4 = min(dev_wb_L_seq_bins_mean_c4(1:100,:))'-90;
dev_wb_L_seq_bins_MAXus_c4 = max(dev_wb_L_seq_bins_mean_c4(50:150,:))'-90;
dev_wb_L_seq_bins_MINus_c4 = min(dev_wb_L_seq_bins_mean_c4(100:end,:))'-90;

Adev_wb_L_seq_bins_DS_c4 = dev_wb_L_seq_bins_MAXds_c4 - dev_wb_L_seq_bins_MINds_c4 ;
Adev_wb_L_seq_bins_US_c4 = dev_wb_L_seq_bins_MAXus_c4 - dev_wb_L_seq_bins_MINus_c4 ;

% R
dev_wb_R_seq_bins_MAXds_c4 = max(dev_wb_R_seq_bins_mean_c4(1:50,:))'-90;
dev_wb_R_seq_bins_MINds_c4 = min(dev_wb_R_seq_bins_mean_c4(1:100,:))'-90;
dev_wb_R_seq_bins_MAXus_c4 = max(dev_wb_R_seq_bins_mean_c4(50:150,:))'-90;
dev_wb_R_seq_bins_MINus_c4 = min(dev_wb_R_seq_bins_mean_c4(100:end,:))'-90;

Adev_wb_R_seq_bins_DS_c4 = dev_wb_R_seq_bins_MAXds_c4 - dev_wb_R_seq_bins_MINds_c4 ;
Adev_wb_R_seq_bins_US_c4 = dev_wb_R_seq_bins_MAXus_c4 - dev_wb_R_seq_bins_MINus_c4 ;

% R+L
Adev_wb_seq_bins_DS_c4 = nanmean([Adev_wb_R_seq_bins_DS_c4 Adev_wb_L_seq_bins_DS_c4]')';
Adev_wb_seq_bins_US_c4 = nanmean([Adev_wb_R_seq_bins_US_c4 Adev_wb_L_seq_bins_US_c4]')';

dev_wb_seq_bins_MAXds_c4 = nanmean([dev_wb_R_seq_bins_MAXds_c4 dev_wb_L_seq_bins_MAXds_c4]')';
dev_wb_seq_bins_MINds_c4 = nanmean([dev_wb_R_seq_bins_MINds_c4 dev_wb_L_seq_bins_MINds_c4 ]')';
dev_wb_seq_bins_MAXus_c4 = nanmean([dev_wb_R_seq_bins_MAXus_c4 dev_wb_L_seq_bins_MAXus_c4]')';
dev_wb_seq_bins_MINus_c4 = nanmean([dev_wb_R_seq_bins_MINus_c4 dev_wb_L_seq_bins_MINus_c4 ]')';

% R-L
dAdev_wb_seq_bins_DS_c4 = Adev_wb_R_seq_bins_DS_c4 - Adev_wb_L_seq_bins_DS_c4;
dAdev_wb_seq_bins_US_c4 = Adev_wb_R_seq_bins_US_c4 - Adev_wb_L_seq_bins_US_c4;

Ddev_wb_seq_bins_MAXds_c4 = dev_wb_R_seq_bins_MAXds_c4 - dev_wb_L_seq_bins_MAXds_c4;
Ddev_wb_seq_bins_MINds_c4 = dev_wb_R_seq_bins_MINds_c4 - dev_wb_L_seq_bins_MINds_c4;
Ddev_wb_seq_bins_MAXus_c4 = dev_wb_R_seq_bins_MAXus_c4 - dev_wb_L_seq_bins_MAXus_c4;
Ddev_wb_seq_bins_MINus_c4 = dev_wb_R_seq_bins_MINus_c4 - dev_wb_L_seq_bins_MINus_c4;
