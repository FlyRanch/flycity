
%% stroke
% L
stroke_wb_L_seq_bins_MAX1_c2 = max(stroke_wb_L_seq_bins_mean_c2(1:100,:))';
stroke_wb_L_seq_bins_MAX2_c2 = max(stroke_wb_L_seq_bins_mean_c2(100:end,:))';
stroke_wb_L_seq_bins_MAX1_c2 = [stroke_wb_L_seq_bins_MAX1_c2;nan];
stroke_wb_L_seq_bins_MAX2_c2 = [nan;stroke_wb_L_seq_bins_MAX2_c2];
stroke_wb_L_seq_bins_MAX_c2 = [stroke_wb_L_seq_bins_MAX1_c2 stroke_wb_L_seq_bins_MAX2_c2];
stroke_wb_L_seq_bins_MAX_c2 = max(stroke_wb_L_seq_bins_MAX_c2(1:end-1,:)')';

stroke_wb_L_seq_bins_MIN_c2 = min(stroke_wb_L_seq_bins_mean_c2)';
Astroke_wb_L_seq_bins_c2 = stroke_wb_L_seq_bins_MAX_c2 - stroke_wb_L_seq_bins_MIN_c2;

% R
stroke_wb_R_seq_bins_MAX1_c2 = max(stroke_wb_R_seq_bins_mean_c2(1:100,:))';
stroke_wb_R_seq_bins_MAX2_c2 = max(stroke_wb_R_seq_bins_mean_c2(100:end,:))';
stroke_wb_R_seq_bins_MAX1_c2 = [stroke_wb_R_seq_bins_MAX1_c2;nan];
stroke_wb_R_seq_bins_MAX2_c2 = [nan;stroke_wb_R_seq_bins_MAX2_c2];
stroke_wb_R_seq_bins_MAX_c2 = [stroke_wb_R_seq_bins_MAX1_c2 stroke_wb_R_seq_bins_MAX2_c2];
stroke_wb_R_seq_bins_MAX_c2 = max(stroke_wb_R_seq_bins_MAX_c2(1:end-1,:)')';

stroke_wb_R_seq_bins_MIN_c2 = min(stroke_wb_R_seq_bins_mean_c2)';
Astroke_wb_R_seq_bins_c2 = stroke_wb_R_seq_bins_MAX_c2 - stroke_wb_R_seq_bins_MIN_c2;

% R+L
Astroke_wb_seq_bins_c2 = nanmean([Astroke_wb_R_seq_bins_c2 Astroke_wb_L_seq_bins_c2]')';
stroke_wb_seq_bins_MAX_c2 = nanmean([stroke_wb_R_seq_bins_MAX_c2 stroke_wb_L_seq_bins_MAX_c2]')';
stroke_wb_seq_bins_MIN_c2 = nanmean([stroke_wb_R_seq_bins_MIN_c2 stroke_wb_L_seq_bins_MIN_c2]')';

% R-L
dAstroke_wb_seq_bins_c2 = Astroke_wb_R_seq_bins_c2 - Astroke_wb_L_seq_bins_c2;
Dstroke_wb_seq_bins_MAX_c2 = stroke_wb_R_seq_bins_MAX_c2 - stroke_wb_L_seq_bins_MAX_c2;
Dstroke_wb_seq_bins_MIN_c2 = stroke_wb_R_seq_bins_MIN_c2 - stroke_wb_L_seq_bins_MIN_c2;

%% wing rotation
% L
pitch_wb_L_seq_bins_MAXmidDS_c2 = max(pitch_wb_L_seq_bins_mean_c2(40:100,:))'-90;
pitch_wb_L_seq_bins_MINmidUS_c2 = min(pitch_wb_L_seq_bins_mean_c2(150:end,:))'-90;
Apitch_wb_L_seq_bins_c2 = pitch_wb_L_seq_bins_MAXmidDS_c2 - pitch_wb_L_seq_bins_MINmidUS_c2;

% R
pitch_wb_R_seq_bins_MAXmidDS_c2 = max(pitch_wb_R_seq_bins_mean_c2(40:100,:))'-90;
pitch_wb_R_seq_bins_MINmidUS_c2 = min(pitch_wb_R_seq_bins_mean_c2(150:end,:))'-90;
Apitch_wb_R_seq_bins_c2 = pitch_wb_R_seq_bins_MAXmidDS_c2 - pitch_wb_R_seq_bins_MINmidUS_c2;

% R+L
Apitch_wb_seq_bins_c2 = nanmean([Apitch_wb_R_seq_bins_c2 Apitch_wb_L_seq_bins_c2]')';
pitch_wb_seq_bins_MAXmidDS_c2 = nanmean([pitch_wb_R_seq_bins_MAXmidDS_c2 pitch_wb_L_seq_bins_MAXmidDS_c2]')';
pitch_wb_seq_bins_MINmidUS_c2 = nanmean([pitch_wb_R_seq_bins_MINmidUS_c2 pitch_wb_L_seq_bins_MINmidUS_c2]')';

% R-L
dApitch_wb_seq_bins_c2 = Apitch_wb_R_seq_bins_c2 - Apitch_wb_L_seq_bins_c2;
Dpitch_wb_seq_bins_MAXmidDS_c2 = pitch_wb_R_seq_bins_MAXmidDS_c2 - pitch_wb_L_seq_bins_MAXmidDS_c2;
Dpitch_wb_seq_bins_MINmidUS_c2 = pitch_wb_R_seq_bins_MINmidUS_c2 - pitch_wb_L_seq_bins_MINmidUS_c2;

%% wing deviation
% L
dev_wb_L_seq_bins_MAXds_c2 = max(dev_wb_L_seq_bins_mean_c2(1:50,:))'-90;
dev_wb_L_seq_bins_MINds_c2 = min(dev_wb_L_seq_bins_mean_c2(1:100,:))'-90;
dev_wb_L_seq_bins_MAXus_c2 = max(dev_wb_L_seq_bins_mean_c2(50:150,:))'-90;
dev_wb_L_seq_bins_MINus_c2 = min(dev_wb_L_seq_bins_mean_c2(100:end,:))'-90;

Adev_wb_L_seq_bins_DS_c2 = dev_wb_L_seq_bins_MAXds_c2 - dev_wb_L_seq_bins_MINds_c2 ;
Adev_wb_L_seq_bins_US_c2 = dev_wb_L_seq_bins_MAXus_c2 - dev_wb_L_seq_bins_MINus_c2 ;

% R
dev_wb_R_seq_bins_MAXds_c2 = max(dev_wb_R_seq_bins_mean_c2(1:50,:))'-90;
dev_wb_R_seq_bins_MINds_c2 = min(dev_wb_R_seq_bins_mean_c2(1:100,:))'-90;
dev_wb_R_seq_bins_MAXus_c2 = max(dev_wb_R_seq_bins_mean_c2(50:150,:))'-90;
dev_wb_R_seq_bins_MINus_c2 = min(dev_wb_R_seq_bins_mean_c2(100:end,:))'-90;

Adev_wb_R_seq_bins_DS_c2 = dev_wb_R_seq_bins_MAXds_c2 - dev_wb_R_seq_bins_MINds_c2 ;
Adev_wb_R_seq_bins_US_c2 = dev_wb_R_seq_bins_MAXus_c2 - dev_wb_R_seq_bins_MINus_c2 ;

% R+L
Adev_wb_seq_bins_DS_c2 = nanmean([Adev_wb_R_seq_bins_DS_c2 Adev_wb_L_seq_bins_DS_c2]')';
Adev_wb_seq_bins_US_c2 = nanmean([Adev_wb_R_seq_bins_US_c2 Adev_wb_L_seq_bins_US_c2]')';

dev_wb_seq_bins_MAXds_c2 = nanmean([dev_wb_R_seq_bins_MAXds_c2 dev_wb_L_seq_bins_MAXds_c2]')';
dev_wb_seq_bins_MINds_c2 = nanmean([dev_wb_R_seq_bins_MINds_c2 dev_wb_L_seq_bins_MINds_c2 ]')';
dev_wb_seq_bins_MAXus_c2 = nanmean([dev_wb_R_seq_bins_MAXus_c2 dev_wb_L_seq_bins_MAXus_c2]')';
dev_wb_seq_bins_MINus_c2 = nanmean([dev_wb_R_seq_bins_MINus_c2 dev_wb_L_seq_bins_MINus_c2 ]')';

% R-L
dAdev_wb_seq_bins_DS_c2 = Adev_wb_R_seq_bins_DS_c2 - Adev_wb_L_seq_bins_DS_c2;
dAdev_wb_seq_bins_US_c2 = Adev_wb_R_seq_bins_US_c2 - Adev_wb_L_seq_bins_US_c2;

Ddev_wb_seq_bins_MAXds_c2 = dev_wb_R_seq_bins_MAXds_c2 - dev_wb_L_seq_bins_MAXds_c2;
Ddev_wb_seq_bins_MINds_c2 = dev_wb_R_seq_bins_MINds_c2 - dev_wb_L_seq_bins_MINds_c2;
Ddev_wb_seq_bins_MAXus_c2 = dev_wb_R_seq_bins_MAXus_c2 - dev_wb_L_seq_bins_MAXus_c2;
Ddev_wb_seq_bins_MINus_c2 = dev_wb_R_seq_bins_MINus_c2 - dev_wb_L_seq_bins_MINus_c2;
