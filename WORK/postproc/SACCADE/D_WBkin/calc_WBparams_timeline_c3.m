
%% stroke
% L
stroke_wb_L_seq_bins_MAX1_c3 = max(stroke_wb_L_seq_bins_mean_c3(1:100,:))';
stroke_wb_L_seq_bins_MAX2_c3 = max(stroke_wb_L_seq_bins_mean_c3(100:end,:))';
stroke_wb_L_seq_bins_MAX1_c3 = [stroke_wb_L_seq_bins_MAX1_c3;nan];
stroke_wb_L_seq_bins_MAX2_c3 = [nan;stroke_wb_L_seq_bins_MAX2_c3];
stroke_wb_L_seq_bins_MAX_c3 = [stroke_wb_L_seq_bins_MAX1_c3 stroke_wb_L_seq_bins_MAX2_c3];
stroke_wb_L_seq_bins_MAX_c3 = max(stroke_wb_L_seq_bins_MAX_c3(1:end-1,:)')';

stroke_wb_L_seq_bins_MIN_c3 = min(stroke_wb_L_seq_bins_mean_c3)';
Astroke_wb_L_seq_bins_c3 = stroke_wb_L_seq_bins_MAX_c3 - stroke_wb_L_seq_bins_MIN_c3;

% R
stroke_wb_R_seq_bins_MAX1_c3 = max(stroke_wb_R_seq_bins_mean_c3(1:100,:))';
stroke_wb_R_seq_bins_MAX2_c3 = max(stroke_wb_R_seq_bins_mean_c3(100:end,:))';
stroke_wb_R_seq_bins_MAX1_c3 = [stroke_wb_R_seq_bins_MAX1_c3;nan];
stroke_wb_R_seq_bins_MAX2_c3 = [nan;stroke_wb_R_seq_bins_MAX2_c3];
stroke_wb_R_seq_bins_MAX_c3 = [stroke_wb_R_seq_bins_MAX1_c3 stroke_wb_R_seq_bins_MAX2_c3];
stroke_wb_R_seq_bins_MAX_c3 = max(stroke_wb_R_seq_bins_MAX_c3(1:end-1,:)')';

stroke_wb_R_seq_bins_MIN_c3 = min(stroke_wb_R_seq_bins_mean_c3)';
Astroke_wb_R_seq_bins_c3 = stroke_wb_R_seq_bins_MAX_c3 - stroke_wb_R_seq_bins_MIN_c3;

% R+L
Astroke_wb_seq_bins_c3 = nanmean([Astroke_wb_R_seq_bins_c3 Astroke_wb_L_seq_bins_c3]')';
stroke_wb_seq_bins_MAX_c3 = nanmean([stroke_wb_R_seq_bins_MAX_c3 stroke_wb_L_seq_bins_MAX_c3]')';
stroke_wb_seq_bins_MIN_c3 = nanmean([stroke_wb_R_seq_bins_MIN_c3 stroke_wb_L_seq_bins_MIN_c3]')';

% R-L
dAstroke_wb_seq_bins_c3 = Astroke_wb_R_seq_bins_c3 - Astroke_wb_L_seq_bins_c3;
Dstroke_wb_seq_bins_MAX_c3 = stroke_wb_R_seq_bins_MAX_c3 - stroke_wb_L_seq_bins_MAX_c3;
Dstroke_wb_seq_bins_MIN_c3 = stroke_wb_R_seq_bins_MIN_c3 - stroke_wb_L_seq_bins_MIN_c3;

%% wing rotation
% L
pitch_wb_L_seq_bins_MAXmidDS_c3 = max(pitch_wb_L_seq_bins_mean_c3(40:100,:))'-90;
pitch_wb_L_seq_bins_MINmidUS_c3 = min(pitch_wb_L_seq_bins_mean_c3(150:end,:))'-90;
Apitch_wb_L_seq_bins_c3 = pitch_wb_L_seq_bins_MAXmidDS_c3 - pitch_wb_L_seq_bins_MINmidUS_c3;

% R
pitch_wb_R_seq_bins_MAXmidDS_c3 = max(pitch_wb_R_seq_bins_mean_c3(40:100,:))'-90;
pitch_wb_R_seq_bins_MINmidUS_c3 = min(pitch_wb_R_seq_bins_mean_c3(150:end,:))'-90;
Apitch_wb_R_seq_bins_c3 = pitch_wb_R_seq_bins_MAXmidDS_c3 - pitch_wb_R_seq_bins_MINmidUS_c3;

% R+L
Apitch_wb_seq_bins_c3 = nanmean([Apitch_wb_R_seq_bins_c3 Apitch_wb_L_seq_bins_c3]')';
pitch_wb_seq_bins_MAXmidDS_c3 = nanmean([pitch_wb_R_seq_bins_MAXmidDS_c3 pitch_wb_L_seq_bins_MAXmidDS_c3]')';
pitch_wb_seq_bins_MINmidUS_c3 = nanmean([pitch_wb_R_seq_bins_MINmidUS_c3 pitch_wb_L_seq_bins_MINmidUS_c3]')';

% R-L
dApitch_wb_seq_bins_c3 = Apitch_wb_R_seq_bins_c3 - Apitch_wb_L_seq_bins_c3;
Dpitch_wb_seq_bins_MAXmidDS_c3 = pitch_wb_R_seq_bins_MAXmidDS_c3 - pitch_wb_L_seq_bins_MAXmidDS_c3;
Dpitch_wb_seq_bins_MINmidUS_c3 = pitch_wb_R_seq_bins_MINmidUS_c3 - pitch_wb_L_seq_bins_MINmidUS_c3;

%% wing deviation
% L
dev_wb_L_seq_bins_MAXds_c3 = max(dev_wb_L_seq_bins_mean_c3(1:50,:))'-90;
dev_wb_L_seq_bins_MINds_c3 = min(dev_wb_L_seq_bins_mean_c3(1:100,:))'-90;
dev_wb_L_seq_bins_MAXus_c3 = max(dev_wb_L_seq_bins_mean_c3(50:150,:))'-90;
dev_wb_L_seq_bins_MINus_c3 = min(dev_wb_L_seq_bins_mean_c3(100:end,:))'-90;

Adev_wb_L_seq_bins_DS_c3 = dev_wb_L_seq_bins_MAXds_c3 - dev_wb_L_seq_bins_MINds_c3 ;
Adev_wb_L_seq_bins_US_c3 = dev_wb_L_seq_bins_MAXus_c3 - dev_wb_L_seq_bins_MINus_c3 ;

% R
dev_wb_R_seq_bins_MAXds_c3 = max(dev_wb_R_seq_bins_mean_c3(1:50,:))'-90;
dev_wb_R_seq_bins_MINds_c3 = min(dev_wb_R_seq_bins_mean_c3(1:100,:))'-90;
dev_wb_R_seq_bins_MAXus_c3 = max(dev_wb_R_seq_bins_mean_c3(50:150,:))'-90;
dev_wb_R_seq_bins_MINus_c3 = min(dev_wb_R_seq_bins_mean_c3(100:end,:))'-90;

Adev_wb_R_seq_bins_DS_c3 = dev_wb_R_seq_bins_MAXds_c3 - dev_wb_R_seq_bins_MINds_c3 ;
Adev_wb_R_seq_bins_US_c3 = dev_wb_R_seq_bins_MAXus_c3 - dev_wb_R_seq_bins_MINus_c3 ;

% R+L
Adev_wb_seq_bins_DS_c3 = nanmean([Adev_wb_R_seq_bins_DS_c3 Adev_wb_L_seq_bins_DS_c3]')';
Adev_wb_seq_bins_US_c3 = nanmean([Adev_wb_R_seq_bins_US_c3 Adev_wb_L_seq_bins_US_c3]')';

dev_wb_seq_bins_MAXds_c3 = nanmean([dev_wb_R_seq_bins_MAXds_c3 dev_wb_L_seq_bins_MAXds_c3]')';
dev_wb_seq_bins_MINds_c3 = nanmean([dev_wb_R_seq_bins_MINds_c3 dev_wb_L_seq_bins_MINds_c3 ]')';
dev_wb_seq_bins_MAXus_c3 = nanmean([dev_wb_R_seq_bins_MAXus_c3 dev_wb_L_seq_bins_MAXus_c3]')';
dev_wb_seq_bins_MINus_c3 = nanmean([dev_wb_R_seq_bins_MINus_c3 dev_wb_L_seq_bins_MINus_c3 ]')';

% R-L
dAdev_wb_seq_bins_DS_c3 = Adev_wb_R_seq_bins_DS_c3 - Adev_wb_L_seq_bins_DS_c3;
dAdev_wb_seq_bins_US_c3 = Adev_wb_R_seq_bins_US_c3 - Adev_wb_L_seq_bins_US_c3;

Ddev_wb_seq_bins_MAXds_c3 = dev_wb_R_seq_bins_MAXds_c3 - dev_wb_L_seq_bins_MAXds_c3;
Ddev_wb_seq_bins_MINds_c3 = dev_wb_R_seq_bins_MINds_c3 - dev_wb_L_seq_bins_MINds_c3;
Ddev_wb_seq_bins_MAXus_c3 = dev_wb_R_seq_bins_MAXus_c3 - dev_wb_L_seq_bins_MAXus_c3;
Ddev_wb_seq_bins_MINus_c3 = dev_wb_R_seq_bins_MINus_c3 - dev_wb_L_seq_bins_MINus_c3;
