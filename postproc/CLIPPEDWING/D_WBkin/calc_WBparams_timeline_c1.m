
%% stroke
% L
stroke_wb_L_seq_bins_MAX1_c1 = max(stroke_wb_L_seq_bins_mean_c1(1:100,:))';
stroke_wb_L_seq_bins_MAX2_c1 = max(stroke_wb_L_seq_bins_mean_c1(100:end,:))';
stroke_wb_L_seq_bins_MAX1_c1 = [stroke_wb_L_seq_bins_MAX1_c1;nan];
stroke_wb_L_seq_bins_MAX2_c1 = [nan;stroke_wb_L_seq_bins_MAX2_c1];
stroke_wb_L_seq_bins_MAX_c1 = [stroke_wb_L_seq_bins_MAX1_c1 stroke_wb_L_seq_bins_MAX2_c1];
stroke_wb_L_seq_bins_MAX_c1 = max(stroke_wb_L_seq_bins_MAX_c1(1:end-1,:)')';

stroke_wb_L_seq_bins_MIN_c1 = min(stroke_wb_L_seq_bins_mean_c1)';
Astroke_wb_L_seq_bins_c1 = stroke_wb_L_seq_bins_MAX_c1 - stroke_wb_L_seq_bins_MIN_c1;

% R
stroke_wb_R_seq_bins_MAX1_c1 = max(stroke_wb_R_seq_bins_mean_c1(1:100,:))';
stroke_wb_R_seq_bins_MAX2_c1 = max(stroke_wb_R_seq_bins_mean_c1(100:end,:))';
stroke_wb_R_seq_bins_MAX1_c1 = [stroke_wb_R_seq_bins_MAX1_c1;nan];
stroke_wb_R_seq_bins_MAX2_c1 = [nan;stroke_wb_R_seq_bins_MAX2_c1];
stroke_wb_R_seq_bins_MAX_c1 = [stroke_wb_R_seq_bins_MAX1_c1 stroke_wb_R_seq_bins_MAX2_c1];
stroke_wb_R_seq_bins_MAX_c1 = max(stroke_wb_R_seq_bins_MAX_c1(1:end-1,:)')';

stroke_wb_R_seq_bins_MIN_c1 = min(stroke_wb_R_seq_bins_mean_c1)';
Astroke_wb_R_seq_bins_c1 = stroke_wb_R_seq_bins_MAX_c1 - stroke_wb_R_seq_bins_MIN_c1;

% R+L
Astroke_wb_seq_bins_c1 = nanmean([Astroke_wb_R_seq_bins_c1 Astroke_wb_L_seq_bins_c1]')';
stroke_wb_seq_bins_MAX_c1 = nanmean([stroke_wb_R_seq_bins_MAX_c1 stroke_wb_L_seq_bins_MAX_c1]')';
stroke_wb_seq_bins_MIN_c1 = nanmean([stroke_wb_R_seq_bins_MIN_c1 stroke_wb_L_seq_bins_MIN_c1]')';

% R-L
dAstroke_wb_seq_bins_c1 = Astroke_wb_R_seq_bins_c1 - Astroke_wb_L_seq_bins_c1;
Dstroke_wb_seq_bins_MAX_c1 = stroke_wb_R_seq_bins_MAX_c1 - stroke_wb_L_seq_bins_MAX_c1;
Dstroke_wb_seq_bins_MIN_c1 = stroke_wb_R_seq_bins_MIN_c1 - stroke_wb_L_seq_bins_MIN_c1;

%% wing rotation
% L
pitch_wb_L_seq_bins_MAXmidDS_c1 = max(pitch_wb_L_seq_bins_mean_c1(40:100,:))'-90;
pitch_wb_L_seq_bins_MINmidUS_c1 = min(pitch_wb_L_seq_bins_mean_c1(150:end,:))'-90;
Apitch_wb_L_seq_bins_c1 = pitch_wb_L_seq_bins_MAXmidDS_c1 - pitch_wb_L_seq_bins_MINmidUS_c1;

% R
pitch_wb_R_seq_bins_MAXmidDS_c1 = max(pitch_wb_R_seq_bins_mean_c1(40:100,:))'-90;
pitch_wb_R_seq_bins_MINmidUS_c1 = min(pitch_wb_R_seq_bins_mean_c1(150:end,:))'-90;
Apitch_wb_R_seq_bins_c1 = pitch_wb_R_seq_bins_MAXmidDS_c1 - pitch_wb_R_seq_bins_MINmidUS_c1;

% R+L
Apitch_wb_seq_bins_c1 = nanmean([Apitch_wb_R_seq_bins_c1 Apitch_wb_L_seq_bins_c1]')';
pitch_wb_seq_bins_MAXmidDS_c1 = nanmean([pitch_wb_R_seq_bins_MAXmidDS_c1 pitch_wb_L_seq_bins_MAXmidDS_c1]')';
pitch_wb_seq_bins_MINmidUS_c1 = nanmean([pitch_wb_R_seq_bins_MINmidUS_c1 pitch_wb_L_seq_bins_MINmidUS_c1]')';

% R-L
dApitch_wb_seq_bins_c1 = Apitch_wb_R_seq_bins_c1 - Apitch_wb_L_seq_bins_c1;
Dpitch_wb_seq_bins_MAXmidDS_c1 = pitch_wb_R_seq_bins_MAXmidDS_c1 - pitch_wb_L_seq_bins_MAXmidDS_c1;
Dpitch_wb_seq_bins_MINmidUS_c1 = pitch_wb_R_seq_bins_MINmidUS_c1 - pitch_wb_L_seq_bins_MINmidUS_c1;

%% wing deviation
% L
dev_wb_L_seq_bins_MAXds_c1 = max(dev_wb_L_seq_bins_mean_c1(1:50,:))'-90;
dev_wb_L_seq_bins_MINds_c1 = min(dev_wb_L_seq_bins_mean_c1(1:100,:))'-90;
dev_wb_L_seq_bins_MAXus_c1 = max(dev_wb_L_seq_bins_mean_c1(50:150,:))'-90;
dev_wb_L_seq_bins_MINus_c1 = min(dev_wb_L_seq_bins_mean_c1(100:end,:))'-90;

Adev_wb_L_seq_bins_DS_c1 = dev_wb_L_seq_bins_MAXds_c1 - dev_wb_L_seq_bins_MINds_c1 ;
Adev_wb_L_seq_bins_US_c1 = dev_wb_L_seq_bins_MAXus_c1 - dev_wb_L_seq_bins_MINus_c1 ;

% R
dev_wb_R_seq_bins_MAXds_c1 = max(dev_wb_R_seq_bins_mean_c1(1:50,:))'-90;
dev_wb_R_seq_bins_MINds_c1 = min(dev_wb_R_seq_bins_mean_c1(1:100,:))'-90;
dev_wb_R_seq_bins_MAXus_c1 = max(dev_wb_R_seq_bins_mean_c1(50:150,:))'-90;
dev_wb_R_seq_bins_MINus_c1 = min(dev_wb_R_seq_bins_mean_c1(100:end,:))'-90;

Adev_wb_R_seq_bins_DS_c1 = dev_wb_R_seq_bins_MAXds_c1 - dev_wb_R_seq_bins_MINds_c1 ;
Adev_wb_R_seq_bins_US_c1 = dev_wb_R_seq_bins_MAXus_c1 - dev_wb_R_seq_bins_MINus_c1 ;

% R+L
Adev_wb_seq_bins_DS_c1 = nanmean([Adev_wb_R_seq_bins_DS_c1 Adev_wb_L_seq_bins_DS_c1]')';
Adev_wb_seq_bins_US_c1 = nanmean([Adev_wb_R_seq_bins_US_c1 Adev_wb_L_seq_bins_US_c1]')';

dev_wb_seq_bins_MAXds_c1 = nanmean([dev_wb_R_seq_bins_MAXds_c1 dev_wb_L_seq_bins_MAXds_c1]')';
dev_wb_seq_bins_MINds_c1 = nanmean([dev_wb_R_seq_bins_MINds_c1 dev_wb_L_seq_bins_MINds_c1 ]')';
dev_wb_seq_bins_MAXus_c1 = nanmean([dev_wb_R_seq_bins_MAXus_c1 dev_wb_L_seq_bins_MAXus_c1]')';
dev_wb_seq_bins_MINus_c1 = nanmean([dev_wb_R_seq_bins_MINus_c1 dev_wb_L_seq_bins_MINus_c1 ]')';

% R-L
dAdev_wb_seq_bins_DS_c1 = Adev_wb_R_seq_bins_DS_c1 - Adev_wb_L_seq_bins_DS_c1;
dAdev_wb_seq_bins_US_c1 = Adev_wb_R_seq_bins_US_c1 - Adev_wb_L_seq_bins_US_c1;

Ddev_wb_seq_bins_MAXds_c1 = dev_wb_R_seq_bins_MAXds_c1 - dev_wb_L_seq_bins_MAXds_c1;
Ddev_wb_seq_bins_MINds_c1 = dev_wb_R_seq_bins_MINds_c1 - dev_wb_L_seq_bins_MINds_c1;
Ddev_wb_seq_bins_MAXus_c1 = dev_wb_R_seq_bins_MAXus_c1 - dev_wb_L_seq_bins_MAXus_c1;
Ddev_wb_seq_bins_MINus_c1 = dev_wb_R_seq_bins_MINus_c1 - dev_wb_L_seq_bins_MINus_c1;
