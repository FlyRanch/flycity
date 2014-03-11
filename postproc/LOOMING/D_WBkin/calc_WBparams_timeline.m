
%% stroke
% L
stroke_wb_L_seq_bins_MAX1 = max(stroke_wb_L_seq_bins_mean(1:100,:))';
stroke_wb_L_seq_bins_MAX2 = max(stroke_wb_L_seq_bins_mean(100:end,:))';
stroke_wb_L_seq_bins_MAX1 = [stroke_wb_L_seq_bins_MAX1;nan];
stroke_wb_L_seq_bins_MAX2 = [nan;stroke_wb_L_seq_bins_MAX2];
stroke_wb_L_seq_bins_MAX = [stroke_wb_L_seq_bins_MAX1 stroke_wb_L_seq_bins_MAX2];
stroke_wb_L_seq_bins_MAX = max(stroke_wb_L_seq_bins_MAX(1:end-1,:)')';

stroke_wb_L_seq_bins_MIN = min(stroke_wb_L_seq_bins_mean)';
Astroke_wb_L_seq_bins = stroke_wb_L_seq_bins_MAX - stroke_wb_L_seq_bins_MIN;

% R
stroke_wb_R_seq_bins_MAX1 = max(stroke_wb_R_seq_bins_mean(1:100,:))';
stroke_wb_R_seq_bins_MAX2 = max(stroke_wb_R_seq_bins_mean(100:end,:))';
stroke_wb_R_seq_bins_MAX1 = [stroke_wb_R_seq_bins_MAX1;nan];
stroke_wb_R_seq_bins_MAX2 = [nan;stroke_wb_R_seq_bins_MAX2];
stroke_wb_R_seq_bins_MAX = [stroke_wb_R_seq_bins_MAX1 stroke_wb_R_seq_bins_MAX2];
stroke_wb_R_seq_bins_MAX = max(stroke_wb_R_seq_bins_MAX(1:end-1,:)')';

stroke_wb_R_seq_bins_MIN = min(stroke_wb_R_seq_bins_mean)';
Astroke_wb_R_seq_bins = stroke_wb_R_seq_bins_MAX - stroke_wb_R_seq_bins_MIN;

% R+L
Astroke_wb_seq_bins = nanmean([Astroke_wb_R_seq_bins Astroke_wb_L_seq_bins]')';
stroke_wb_seq_bins_MAX = nanmean([stroke_wb_R_seq_bins_MAX stroke_wb_L_seq_bins_MAX]')';
stroke_wb_seq_bins_MIN = nanmean([stroke_wb_R_seq_bins_MIN stroke_wb_L_seq_bins_MIN]')';

% R-L
dAstroke_wb_seq_bins = Astroke_wb_R_seq_bins - Astroke_wb_L_seq_bins;
Dstroke_wb_seq_bins_MAX = stroke_wb_R_seq_bins_MAX - stroke_wb_L_seq_bins_MAX;
Dstroke_wb_seq_bins_MIN = stroke_wb_R_seq_bins_MIN - stroke_wb_L_seq_bins_MIN;

%% wing rotation
% L
pitch_wb_L_seq_bins_MAXmidDS = max(pitch_wb_L_seq_bins_mean(40:100,:))'-90;
pitch_wb_L_seq_bins_MINmidUS = min(pitch_wb_L_seq_bins_mean(150:end,:))'-90;
Apitch_wb_L_seq_bins = pitch_wb_L_seq_bins_MAXmidDS - pitch_wb_L_seq_bins_MINmidUS;

% R
pitch_wb_R_seq_bins_MAXmidDS = max(pitch_wb_R_seq_bins_mean(40:100,:))'-90;
pitch_wb_R_seq_bins_MINmidUS = min(pitch_wb_R_seq_bins_mean(150:end,:))'-90;
Apitch_wb_R_seq_bins = pitch_wb_R_seq_bins_MAXmidDS - pitch_wb_R_seq_bins_MINmidUS;

% R+L
Apitch_wb_seq_bins = nanmean([Apitch_wb_R_seq_bins Apitch_wb_L_seq_bins]')';
pitch_wb_seq_bins_MAXmidDS = nanmean([pitch_wb_R_seq_bins_MAXmidDS pitch_wb_L_seq_bins_MAXmidDS]')';
pitch_wb_seq_bins_MINmidUS = nanmean([pitch_wb_R_seq_bins_MINmidUS pitch_wb_L_seq_bins_MINmidUS]')';

% R-L
dApitch_wb_seq_bins = Apitch_wb_R_seq_bins - Apitch_wb_L_seq_bins;
Dpitch_wb_seq_bins_MAXmidDS = pitch_wb_R_seq_bins_MAXmidDS - pitch_wb_L_seq_bins_MAXmidDS;
Dpitch_wb_seq_bins_MINmidUS = pitch_wb_R_seq_bins_MINmidUS - pitch_wb_L_seq_bins_MINmidUS;

%% wing deviation
% L
dev_wb_L_seq_bins_MAXds = max(dev_wb_L_seq_bins_mean(1:50,:))'-90;
dev_wb_L_seq_bins_MINds = min(dev_wb_L_seq_bins_mean(1:100,:))'-90;
dev_wb_L_seq_bins_MAXus = max(dev_wb_L_seq_bins_mean(50:150,:))'-90;
dev_wb_L_seq_bins_MINus = min(dev_wb_L_seq_bins_mean(100:end,:))'-90;

Adev_wb_L_seq_bins_DS = dev_wb_L_seq_bins_MAXds - dev_wb_L_seq_bins_MINds;
Adev_wb_L_seq_bins_US = dev_wb_L_seq_bins_MAXus - dev_wb_L_seq_bins_MINus;

% R
dev_wb_R_seq_bins_MAXds = max(dev_wb_R_seq_bins_mean(1:50,:))'-90;
dev_wb_R_seq_bins_MINds = min(dev_wb_R_seq_bins_mean(1:100,:))'-90;
dev_wb_R_seq_bins_MAXus = max(dev_wb_R_seq_bins_mean(50:150,:))'-90;
dev_wb_R_seq_bins_MINus = min(dev_wb_R_seq_bins_mean(100:end,:))'-90;

Adev_wb_R_seq_bins_DS = dev_wb_R_seq_bins_MAXds - dev_wb_R_seq_bins_MINds;
Adev_wb_R_seq_bins_US = dev_wb_R_seq_bins_MAXus - dev_wb_R_seq_bins_MINus;

% R+L
Adev_wb_seq_bins_DS = nanmean([Adev_wb_R_seq_bins_DS Adev_wb_L_seq_bins_DS]')';
Adev_wb_seq_bins_US = nanmean([Adev_wb_R_seq_bins_US Adev_wb_L_seq_bins_US]')';

dev_wb_seq_bins_MAXds = nanmean([dev_wb_R_seq_bins_MAXds dev_wb_L_seq_bins_MAXds]')';
dev_wb_seq_bins_MINds = nanmean([dev_wb_R_seq_bins_MINds dev_wb_L_seq_bins_MINds]')';
dev_wb_seq_bins_MAXus = nanmean([dev_wb_R_seq_bins_MAXus dev_wb_L_seq_bins_MAXus]')';
dev_wb_seq_bins_MINus = nanmean([dev_wb_R_seq_bins_MINus dev_wb_L_seq_bins_MINus]')';

% R-L
dAdev_wb_seq_bins_DS = Adev_wb_R_seq_bins_DS - Adev_wb_L_seq_bins_DS;
dAdev_wb_seq_bins_US = Adev_wb_R_seq_bins_US - Adev_wb_L_seq_bins_US;

Ddev_wb_seq_bins_MAXds = dev_wb_R_seq_bins_MAXds - dev_wb_L_seq_bins_MAXds;
Ddev_wb_seq_bins_MINds = dev_wb_R_seq_bins_MINds - dev_wb_L_seq_bins_MINds;
Ddev_wb_seq_bins_MAXus = dev_wb_R_seq_bins_MAXus - dev_wb_L_seq_bins_MAXus;
Ddev_wb_seq_bins_MINus = dev_wb_R_seq_bins_MINus - dev_wb_L_seq_bins_MINus;
