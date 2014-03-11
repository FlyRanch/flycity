t_wb_bin = [0:1/(n-1):1]';
t_ds_bin = [-1:1/(n-1):0]';
t_us_bin = [0:1/(n-1):1]';

% seq data
n_wb_L = nan(size(n_down_start_L));
t_wb_L = nan(size(n_down_start_L));
IDX_wb_L = nan(size(n_down_start_L));
n_wb_R = nan(size(n_down_start_L));
t_wb_R = nan(size(n_down_start_L));
IDX_wb_R = nan(size(n_down_start_L));

% body kin data
V_mean_wb = nan(size(n_down_start_L));
pitch_global_mean_wb = nan(size(n_down_start_L));

% freq data
f_wb_L = nan(size(n_down_start_L));
f_wb_R = nan(size(n_down_start_L));
f_ds_L = nan(size(n_down_start_L));
f_ds_R = nan(size(n_down_start_L));
f_us_L = nan(size(n_down_start_L));
f_us_R = nan(size(n_down_start_L));

dt_wb_L = nan(size(n_down_start_L));
dt_wb_R = nan(size(n_down_start_L));
dt_ds_L = nan(size(n_down_start_L));
dt_ds_R = nan(size(n_down_start_L));
dt_us_L = nan(size(n_down_start_L));
dt_us_R = nan(size(n_down_start_L));

df_wb = nan(size(n_down_start_L));
df_ds = nan(size(n_down_start_L));
df_us = nan(size(n_down_start_L));

ddt_wb = nan(size(n_down_start_L));
ddt_ds = nan(size(n_down_start_L));
ddt_us = nan(size(n_down_start_L));

% stroke data
stroke_mean_wb_L = nan(size(n_down_start_L));
stroke_mean_ds_L = nan(size(n_down_start_L));
stroke_mean_us_L = nan(size(n_down_start_L));
stroke_mean_wb_R = nan(size(n_down_start_L));
stroke_mean_ds_R = nan(size(n_down_start_L));
stroke_mean_us_R = nan(size(n_down_start_L));

stroke_max_wb_L = nan(size(n_down_start_L));
stroke_min_wb_L = nan(size(n_down_start_L));
stroke_max_wb_R = nan(size(n_down_start_L));
stroke_min_wb_R = nan(size(n_down_start_L));

stroke_max_udsPREV_L = nan(size(n_down_start_L));
stroke_max_udsNEXT_L = nan(size(n_down_start_L));
stroke_max_udsPREV_R = nan(size(n_down_start_L));
stroke_max_udsNEXT_R = nan(size(n_down_start_L));

Astroke_wb_L = nan(size(n_down_start_L));
Astroke_wb_R = nan(size(n_down_start_L));
dAstroke_wb = nan(size(n_down_start_L));

Astroke_ds_L = nan(size(n_down_start_L));
Astroke_ds_R = nan(size(n_down_start_L));
dAstroke_ds = nan(size(n_down_start_L));

Astroke_us_L = nan(size(n_down_start_L));
Astroke_us_R = nan(size(n_down_start_L));
dAstroke_us = nan(size(n_down_start_L));

dstroke_mean_wb = nan(size(n_down_start_L));
dstroke_mean_ds = nan(size(n_down_start_L));
dstroke_mean_us = nan(size(n_down_start_L));

dstroke_max_wb = nan(size(n_down_start_L));
dstroke_min_wb = nan(size(n_down_start_L));
dstroke_max_udsPREV = nan(size(n_down_start_L));
dstroke_max_udsNEXT = nan(size(n_down_start_L));

% pitch data
pitch_mean_wb_L = nan(size(n_down_start_L));
pitch_mean_wb_R = nan(size(n_down_start_L));
dpitch_mean_wb = nan(size(n_down_start_L));

pitch_max_wb_L = nan(size(n_down_start_L));
pitch_max_wb_R = nan(size(n_down_start_L));
dpitch_max_wb = nan(size(n_down_start_L));

pitch_min_wb_L = nan(size(n_down_start_L));
pitch_min_wb_R = nan(size(n_down_start_L));
dpitch_min_wb = nan(size(n_down_start_L));

pitch_mean_ds_L = nan(size(n_down_start_L));
pitch_mean_ds_R = nan(size(n_down_start_L));
dpitch_mean_ds = nan(size(n_down_start_L));

pitch_mean_us_L = nan(size(n_down_start_L));
pitch_mean_us_R = nan(size(n_down_start_L));
dpitch_mean_us = nan(size(n_down_start_L));

pitch_max_ds_L = nan(size(n_down_start_L));
pitch_max_ds_R = nan(size(n_down_start_L));
dpitch_max_ds = nan(size(n_down_start_L));

pitch_min_us_L = nan(size(n_down_start_L));
pitch_min_us_R = nan(size(n_down_start_L));
dpitch_min_us = nan(size(n_down_start_L));

Apitch_wb_L = nan(size(n_down_start_L));
Apitch_wb_R = nan(size(n_down_start_L));
dApitch_wb = nan(size(n_down_start_L));

Apitch_ds_L = nan(size(n_down_start_L));
Apitch_ds_R = nan(size(n_down_start_L));
dApitch_ds = nan(size(n_down_start_L));

Apitch_us_L = nan(size(n_down_start_L));
Apitch_us_R = nan(size(n_down_start_L));
dApitch_us = nan(size(n_down_start_L));

pitch_mid_ds_L = nan(size(n_down_start_L));
pitch_mid_ds_R = nan(size(n_down_start_L));

pitch_mid_us_L = nan(size(n_down_start_L));
pitch_mid_us_R = nan(size(n_down_start_L));

Apitch_mid_L = nan(size(n_down_start_L));
Apitch_mid_R = nan(size(n_down_start_L));

dpitch_mid_ds = nan(size(n_down_start_L));
dpitch_mid_us = nan(size(n_down_start_L));
dApitch_mid = nan(size(n_down_start_L));

% dev data
dev_mean_wb_L = nan(size(n_down_start_L));
dev_mean_ds_L = nan(size(n_down_start_L));
dev_mean_us_L = nan(size(n_down_start_L));
dev_mean_dus_L = nan(size(n_down_start_L));
dev_mean_udsPREV_L = nan(size(n_down_start_L));
dev_mean_udsNEXT_L = nan(size(n_down_start_L));

dev_mean_wb_R = nan(size(n_down_start_L));
dev_mean_ds_R = nan(size(n_down_start_L));
dev_mean_us_R = nan(size(n_down_start_L));
dev_mean_dus_R = nan(size(n_down_start_L));
dev_mean_udsPREV_R = nan(size(n_down_start_L));
dev_mean_udsNEXT_R = nan(size(n_down_start_L));

dev_max_dus_L = nan(size(n_down_start_L));
dev_max_udsPREV_L = nan(size(n_down_start_L));
dev_max_udsNEXT_L = nan(size(n_down_start_L));
dev_min_ds_L = nan(size(n_down_start_L));
dev_min_us_L = nan(size(n_down_start_L));

dev_max_dus_R = nan(size(n_down_start_L));
dev_max_udsPREV_R = nan(size(n_down_start_L));
dev_max_udsNEXT_R = nan(size(n_down_start_L));
dev_min_ds_R = nan(size(n_down_start_L));
dev_min_us_R = nan(size(n_down_start_L));

Adev_ds_L = nan(size(n_down_start_L));
Adev_us_L = nan(size(n_down_start_L));
Adev_ds_R = nan(size(n_down_start_L));
Adev_us_R = nan(size(n_down_start_L));
dAdev_ds = nan(size(n_down_start_L));
dAdev_us = nan(size(n_down_start_L));

ddev_mean_wb = nan(size(n_down_start_L));
ddev_mean_ds = nan(size(n_down_start_L));
ddev_mean_us = nan(size(n_down_start_L));
ddev_mean_dus = nan(size(n_down_start_L));
ddev_mean_udsPREV = nan(size(n_down_start_L));
ddev_mean_udsNEXT = nan(size(n_down_start_L));

ddev_max_dus = nan(size(n_down_start_L));
ddev_max_udsPREV = nan(size(n_down_start_L));
ddev_max_udsNEXT = nan(size(n_down_start_L));
ddev_min_ds = nan(size(n_down_start_L));
ddev_min_us = nan(size(n_down_start_L));

% Uwing data
U_mean_wb_L = nan(size(n_down_start_L));
U_mean_ds_L = nan(size(n_down_start_L));
U_mean_us_L = nan(size(n_down_start_L));
U_mean_dus_L = nan(size(n_down_start_L));
U_mean_udsPREV_L = nan(size(n_down_start_L));
U_mean_udsNEXT_L = nan(size(n_down_start_L));

U_mean_wb_R = nan(size(n_down_start_L));
U_mean_ds_R = nan(size(n_down_start_L));
U_mean_us_R = nan(size(n_down_start_L));
U_mean_dus_R = nan(size(n_down_start_L));
U_mean_udsPREV_R = nan(size(n_down_start_L));
U_mean_udsNEXT_R = nan(size(n_down_start_L));

U_max_dus_L = nan(size(n_down_start_L));
U_max_udsPREV_L = nan(size(n_down_start_L));
U_max_udsNEXT_L = nan(size(n_down_start_L));
U_min_ds_L = nan(size(n_down_start_L));
U_min_us_L = nan(size(n_down_start_L));

U_max_dus_R = nan(size(n_down_start_L));
U_max_udsPREV_R = nan(size(n_down_start_L));
U_max_udsNEXT_R = nan(size(n_down_start_L));
U_min_ds_R = nan(size(n_down_start_L));
U_min_us_R = nan(size(n_down_start_L));

AU_ds_L = nan(size(n_down_start_L));
AU_us_L = nan(size(n_down_start_L));
AU_ds_R = nan(size(n_down_start_L));
AU_us_R = nan(size(n_down_start_L));
dAU_ds = nan(size(n_down_start_L));
dAU_us = nan(size(n_down_start_L));

dU_mean_wb = nan(size(n_down_start_L));
dU_mean_ds = nan(size(n_down_start_L));
dU_mean_us = nan(size(n_down_start_L));
dU_mean_dus = nan(size(n_down_start_L));
dU_mean_udsPREV = nan(size(n_down_start_L));
dU_mean_udsNEXT = nan(size(n_down_start_L));

dU_max_dus = nan(size(n_down_start_L));
dU_max_udsPREV = nan(size(n_down_start_L));
dU_max_udsNEXT = nan(size(n_down_start_L));
dU_min_ds = nan(size(n_down_start_L));
dU_min_us = nan(size(n_down_start_L));

% aoa wing data
aoa_mean_wb_L = nan(size(n_down_start_L));
aoa_mean_ds_L = nan(size(n_down_start_L));
aoa_mean_us_L = nan(size(n_down_start_L));
aoa_mean_dus_L = nan(size(n_down_start_L));
aoa_mean_udsPREV_L = nan(size(n_down_start_L));
aoa_mean_udsNEXT_L = nan(size(n_down_start_L));

aoa_mean_wb_R = nan(size(n_down_start_L));
aoa_mean_ds_R = nan(size(n_down_start_L));
aoa_mean_us_R = nan(size(n_down_start_L));
aoa_mean_dus_R = nan(size(n_down_start_L));
aoa_mean_udsPREV_R = nan(size(n_down_start_L));
aoa_mean_udsNEXT_R = nan(size(n_down_start_L));

aoa_max_dus_L = nan(size(n_down_start_L));
aoa_max_udsPREV_L = nan(size(n_down_start_L));
aoa_max_udsNEXT_L = nan(size(n_down_start_L));
aoa_min_ds_L = nan(size(n_down_start_L));
aoa_min_us_L = nan(size(n_down_start_L));

aoa_max_dus_R = nan(size(n_down_start_L));
aoa_max_udsPREV_R = nan(size(n_down_start_L));
aoa_max_udsNEXT_R = nan(size(n_down_start_L));
aoa_min_ds_R = nan(size(n_down_start_L));
aoa_min_us_R = nan(size(n_down_start_L));

Aaoa_ds_L = nan(size(n_down_start_L));
Aaoa_us_L = nan(size(n_down_start_L));
Aaoa_ds_R = nan(size(n_down_start_L));
Aaoa_us_R = nan(size(n_down_start_L));
dAaoa_ds = nan(size(n_down_start_L));
dAaoa_us = nan(size(n_down_start_L));

daoa_mean_wb = nan(size(n_down_start_L));
daoa_mean_ds = nan(size(n_down_start_L));
daoa_mean_us = nan(size(n_down_start_L));
daoa_mean_dus = nan(size(n_down_start_L));
daoa_mean_udsPREV = nan(size(n_down_start_L));
daoa_mean_udsNEXT = nan(size(n_down_start_L));

daoa_max_dus = nan(size(n_down_start_L));
daoa_max_udsPREV = nan(size(n_down_start_L));
daoa_max_udsNEXT = nan(size(n_down_start_L));
daoa_min_ds = nan(size(n_down_start_L));
daoa_min_us = nan(size(n_down_start_L));

% body kin data
roll_mean_wb = nan(size(n_down_start_L));
roll_mean_ds = nan(size(n_down_start_L));
roll_mean_us = nan(size(n_down_start_L));
pitch_mean_wb = nan(size(n_down_start_L));
pitch_mean_ds = nan(size(n_down_start_L));
pitch_mean_us = nan(size(n_down_start_L));
yaw_mean_wb = nan(size(n_down_start_L));
yaw_mean_ds = nan(size(n_down_start_L));
yaw_mean_us = nan(size(n_down_start_L));

roll_dot_mean_wb = nan(size(n_down_start_L));
roll_dot_mean_ds = nan(size(n_down_start_L));
roll_dot_mean_us = nan(size(n_down_start_L));
pitch_dot_mean_wb = nan(size(n_down_start_L));
pitch_dot_mean_ds = nan(size(n_down_start_L));
pitch_dot_mean_us = nan(size(n_down_start_L));
yaw_dot_mean_wb = nan(size(n_down_start_L));
yaw_dot_mean_ds = nan(size(n_down_start_L));
yaw_dot_mean_us = nan(size(n_down_start_L));

roll_dot_dot_mean_wb = nan(size(n_down_start_L));
roll_dot_dot_mean_ds = nan(size(n_down_start_L));
roll_dot_dot_mean_us = nan(size(n_down_start_L));
pitch_dot_dot_mean_wb = nan(size(n_down_start_L));
pitch_dot_dot_mean_ds = nan(size(n_down_start_L));
pitch_dot_dot_mean_us = nan(size(n_down_start_L));
yaw_dot_dot_mean_wb = nan(size(n_down_start_L));
yaw_dot_dot_mean_ds = nan(size(n_down_start_L));
yaw_dot_dot_mean_us = nan(size(n_down_start_L));

F_mean_wb = nan(size(n_down_start_L));
F_mean_ds = nan(size(n_down_start_L));
F_mean_us = nan(size(n_down_start_L));
Fsp_roll_mean_wb = nan(size(n_down_start_L));
Fsp_roll_mean_ds = nan(size(n_down_start_L));
Fsp_roll_mean_us = nan(size(n_down_start_L));
Fsp_pitch_mean_wb = nan(size(n_down_start_L));
Fsp_pitch_mean_ds = nan(size(n_down_start_L));
Fsp_pitch_mean_us = nan(size(n_down_start_L));

t_wb_L = [];
t_wb_R = [];

t_ds_L = [];
t_ds_R = [];

t_us_L = [];
t_us_R = [];

stroke_wb_L = [];
stroke_wb_R = [];
pitch_wb_L = [];
pitch_wb_R = [];
dev_wb_L = [];
dev_wb_R = [];
aoa_wb_L = [];
aoa_wb_R = [];
U_wb_L = [];
U_wb_R = [];
Dstroke_wb = [];
Dpitch_wb = [];
Ddev_wb = [];
Daoa_wb = [];
DU_wb = [];

stroke_ds_L = [];
stroke_ds_R = [];
pitch_ds_L = [];
pitch_ds_R = [];
dev_ds_L = [];
dev_ds_R = [];
aoa_ds_L = [];
aoa_ds_R = [];
U_ds_L = [];
U_ds_R = [];
Dstroke_ds = [];
Dpitch_ds = [];
Ddev_ds = [];
Daoa_ds = [];
DU_ds = [];

stroke_us_L = [];
stroke_us_R = [];
pitch_us_L = [];
pitch_us_R = [];
dev_us_L = [];
dev_us_R = [];
aoa_us_L = [];
aoa_us_R = [];
U_us_L = [];
U_us_R = [];
Dstroke_us = [];
Dpitch_us = [];
Ddev_us = [];
Daoa_us = [];
DU_us = [];

t_wb_bins = [];
t_ds_bins = [];
t_us_bins = [];

stroke_wb_L_bins = [];
stroke_wb_R_bins = [];
pitch_wb_L_bins = [];
pitch_wb_R_bins = [];
dev_wb_L_bins = [];
dev_wb_R_bins = [];
aoa_wb_L_bins = [];
aoa_wb_R_bins = [];
U_wb_L_bins = [];
U_wb_R_bins = [];
Dstroke_wb_bins = [];
Dpitch_wb_bins = [];
Ddev_wb_bins = [];
Daoa_wb_bins = [];
DU_wb_bins = [];

stroke_ds_L_bins = [];
stroke_ds_R_bins = [];
pitch_ds_L_bins = [];
pitch_ds_R_bins = [];
dev_ds_L_bins = [];
dev_ds_R_bins = [];
aoa_ds_L_bins = [];
aoa_ds_R_bins = [];
U_ds_L_bins = [];
U_ds_R_bins = [];
Dstroke_ds_bins = [];
Dpitch_ds_bins = [];
Ddev_ds_bins = [];
Daoa_ds_bins = [];
DU_ds_bins = [];

stroke_us_L_bins = [];
stroke_us_R_bins = [];
pitch_us_L_bins = [];
pitch_us_R_bins = [];
dev_us_L_bins = [];
dev_us_R_bins = [];
aoa_us_L_bins = [];
aoa_us_R_bins = [];
U_us_L_bins = [];
U_us_R_bins = [];
Dstroke_us_bins = [];
Dpitch_us_bins = [];
Ddev_us_bins = [];
Daoa_us_bins = [];
DU_us_bins = [];

dt_ds = [];
dt_us = [];
dt_wb = [];
f_wb = [];
Rds = [];

n_now = 0;
% force & accel data
for seq = 1:size(n_down_start_L,2)
    for wb = 2:size(n_down_start_L,1)-1
        
%         counter = size(n_down_start_L,2)*100+size(n_down_start_L,1)-seq*100-wb
        
        n_down_start_L_now = n_down_start_L(wb,seq);
        n_down_stop_L_now = n_up_start_L(wb,seq);
        n_up_start_L_now = n_down_stop_L_now;
        n_up_stop_L_now = n_down_start_L(wb+1,seq);
        n_NEXTdown_start_L_now = n_up_stop_L_now;
        n_NEXTdown_stop_L_now = n_up_start_L(wb+1,seq);
        n_PREVup_start_L_now = n_up_start_L(wb-1,seq);
        n_PREVup_stop_L_now = n_down_start_L_now;
        
        n_down_start_R_now = n_down_start_R(wb,seq);
        n_down_stop_R_now = n_up_start_R(wb,seq);
        n_up_start_R_now = n_up_start_R(wb,seq);
        n_up_stop_R_now = n_down_start_R(wb+1,seq);
        n_NEXTdown_start_R_now = n_up_stop_R_now;
        n_NEXTdown_stop_R_now = n_up_start_R(wb+1,seq);
        n_PREVup_start_R_now = n_up_start_R(wb-1,seq);
        n_PREVup_stop_R_now = n_down_start_R_now;

        n_down_mid_L_now = n_down_start_L_now + round((n_down_stop_L_now - n_down_start_L_now)/2);
        n_down_mid_R_now = n_down_start_R_now + round((n_down_stop_R_now - n_down_start_R_now)/2);
        n_up_mid_L_now = n_up_start_L_now + round((n_up_stop_L_now - n_up_start_L_now)/2);
        n_up_mid_R_now = n_up_start_R_now + round((n_up_stop_R_now - n_up_start_R_now)/2);
        n_NEXTdown_mid_L_now = n_NEXTdown_start_L_now + round((n_NEXTdown_stop_L_now - n_NEXTdown_start_L_now)/2);
        n_NEXTdown_mid_R_now = n_NEXTdown_start_R_now + round((n_NEXTdown_stop_R_now - n_NEXTdown_start_R_now)/2);
        n_PREVup_mid_L_now = n_PREVup_start_L_now + round((n_PREVup_stop_L_now - n_PREVup_start_L_now)/2);
        n_PREVup_mid_R_now = n_PREVup_start_R_now + round((n_PREVup_stop_R_now - n_PREVup_start_R_now)/2);
                
        if isnan(n_down_start_L_now)==0 && isnan(n_up_stop_L_now)==0

            n_now = n_now + 1;
            
            stroke_wb_L_now = stroke_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            stroke_wb_R_now = stroke_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            stroke_ds_L_now = stroke_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            stroke_ds_R_now = stroke_R_mirror(n_down_start_R_now:n_down_stop_R_now,seq);
            stroke_us_L_now = stroke_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            stroke_us_R_now = stroke_R_mirror(n_up_start_R_now:n_up_stop_R_now,seq);
            
            stroke_dus_L_now = stroke_L_mirror(n_down_mid_L_now:n_up_mid_L_now,seq);
            stroke_uds_L_next = stroke_L_mirror(n_up_mid_L_now:n_NEXTdown_mid_L_now,seq);
            stroke_uds_L_prev = stroke_L_mirror(n_PREVup_mid_L_now:n_down_mid_L_now,seq);
            stroke_dus_R_now = stroke_R_mirror(n_down_mid_R_now:n_up_mid_R_now,seq);
            stroke_uds_R_next = stroke_R_mirror(n_up_mid_R_now:n_NEXTdown_mid_R_now,seq);
            stroke_uds_R_prev = stroke_R_mirror(n_PREVup_mid_R_now:n_down_mid_R_now,seq);

            pitch_wb_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq));
            pitch_wb_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq));
            pitch_ds_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq));
            pitch_ds_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_down_stop_R_now,seq));
            pitch_us_L_now = unwrap(pitch_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq));
            pitch_us_R_now = unwrap(pitch_R_mirror(n_up_start_R_now:n_up_stop_R_now,seq));
            
            pitch_dus_L_now = unwrap(pitch_L_mirror(n_down_mid_L_now:n_up_mid_L_now,seq));
            pitch_uds_L_next = unwrap(pitch_L_mirror(n_up_mid_L_now:n_NEXTdown_mid_L_now,seq));
            pitch_uds_L_prev = unwrap(pitch_L_mirror(n_PREVup_mid_L_now:n_down_mid_L_now,seq));
            pitch_dus_R_now = unwrap(pitch_R_mirror(n_down_mid_R_now:n_up_mid_R_now,seq));
            pitch_uds_R_next = unwrap(pitch_R_mirror(n_up_mid_R_now:n_NEXTdown_mid_R_now,seq));
            pitch_uds_R_prev = unwrap(pitch_R_mirror(n_PREVup_mid_R_now:n_down_mid_R_now,seq));

            dev_wb_L_now = dev_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            dev_wb_R_now = dev_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            dev_ds_L_now = dev_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            dev_ds_R_now = dev_R_mirror(n_down_start_R_now:n_down_stop_R_now,seq);
            dev_us_L_now = dev_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            dev_us_R_now = dev_R_mirror(n_up_start_R_now:n_up_stop_R_now,seq);
            
            dev_dus_L_now = dev_L_mirror(n_down_mid_L_now:n_up_mid_L_now,seq);
            dev_uds_L_next = dev_L_mirror(n_up_mid_L_now:n_NEXTdown_mid_L_now,seq);
            dev_uds_L_prev = dev_L_mirror(n_PREVup_mid_L_now:n_down_mid_L_now,seq);
            dev_dus_R_now = dev_R_mirror(n_down_mid_R_now:n_up_mid_R_now,seq);
            dev_uds_R_next = dev_R_mirror(n_up_mid_R_now:n_NEXTdown_mid_R_now,seq);
            dev_uds_R_prev = dev_R_mirror(n_PREVup_mid_R_now:n_down_mid_R_now,seq);
            
            U_wb_L_now = U_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            U_wb_R_now = U_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            U_ds_L_now = U_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            U_ds_R_now = U_R_mirror(n_down_start_R_now:n_down_stop_R_now,seq);
            U_us_L_now = U_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            U_us_R_now = U_R_mirror(n_up_start_R_now:n_up_stop_R_now,seq);
            
            U_dus_L_now = U_L_mirror(n_down_mid_L_now:n_up_mid_L_now,seq);
            U_uds_L_next = U_L_mirror(n_up_mid_L_now:n_NEXTdown_mid_L_now,seq);
            U_uds_L_prev = U_L_mirror(n_PREVup_mid_L_now:n_down_mid_L_now,seq);
            U_dus_R_now = U_R_mirror(n_down_mid_R_now:n_up_mid_R_now,seq);
            U_uds_R_next = U_R_mirror(n_up_mid_R_now:n_NEXTdown_mid_R_now,seq);
            U_uds_R_prev = U_R_mirror(n_PREVup_mid_R_now:n_down_mid_R_now,seq);
            
            aoa_wb_L_now = aoa_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            aoa_wb_R_now = aoa_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            aoa_ds_L_now = aoa_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            aoa_ds_R_now = aoa_R_mirror(n_down_start_R_now:n_down_stop_R_now,seq);
            aoa_us_L_now = aoa_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            aoa_us_R_now = aoa_R_mirror(n_up_start_R_now:n_up_stop_R_now,seq);
            
            aoa_dus_L_now = aoa_L_mirror(n_down_mid_L_now:n_up_mid_L_now,seq);
            aoa_uds_L_next = aoa_L_mirror(n_up_mid_L_now:n_NEXTdown_mid_L_now,seq);
            aoa_uds_L_prev = aoa_L_mirror(n_PREVup_mid_L_now:n_down_mid_L_now,seq);
            aoa_dus_R_now = aoa_R_mirror(n_down_mid_R_now:n_up_mid_R_now,seq);
            aoa_uds_R_next = aoa_R_mirror(n_up_mid_R_now:n_NEXTdown_mid_R_now,seq);
            aoa_uds_R_prev = aoa_R_mirror(n_PREVup_mid_R_now:n_down_mid_R_now,seq);
            
            % body kin
            V_wb_now = nanmean(V(n_down_start_L_now:n_up_stop_L_now,seq));
            pitch_global_wb_now = nanmean(pitch_global(n_down_start_L_now:n_up_stop_L_now,seq));
            
            roll_wb_now = roll_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            roll_ds_now = roll_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            roll_us_now = roll_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            pitch_wb_now = pitch(n_down_start_L_now:n_up_stop_L_now,seq);
            pitch_ds_now = pitch(n_down_start_L_now:n_down_stop_L_now,seq);
            pitch_us_now = pitch(n_up_start_L_now:n_up_stop_L_now,seq);
            yaw_wb_now = yaw_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            yaw_ds_now = yaw_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            yaw_us_now = yaw_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            roll_dot_wb_now = roll_dot_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            roll_dot_ds_now = roll_dot_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            roll_dot_us_now = roll_dot_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            pitch_dot_wb_now = pitch_dot(n_down_start_L_now:n_up_stop_L_now,seq);
            pitch_dot_ds_now = pitch_dot(n_down_start_L_now:n_down_stop_L_now,seq);
            pitch_dot_us_now = pitch_dot(n_up_start_L_now:n_up_stop_L_now,seq);
            yaw_dot_wb_now = yaw_dot_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            yaw_dot_ds_now = yaw_dot_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            yaw_dot_us_now = yaw_dot_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            roll_dot_dot_wb_now = roll_dot_dot_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            roll_dot_dot_ds_now = roll_dot_dot_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            roll_dot_dot_us_now = roll_dot_dot_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            pitch_dot_dot_wb_now = pitch_dot_dot(n_down_start_L_now:n_up_stop_L_now,seq);
            pitch_dot_dot_ds_now = pitch_dot_dot(n_down_start_L_now:n_down_stop_L_now,seq);
            pitch_dot_dot_us_now = pitch_dot_dot(n_up_start_L_now:n_up_stop_L_now,seq);
            yaw_dot_dot_wb_now = yaw_dot_dot_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            yaw_dot_dot_ds_now = yaw_dot_dot_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            yaw_dot_dot_us_now = yaw_dot_dot_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            F_wb_now = F(n_down_start_L_now:n_up_stop_L_now,seq);
            F_ds_now = F(n_down_start_L_now:n_down_stop_L_now,seq);
            F_us_now = F(n_up_start_L_now:n_up_stop_L_now,seq);
            
            Fsp_roll_wb_now = Fsp_roll_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            Fsp_roll_ds_now = Fsp_roll_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            Fsp_roll_us_now = Fsp_roll_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            Fsp_pitch_wb_now = Fsp_pitch(n_down_start_L_now:n_up_stop_L_now,seq);
            Fsp_pitch_ds_now = Fsp_pitch(n_down_start_L_now:n_down_stop_L_now,seq);
            Fsp_pitch_us_now = Fsp_pitch(n_up_start_L_now:n_up_stop_L_now,seq);
% 
%             % !!!!! WB means inc NEXT ds mid !!!!!!!!!!
%             roll_wb_now = roll_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             roll_ds_now = roll_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             roll_us_now = roll_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_wb_now = pitch(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_ds_now = pitch(n_down_start_L_now:n_up_mid_L_now,seq);
%             pitch_us_now = pitch(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_wb_now = yaw_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_ds_now = yaw_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             yaw_us_now = yaw_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             
%             roll_dot_wb_now = roll_dot_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             roll_dot_ds_now = roll_dot_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             roll_dot_us_now = roll_dot_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_dot_wb_now = pitch_dot(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_dot_ds_now = pitch_dot(n_down_start_L_now:n_up_mid_L_now,seq);
%             pitch_dot_us_now = pitch_dot(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_dot_wb_now = yaw_dot_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_dot_ds_now = yaw_dot_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             yaw_dot_us_now = yaw_dot_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             
%             roll_dot_dot_wb_now = roll_dot_dot_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             roll_dot_dot_ds_now = roll_dot_dot_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             roll_dot_dot_us_now = roll_dot_dot_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_dot_dot_wb_now = pitch_dot_dot(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             pitch_dot_dot_ds_now = pitch_dot_dot(n_down_start_L_now:n_up_mid_L_now,seq);
%             pitch_dot_dot_us_now = pitch_dot_dot(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_dot_dot_wb_now = yaw_dot_dot_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             yaw_dot_dot_ds_now = yaw_dot_dot_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             yaw_dot_dot_us_now = yaw_dot_dot_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             
%             F_wb_now = F(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             F_ds_now = F(n_down_start_L_now:n_up_mid_L_now,seq);
%             F_us_now = F(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             
%             Fsp_roll_wb_now = Fsp_roll_mirror(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             Fsp_roll_ds_now = Fsp_roll_mirror(n_down_start_L_now:n_up_mid_L_now,seq);
%             Fsp_roll_us_now = Fsp_roll_mirror(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);
%             
%             Fsp_pitch_wb_now = Fsp_pitch(n_down_start_L_now:n_NEXTdown_mid_L_now,seq);
%             Fsp_pitch_ds_now = Fsp_pitch(n_down_start_L_now:n_up_mid_L_now,seq);
%             Fsp_pitch_us_now = Fsp_pitch(n_up_start_L_now:n_NEXTdown_mid_L_now,seq);


            %% wingbeat distributions

            % time steps
            t_wb_L_now = [0:1/(n_up_stop_L_now-n_down_start_L_now):1]';
            t_wb_R_now = [0:1/(n_up_stop_R_now-n_down_start_R_now):1]';
            ni_wb_L = [1:(length(stroke_wb_L_now)-1)/(n-1):length(stroke_wb_L_now)]';
            ni_wb_R = [1:(length(stroke_wb_R_now)-1)/(n-1):length(stroke_wb_R_now)]';

            t_ds_L_now = [-1:1/(n_down_stop_L_now-n_down_start_L_now):0]';
            t_ds_R_now = [-1:1/(n_down_stop_R_now-n_down_start_R_now):0]';
            ni_ds_L = [1:(length(stroke_ds_L_now)-1)/(n-1):length(stroke_ds_L_now)]';
            ni_ds_R = [1:(length(stroke_ds_R_now)-1)/(n-1):length(stroke_ds_R_now)]';

            t_us_L_now = [0:1/(n_up_stop_L_now-n_down_stop_L_now):1]';
            t_us_R_now = [0:1/(n_up_stop_R_now-n_down_stop_R_now):1]';
            ni_us_L = [1:(length(stroke_us_L_now)-1)/(n-1)/2:length(stroke_us_L_now)]';
            ni_us_R = [1:(length(stroke_us_R_now)-1)/(n-1)/2:length(stroke_us_R_now)]';

            % interpolate
            stroke_wb_L_interp = rad2deg(interp1(stroke_wb_L_now,ni_wb_L));
            stroke_wb_R_interp = rad2deg(interp1(stroke_wb_R_now,ni_wb_R));
            pitch_wb_L_interp = rad2deg(interp1(pitch_wb_L_now,ni_wb_L));
            pitch_wb_R_interp = rad2deg(interp1(pitch_wb_R_now,ni_wb_R));
            dev_wb_L_interp = rad2deg(interp1(dev_wb_L_now,ni_wb_L));
            dev_wb_R_interp = rad2deg(interp1(dev_wb_R_now,ni_wb_R));
            aoa_wb_L_interp = rad2deg(interp1(aoa_wb_L_now,ni_wb_L));
            aoa_wb_R_interp = rad2deg(interp1(aoa_wb_R_now,ni_wb_R));
            U_wb_L_interp = interp1(U_wb_L_now,ni_wb_L)/1000;
            U_wb_R_interp = interp1(U_wb_R_now,ni_wb_R)/1000;
            
            stroke_ds_L_interp = rad2deg(interp1(stroke_ds_L_now,ni_ds_L));
            stroke_ds_R_interp = rad2deg(interp1(stroke_ds_R_now,ni_ds_R));
            pitch_ds_L_interp = rad2deg(interp1(pitch_ds_L_now,ni_ds_L));
            pitch_ds_R_interp = rad2deg(interp1(pitch_ds_R_now,ni_ds_R));
            dev_ds_L_interp = rad2deg(interp1(dev_ds_L_now,ni_ds_L));
            dev_ds_R_interp = rad2deg(interp1(dev_ds_R_now,ni_ds_R));
            aoa_ds_L_interp = rad2deg(interp1(aoa_ds_L_now,ni_ds_L));
            aoa_ds_R_interp = rad2deg(interp1(aoa_ds_R_now,ni_ds_R));
            U_ds_L_interp = interp1(U_ds_L_now,ni_ds_L)/1000;
            U_ds_R_interp = interp1(U_ds_R_now,ni_ds_R)/1000;
            
            stroke_us_L_interp = rad2deg(interp1(stroke_us_L_now,ni_us_L));
            stroke_us_R_interp = rad2deg(interp1(stroke_us_R_now,ni_us_R));
            pitch_us_L_interp = rad2deg(interp1(pitch_us_L_now,ni_us_L));
            pitch_us_R_interp = rad2deg(interp1(pitch_us_R_now,ni_us_R));
            dev_us_L_interp = rad2deg(interp1(dev_us_L_now,ni_us_L));
            dev_us_R_interp = rad2deg(interp1(dev_us_R_now,ni_us_R));
            aoa_us_L_interp = rad2deg(interp1(aoa_us_L_now,ni_us_L));
            aoa_us_R_interp = rad2deg(interp1(aoa_us_R_now,ni_us_R));
            U_us_L_interp = interp1(U_us_L_now,ni_us_L)/1000;
            U_us_R_interp = interp1(U_us_R_now,ni_us_R)/1000;
            
            Dstroke_wb_interp = stroke_wb_L_interp - stroke_wb_R_interp;
            Dpitch_wb_interp = pitch_wb_L_interp - pitch_wb_R_interp;
            Ddev_wb_interp = dev_wb_L_interp - dev_wb_R_interp;
            Daoa_wb_interp = aoa_wb_L_interp - aoa_wb_R_interp;
            DU_wb_interp = U_wb_L_interp - U_wb_R_interp;
            
            Dstroke_ds_interp = stroke_ds_L_interp - stroke_ds_R_interp;
            Dpitch_ds_interp = pitch_ds_L_interp - pitch_ds_R_interp;
            Ddev_ds_interp = dev_ds_L_interp - dev_ds_R_interp;
            Daoa_ds_interp = aoa_ds_L_interp - aoa_ds_R_interp;
            DU_ds_interp = U_ds_L_interp - U_ds_R_interp;
            
            Dstroke_us_interp = stroke_us_L_interp - stroke_us_R_interp;
            Dpitch_us_interp = pitch_us_L_interp - pitch_us_R_interp;
            Ddev_us_interp = dev_us_L_interp - dev_us_R_interp;
            Daoa_us_interp = aoa_us_L_interp - aoa_us_R_interp;
            DU_us_interp = U_us_L_interp - U_us_R_interp;
            
            % store data
            t_wb_L(end+1:end+length(t_wb_L_now),1) = t_wb_L_now;
            t_wb_R(end+1:end+length(t_wb_R_now),1) = t_wb_R_now;
            stroke_wb_L(end+1:end+length(t_wb_L_now),1) = rad2deg(stroke_wb_L_now);
            stroke_wb_R(end+1:end+length(t_wb_R_now),1) = rad2deg(stroke_wb_R_now);
            pitch_wb_L(end+1:end+length(t_wb_L_now),1) = rad2deg(pitch_wb_L_now);
            pitch_wb_R(end+1:end+length(t_wb_R_now),1) = rad2deg(pitch_wb_R_now);
            dev_wb_L(end+1:end+length(t_wb_L_now),1) = rad2deg(dev_wb_L_now);
            dev_wb_R(end+1:end+length(t_wb_R_now),1) = rad2deg(dev_wb_R_now);
            aoa_wb_L(end+1:end+length(t_wb_L_now),1) = rad2deg(aoa_wb_L_now);
            aoa_wb_R(end+1:end+length(t_wb_R_now),1) = rad2deg(aoa_wb_R_now);
            U_wb_L(end+1:end+length(t_wb_L_now),1) = U_wb_L_now/1000;
            U_wb_R(end+1:end+length(t_wb_R_now),1) = U_wb_R_now/1000;
            
            % store interp binned data separate rows
            t_wb_bins(:,n_now) = t_wb_bin;
            stroke_wb_L_bins(:,n_now) = stroke_wb_L_interp;
            stroke_wb_R_bins(:,n_now) = stroke_wb_R_interp;
            pitch_wb_L_bins(:,n_now) = pitch_wb_L_interp;
            pitch_wb_R_bins(:,n_now) = pitch_wb_R_interp;
            dev_wb_L_bins(:,n_now) = dev_wb_L_interp;
            dev_wb_R_bins(:,n_now) = dev_wb_R_interp;
            aoa_wb_L_bins(:,n_now) = aoa_wb_L_interp;
            aoa_wb_R_bins(:,n_now) = aoa_wb_R_interp;
            U_wb_L_bins(:,n_now) = U_wb_L_interp;
            U_wb_R_bins(:,n_now) = U_wb_R_interp;
            Dstroke_wb_bins(:,n_now) = Dstroke_wb_interp;
            Dpitch_wb_bins(:,n_now) = Dpitch_wb_interp;
            Ddev_wb_bins(:,n_now) = Ddev_wb_interp;
            Daoa_wb_bins(:,n_now) = Daoa_wb_interp;
            DU_wb_bins(:,n_now) = DU_wb_interp;
            
            t_ds_bins(:,n_now) = t_ds_bin;
            stroke_ds_L_bins(:,n_now) = stroke_ds_L_interp;
            stroke_ds_R_bins(:,n_now) = stroke_ds_R_interp;
            pitch_ds_L_bins(:,n_now) = pitch_ds_L_interp;
            pitch_ds_R_bins(:,n_now) = pitch_ds_R_interp;
            dev_ds_L_bins(:,n_now) = dev_ds_L_interp;
            dev_ds_R_bins(:,n_now) = dev_ds_R_interp;
            aoa_ds_L_bins(:,n_now) = aoa_ds_L_interp;
            aoa_ds_R_bins(:,n_now) = aoa_ds_R_interp;
            U_ds_L_bins(:,n_now) = U_ds_L_interp;
            U_ds_R_bins(:,n_now) = U_ds_R_interp;
            Dstroke_ds_bins(:,n_now) = Dstroke_ds_interp;
            Dpitch_ds_bins(:,n_now) = Dpitch_ds_interp;
            Ddev_ds_bins(:,n_now) = Ddev_ds_interp;
            Daoa_ds_bins(:,n_now) = Daoa_ds_interp;
            DU_ds_bins(:,n_now) = DU_ds_interp;
            
            t_us_bins(:,n_now) = t_us_bin;
            stroke_us_L_bins(:,n_now) = stroke_us_L_interp;
            stroke_us_R_bins(:,n_now) = stroke_us_R_interp;
            pitch_us_L_bins(:,n_now) = pitch_us_L_interp;
            pitch_us_R_bins(:,n_now) = pitch_us_R_interp;
            dev_us_L_bins(:,n_now) = dev_us_L_interp;
            dev_us_R_bins(:,n_now) = dev_us_R_interp;
            aoa_us_L_bins(:,n_now) = aoa_us_L_interp;
            aoa_us_R_bins(:,n_now) = aoa_us_R_interp;
            U_us_L_bins(:,n_now) = U_us_L_interp;
            U_us_R_bins(:,n_now) = U_us_R_interp;
            Dstroke_us_bins(:,n_now) = Dstroke_us_interp;
            Dpitch_us_bins(:,n_now) = Dpitch_us_interp;
            Ddev_us_bins(:,n_now) = Ddev_us_interp;
            Daoa_us_bins(:,n_now) = Daoa_us_interp;
            DU_us_bins(:,n_now) = DU_us_interp;
                        
            % freq data
            dt_wb_L(wb,seq) = (n_up_stop_L_now-n_down_start_L_now)/fps;
            dt_wb_R(wb,seq) = (n_up_stop_R_now-n_down_start_R_now)/fps;
            dt_ds_L(wb,seq) = (n_down_stop_L_now-n_down_start_L_now)/fps;
            dt_ds_R(wb,seq) = (n_down_stop_R_now-n_down_start_R_now)/fps;
            dt_us_L(wb,seq) = (n_up_stop_L_now-n_up_start_L_now)/fps;
            dt_us_R(wb,seq) = (n_up_stop_R_now-n_up_start_R_now)/fps;
            
            f_wb_L(wb,seq) = fps/(n_up_stop_L_now-n_down_start_L_now);
            f_wb_R(wb,seq) = fps/(n_up_stop_R_now-n_down_start_R_now);
            f_ds_L(wb,seq) = fps/(n_down_stop_L_now-n_down_start_L_now);
            f_ds_R(wb,seq) = fps/(n_down_stop_R_now-n_down_start_R_now);
            f_us_L(wb,seq) = fps/(n_up_stop_L_now-n_up_start_L_now);
            f_us_R(wb,seq) = fps/(n_up_stop_R_now-n_up_start_R_now);
            
            ddt_wb(wb,seq) = dt_wb_L(wb,seq) - dt_wb_R(wb,seq);
            ddt_ds(wb,seq) = dt_ds_L(wb,seq) - dt_ds_R(wb,seq);
            ddt_us(wb,seq) = dt_us_L(wb,seq) - dt_us_R(wb,seq);
            
            df_wb(wb,seq) = f_wb_L(wb,seq) - f_wb_R(wb,seq);
            df_ds(wb,seq) = f_ds_L(wb,seq) - f_ds_R(wb,seq);
            df_us(wb,seq) = f_us_L(wb,seq) - f_us_R(wb,seq);
            
            dt_wb_now = ([dt_wb_L(wb,seq),dt_wb_R(wb,seq)]);
            dt_ds_now = ([dt_ds_L(wb,seq),dt_ds_R(wb,seq)]);
            dt_us_now = ([dt_us_L(wb,seq),dt_us_R(wb,seq)]);
            f_wb_now = ([f_wb_L(wb,seq),f_wb_R(wb,seq)]);
            Rds_now = dt_ds_now ./ (dt_ds_now + dt_us_now);
            
            dt_wb(end+1:end+2,1) = dt_wb_now;
            dt_ds(end+1:end+2,1) = dt_ds_now;
            dt_us(end+1:end+2,1) = dt_us_now;
            f_wb(end+1:end+2,1) = f_wb_now;
            Rds(end+1:end+2,1) = Rds_now;
            
            %% MEAN VALUES

            % stroke data
            stroke_mean_wb_L(wb,seq) = circ_mean(stroke_wb_L_now);
            stroke_mean_ds_L(wb,seq) = circ_mean(stroke_ds_L_now);
            stroke_mean_us_L(wb,seq) = circ_mean(stroke_us_L_now);
            stroke_mean_wb_R(wb,seq) = circ_mean(stroke_wb_R_now);
            stroke_mean_ds_R(wb,seq) = circ_mean(stroke_ds_R_now);
            stroke_mean_us_R(wb,seq) = circ_mean(stroke_us_R_now);
            
            stroke_max_wb_L(wb,seq) = max(stroke_wb_L_now);
            stroke_min_wb_L(wb,seq) = min(stroke_wb_L_now);
            stroke_max_wb_R(wb,seq) = max(stroke_wb_R_now);
            stroke_min_wb_R(wb,seq) = min(stroke_wb_R_now);
            
            stroke_max_udsPREV_L(wb,seq) = max(stroke_uds_L_prev);
            stroke_max_udsNEXT_L(wb,seq) = max(stroke_uds_L_next);
            stroke_max_udsPREV_R(wb,seq) = max(stroke_uds_R_prev);
            stroke_max_udsNEXT_R(wb,seq) = max(stroke_uds_R_next);
            
            Astroke_wb_L(wb,seq) = stroke_max_wb_L(wb,seq) - stroke_min_wb_L(wb,seq);
            Astroke_wb_R(wb,seq) = stroke_max_wb_R(wb,seq) - stroke_min_wb_R(wb,seq);
            dAstroke_wb(wb,seq) = Astroke_wb_L(wb,seq) - Astroke_wb_R(wb,seq);

            Astroke_ds_L(wb,seq) = stroke_max_udsPREV_L(wb,seq) - stroke_min_wb_L(wb,seq);
            Astroke_ds_R(wb,seq) = stroke_max_udsPREV_R(wb,seq) - stroke_min_wb_R(wb,seq);
            dAstroke_ds(wb,seq) = Astroke_ds_L(wb,seq) - Astroke_ds_R(wb,seq);
            
            Astroke_us_L(wb,seq) = stroke_max_udsNEXT_L(wb,seq) - stroke_min_wb_L(wb,seq);
            Astroke_us_R(wb,seq) = stroke_max_udsNEXT_R(wb,seq) - stroke_min_wb_R(wb,seq);
            dAstroke_us(wb,seq) = Astroke_us_L(wb,seq) - Astroke_us_R(wb,seq);
            
            dstroke_mean_wb(wb,seq) = stroke_mean_wb_L(wb,seq) - stroke_mean_wb_R(wb,seq);
            dstroke_mean_ds(wb,seq) = stroke_mean_ds_L(wb,seq) - stroke_mean_ds_R(wb,seq);
            dstroke_mean_us(wb,seq) = stroke_mean_us_L(wb,seq) - stroke_mean_us_R(wb,seq);
            
            dstroke_max_wb(wb,seq) = stroke_max_wb_L(wb,seq) - stroke_max_wb_R(wb,seq);
            dstroke_min_wb(wb,seq) = stroke_min_wb_L(wb,seq) - stroke_min_wb_R(wb,seq);
            dstroke_max_udsPREV(wb,seq) = stroke_max_udsPREV_L(wb,seq) - stroke_max_udsPREV_R(wb,seq);
            dstroke_max_udsNEXT(wb,seq) = stroke_max_udsNEXT_L(wb,seq) - stroke_max_udsNEXT_R(wb,seq);
            
            % pitch data
            pitch_mean_wb_L(wb,seq) = circ_mean(pitch_wb_L_now);
            pitch_mean_wb_R(wb,seq) = circ_mean(pitch_wb_R_now);
            dpitch_mean_wb(wb,seq) = pitch_mean_wb_L(wb,seq) - pitch_mean_wb_R(wb,seq);

            pitch_max_wb_L(wb,seq) = max(pitch_wb_L_now);
            pitch_max_wb_R(wb,seq) = max(pitch_wb_R_now);
            dpitch_max_wb(wb,seq) = max(pitch_wb_L_now) - max(pitch_wb_R_now);
            
            pitch_min_wb_L(wb,seq) = min(pitch_wb_L_now);
            pitch_min_wb_R(wb,seq) = min(pitch_wb_R_now);
            dpitch_min_wb(wb,seq) = min(pitch_wb_L_now) - min(pitch_wb_R_now);
            
            pitch_mean_ds_L(wb,seq) = circ_mean(pitch_ds_L_now);
            pitch_mean_ds_R(wb,seq) = circ_mean(pitch_ds_R_now);
            dpitch_mean_ds(wb,seq) = pitch_mean_ds_L(wb,seq) - pitch_mean_ds_R(wb,seq);
            
            pitch_mean_us_L(wb,seq) = circ_mean(pitch_us_L_now);
            pitch_mean_us_R(wb,seq) = circ_mean(pitch_us_R_now);
            dpitch_mean_us(wb,seq) = pitch_mean_us_L(wb,seq) - pitch_mean_us_R(wb,seq);

            pitch_max_ds_L(wb,seq) = max(pitch_ds_L_now);
            pitch_max_ds_R(wb,seq) = max(pitch_ds_R_now);
            dpitch_max_ds(wb,seq) = max(pitch_ds_L_now) - max(pitch_ds_R_now);
            
            pitch_min_us_L(wb,seq) = min(pitch_us_L_now);
            pitch_min_us_R(wb,seq) = min(pitch_us_R_now);
            dpitch_min_us(wb,seq) = min(pitch_us_L_now) - min(pitch_us_R_now);
            
            Apitch_wb_L(wb,seq) = max(pitch_wb_L_now)-min(pitch_wb_L_now);
            Apitch_wb_R(wb,seq) = max(pitch_wb_R_now)-min(pitch_wb_R_now);
            dApitch_wb(wb,seq) = Apitch_wb_L(wb,seq) - Apitch_wb_R(wb,seq);
            
            Apitch_ds_L(wb,seq) = max(pitch_ds_L_now)-min(pitch_ds_L_now);
            Apitch_ds_R(wb,seq) = max(pitch_ds_R_now)-min(pitch_ds_R_now);
            dApitch_ds(wb,seq) = Apitch_ds_L(wb,seq) - Apitch_ds_R(wb,seq);
            
            Apitch_us_L(wb,seq) = max(pitch_us_L_now)-min(pitch_us_L_now);
            Apitch_us_R(wb,seq) = max(pitch_us_R_now)-min(pitch_us_R_now);
            dApitch_us(wb,seq) = Apitch_us_L(wb,seq) - Apitch_us_R(wb,seq);
            
            pitch_mid_ds_L(wb,seq) = pitch_ds_L_now(round(length(pitch_ds_L_now)/2));
            pitch_mid_ds_R(wb,seq) = pitch_ds_R_now(round(length(pitch_ds_R_now)/2));
            dpitch_mid_ds(wb,seq) = pitch_mid_ds_L(wb,seq) - pitch_mid_ds_R(wb,seq);

            pitch_mid_us_L(wb,seq) = pitch_us_L_now(round(length(pitch_us_L_now)/2));
            pitch_mid_us_R(wb,seq) = pitch_us_R_now(round(length(pitch_us_R_now)/2));
            dpitch_mid_us(wb,seq) = pitch_mid_us_L(wb,seq) - pitch_mid_us_R(wb,seq);

            Apitch_mid_L(wb,seq) = pitch_mid_ds_L(wb,seq) - pitch_mid_us_L(wb,seq);
            Apitch_mid_R(wb,seq) = pitch_mid_ds_R(wb,seq) - pitch_mid_us_R(wb,seq);
            dApitch_mid(wb,seq) = Apitch_mid_L(wb,seq) - Apitch_mid_R(wb,seq);
            
            % dev data
            dev_mean_wb_L(wb,seq) = circ_mean(dev_wb_L_now);
            dev_mean_ds_L(wb,seq) = circ_mean(dev_ds_L_now);
            dev_mean_us_L(wb,seq) = circ_mean(dev_us_L_now);
            dev_mean_dus_L(wb,seq) = circ_mean(dev_dus_L_now);
            dev_mean_udsPREV_L(wb,seq) = circ_mean(dev_uds_L_prev);
            dev_mean_udsNEXT_L(wb,seq) = circ_mean(dev_uds_L_next);
            
            dev_mean_wb_R(wb,seq) = circ_mean(dev_wb_R_now);
            dev_mean_ds_R(wb,seq) = circ_mean(dev_ds_R_now);
            dev_mean_us_R(wb,seq) = circ_mean(dev_us_R_now);
            dev_mean_dus_R(wb,seq) = circ_mean(dev_dus_R_now);
            dev_mean_udsPREV_R(wb,seq) = circ_mean(dev_uds_R_prev);
            dev_mean_udsNEXT_R(wb,seq) = circ_mean(dev_uds_R_next);
           
            dev_max_dus_L(wb,seq) = max(dev_dus_L_now);
            dev_max_udsPREV_L(wb,seq) = max(dev_uds_L_prev);
            dev_max_udsNEXT_L(wb,seq) = max(dev_uds_L_next);
            dev_min_ds_L(wb,seq) = min(dev_ds_L_now);
            dev_min_us_L(wb,seq) = min(dev_us_L_now);
            
            dev_max_dus_R(wb,seq) = max(dev_dus_R_now);
            dev_max_udsPREV_R(wb,seq) = max(dev_uds_R_prev);
            dev_max_udsNEXT_R(wb,seq) = max(dev_uds_R_next);
            dev_min_ds_R(wb,seq) = min(dev_ds_R_now);
            dev_min_us_R(wb,seq) = min(dev_us_R_now);
            
            Adev_ds_L(wb,seq) = dev_max_udsPREV_L(wb,seq) - dev_min_ds_L(wb,seq);
            Adev_us_L(wb,seq) = dev_max_dus_L(wb,seq) - dev_min_us_L(wb,seq);
            Adev_ds_R(wb,seq) = dev_max_udsPREV_R(wb,seq) - dev_min_ds_R(wb,seq);
            Adev_us_R(wb,seq) = dev_max_dus_R(wb,seq) - dev_min_us_R(wb,seq);
            dAdev_ds(wb,seq) = Adev_ds_L(wb,seq) - Adev_ds_R(wb,seq);
            dAdev_us(wb,seq) = Adev_us_L(wb,seq) - Adev_us_R(wb,seq);

            ddev_mean_wb(wb,seq) = dev_mean_wb_L(wb,seq) - dev_mean_wb_R(wb,seq);
            ddev_mean_ds(wb,seq) = dev_mean_ds_L(wb,seq) - dev_mean_ds_R(wb,seq);
            ddev_mean_us(wb,seq) = dev_mean_us_L(wb,seq) - dev_mean_us_R(wb,seq);
            ddev_mean_dus(wb,seq) = dev_mean_dus_L(wb,seq) - dev_mean_dus_R(wb,seq);
            ddev_mean_udsPREV(wb,seq) = dev_mean_udsPREV_L(wb,seq) - dev_mean_udsPREV_R(wb,seq);
            ddev_mean_udsNEXT(wb,seq) = dev_mean_udsNEXT_L(wb,seq) - dev_mean_udsNEXT_R(wb,seq);
            
            ddev_max_dus(wb,seq) = dev_max_dus_L(wb,seq) - dev_max_dus_R(wb,seq);
            ddev_max_udsPREV(wb,seq) = dev_max_udsPREV_L(wb,seq) - dev_max_udsPREV_R(wb,seq);
            ddev_max_udsNEXT(wb,seq) = dev_max_udsNEXT_L(wb,seq) - dev_max_udsNEXT_R(wb,seq);
            ddev_min_ds(wb,seq) = dev_min_ds_L(wb,seq) - dev_min_ds_R(wb,seq);
            ddev_min_us(wb,seq) = dev_min_us_L(wb,seq) - dev_min_us_R(wb,seq);
            
            
            % U data
            U_mean_wb_L(wb,seq) = circ_mean(U_wb_L_now);
            U_mean_ds_L(wb,seq) = circ_mean(U_ds_L_now);
            U_mean_us_L(wb,seq) = circ_mean(U_us_L_now);
            U_mean_dus_L(wb,seq) = circ_mean(U_dus_L_now);
            U_mean_udsPREV_L(wb,seq) = circ_mean(U_uds_L_prev);
            U_mean_udsNEXT_L(wb,seq) = circ_mean(U_uds_L_next);
            
            U_mean_wb_R(wb,seq) = circ_mean(U_wb_R_now);
            U_mean_ds_R(wb,seq) = circ_mean(U_ds_R_now);
            U_mean_us_R(wb,seq) = circ_mean(U_us_R_now);
            U_mean_dus_R(wb,seq) = circ_mean(U_dus_R_now);
            U_mean_udsPREV_R(wb,seq) = circ_mean(U_uds_R_prev);
            U_mean_udsNEXT_R(wb,seq) = circ_mean(U_uds_R_next);
           
            U_max_dus_L(wb,seq) = max(U_dus_L_now);
            U_max_udsPREV_L(wb,seq) = max(U_uds_L_prev);
            U_max_udsNEXT_L(wb,seq) = max(U_uds_L_next);
            U_min_ds_L(wb,seq) = min(U_ds_L_now);
            U_min_us_L(wb,seq) = min(U_us_L_now);
            
            U_max_dus_R(wb,seq) = max(U_dus_R_now);
            U_max_udsPREV_R(wb,seq) = max(U_uds_R_prev);
            U_max_udsNEXT_R(wb,seq) = max(U_uds_R_next);
            U_min_ds_R(wb,seq) = min(U_ds_R_now);
            U_min_us_R(wb,seq) = min(U_us_R_now);
            
            AU_ds_L(wb,seq) = U_max_udsPREV_L(wb,seq) - U_min_ds_L(wb,seq);
            AU_us_L(wb,seq) = U_max_dus_L(wb,seq) - U_min_us_L(wb,seq);
            AU_ds_R(wb,seq) = U_max_udsPREV_R(wb,seq) - U_min_ds_R(wb,seq);
            AU_us_R(wb,seq) = U_max_dus_R(wb,seq) - U_min_us_R(wb,seq);
            dAU_ds(wb,seq) = AU_ds_L(wb,seq) - AU_ds_R(wb,seq);
            dAU_us(wb,seq) = AU_us_L(wb,seq) - AU_us_R(wb,seq);

            dU_mean_wb(wb,seq) = U_mean_wb_L(wb,seq) - U_mean_wb_R(wb,seq);
            dU_mean_ds(wb,seq) = U_mean_ds_L(wb,seq) - U_mean_ds_R(wb,seq);
            dU_mean_us(wb,seq) = U_mean_us_L(wb,seq) - U_mean_us_R(wb,seq);
            dU_mean_dus(wb,seq) = U_mean_dus_L(wb,seq) - U_mean_dus_R(wb,seq);
            dU_mean_udsPREV(wb,seq) = U_mean_udsPREV_L(wb,seq) - U_mean_udsPREV_R(wb,seq);
            dU_mean_udsNEXT(wb,seq) = U_mean_udsNEXT_L(wb,seq) - U_mean_udsNEXT_R(wb,seq);
            
            dU_max_dus(wb,seq) = U_max_dus_L(wb,seq) - U_max_dus_R(wb,seq);
            dU_max_udsPREV(wb,seq) = U_max_udsPREV_L(wb,seq) - U_max_udsPREV_R(wb,seq);
            dU_max_udsNEXT(wb,seq) = U_max_udsNEXT_L(wb,seq) - U_max_udsNEXT_R(wb,seq);
            dU_min_ds(wb,seq) = U_min_ds_L(wb,seq) - U_min_ds_R(wb,seq);
            dU_min_us(wb,seq) = U_min_us_L(wb,seq) - U_min_us_R(wb,seq);
            
            
            % aoa data
            aoa_mean_wb_L(wb,seq) = circ_mean(aoa_wb_L_now);
            aoa_mean_ds_L(wb,seq) = circ_mean(aoa_ds_L_now);
            aoa_mean_us_L(wb,seq) = circ_mean(aoa_us_L_now);
            aoa_mean_dus_L(wb,seq) = circ_mean(aoa_dus_L_now);
            aoa_mean_udsPREV_L(wb,seq) = circ_mean(aoa_uds_L_prev);
            aoa_mean_udsNEXT_L(wb,seq) = circ_mean(aoa_uds_L_next);
            
            aoa_mean_wb_R(wb,seq) = circ_mean(aoa_wb_R_now);
            aoa_mean_ds_R(wb,seq) = circ_mean(aoa_ds_R_now);
            aoa_mean_us_R(wb,seq) = circ_mean(aoa_us_R_now);
            aoa_mean_dus_R(wb,seq) = circ_mean(aoa_dus_R_now);
            aoa_mean_udsPREV_R(wb,seq) = circ_mean(aoa_uds_R_prev);
            aoa_mean_udsNEXT_R(wb,seq) = circ_mean(aoa_uds_R_next);
           
            aoa_max_dus_L(wb,seq) = max(aoa_dus_L_now);
            aoa_max_udsPREV_L(wb,seq) = max(aoa_uds_L_prev);
            aoa_max_udsNEXT_L(wb,seq) = max(aoa_uds_L_next);
            aoa_min_ds_L(wb,seq) = min(aoa_ds_L_now);
            aoa_min_us_L(wb,seq) = min(aoa_us_L_now);
            
            aoa_max_dus_R(wb,seq) = max(aoa_dus_R_now);
            aoa_max_udsPREV_R(wb,seq) = max(aoa_uds_R_prev);
            aoa_max_udsNEXT_R(wb,seq) = max(aoa_uds_R_next);
            aoa_min_ds_R(wb,seq) = min(aoa_ds_R_now);
            aoa_min_us_R(wb,seq) = min(aoa_us_R_now);
            
            Aaoa_ds_L(wb,seq) = aoa_max_udsPREV_L(wb,seq) - aoa_min_ds_L(wb,seq);
            Aaoa_us_L(wb,seq) = aoa_max_dus_L(wb,seq) - aoa_min_us_L(wb,seq);
            Aaoa_ds_R(wb,seq) = aoa_max_udsPREV_R(wb,seq) - aoa_min_ds_R(wb,seq);
            Aaoa_us_R(wb,seq) = aoa_max_dus_R(wb,seq) - aoa_min_us_R(wb,seq);
            dAaoa_ds(wb,seq) = Aaoa_ds_L(wb,seq) - Aaoa_ds_R(wb,seq);
            dAaoa_us(wb,seq) = Aaoa_us_L(wb,seq) - Aaoa_us_R(wb,seq);

            daoa_mean_wb(wb,seq) = aoa_mean_wb_L(wb,seq) - aoa_mean_wb_R(wb,seq);
            daoa_mean_ds(wb,seq) = aoa_mean_ds_L(wb,seq) - aoa_mean_ds_R(wb,seq);
            daoa_mean_us(wb,seq) = aoa_mean_us_L(wb,seq) - aoa_mean_us_R(wb,seq);
            daoa_mean_dus(wb,seq) = aoa_mean_dus_L(wb,seq) - aoa_mean_dus_R(wb,seq);
            daoa_mean_udsPREV(wb,seq) = aoa_mean_udsPREV_L(wb,seq) - aoa_mean_udsPREV_R(wb,seq);
            daoa_mean_udsNEXT(wb,seq) = aoa_mean_udsNEXT_L(wb,seq) - aoa_mean_udsNEXT_R(wb,seq);
            
            daoa_max_dus(wb,seq) = aoa_max_dus_L(wb,seq) - aoa_max_dus_R(wb,seq);
            daoa_max_udsPREV(wb,seq) = aoa_max_udsPREV_L(wb,seq) - aoa_max_udsPREV_R(wb,seq);
            daoa_max_udsNEXT(wb,seq) = aoa_max_udsNEXT_L(wb,seq) - aoa_max_udsNEXT_R(wb,seq);
            daoa_min_ds(wb,seq) = aoa_min_ds_L(wb,seq) - aoa_min_ds_R(wb,seq);
            daoa_min_us(wb,seq) = aoa_min_us_L(wb,seq) - aoa_min_us_R(wb,seq);
            
            % body kin data
            V_mean_wb(wb,seq) = nanmean(V_wb_now);
            pitch_global_mean_wb(wb,seq) = nanmean(pitch_global_wb_now);
                        
            roll_mean_wb(wb,seq) = nanmean(roll_wb_now);
            roll_mean_ds(wb,seq) = nanmean(roll_ds_now);
            roll_mean_us(wb,seq) = nanmean(roll_us_now);
            pitch_mean_wb(wb,seq) = nanmean(pitch_wb_now);
            pitch_mean_ds(wb,seq) = nanmean(pitch_ds_now);
            pitch_mean_us(wb,seq) = nanmean(pitch_us_now);
            yaw_mean_wb(wb,seq) = nanmean(yaw_wb_now);
            yaw_mean_ds(wb,seq) = nanmean(yaw_ds_now);
            yaw_mean_us(wb,seq) = nanmean(yaw_us_now);
            
            roll_dot_mean_wb(wb,seq) = nanmean(roll_dot_wb_now);
            roll_dot_mean_ds(wb,seq) = nanmean(roll_dot_ds_now);
            roll_dot_mean_us(wb,seq) = nanmean(roll_dot_us_now);
            pitch_dot_mean_wb(wb,seq) = nanmean(pitch_dot_wb_now);
            pitch_dot_mean_ds(wb,seq) = nanmean(pitch_dot_ds_now);
            pitch_dot_mean_us(wb,seq) = nanmean(pitch_dot_us_now);
            yaw_dot_mean_wb(wb,seq) = nanmean(yaw_dot_wb_now);
            yaw_dot_mean_ds(wb,seq) = nanmean(yaw_dot_ds_now);
            yaw_dot_mean_us(wb,seq) = nanmean(yaw_dot_us_now);
            
            roll_dot_dot_mean_wb(wb,seq) = nanmean(roll_dot_dot_wb_now);
            roll_dot_dot_mean_ds(wb,seq) = nanmean(roll_dot_dot_ds_now);
            roll_dot_dot_mean_us(wb,seq) = nanmean(roll_dot_dot_us_now);
            pitch_dot_dot_mean_wb(wb,seq) = nanmean(pitch_dot_dot_wb_now);
            pitch_dot_dot_mean_ds(wb,seq) = nanmean(pitch_dot_dot_ds_now);
            pitch_dot_dot_mean_us(wb,seq) = nanmean(pitch_dot_dot_us_now);
            yaw_dot_dot_mean_wb(wb,seq) = nanmean(yaw_dot_dot_wb_now);
            yaw_dot_dot_mean_ds(wb,seq) = nanmean(yaw_dot_dot_ds_now);
            yaw_dot_dot_mean_us(wb,seq) = nanmean(yaw_dot_dot_us_now);
            
            F_mean_wb(wb,seq) = nanmean(F_wb_now);
            F_mean_ds(wb,seq) = nanmean(F_ds_now);
            F_mean_us(wb,seq) = nanmean(F_us_now);
            Fsp_roll_mean_wb(wb,seq) = circ_mean(Fsp_roll_wb_now);
            Fsp_roll_mean_ds(wb,seq) = circ_mean(Fsp_roll_ds_now);
            Fsp_roll_mean_us(wb,seq) = circ_mean(Fsp_roll_us_now);
            Fsp_pitch_mean_wb(wb,seq) = circ_mean(Fsp_pitch_wb_now);
            Fsp_pitch_mean_ds(wb,seq) = circ_mean(Fsp_pitch_ds_now);
            Fsp_pitch_mean_us(wb,seq) = circ_mean(Fsp_pitch_us_now);
        end
    end
end

