t_wb_bin = [0:1/(n-1):1]';
t_ds_bin = [-1:1/(n-1):0]';
t_us_bin = [0:1/(n-1):1]';

% seq data
seq_nr = [];
wb_nr = [];

% body kin data
V_mean_wb = [];
pitch_global_mean_wb = [];
l_wing_mean_wb = [];

% freq data
f_wb_L = [];
f_wb_R = [];
f_ds_L = [];
f_ds_R = [];
f_us_L = [];
f_us_R = [];

dt_wb_L = [];
dt_wb_R = [];
dt_ds_L = [];
dt_ds_R = [];
dt_us_L = [];
dt_us_R = [];

df_wb = [];
df_ds = [];
df_us = [];

ddt_wb = [];
ddt_ds = [];
ddt_us = [];

% stroke data
stroke_mean_wb_L = [];
stroke_mean_ds_L = [];
stroke_mean_us_L = [];
stroke_mean_wb_R = [];
stroke_mean_ds_R = [];
stroke_mean_us_R = [];

stroke_max_wb_L = [];
stroke_min_wb_L = [];
stroke_max_wb_R = [];
stroke_min_wb_R = [];

stroke_max_udsPREV_L = [];
stroke_max_udsNEXT_L = [];
stroke_max_udsPREV_R = [];
stroke_max_udsNEXT_R = [];

Astroke_wb_L = [];
Astroke_wb_R = [];
dAstroke_wb = [];

Astroke_ds_L = [];
Astroke_ds_R = [];
dAstroke_ds = [];

Astroke_us_L = [];
Astroke_us_R = [];
dAstroke_us = [];

dstroke_mean_wb = [];
dstroke_mean_ds = [];
dstroke_mean_us = [];

dstroke_max_wb = [];
dstroke_min_wb = [];
dstroke_max_udsPREV = [];
dstroke_max_udsNEXT = [];

% pitch data
pitch_mean_wb_L = [];
pitch_mean_wb_R = [];
dpitch_mean_wb = [];

pitch_max_wb_L = [];
pitch_max_wb_R = [];
dpitch_max_wb = [];

pitch_min_wb_L = [];
pitch_min_wb_R = [];
dpitch_min_wb = [];

pitch_mean_ds_L = [];
pitch_mean_ds_R = [];
dpitch_mean_ds = [];

pitch_mean_us_L = [];
pitch_mean_us_R = [];
dpitch_mean_us = [];

pitch_max_ds_L = [];
pitch_max_ds_R = [];
dpitch_max_ds = [];

pitch_min_us_L = [];
pitch_min_us_R = [];
dpitch_min_us = [];

Apitch_wb_L = [];
Apitch_wb_R = [];
dApitch_wb = [];

Apitch_ds_L = [];
Apitch_ds_R = [];
dApitch_ds = [];

Apitch_us_L = [];
Apitch_us_R = [];
dApitch_us = [];

pitch_mid_ds_L = [];
pitch_mid_ds_R = [];

pitch_mid_us_L = [];
pitch_mid_us_R = [];

Apitch_mid_L = [];
Apitch_mid_R = [];

dpitch_mid_ds = [];
dpitch_mid_us = [];
dApitch_mid = [];

% dev data
dev_mean_wb_L = [];
dev_mean_ds_L = [];
dev_mean_us_L = [];
dev_mean_dus_L = [];
dev_mean_udsPREV_L = [];
dev_mean_udsNEXT_L = [];

dev_mean_wb_R = [];
dev_mean_ds_R = [];
dev_mean_us_R = [];
dev_mean_dus_R = [];
dev_mean_udsPREV_R = [];
dev_mean_udsNEXT_R = [];

dev_max_dus_L = [];
dev_max_udsPREV_L = [];
dev_max_udsNEXT_L = [];
dev_min_ds_L = [];
dev_min_us_L = [];

dev_max_dus_R = [];
dev_max_udsPREV_R = [];
dev_max_udsNEXT_R = [];
dev_min_ds_R = [];
dev_min_us_R = [];

Adev_ds_L = [];
Adev_us_L = [];
Adev_ds_R = [];
Adev_us_R = [];
dAdev_ds = [];
dAdev_us = [];

ddev_mean_wb = [];
ddev_mean_ds = [];
ddev_mean_us = [];
ddev_mean_dus = [];
ddev_mean_udsPREV = [];
ddev_mean_udsNEXT = [];

ddev_max_dus = [];
ddev_max_udsPREV = [];
ddev_max_udsNEXT = [];
ddev_min_ds = [];
ddev_min_us = [];

% Uwing data
U_mean_wb_L = [];
U_mean_ds_L = [];
U_mean_us_L = [];
U_mean_dus_L = [];
U_mean_udsPREV_L = [];
U_mean_udsNEXT_L = [];

U_mean_wb_R = [];
U_mean_ds_R = [];
U_mean_us_R = [];
U_mean_dus_R = [];
U_mean_udsPREV_R = [];
U_mean_udsNEXT_R = [];

U_max_dus_L = [];
U_max_udsPREV_L = [];
U_max_udsNEXT_L = [];
U_min_ds_L = [];
U_min_us_L = [];

U_max_dus_R = [];
U_max_udsPREV_R = [];
U_max_udsNEXT_R = [];
U_min_ds_R = [];
U_min_us_R = [];

AU_ds_L = [];
AU_us_L = [];
AU_ds_R = [];
AU_us_R = [];
dAU_ds = [];
dAU_us = [];

dU_mean_wb = [];
dU_mean_ds = [];
dU_mean_us = [];
dU_mean_dus = [];
dU_mean_udsPREV = [];
dU_mean_udsNEXT = [];

dU_max_dus = [];
dU_max_udsPREV = [];
dU_max_udsNEXT = [];
dU_min_ds = [];
dU_min_us = [];

% aoa wing data
aoa_mean_wb_L = [];
aoa_mean_ds_L = [];
aoa_mean_us_L = [];
aoa_mean_dus_L = [];
aoa_mean_udsPREV_L = [];
aoa_mean_udsNEXT_L = [];

aoa_mean_wb_R = [];
aoa_mean_ds_R = [];
aoa_mean_us_R = [];
aoa_mean_dus_R = [];
aoa_mean_udsPREV_R = [];
aoa_mean_udsNEXT_R = [];

aoa_max_dus_L = [];
aoa_max_udsPREV_L = [];
aoa_max_udsNEXT_L = [];
aoa_min_ds_L = [];
aoa_min_us_L = [];

aoa_max_dus_R = [];
aoa_max_udsPREV_R = [];
aoa_max_udsNEXT_R = [];
aoa_min_ds_R = [];
aoa_min_us_R = [];

Aaoa_ds_L = [];
Aaoa_us_L = [];
Aaoa_ds_R = [];
Aaoa_us_R = [];
dAaoa_ds = [];
dAaoa_us = [];

daoa_mean_wb = [];
daoa_mean_ds = [];
daoa_mean_us = [];
daoa_mean_dus = [];
daoa_mean_udsPREV = [];
daoa_mean_udsNEXT = [];

daoa_max_dus = [];
daoa_max_udsPREV = [];
daoa_max_udsNEXT = [];
daoa_min_ds = [];
daoa_min_us = [];

% body kin data
roll_mean_wb = [];
roll_mean_ds = [];
roll_mean_us = [];
pitch_mean_wb = [];
pitch_mean_ds = [];
pitch_mean_us = [];
yaw_mean_wb = [];
yaw_mean_ds = [];
yaw_mean_us = [];

roll_dot_mean_wb = [];
roll_dot_mean_ds = [];
roll_dot_mean_us = [];
pitch_dot_mean_wb = [];
pitch_dot_mean_ds = [];
pitch_dot_mean_us = [];
yaw_dot_mean_wb = [];
yaw_dot_mean_ds = [];
yaw_dot_mean_us = [];

roll_dot_dot_mean_wb = [];
roll_dot_dot_mean_ds = [];
roll_dot_dot_mean_us = [];
pitch_dot_dot_mean_wb = [];
pitch_dot_dot_mean_ds = [];
pitch_dot_dot_mean_us = [];
yaw_dot_dot_mean_wb = [];
yaw_dot_dot_mean_ds = [];
yaw_dot_dot_mean_us = [];

drot_L_mean_wb = [];
drot_L_mean_ds = [];
drot_L_mean_us = [];

drot_R_mean_wb = [];
drot_R_mean_ds = [];
drot_R_mean_us = [];

rot_dot_L_mean_wb = [];
rot_dot_L_mean_ds = [];
rot_dot_L_mean_us = [];

rot_dot_R_mean_wb = [];
rot_dot_R_mean_ds = [];
rot_dot_R_mean_us = [];

rot_dot_dot_L_mean_wb = [];
rot_dot_dot_L_mean_ds = [];
rot_dot_dot_L_mean_us = [];

rot_dot_dot_R_mean_wb = [];
rot_dot_dot_R_mean_ds = [];
rot_dot_dot_R_mean_us = [];

F_mean_wb = [];
F_mean_ds = [];
F_mean_us = [];
Fsp_roll_mean_wb = [];
Fsp_roll_mean_ds = [];
Fsp_roll_mean_us = [];
Fsp_pitch_mean_wb = [];
Fsp_pitch_mean_ds = [];
Fsp_pitch_mean_us = [];

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
            l_wing_wb_now = wing_length(seq);
            
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
            
            drot_L_wb_now = drot_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            drot_L_ds_now = drot_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            drot_L_us_now = drot_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            drot_R_wb_now = drot_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            drot_R_ds_now = drot_R_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            drot_R_us_now = drot_R_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            rot_dot_L_wb_now = rot_dot_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            rot_dot_L_ds_now = rot_dot_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            rot_dot_L_us_now = rot_dot_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            rot_dot_R_wb_now = rot_dot_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            rot_dot_R_ds_now = rot_dot_R_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            rot_dot_R_us_now = rot_dot_R_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            rot_dot_dot_L_wb_now = rot_dot_dot_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            rot_dot_dot_L_ds_now = rot_dot_dot_L_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            rot_dot_dot_L_us_now = rot_dot_dot_L_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
            rot_dot_dot_R_wb_now = rot_dot_dot_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            rot_dot_dot_R_ds_now = rot_dot_dot_R_mirror(n_down_start_L_now:n_down_stop_L_now,seq);
            rot_dot_dot_R_us_now = rot_dot_dot_R_mirror(n_up_start_L_now:n_up_stop_L_now,seq);
            
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
            ni_us_L = [1:(length(stroke_us_L_now)-1)/(n-1):length(stroke_us_L_now)]';
            ni_us_R = [1:(length(stroke_us_R_now)-1)/(n-1):length(stroke_us_R_now)]';

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
            t_wb_L(1:fr_max,n_now) = nan;
            t_wb_R(1:fr_max,n_now) = nan;
            stroke_wb_L(1:fr_max,n_now) = nan;
            stroke_wb_R(1:fr_max,n_now) = nan;
            pitch_wb_L(1:fr_max,n_now) = nan;
            pitch_wb_R(1:fr_max,n_now) = nan;
            dev_wb_L(1:fr_max,n_now) = nan;
            dev_wb_R(1:fr_max,n_now) = nan;
            aoa_wb_L(1:fr_max,n_now) = nan;
            aoa_wb_R(1:fr_max,n_now) = nan;
            U_wb_L(1:fr_max,n_now) = nan;
            U_wb_R(1:fr_max,n_now) = nan;
            
            t_ds_L(1:fr_max,n_now) = nan;
            t_ds_R(1:fr_max,n_now) = nan;
            stroke_ds_L(1:fr_max,n_now) = nan;
            stroke_ds_R(1:fr_max,n_now) = nan;
            pitch_ds_L(1:fr_max,n_now) = nan;
            pitch_ds_R(1:fr_max,n_now) = nan;
            dev_ds_L(1:fr_max,n_now) = nan;
            dev_ds_R(1:fr_max,n_now) = nan;
            aoa_ds_L(1:fr_max,n_now) = nan;
            aoa_ds_R(1:fr_max,n_now) = nan;
            U_ds_L(1:fr_max,n_now) = nan;
            U_ds_R(1:fr_max,n_now) = nan;
            
            t_us_L(1:fr_max,n_now) = nan;
            t_us_R(1:fr_max,n_now) = nan;
            stroke_us_L(1:fr_max,n_now) = nan;
            stroke_us_R(1:fr_max,n_now) = nan;
            pitch_us_L(1:fr_max,n_now) = nan;
            pitch_us_R(1:fr_max,n_now) = nan;
            dev_us_L(1:fr_max,n_now) = nan;
            dev_us_R(1:fr_max,n_now) = nan;
            aoa_us_L(1:fr_max,n_now) = nan;
            aoa_us_R(1:fr_max,n_now) = nan;
            U_us_L(1:fr_max,n_now) = nan;
            U_us_R(1:fr_max,n_now) = nan;
            
            t_wb_L(1:length(t_wb_L_now),n_now) = t_wb_L_now;
            t_wb_R(1:length(t_wb_R_now),n_now) = t_wb_R_now;
            stroke_wb_L(1:length(t_wb_L_now),n_now) = rad2deg(stroke_wb_L_now);
            stroke_wb_R(1:length(t_wb_R_now),n_now) = rad2deg(stroke_wb_R_now);
            pitch_wb_L(1:length(t_wb_L_now),n_now) = rad2deg(pitch_wb_L_now);
            pitch_wb_R(1:length(t_wb_R_now),n_now) = rad2deg(pitch_wb_R_now);
            dev_wb_L(1:length(t_wb_L_now),n_now) = rad2deg(dev_wb_L_now);
            dev_wb_R(1:length(t_wb_R_now),n_now) = rad2deg(dev_wb_R_now);
            aoa_wb_L(1:length(t_wb_L_now),n_now) = rad2deg(aoa_wb_L_now);
            aoa_wb_R(1:length(t_wb_R_now),n_now) = rad2deg(aoa_wb_R_now);
            U_wb_L(1:length(t_wb_L_now),n_now) = U_wb_L_now/1000;
            U_wb_R(1:length(t_wb_R_now),n_now) = U_wb_R_now/1000;
            
            t_ds_L(1:length(t_ds_L_now),n_now) = t_ds_L_now;
            t_ds_R(1:length(t_ds_R_now),n_now) = t_ds_R_now;
            stroke_ds_L(1:length(t_ds_L_now),n_now) = rad2deg(stroke_ds_L_now);
            stroke_ds_R(1:length(t_ds_R_now),n_now) = rad2deg(stroke_ds_R_now);
            pitch_ds_L(1:length(t_ds_L_now),n_now) = rad2deg(pitch_ds_L_now);
            pitch_ds_R(1:length(t_ds_R_now),n_now) = rad2deg(pitch_ds_R_now);
            dev_ds_L(1:length(t_ds_L_now),n_now) = rad2deg(dev_ds_L_now);
            dev_ds_R(1:length(t_ds_R_now),n_now) = rad2deg(dev_ds_R_now);
            aoa_ds_L(1:length(t_ds_L_now),n_now) = rad2deg(aoa_ds_L_now);
            aoa_ds_R(1:length(t_ds_R_now),n_now) = rad2deg(aoa_ds_R_now);
            U_ds_L(1:length(t_ds_L_now),n_now) = U_ds_L_now/1000;
            U_ds_R(1:length(t_ds_R_now),n_now) = U_ds_R_now/1000;
            
            t_us_L(1:length(t_us_L_now),n_now) = t_us_L_now;
            t_us_R(1:length(t_us_R_now),n_now) = t_us_R_now;
            stroke_us_L(1:length(t_us_L_now),n_now) = rad2deg(stroke_us_L_now);
            stroke_us_R(1:length(t_us_R_now),n_now) = rad2deg(stroke_us_R_now);
            pitch_us_L(1:length(t_us_L_now),n_now) = rad2deg(pitch_us_L_now);
            pitch_us_R(1:length(t_us_R_now),n_now) = rad2deg(pitch_us_R_now);
            dev_us_L(1:length(t_us_L_now),n_now) = rad2deg(dev_us_L_now);
            dev_us_R(1:length(t_us_R_now),n_now) = rad2deg(dev_us_R_now);
            aoa_us_L(1:length(t_us_L_now),n_now) = rad2deg(aoa_us_L_now);
            aoa_us_R(1:length(t_us_R_now),n_now) = rad2deg(aoa_us_R_now);
            U_us_L(1:length(t_us_L_now),n_now) = U_us_L_now/1000;
            U_us_R(1:length(t_us_R_now),n_now) = U_us_R_now/1000;
            
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
            dt_wb_L(n_now,1) = (n_up_stop_L_now-n_down_start_L_now)/fps;
            dt_wb_R(n_now,1) = (n_up_stop_R_now-n_down_start_R_now)/fps;
            dt_ds_L(n_now,1) = (n_down_stop_L_now-n_down_start_L_now)/fps;
            dt_ds_R(n_now,1) = (n_down_stop_R_now-n_down_start_R_now)/fps;
            dt_us_L(n_now,1) = (n_up_stop_L_now-n_up_start_L_now)/fps;
            dt_us_R(n_now,1) = (n_up_stop_R_now-n_up_start_R_now)/fps;
            
            f_wb_L(n_now,1) = fps/(n_up_stop_L_now-n_down_start_L_now);
            f_wb_R(n_now,1) = fps/(n_up_stop_R_now-n_down_start_R_now);
            f_ds_L(n_now,1) = fps/(n_down_stop_L_now-n_down_start_L_now);
            f_ds_R(n_now,1) = fps/(n_down_stop_R_now-n_down_start_R_now);
            f_us_L(n_now,1) = fps/(n_up_stop_L_now-n_up_start_L_now);
            f_us_R(n_now,1) = fps/(n_up_stop_R_now-n_up_start_R_now);
            
            ddt_wb(n_now,1) = dt_wb_L(n_now,1) - dt_wb_R(n_now,1);
            ddt_ds(n_now,1) = dt_ds_L(n_now,1) - dt_ds_R(n_now,1);
            ddt_us(n_now,1) = dt_us_L(n_now,1) - dt_us_R(n_now,1);
            
            df_wb(n_now,1) = f_wb_L(n_now,1) - f_wb_R(n_now,1);
            df_ds(n_now,1) = f_ds_L(n_now,1) - f_ds_R(n_now,1);
            df_us(n_now,1) = f_us_L(n_now,1) - f_us_R(n_now,1);
            
            seq_nr(n_now,1) = seq;
            wb_nr(n_now,1) = wb;
            
            %% MEAN VALUES

            % stroke data
            stroke_mean_wb_L(n_now,1) = circ_mean(stroke_wb_L_now);
            stroke_mean_ds_L(n_now,1) = circ_mean(stroke_ds_L_now);
            stroke_mean_us_L(n_now,1) = circ_mean(stroke_us_L_now);
            stroke_mean_wb_R(n_now,1) = circ_mean(stroke_wb_R_now);
            stroke_mean_ds_R(n_now,1) = circ_mean(stroke_ds_R_now);
            stroke_mean_us_R(n_now,1) = circ_mean(stroke_us_R_now);
            
            stroke_max_wb_L(n_now,1) = max(stroke_wb_L_now);
            stroke_min_wb_L(n_now,1) = min(stroke_wb_L_now);
            stroke_max_wb_R(n_now,1) = max(stroke_wb_R_now);
            stroke_min_wb_R(n_now,1) = min(stroke_wb_R_now);
            
            stroke_max_udsPREV_L(n_now,1) = max(stroke_uds_L_prev);
            stroke_max_udsNEXT_L(n_now,1) = max(stroke_uds_L_next);
            stroke_max_udsPREV_R(n_now,1) = max(stroke_uds_R_prev);
            stroke_max_udsNEXT_R(n_now,1) = max(stroke_uds_R_next);
            
            Astroke_wb_L(n_now,1) = stroke_max_wb_L(n_now,1) - stroke_min_wb_L(n_now,1);
            Astroke_wb_R(n_now,1) = stroke_max_wb_R(n_now,1) - stroke_min_wb_R(n_now,1);
            dAstroke_wb(n_now,1) = Astroke_wb_L(n_now,1) - Astroke_wb_R(n_now,1);

            Astroke_ds_L(n_now,1) = stroke_max_udsPREV_L(n_now,1) - stroke_min_wb_L(n_now,1);
            Astroke_ds_R(n_now,1) = stroke_max_udsPREV_R(n_now,1) - stroke_min_wb_R(n_now,1);
            dAstroke_ds(n_now,1) = Astroke_ds_L(n_now,1) - Astroke_ds_R(n_now,1);
            
            Astroke_us_L(n_now,1) = stroke_max_udsNEXT_L(n_now,1) - stroke_min_wb_L(n_now,1);
            Astroke_us_R(n_now,1) = stroke_max_udsNEXT_R(n_now,1) - stroke_min_wb_R(n_now,1);
            dAstroke_us(n_now,1) = Astroke_us_L(n_now,1) - Astroke_us_R(n_now,1);
            
            dstroke_mean_wb(n_now,1) = stroke_mean_wb_L(n_now,1) - stroke_mean_wb_R(n_now,1);
            dstroke_mean_ds(n_now,1) = stroke_mean_ds_L(n_now,1) - stroke_mean_ds_R(n_now,1);
            dstroke_mean_us(n_now,1) = stroke_mean_us_L(n_now,1) - stroke_mean_us_R(n_now,1);
            
            dstroke_max_wb(n_now,1) = stroke_max_wb_L(n_now,1) - stroke_max_wb_R(n_now,1);
            dstroke_min_wb(n_now,1) = stroke_min_wb_L(n_now,1) - stroke_min_wb_R(n_now,1);
            dstroke_max_udsPREV(n_now,1) = stroke_max_udsPREV_L(n_now,1) - stroke_max_udsPREV_R(n_now,1);
            dstroke_max_udsNEXT(n_now,1) = stroke_max_udsNEXT_L(n_now,1) - stroke_max_udsNEXT_R(n_now,1);
            
            % pitch data
            pitch_mean_wb_L(n_now,1) = circ_mean(pitch_wb_L_now);
            pitch_mean_wb_R(n_now,1) = circ_mean(pitch_wb_R_now);
            dpitch_mean_wb(n_now,1) = pitch_mean_wb_L(n_now,1) - pitch_mean_wb_R(n_now,1);

            pitch_max_wb_L(n_now,1) = max(pitch_wb_L_now);
            pitch_max_wb_R(n_now,1) = max(pitch_wb_R_now);
            dpitch_max_wb(n_now,1) = max(pitch_wb_L_now) - max(pitch_wb_R_now);
            
            pitch_min_wb_L(n_now,1) = min(pitch_wb_L_now);
            pitch_min_wb_R(n_now,1) = min(pitch_wb_R_now);
            dpitch_min_wb(n_now,1) = min(pitch_wb_L_now) - min(pitch_wb_R_now);
            
            pitch_mean_ds_L(n_now,1) = circ_mean(pitch_ds_L_now);
            pitch_mean_ds_R(n_now,1) = circ_mean(pitch_ds_R_now);
            dpitch_mean_ds(n_now,1) = pitch_mean_ds_L(n_now,1) - pitch_mean_ds_R(n_now,1);
            
            pitch_mean_us_L(n_now,1) = circ_mean(pitch_us_L_now);
            pitch_mean_us_R(n_now,1) = circ_mean(pitch_us_R_now);
            dpitch_mean_us(n_now,1) = pitch_mean_us_L(n_now,1) - pitch_mean_us_R(n_now,1);

            pitch_max_ds_L(n_now,1) = max(pitch_ds_L_now);
            pitch_max_ds_R(n_now,1) = max(pitch_ds_R_now);
            dpitch_max_ds(n_now,1) = max(pitch_ds_L_now) - max(pitch_ds_R_now);
            
            pitch_min_us_L(n_now,1) = min(pitch_us_L_now);
            pitch_min_us_R(n_now,1) = min(pitch_us_R_now);
            dpitch_min_us(n_now,1) = min(pitch_us_L_now) - min(pitch_us_R_now);
            
            Apitch_wb_L(n_now,1) = max(pitch_wb_L_now)-min(pitch_wb_L_now);
            Apitch_wb_R(n_now,1) = max(pitch_wb_R_now)-min(pitch_wb_R_now);
            dApitch_wb(n_now,1) = Apitch_wb_L(n_now,1) - Apitch_wb_R(n_now,1);
            
            Apitch_ds_L(n_now,1) = max(pitch_ds_L_now)-min(pitch_ds_L_now);
            Apitch_ds_R(n_now,1) = max(pitch_ds_R_now)-min(pitch_ds_R_now);
            dApitch_ds(n_now,1) = Apitch_ds_L(n_now,1) - Apitch_ds_R(n_now,1);
            
            Apitch_us_L(n_now,1) = max(pitch_us_L_now)-min(pitch_us_L_now);
            Apitch_us_R(n_now,1) = max(pitch_us_R_now)-min(pitch_us_R_now);
            dApitch_us(n_now,1) = Apitch_us_L(n_now,1) - Apitch_us_R(n_now,1);
            
            pitch_mid_ds_L(n_now,1) = pitch_ds_L_now(round(length(pitch_ds_L_now)/2));
            pitch_mid_ds_R(n_now,1) = pitch_ds_R_now(round(length(pitch_ds_R_now)/2));
            dpitch_mid_ds(n_now,1) = pitch_mid_ds_L(n_now,1) - pitch_mid_ds_R(n_now,1);

            pitch_mid_us_L(n_now,1) = pitch_us_L_now(round(length(pitch_us_L_now)/2));
            pitch_mid_us_R(n_now,1) = pitch_us_R_now(round(length(pitch_us_R_now)/2));
            dpitch_mid_us(n_now,1) = pitch_mid_us_L(n_now,1) - pitch_mid_us_R(n_now,1);

            Apitch_mid_L(n_now,1) = pitch_mid_ds_L(n_now,1) - pitch_mid_us_L(n_now,1);
            Apitch_mid_R(n_now,1) = pitch_mid_ds_R(n_now,1) - pitch_mid_us_R(n_now,1);
            dApitch_mid(n_now,1) = Apitch_mid_L(n_now,1) - Apitch_mid_R(n_now,1);
            
            % dev data
            dev_mean_wb_L(n_now,1) = circ_mean(dev_wb_L_now);
            dev_mean_ds_L(n_now,1) = circ_mean(dev_ds_L_now);
            dev_mean_us_L(n_now,1) = circ_mean(dev_us_L_now);
            dev_mean_dus_L(n_now,1) = circ_mean(dev_dus_L_now);
            dev_mean_udsPREV_L(n_now,1) = circ_mean(dev_uds_L_prev);
            dev_mean_udsNEXT_L(n_now,1) = circ_mean(dev_uds_L_next);
            
            dev_mean_wb_R(n_now,1) = circ_mean(dev_wb_R_now);
            dev_mean_ds_R(n_now,1) = circ_mean(dev_ds_R_now);
            dev_mean_us_R(n_now,1) = circ_mean(dev_us_R_now);
            dev_mean_dus_R(n_now,1) = circ_mean(dev_dus_R_now);
            dev_mean_udsPREV_R(n_now,1) = circ_mean(dev_uds_R_prev);
            dev_mean_udsNEXT_R(n_now,1) = circ_mean(dev_uds_R_next);
           
            dev_max_dus_L(n_now,1) = max(dev_dus_L_now);
            dev_max_udsPREV_L(n_now,1) = max(dev_uds_L_prev);
            dev_max_udsNEXT_L(n_now,1) = max(dev_uds_L_next);
            dev_min_ds_L(n_now,1) = min(dev_ds_L_now);
            dev_min_us_L(n_now,1) = min(dev_us_L_now);
            
            dev_max_dus_R(n_now,1) = max(dev_dus_R_now);
            dev_max_udsPREV_R(n_now,1) = max(dev_uds_R_prev);
            dev_max_udsNEXT_R(n_now,1) = max(dev_uds_R_next);
            dev_min_ds_R(n_now,1) = min(dev_ds_R_now);
            dev_min_us_R(n_now,1) = min(dev_us_R_now);
            
            Adev_ds_L(n_now,1) = dev_max_udsPREV_L(n_now,1) - dev_min_ds_L(n_now,1);
            Adev_us_L(n_now,1) = dev_max_dus_L(n_now,1) - dev_min_us_L(n_now,1);
            Adev_ds_R(n_now,1) = dev_max_udsPREV_R(n_now,1) - dev_min_ds_R(n_now,1);
            Adev_us_R(n_now,1) = dev_max_dus_R(n_now,1) - dev_min_us_R(n_now,1);
            dAdev_ds(n_now,1) = Adev_ds_L(n_now,1) - Adev_ds_R(n_now,1);
            dAdev_us(n_now,1) = Adev_us_L(n_now,1) - Adev_us_R(n_now,1);

            ddev_mean_wb(n_now,1) = dev_mean_wb_L(n_now,1) - dev_mean_wb_R(n_now,1);
            ddev_mean_ds(n_now,1) = dev_mean_ds_L(n_now,1) - dev_mean_ds_R(n_now,1);
            ddev_mean_us(n_now,1) = dev_mean_us_L(n_now,1) - dev_mean_us_R(n_now,1);
            ddev_mean_dus(n_now,1) = dev_mean_dus_L(n_now,1) - dev_mean_dus_R(n_now,1);
            ddev_mean_udsPREV(n_now,1) = dev_mean_udsPREV_L(n_now,1) - dev_mean_udsPREV_R(n_now,1);
            ddev_mean_udsNEXT(n_now,1) = dev_mean_udsNEXT_L(n_now,1) - dev_mean_udsNEXT_R(n_now,1);
            
            ddev_max_dus(n_now,1) = dev_max_dus_L(n_now,1) - dev_max_dus_R(n_now,1);
            ddev_max_udsPREV(n_now,1) = dev_max_udsPREV_L(n_now,1) - dev_max_udsPREV_R(n_now,1);
            ddev_max_udsNEXT(n_now,1) = dev_max_udsNEXT_L(n_now,1) - dev_max_udsNEXT_R(n_now,1);
            ddev_min_ds(n_now,1) = dev_min_ds_L(n_now,1) - dev_min_ds_R(n_now,1);
            ddev_min_us(n_now,1) = dev_min_us_L(n_now,1) - dev_min_us_R(n_now,1);
            
            
            % U data
            U_mean_wb_L(n_now,1) = circ_mean(U_wb_L_now);
            U_mean_ds_L(n_now,1) = circ_mean(U_ds_L_now);
            U_mean_us_L(n_now,1) = circ_mean(U_us_L_now);
            U_mean_dus_L(n_now,1) = circ_mean(U_dus_L_now);
            U_mean_udsPREV_L(n_now,1) = circ_mean(U_uds_L_prev);
            U_mean_udsNEXT_L(n_now,1) = circ_mean(U_uds_L_next);
            
            U_mean_wb_R(n_now,1) = circ_mean(U_wb_R_now);
            U_mean_ds_R(n_now,1) = circ_mean(U_ds_R_now);
            U_mean_us_R(n_now,1) = circ_mean(U_us_R_now);
            U_mean_dus_R(n_now,1) = circ_mean(U_dus_R_now);
            U_mean_udsPREV_R(n_now,1) = circ_mean(U_uds_R_prev);
            U_mean_udsNEXT_R(n_now,1) = circ_mean(U_uds_R_next);
           
            U_max_dus_L(n_now,1) = max(U_dus_L_now);
            U_max_udsPREV_L(n_now,1) = max(U_uds_L_prev);
            U_max_udsNEXT_L(n_now,1) = max(U_uds_L_next);
            U_min_ds_L(n_now,1) = min(U_ds_L_now);
            U_min_us_L(n_now,1) = min(U_us_L_now);
            
            U_max_dus_R(n_now,1) = max(U_dus_R_now);
            U_max_udsPREV_R(n_now,1) = max(U_uds_R_prev);
            U_max_udsNEXT_R(n_now,1) = max(U_uds_R_next);
            U_min_ds_R(n_now,1) = min(U_ds_R_now);
            U_min_us_R(n_now,1) = min(U_us_R_now);
            
            AU_ds_L(n_now,1) = U_max_udsPREV_L(n_now,1) - U_min_ds_L(n_now,1);
            AU_us_L(n_now,1) = U_max_dus_L(n_now,1) - U_min_us_L(n_now,1);
            AU_ds_R(n_now,1) = U_max_udsPREV_R(n_now,1) - U_min_ds_R(n_now,1);
            AU_us_R(n_now,1) = U_max_dus_R(n_now,1) - U_min_us_R(n_now,1);
            dAU_ds(n_now,1) = AU_ds_L(n_now,1) - AU_ds_R(n_now,1);
            dAU_us(n_now,1) = AU_us_L(n_now,1) - AU_us_R(n_now,1);

            dU_mean_wb(n_now,1) = U_mean_wb_L(n_now,1) - U_mean_wb_R(n_now,1);
            dU_mean_ds(n_now,1) = U_mean_ds_L(n_now,1) - U_mean_ds_R(n_now,1);
            dU_mean_us(n_now,1) = U_mean_us_L(n_now,1) - U_mean_us_R(n_now,1);
            dU_mean_dus(n_now,1) = U_mean_dus_L(n_now,1) - U_mean_dus_R(n_now,1);
            dU_mean_udsPREV(n_now,1) = U_mean_udsPREV_L(n_now,1) - U_mean_udsPREV_R(n_now,1);
            dU_mean_udsNEXT(n_now,1) = U_mean_udsNEXT_L(n_now,1) - U_mean_udsNEXT_R(n_now,1);
            
            dU_max_dus(n_now,1) = U_max_dus_L(n_now,1) - U_max_dus_R(n_now,1);
            dU_max_udsPREV(n_now,1) = U_max_udsPREV_L(n_now,1) - U_max_udsPREV_R(n_now,1);
            dU_max_udsNEXT(n_now,1) = U_max_udsNEXT_L(n_now,1) - U_max_udsNEXT_R(n_now,1);
            dU_min_ds(n_now,1) = U_min_ds_L(n_now,1) - U_min_ds_R(n_now,1);
            dU_min_us(n_now,1) = U_min_us_L(n_now,1) - U_min_us_R(n_now,1);
            
            
            % aoa data
            aoa_mean_wb_L(n_now,1) = circ_mean(aoa_wb_L_now);
            aoa_mean_ds_L(n_now,1) = circ_mean(aoa_ds_L_now);
            aoa_mean_us_L(n_now,1) = circ_mean(aoa_us_L_now);
            aoa_mean_dus_L(n_now,1) = circ_mean(aoa_dus_L_now);
            aoa_mean_udsPREV_L(n_now,1) = circ_mean(aoa_uds_L_prev);
            aoa_mean_udsNEXT_L(n_now,1) = circ_mean(aoa_uds_L_next);
            
            aoa_mean_wb_R(n_now,1) = circ_mean(aoa_wb_R_now);
            aoa_mean_ds_R(n_now,1) = circ_mean(aoa_ds_R_now);
            aoa_mean_us_R(n_now,1) = circ_mean(aoa_us_R_now);
            aoa_mean_dus_R(n_now,1) = circ_mean(aoa_dus_R_now);
            aoa_mean_udsPREV_R(n_now,1) = circ_mean(aoa_uds_R_prev);
            aoa_mean_udsNEXT_R(n_now,1) = circ_mean(aoa_uds_R_next);
           
            aoa_max_dus_L(n_now,1) = max(aoa_dus_L_now);
            aoa_max_udsPREV_L(n_now,1) = max(aoa_uds_L_prev);
            aoa_max_udsNEXT_L(n_now,1) = max(aoa_uds_L_next);
            aoa_min_ds_L(n_now,1) = min(aoa_ds_L_now);
            aoa_min_us_L(n_now,1) = min(aoa_us_L_now);
            
            aoa_max_dus_R(n_now,1) = max(aoa_dus_R_now);
            aoa_max_udsPREV_R(n_now,1) = max(aoa_uds_R_prev);
            aoa_max_udsNEXT_R(n_now,1) = max(aoa_uds_R_next);
            aoa_min_ds_R(n_now,1) = min(aoa_ds_R_now);
            aoa_min_us_R(n_now,1) = min(aoa_us_R_now);
            
            Aaoa_ds_L(n_now,1) = aoa_max_udsPREV_L(n_now,1) - aoa_min_ds_L(n_now,1);
            Aaoa_us_L(n_now,1) = aoa_max_dus_L(n_now,1) - aoa_min_us_L(n_now,1);
            Aaoa_ds_R(n_now,1) = aoa_max_udsPREV_R(n_now,1) - aoa_min_ds_R(n_now,1);
            Aaoa_us_R(n_now,1) = aoa_max_dus_R(n_now,1) - aoa_min_us_R(n_now,1);
            dAaoa_ds(n_now,1) = Aaoa_ds_L(n_now,1) - Aaoa_ds_R(n_now,1);
            dAaoa_us(n_now,1) = Aaoa_us_L(n_now,1) - Aaoa_us_R(n_now,1);

            daoa_mean_wb(n_now,1) = aoa_mean_wb_L(n_now,1) - aoa_mean_wb_R(n_now,1);
            daoa_mean_ds(n_now,1) = aoa_mean_ds_L(n_now,1) - aoa_mean_ds_R(n_now,1);
            daoa_mean_us(n_now,1) = aoa_mean_us_L(n_now,1) - aoa_mean_us_R(n_now,1);
            daoa_mean_dus(n_now,1) = aoa_mean_dus_L(n_now,1) - aoa_mean_dus_R(n_now,1);
            daoa_mean_udsPREV(n_now,1) = aoa_mean_udsPREV_L(n_now,1) - aoa_mean_udsPREV_R(n_now,1);
            daoa_mean_udsNEXT(n_now,1) = aoa_mean_udsNEXT_L(n_now,1) - aoa_mean_udsNEXT_R(n_now,1);
            
            daoa_max_dus(n_now,1) = aoa_max_dus_L(n_now,1) - aoa_max_dus_R(n_now,1);
            daoa_max_udsPREV(n_now,1) = aoa_max_udsPREV_L(n_now,1) - aoa_max_udsPREV_R(n_now,1);
            daoa_max_udsNEXT(n_now,1) = aoa_max_udsNEXT_L(n_now,1) - aoa_max_udsNEXT_R(n_now,1);
            daoa_min_ds(n_now,1) = aoa_min_ds_L(n_now,1) - aoa_min_ds_R(n_now,1);
            daoa_min_us(n_now,1) = aoa_min_us_L(n_now,1) - aoa_min_us_R(n_now,1);
            
            % body kin data
            V_mean_wb(n_now,1) = nanmean(V_wb_now);
            pitch_global_mean_wb(n_now,1) = nanmean(pitch_global_wb_now);
            l_wing_mean_wb(n_now,1) = nanmean(l_wing_wb_now);

            roll_mean_wb(n_now,1) = nanmean(roll_wb_now);
            roll_mean_ds(n_now,1) = nanmean(roll_ds_now);
            roll_mean_us(n_now,1) = nanmean(roll_us_now);
            pitch_mean_wb(n_now,1) = nanmean(pitch_wb_now);
            pitch_mean_ds(n_now,1) = nanmean(pitch_ds_now);
            pitch_mean_us(n_now,1) = nanmean(pitch_us_now);
            yaw_mean_wb(n_now,1) = nanmean(yaw_wb_now);
            yaw_mean_ds(n_now,1) = nanmean(yaw_ds_now);
            yaw_mean_us(n_now,1) = nanmean(yaw_us_now);
            
            roll_dot_mean_wb(n_now,1) = nanmean(roll_dot_wb_now);
            roll_dot_mean_ds(n_now,1) = nanmean(roll_dot_ds_now);
            roll_dot_mean_us(n_now,1) = nanmean(roll_dot_us_now);
            pitch_dot_mean_wb(n_now,1) = nanmean(pitch_dot_wb_now);
            pitch_dot_mean_ds(n_now,1) = nanmean(pitch_dot_ds_now);
            pitch_dot_mean_us(n_now,1) = nanmean(pitch_dot_us_now);
            yaw_dot_mean_wb(n_now,1) = nanmean(yaw_dot_wb_now);
            yaw_dot_mean_ds(n_now,1) = nanmean(yaw_dot_ds_now);
            yaw_dot_mean_us(n_now,1) = nanmean(yaw_dot_us_now);
            
            roll_dot_dot_mean_wb(n_now,1) = nanmean(roll_dot_dot_wb_now);
            roll_dot_dot_mean_ds(n_now,1) = nanmean(roll_dot_dot_ds_now);
            roll_dot_dot_mean_us(n_now,1) = nanmean(roll_dot_dot_us_now);
            pitch_dot_dot_mean_wb(n_now,1) = nanmean(pitch_dot_dot_wb_now);
            pitch_dot_dot_mean_ds(n_now,1) = nanmean(pitch_dot_dot_ds_now);
            pitch_dot_dot_mean_us(n_now,1) = nanmean(pitch_dot_dot_us_now);
            yaw_dot_dot_mean_wb(n_now,1) = nanmean(yaw_dot_dot_wb_now);
            yaw_dot_dot_mean_ds(n_now,1) = nanmean(yaw_dot_dot_ds_now);
            yaw_dot_dot_mean_us(n_now,1) = nanmean(yaw_dot_dot_us_now);
            
            drot_L_mean_wb(n_now,1) = nanmean(drot_L_wb_now);
            drot_L_mean_ds(n_now,1) = nanmean(drot_L_ds_now);
            drot_L_mean_us(n_now,1) = nanmean(drot_L_us_now);

            drot_R_mean_wb(n_now,1) = nanmean(drot_R_wb_now);
            drot_R_mean_ds(n_now,1) = nanmean(drot_R_ds_now);
            drot_R_mean_us(n_now,1) = nanmean(drot_R_us_now);

            rot_dot_L_mean_wb(n_now,1) = nanmean(rot_dot_L_wb_now);
            rot_dot_L_mean_ds(n_now,1) = nanmean(rot_dot_L_ds_now);
            rot_dot_L_mean_us(n_now,1) = nanmean(rot_dot_L_us_now);

            rot_dot_R_mean_wb(n_now,1) = nanmean(rot_dot_R_wb_now);
            rot_dot_R_mean_ds(n_now,1) = nanmean(rot_dot_R_ds_now);
            rot_dot_R_mean_us(n_now,1) = nanmean(rot_dot_R_us_now);

            rot_dot_dot_L_mean_wb(n_now,1) = nanmean(rot_dot_dot_L_wb_now);
            rot_dot_dot_L_mean_ds(n_now,1) = nanmean(rot_dot_dot_L_ds_now);
            rot_dot_dot_L_mean_us(n_now,1) = nanmean(rot_dot_dot_L_us_now);

            rot_dot_dot_R_mean_wb(n_now,1) = nanmean(rot_dot_dot_R_wb_now);
            rot_dot_dot_R_mean_ds(n_now,1) = nanmean(rot_dot_dot_R_ds_now);
            rot_dot_dot_R_mean_us(n_now,1) = nanmean(rot_dot_dot_R_us_now);

            F_mean_wb(n_now,1) = nanmean(F_wb_now);
            F_mean_ds(n_now,1) = nanmean(F_ds_now);
            F_mean_us(n_now,1) = nanmean(F_us_now);
            Fsp_roll_mean_wb(n_now,1) = circ_mean(Fsp_roll_wb_now);
            Fsp_roll_mean_ds(n_now,1) = circ_mean(Fsp_roll_ds_now);
            Fsp_roll_mean_us(n_now,1) = circ_mean(Fsp_roll_us_now);
            Fsp_pitch_mean_wb(n_now,1) = circ_mean(Fsp_pitch_wb_now);
            Fsp_pitch_mean_ds(n_now,1) = circ_mean(Fsp_pitch_ds_now);
            Fsp_pitch_mean_us(n_now,1) = circ_mean(Fsp_pitch_us_now);
        end
    end
end
n_wb = n_now

% Torques
Mroll_mean_wb_accel  = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Iroll  .* roll_dot_dot_mean_wb;
Mpitch_mean_wb_accel = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Ipitch .* pitch_dot_dot_mean_wb;
Myaw_mean_wb_accel   = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Iyaw   .* yaw_dot_dot_mean_wb;

Mroll_mean_wb_damp  = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Croll  .* roll_dot_mean_wb;
Mpitch_mean_wb_damp = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Cpitch .* pitch_dot_mean_wb;
Myaw_mean_wb_damp   = l_wing_mean_wb.^5 .* f_wb_L.^2 .* Cyaw   .* yaw_dot_mean_wb;

Mroll_mean_wb  = Mroll_mean_wb_accel  + Mroll_mean_wb_damp;
Mpitch_mean_wb = Mpitch_mean_wb_accel + Mpitch_mean_wb_damp;
Myaw_mean_wb   = Myaw_mean_wb_accel   + Myaw_mean_wb_damp;

% Normalized accel
roll_dot_dot_mean_wb_NOnorm  =  roll_dot_dot_mean_wb;
pitch_dot_dot_mean_wb_NOnorm  =  pitch_dot_dot_mean_wb;
yaw_dot_dot_mean_wb_NOnorm  =  yaw_dot_dot_mean_wb;

roll_dot_dot_mean_wb  =  roll_dot_dot_mean_wb ./ f_wb_L.^2;
pitch_dot_dot_mean_wb  =  pitch_dot_dot_mean_wb ./ f_wb_L.^2;
yaw_dot_dot_mean_wb  =  yaw_dot_dot_mean_wb ./ f_wb_L.^2;

roll_dot_dot_mean_wb_Norm  =  roll_dot_dot_mean_wb;
pitch_dot_dot_mean_wb_Norm  =  pitch_dot_dot_mean_wb;
yaw_dot_dot_mean_wb_Norm  =  yaw_dot_dot_mean_wb;

rot_dot_dot_L_mean_wb_NOnorm  =  rot_dot_dot_L_mean_wb;
rot_dot_dot_R_mean_wb_NOnorm  =  rot_dot_dot_R_mean_wb;

rot_dot_dot_L_mean_wb  =  rot_dot_dot_L_mean_wb ./ f_wb_L.^2;
rot_dot_dot_R_mean_wb  =  rot_dot_dot_R_mean_wb ./ f_wb_L.^2;

rot_dot_dot_L_mean_wb_Norm  =  rot_dot_dot_L_mean_wb;
rot_dot_dot_R_mean_wb_Norm  =  rot_dot_dot_R_mean_wb;

% save wb
save(['WBdataset_all_',num2str(n_wb),'WBs_accelNorm.mat'],...
...
't_wb_bin',...
't_ds_bin',...
't_us_bin',...
...
'seq_nr',...
'wb_nr',...
...
'V_mean_wb',...
'pitch_global_mean_wb',...
'l_wing_mean_wb',...
...
'Iroll',...
'Ipitch',...
'Iyaw',...
...
'Croll',...
'Cpitch',...
'Cyaw',...
...
'Mroll_mean_wb',...
'Mpitch_mean_wb',...
'Myaw_mean_wb',...
...
'Mroll_mean_wb_accel',...
'Mpitch_mean_wb_accel',...
'Myaw_mean_wb_accel',...
...
'Mroll_mean_wb_damp',...
'Mpitch_mean_wb_damp',...
'Myaw_mean_wb_damp',...
...
'roll_dot_dot_mean_wb_NOnorm',...
'pitch_dot_dot_mean_wb_NOnorm',...
'yaw_dot_dot_mean_wb_NOnorm',...
...
'roll_dot_dot_mean_wb_Norm',...
'pitch_dot_dot_mean_wb_Norm',...
'yaw_dot_dot_mean_wb_Norm',...
...
'rot_dot_dot_L_mean_wb_NOnorm',...
'rot_dot_dot_R_mean_wb_NOnorm',...
...
'rot_dot_dot_L_mean_wb_Norm',...
'rot_dot_dot_R_mean_wb_Norm',...
...
'f_wb_L',...
'f_wb_R',...
'f_ds_L',...
'f_ds_R',...
'f_us_L',...
'f_us_R',...
...
'dt_wb_L',...
'dt_wb_R',...
'dt_ds_L',...
'dt_ds_R',...
'dt_us_L',...
'dt_us_R',...
...
'df_wb',...
'df_ds',...
'df_us',...
...
'ddt_wb',...
'ddt_ds',...
'ddt_us',...
...
'stroke_mean_wb_L',...
'stroke_mean_ds_L',...
'stroke_mean_us_L',...
'stroke_mean_wb_R',...
'stroke_mean_ds_R',...
'stroke_mean_us_R',...
...
'stroke_max_wb_L',...
'stroke_min_wb_L',...
'stroke_max_wb_R',...
'stroke_min_wb_R',...
...
'stroke_max_udsPREV_L',...
'stroke_max_udsNEXT_L',...
'stroke_max_udsPREV_R',...
'stroke_max_udsNEXT_R',...
...
'Astroke_wb_L',...
'Astroke_wb_R',...
'dAstroke_wb',...
...
'Astroke_ds_L',...
'Astroke_ds_R',...
'dAstroke_ds',...
...
'Astroke_us_L',...
'Astroke_us_R',...
'dAstroke_us',...
...
'dstroke_mean_wb',...
'dstroke_mean_ds',...
'dstroke_mean_us',...
...
'dstroke_max_wb',...
'dstroke_min_wb',...
'dstroke_max_udsPREV',...
'dstroke_max_udsNEXT',...
...
'pitch_mean_wb_L',...
'pitch_mean_wb_R',...
'dpitch_mean_wb',...
...
'pitch_max_wb_L',...
'pitch_max_wb_R',...
'dpitch_max_wb',...
...
'pitch_min_wb_L',...
'pitch_min_wb_R',...
'dpitch_min_wb',...
...
'pitch_mean_ds_L',...
'pitch_mean_ds_R',...
'dpitch_mean_ds',...
...
'pitch_mean_us_L',...
'pitch_mean_us_R',...
'dpitch_mean_us',...
...
'pitch_max_ds_L',...
'pitch_max_ds_R',...
'dpitch_max_ds',...
...
'pitch_min_us_L',...
'pitch_min_us_R',...
'dpitch_min_us',...
...
'Apitch_wb_L',...
'Apitch_wb_R',...
'dApitch_wb',...
...
'Apitch_ds_L',...
'Apitch_ds_R',...
'dApitch_ds',...
...
'Apitch_us_L',...
'Apitch_us_R',...
'dApitch_us',...
...
'pitch_mid_ds_L',...
'pitch_mid_ds_R',...
...
'pitch_mid_us_L',...
'pitch_mid_us_R',...
...
'Apitch_mid_L',...
'Apitch_mid_R',...
...
'dpitch_mid_ds',...
'dpitch_mid_us',...
'dApitch_mid',...
...
'dev_mean_wb_L',...
'dev_mean_ds_L',...
'dev_mean_us_L',...
'dev_mean_dus_L',...
'dev_mean_udsPREV_L',...
'dev_mean_udsNEXT_L',...
...
'dev_mean_wb_R',...
'dev_mean_ds_R',...
'dev_mean_us_R',...
'dev_mean_dus_R',...
'dev_mean_udsPREV_R',...
'dev_mean_udsNEXT_R',...
...
'dev_max_dus_L',...
'dev_max_udsPREV_L',...
'dev_max_udsNEXT_L',...
'dev_min_ds_L',...
'dev_min_us_L',...
...
'dev_max_dus_R',...
'dev_max_udsPREV_R',...
'dev_max_udsNEXT_R',...
'dev_min_ds_R',...
'dev_min_us_R',...
...
'Adev_ds_L',...
'Adev_us_L',...
'Adev_ds_R',...
'Adev_us_R',...
'dAdev_ds',...
'dAdev_us',...
...
'ddev_mean_wb',...
'ddev_mean_ds',...
'ddev_mean_us',...
'ddev_mean_dus',...
'ddev_mean_udsPREV',...
'ddev_mean_udsNEXT',...
...
'ddev_max_dus',...
'ddev_max_udsPREV',...
'ddev_max_udsNEXT',...
'ddev_min_ds',...
'ddev_min_us',...
...
'U_mean_wb_L',...
'U_mean_ds_L',...
'U_mean_us_L',...
'U_mean_dus_L',...
'U_mean_udsPREV_L',...
'U_mean_udsNEXT_L',...
...
'U_mean_wb_R',...
'U_mean_ds_R',...
'U_mean_us_R',...
'U_mean_dus_R',...
'U_mean_udsPREV_R',...
'U_mean_udsNEXT_R',...
...
'U_max_dus_L',...
'U_max_udsPREV_L',...
'U_max_udsNEXT_L',...
'U_min_ds_L',...
'U_min_us_L',...
...
'U_max_dus_R',...
'U_max_udsPREV_R',...
'U_max_udsNEXT_R',...
'U_min_ds_R',...
'U_min_us_R',...
...
'AU_ds_L',...
'AU_us_L',...
'AU_ds_R',...
'AU_us_R',...
'dAU_ds',...
'dAU_us',...
...
'dU_mean_wb',...
'dU_mean_ds',...
'dU_mean_us',...
'dU_mean_dus',...
'dU_mean_udsPREV',...
'dU_mean_udsNEXT',...
...
'dU_max_dus',...
'dU_max_udsPREV',...
'dU_max_udsNEXT',...
'dU_min_ds',...
'dU_min_us',...
...
'aoa_mean_wb_L',...
'aoa_mean_ds_L',...
'aoa_mean_us_L',...
'aoa_mean_dus_L',...
'aoa_mean_udsPREV_L',...
'aoa_mean_udsNEXT_L',...
...
'aoa_mean_wb_R',...
'aoa_mean_ds_R',...
'aoa_mean_us_R',...
'aoa_mean_dus_R',...
'aoa_mean_udsPREV_R',...
'aoa_mean_udsNEXT_R',...
...
'aoa_max_dus_L',...
'aoa_max_udsPREV_L',...
'aoa_max_udsNEXT_L',...
'aoa_min_ds_L',...
'aoa_min_us_L',...
...
'aoa_max_dus_R',...
'aoa_max_udsPREV_R',...
'aoa_max_udsNEXT_R',...
'aoa_min_ds_R',...
'aoa_min_us_R',...
...
'Aaoa_ds_L',...
'Aaoa_us_L',...
'Aaoa_ds_R',...
'Aaoa_us_R',...
'dAaoa_ds',...
'dAaoa_us',...
...
'daoa_mean_wb',...
'daoa_mean_ds',...
'daoa_mean_us',...
'daoa_mean_dus',...
'daoa_mean_udsPREV',...
'daoa_mean_udsNEXT',...
...
'daoa_max_dus',...
'daoa_max_udsPREV',...
'daoa_max_udsNEXT',...
'daoa_min_ds',...
'daoa_min_us',...
...
'roll_mean_wb',...
'roll_mean_ds',...
'roll_mean_us',...
'pitch_mean_wb',...
'pitch_mean_ds',...
'pitch_mean_us',...
'yaw_mean_wb',...
'yaw_mean_ds',...
'yaw_mean_us',...
...
'roll_dot_mean_wb',...
'roll_dot_mean_ds',...
'roll_dot_mean_us',...
'pitch_dot_mean_wb',...
'pitch_dot_mean_ds',...
'pitch_dot_mean_us',...
'yaw_dot_mean_wb',...
'yaw_dot_mean_ds',...
'yaw_dot_mean_us',...
...
'roll_dot_dot_mean_wb',...
'roll_dot_dot_mean_ds',...
'roll_dot_dot_mean_us',...
'pitch_dot_dot_mean_wb',...
'pitch_dot_dot_mean_ds',...
'pitch_dot_dot_mean_us',...
'yaw_dot_dot_mean_wb',...
'yaw_dot_dot_mean_ds',...
'yaw_dot_dot_mean_us',...
...
'drot_L_mean_wb',...
'drot_L_mean_ds',...
'drot_L_mean_us',...
...
'drot_R_mean_wb',...
'drot_R_mean_ds',...
'drot_R_mean_us',...
...
'rot_dot_L_mean_wb',...
'rot_dot_L_mean_ds',...
'rot_dot_L_mean_us',...
...
'rot_dot_R_mean_wb',...
'rot_dot_R_mean_ds',...
'rot_dot_R_mean_us',...
...
'rot_dot_dot_L_mean_wb',...
'rot_dot_dot_L_mean_ds',...
'rot_dot_dot_L_mean_us',...
...
'rot_dot_dot_R_mean_wb',...
'rot_dot_dot_R_mean_ds',...
'rot_dot_dot_R_mean_us',...
...
'F_mean_wb',...
'F_mean_ds',...
'F_mean_us',...
'Fsp_roll_mean_wb',...
'Fsp_roll_mean_ds',...
'Fsp_roll_mean_us',...
'Fsp_pitch_mean_wb',...
'Fsp_pitch_mean_ds',...
'Fsp_pitch_mean_us',...
...
't_wb_L',...
't_wb_R',...
...
't_ds_L',...
't_ds_R',...
...
't_us_L',...
't_us_R',...
...
'stroke_wb_L',...
'stroke_wb_R',...
'pitch_wb_L',...
'pitch_wb_R',...
'dev_wb_L',...
'dev_wb_R',...
'aoa_wb_L',...
'aoa_wb_R',...
'U_wb_L',...
'U_wb_R',...
'Dstroke_wb',...
'Dpitch_wb',...
'Ddev_wb',...
'Daoa_wb',...
'DU_wb',...
...
'stroke_ds_L',...
'stroke_ds_R',...
'pitch_ds_L',...
'pitch_ds_R',...
'dev_ds_L',...
'dev_ds_R',...
'aoa_ds_L',...
'aoa_ds_R',...
'U_ds_L',...
'U_ds_R',...
'Dstroke_ds',...
'Dpitch_ds',...
'Ddev_ds',...
'Daoa_ds',...
'DU_ds',...
...
'stroke_us_L',...
'stroke_us_R',...
'pitch_us_L',...
'pitch_us_R',...
'dev_us_L',...
'dev_us_R',...
'aoa_us_L',...
'aoa_us_R',...
'U_us_L',...
'U_us_R',...
'Dstroke_us',...
'Dpitch_us',...
'Ddev_us',...
'Daoa_us',...
'DU_us',...
...
't_wb_bins',...
't_ds_bins',...
't_us_bins',...
...
'stroke_wb_L_bins',...
'stroke_wb_R_bins',...
'pitch_wb_L_bins',...
'pitch_wb_R_bins',...
'dev_wb_L_bins',...
'dev_wb_R_bins',...
'aoa_wb_L_bins',...
'aoa_wb_R_bins',...
'U_wb_L_bins',...
'U_wb_R_bins',...
'Dstroke_wb_bins',...
'Dpitch_wb_bins',...
'Ddev_wb_bins',...
'Daoa_wb_bins',...
'DU_wb_bins',...
...
'stroke_ds_L_bins',...
'stroke_ds_R_bins',...
'pitch_ds_L_bins',...
'pitch_ds_R_bins',...
'dev_ds_L_bins',...
'dev_ds_R_bins',...
'aoa_ds_L_bins',...
'aoa_ds_R_bins',...
'U_ds_L_bins',...
'U_ds_R_bins',...
'Dstroke_ds_bins',...
'Dpitch_ds_bins',...
'Ddev_ds_bins',...
'Daoa_ds_bins',...
'DU_ds_bins',...
...
'stroke_us_L_bins',...
'stroke_us_R_bins',...
'pitch_us_L_bins',...
'pitch_us_R_bins',...
'dev_us_L_bins',...
'dev_us_R_bins',...
'aoa_us_L_bins',...
'aoa_us_R_bins',...
'U_us_L_bins',...
'U_us_R_bins',...
'Dstroke_us_bins',...
'Dpitch_us_bins',...
'Ddev_us_bins',...
'Daoa_us_bins',...
'DU_us_bins',...
...
'dt_ds',...
'dt_us',...
'dt_wb',...
'f_wb');
