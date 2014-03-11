%% accelerations
% A
[A_max n_Amax] = calc_max_value(A,n_pre,n_post);
[An_max n_Anmax] = calc_max_value(An,n_pre,n_post);
[At_max n_Atmax] = calc_max_value(At,n_pre,n_post);

% A_hor
[A_hor_max n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

[An_hor_max n_Anhormax] = calc_max_abs_value(An_hor,n_pre,n_post);
[At_hor_max n_Athormax] = calc_max_value(At_hor,n_pre,n_post);
[At_hor_min n_Athormin] = calc_min_value(At_hor,n_pre,n_post);

for i = 1:length(n_turn_max)
    if isnan(n_turn_max(i))==0
        An_hor_max(i) = An_hor(n_turn_max(i),i);
        n_Anhormax(i) = n_turn_max(i);
    end
    
    if isnan(n_accel_max(i))==0
        At_hor_max(i) = At_hor(n_accel_max(i),i);
        n_Athormax(i) = n_accel_max(i);
    end
    
    if isnan(n_decel_min(i))==0
        At_hor_min(i) = At_hor(n_decel_min(i),i);
        n_Athormin(i) = n_decel_min(i);
    end
end

An_hor_mirror_max = calc_value(An_hor_mirror,n_Anhormax);

%% force
% F
[F_max n_Fmax] = calc_max_value(F,n_pre,n_post);
% [F_WBmean_max n_F_WBmeanmax] = calc_max_value(F_WBmean,n_pre,n_post);

% F_hor
[F_hor_max n_Fhormax] = calc_max_value(F_hor,n_pre,n_post);

% F_ver
[F_ver_min n_Fvermin] = calc_min_value(F_ver,n_pre,n_post);

% % F_norm
% [F_norm_max n_F_normmax] = calc_max_value(Fenhance_Norm,n_pre,n_post);
% % [F_WBmean_norm_max n_F_WBmean_normmax] = calc_max_value(F_WBmean_norm,n_pre,n_post);

%% ATTITUDE
[yaw_absmax n_yawabsmax] = calc_max_abs_value(yaw,n_pre,n_post);

[pitch_max n_pitchmax] = calc_max_value(pitch,n_pre,n_post);
[pitch_min n_pitchmin] = calc_min_value(pitch,n_pre,n_post);
[pitch_absmax n_pitchabsmax] = calc_max_abs_value(pitch,n_pre,n_post);

% MIRROR
[roll_mirror_max n_rollmirrormax] = calc_max_value(roll_mirror,n_pre,n_post);
[yaw_mirror_max n_yawmirrormax] = calc_max_value(yaw_mirror,n_pre,n_post);

[roll_mirror_min n_rollmirrormin] = calc_min_value(roll_mirror,n_pre,n_post);
[yaw_mirror_min n_yawmirrormin] = calc_min_value(yaw_mirror,n_pre,n_post);

%% ATTITUDE DOT
[pitch_dot_max n_pitchdotmax] = calc_max_value(pitch_dot,n_pre,n_post);
[pitch_dot_min n_pitchdotmin] = calc_min_value(pitch_dot,n_pre,n_post);
[pitch_dot_absmax n_pitchdotabsmax] = calc_max_abs_value(pitch_dot,n_pre,n_post);

%MIRROR
[roll_dot_mirror_max n_rolldotmirrormax] = calc_max_value(roll_dot_mirror,n_pre,n_post);
[yaw_dot_mirror_max n_yawdotmirrormax] = calc_max_value(yaw_dot_mirror,n_pre,n_post);

[roll_dot_mirror_min n_rolldotmirrormin] = calc_min_value(roll_dot_mirror,n_pre,n_post);
[yaw_dot_mirror_min n_yawdotmirrormin] = calc_min_value(yaw_dot_mirror,n_pre,n_post);

%% ATTITUDE DOT WBmean
% [pitch_dot_WBmean_max n_pitchdot_WBmeanmax] = calc_max_value(pitch_dot_WBmean,n_pre,n_post);
% [pitch_dot_WBmean_min n_pitchdot_WBmeanmin] = calc_min_value(pitch_dot_WBmean,n_pre,n_post);
% [pitch_dot_WBmean_absmax n_pitchdot_WBmeanabsmax] = calc_max_abs_value(pitch_dot_WBmean,n_pre,n_post);
% 
% %MIRROR
% [roll_dot_WBmean_mirror_max n_rolldot_WBmeanmirrormax] = calc_max_value(roll_dot_WBmean_mirror,n_pre,n_post);
% [yaw_dot_WBmean_mirror_max n_yawdot_WBmeanmirrormax] = calc_max_value(yaw_dot_WBmean_mirror,n_pre,n_post);
% 
% [roll_dot_WBmean_mirror_min n_rolldot_WBmeanmirrormin] = calc_min_value(roll_dot_WBmean_mirror,n_pre,n_post);
% [yaw_dot_WBmean_mirror_min n_yawdot_WBmeanmirrormin] = calc_min_value(yaw_dot_WBmean_mirror,n_pre,n_post);

%% ATTITUDE DOT DOT
[pitch_dot_dot_max n_pitchdot_dotmax] = calc_max_value(pitch_dot_dot,n_pre,n_post);
[pitch_dot_dot_min n_pitchdot_dotmin] = calc_min_value(pitch_dot_dot,n_pre,n_post);
[pitch_dot_dot_absmax n_pitchdot_dotabsmax] = calc_max_abs_value(pitch_dot_dot,n_pre,n_post);

%MIRROR
[roll_dot_dot_mirror_max n_rolldot_dotmirrormax] = calc_max_value(roll_dot_dot_mirror,n_pre,n_post);
[yaw_dot_dot_mirror_max n_yawdot_dotmirrormax] = calc_max_value(yaw_dot_dot_mirror,n_pre,n_post);

[roll_dot_dot_mirror_min n_rolldot_dotmirrormin] = calc_min_value(roll_dot_dot_mirror,n_pre,n_post);
[yaw_dot_dot_mirror_min n_yawdot_dotmirrormin] = calc_min_value(yaw_dot_dot_mirror,n_pre,n_post);

%% ATTITUDE DOT DOT NORM
% [pitch_dot_dot_norm_max n_pitchdot_dot_normmax] = calc_max_value(pitch_dot_dot_norm,n_pre,n_post);
% [pitch_dot_dot_norm_min n_pitchdot_dot_normmin] = calc_min_value(pitch_dot_dot_norm,n_pre,n_post);
% [pitch_dot_dot_norm_absmax n_pitchdot_dot_normabsmax] = calc_max_abs_value(pitch_dot_dot_norm,n_pre,n_post);
% 
% %MIRROR
% [roll_dot_dot_norm_mirror_max n_rolldot_dot_normmirrormax] = calc_max_value(roll_dot_dot_norm_mirror,n_pre,n_post);
% [yaw_dot_dot_norm_mirror_max n_yawdot_dot_normmirrormax] = calc_max_value(yaw_dot_dot_norm_mirror,n_pre,n_post);
% 
% [roll_dot_dot_norm_mirror_min n_rolldot_dot_normmirrormin] = calc_min_value(roll_dot_dot_norm_mirror,n_pre,n_post);
% [yaw_dot_dot_norm_mirror_min n_yawdot_dot_normmirrormin] = calc_min_value(yaw_dot_dot_norm_mirror,n_pre,n_post);

%% ATTITUDE DOT DOT WBmean
% [pitch_dot_dot_WBmean_max n_pitchdot_dot_WBmeanmax] = calc_max_value(pitch_dot_dot_WBmean,n_pre,n_post);
% [pitch_dot_dot_WBmean_min n_pitchdot_dot_WBmeanmin] = calc_min_value(pitch_dot_dot_WBmean,n_pre,n_post);
% [pitch_dot_dot_WBmean_absmax n_pitchdot_dot_WBmeanabsmax] = calc_max_abs_value(pitch_dot_dot_WBmean,n_pre,n_post);
% 
% %MIRROR
% [roll_dot_dot_WBmean_mirror_max n_rolldot_dot_WBmeanmirrormax] = calc_max_value(roll_dot_dot_WBmean_mirror,n_pre,n_post);
% [yaw_dot_dot_WBmean_mirror_max n_yawdot_dot_WBmeanmirrormax] = calc_max_value(yaw_dot_dot_WBmean_mirror,n_pre,n_post);
% 
% [roll_dot_dot_WBmean_mirror_min n_rolldot_dot_WBmeanmirrormin] = calc_min_value(roll_dot_dot_WBmean_mirror,n_pre,n_post);
% [yaw_dot_dot_WBmean_mirror_min n_yawdot_dot_WBmeanmirrormin] = calc_min_value(yaw_dot_dot_WBmean_mirror,n_pre,n_post);

%% ATTITUDE DOT DOT WBmean NORM
% [pitch_dot_dot_WBmean_norm_max n_pitchdot_dot_WBmean_normmax] = calc_max_value(pitch_dot_dot_WBmean_norm,n_pre,n_post);
% [pitch_dot_dot_WBmean_norm_min n_pitchdot_dot_WBmean_normmin] = calc_min_value(pitch_dot_dot_WBmean_norm,n_pre,n_post);
% [pitch_dot_dot_WBmean_norm_absmax n_pitchdot_dot_WBmean_normabsmax] = calc_max_abs_value(pitch_dot_dot_WBmean_norm,n_pre,n_post);
% 
% %MIRROR
% [roll_dot_dot_WBmean_norm_mirror_max n_rolldot_dot_WBmean_normmirrormax] = calc_max_value(roll_dot_dot_WBmean_norm_mirror,n_pre,n_post);
% [yaw_dot_dot_WBmean_norm_mirror_max n_yawdot_dot_WBmean_normmirrormax] = calc_max_value(yaw_dot_dot_WBmean_norm_mirror,n_pre,n_post);
% 
% [roll_dot_dot_WBmean_norm_mirror_min n_rolldot_dot_WBmean_normmirrormin] = calc_min_value(roll_dot_dot_WBmean_norm_mirror,n_pre,n_post);
% [yaw_dot_dot_WBmean_norm_mirror_min n_yawdot_dot_WBmean_normmirrormin] = calc_min_value(yaw_dot_dot_WBmean_norm_mirror,n_pre,n_post);

%% TORQUE
% [Mpitch_max n_Mpitch_max] = calc_max_value(Mpitch,n_pre,n_post);
% [Mpitch_min n_Mpitch_min] = calc_min_value(Mpitch,n_pre,n_post);
% [Mpitch_absmax n_Mpitch_absmax] = calc_max_abs_value(Mpitch,n_pre,n_post);
% 
% [Mroll_max n_Mroll_max] = calc_max_value(Mroll,n_pre,n_post);
% [Myaw_max n_Myaw_max] = calc_max_value(Myaw,n_pre,n_post);
% 
% [Mroll_min n_Mroll_min] = calc_min_value(Mroll,n_pre,n_post);
% [Myaw_min n_Myaw_min] = calc_min_value(Myaw,n_pre,n_post);

%% TORQUE NORM
% [Mpitch_norm_max n_Mpitch_norm_max] = calc_max_value(Mpitch_Norm,n_pre,n_post);
% [Mpitch_norm_min n_Mpitch_norm_min] = calc_min_value(Mpitch_Norm,n_pre,n_post);
% [Mpitch_norm_absmax n_Mpitch_norm_absmax] = calc_max_abs_value(Mpitch_Norm,n_pre,n_post);
% 
% [Mroll_norm_max n_Mroll_norm_max] = calc_max_value(Mroll_Norm,n_pre,n_post);
% [Myaw_norm_max n_Myaw_norm_max] = calc_max_value(Myaw_Norm,n_pre,n_post);
% 
% [Mroll_norm_min n_Mroll_norm_min] = calc_min_value(Mroll_Norm,n_pre,n_post);
% [Myaw_norm_min n_Myaw_norm_min] = calc_min_value(Myaw_Norm,n_pre,n_post);

%% TORQUE wingbeat mean
% [Mpitch_WBmean_max n_Mpitch_WBmean_max] = calc_max_value(Mpitch_mean_wb,n_pre,n_post);
% [Mpitch_WBmean_min n_Mpitch_WBmean_min] = calc_min_value(Mpitch_mean_wb,n_pre,n_post);
% [Mpitch_WBmean_absmax n_Mpitch_WBmean_absmax] = calc_max_abs_value(Mpitch_mean_wb,n_pre,n_post);
% 
% [Mroll_WBmean_max n_Mroll_WBmean_max] = calc_max_value(Mroll_mean_wb,n_pre,n_post);
% [Myaw_WBmean_max n_Myaw_WBmean_max] = calc_max_value(Myaw_mean_wb,n_pre,n_post);
% 
% [Mroll_WBmean_min n_Mroll_WBmean_min] = calc_min_value(Mroll_mean_wb,n_pre,n_post);
% [Myaw_WBmean_min n_Myaw_WBmean_min] = calc_min_value(Myaw_mean_wb,n_pre,n_post);

%% TORQUE wingbeat mean NORM
% [Mpitch_WBmean_norm_max n_Mpitch_WBmean_norm_max] = calc_max_value(Mpitch_WBmean_Norm,n_pre,n_post);
% [Mpitch_WBmean_norm_min n_Mpitch_WBmean_norm_min] = calc_min_value(Mpitch_WBmean_Norm,n_pre,n_post);
% [Mpitch_WBmean_norm_absmax n_Mpitch_WBmean_norm_absmax] = calc_max_abs_value(Mpitch_WBmean_Norm,n_pre,n_post);
% 
% [Mroll_WBmean_norm_max n_Mroll_WBmean_norm_max] = calc_max_value(Mroll_WBmean_Norm,n_pre,n_post);
% [Myaw_WBmean_norm_max n_Myaw_WBmean_norm_max] = calc_max_value(Myaw_WBmean_Norm,n_pre,n_post);
% 
% [Mroll_WBmean_norm_min n_Mroll_WBmean_norm_min] = calc_min_value(Mroll_WBmean_Norm,n_pre,n_post);
% [Myaw_WBmean_norm_min n_Myaw_WBmean_norm_min] = calc_min_value(Myaw_WBmean_Norm,n_pre,n_post);


