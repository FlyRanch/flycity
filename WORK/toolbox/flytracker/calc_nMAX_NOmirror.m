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

% An_hor_mirror_max = calc_value(An_hor_mirror,n_Anhormax);
% 
% %% force
% % F
% [F_max n_Fmax] = calc_max_value(F,n_pre,n_post);
% 
% % F_hor
% [F_hor_max n_Fhormax] = calc_max_value(F_hor,n_pre,n_post);
% 
% % F_ver
% [F_ver_min n_Fvermin] = calc_min_value(F_ver,n_pre,n_post);
% 
% %% ATTITUDE
% [yaw_absmax n_yawabsmax] = calc_max_abs_value(yaw,n_pre,n_post);
% 
% [pitch_max n_pitchmax] = calc_max_value(pitch,n_pre,n_post);
% [pitch_min n_pitchmin] = calc_min_value(pitch,n_pre,n_post);
% [pitch_absmax n_pitchabsmax] = calc_max_abs_value(pitch,n_pre,n_post);
% 
% % MIRROR
% [roll_mirror_max n_rollmirrormax] = calc_max_value(roll_mirror,n_pre,n_post);
% [yaw_mirror_max n_yawmirrormax] = calc_max_value(yaw_mirror,n_pre,n_post);
% 
% [roll_mirror_min n_rollmirrormin] = calc_min_value(roll_mirror,n_pre,n_post);
% [yaw_mirror_min n_yawmirrormin] = calc_min_value(yaw_mirror,n_pre,n_post);
% 
% %% ATTITUDE DOT
% [pitch_dot_max n_pitchdotmax] = calc_max_value(pitch_dot,n_pre,n_post);
% [pitch_dot_min n_pitchdotmin] = calc_min_value(pitch_dot,n_pre,n_post);
% [pitch_dot_absmax n_pitchdotabsmax] = calc_max_abs_value(pitch_dot,n_pre,n_post);
% 
% %MIRROR
% [roll_dot_mirror_max n_rolldotmirrormax] = calc_max_value(roll_dot_mirror,n_pre,n_post);
% [yaw_dot_mirror_max n_yawdotmirrormax] = calc_max_value(yaw_dot_mirror,n_pre,n_post);
% 
% [roll_dot_mirror_min n_rolldotmirrormin] = calc_min_value(roll_dot_mirror,n_pre,n_post);
% [yaw_dot_mirror_min n_yawdotmirrormin] = calc_min_value(yaw_dot_mirror,n_pre,n_post);
