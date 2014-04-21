% n_pre & n_post
n_pre = min([n_turn_start n_accel_start n_decel_start]')';
n_post = max([n_turn_stop n_accel_stop n_decel_stop]')'-1;
for i = 1:length(n_post)
    if isnan(n_post(i)) == 1
        V_temp = V(:,i);
        n_post(i,1) = find(isnan(V_temp)==0, 1, 'last' );
    end
end


%% accelerations

% A_hor
[A_hor_max n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

n_Anhormax = n_turn_max;
n_Athormax = n_accel_max;
n_Athormin = n_accel_min;

An_hor_max = calc_value(An_hor,n_turn_max);
At_hor_max = calc_value(At_hor,n_accel_max);
At_hor_min = calc_value(At_hor,n_decel_min);

% A
[A_max n_Amax] = calc_max_value(A,n_pre,n_post);
[An_max n_Anmax] = calc_max_value(An,n_pre,n_post);
[At_max n_Atmax] = calc_max_value(At,n_pre,n_post);

%% force
% F
[F_max n_Fmax] = calc_max_value(F,n_pre,n_post);

% F_hor
[F_hor_max n_Fhormax] = calc_max_value(F_hor,n_pre,n_post);

% F_ver
[F_ver_min n_Fvermin] = calc_min_value(F_ver,n_pre,n_post);

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
