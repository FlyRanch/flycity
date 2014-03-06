% n_pre & n_post
n_pre = min([n_turn_start n_accel_start n_decel_start]')';
n_post = max([n_turn_stop n_accel_stop n_decel_stop]')'-1;
for i = 1:length(n_post)
    if isnan(n_post(i)) == 1
        V_temp = V(:,i);
        n_post(i,1) = find(isnan(V_temp)==0, 1, 'last' );
    end
end

%% time
% t_pre
t_pre = calc_t_value(t,n_pre);
t_post = calc_t_value(t,n_post);

t_turn_start = calc_t_value(t,n_turn_start);
t_turn_stop = calc_t_value(t,n_turn_stop);
t_turn_max = calc_t_value(t,n_turn_max);

t_accel_start = calc_t_value(t,n_accel_start);
t_accel_stop = calc_t_value(t,n_accel_stop);
t_accel_max = calc_t_value(t,n_accel_max);

t_decel_start = calc_t_value(t,n_decel_start);
t_decel_stop = calc_t_value(t,n_decel_stop);
t_decel_min = calc_t_value(t,n_decel_min);

%% V
V_pre = calc_value(V,n_pre);
V_post = calc_value(V,n_post);
V_mean = calc_mean_value(V,n_pre,n_post);
V_trig2resp = calc_mean_value(V,n_first,n_resp);

%% accelerations
An_hor_max = calc_value(An_hor,n_turn_max);
At_hor_max = calc_value(At_hor,n_accel_max);
At_hor_min = calc_value(At_hor,n_decel_min);

% An_hor
An_hor_pre = calc_value(An_hor,n_pre);
An_hor_post = calc_value(An_hor,n_post);
An_hor_mean = calc_mean_value(An_hor,n_pre,n_post);

% At_hor
At_hor_pre = calc_value(At_hor,n_pre);
At_hor_post = calc_value(At_hor,n_post);
At_hor_mean = calc_mean_value(At_hor,n_pre,n_post);

% A_hor norm
A_hor_pre = calc_value(A_hor,n_pre);
A_hor_post = calc_value(A_hor,n_post);
A_hor_mean = calc_mean_value(A_hor,n_pre,n_post);

% A norm
A_pre = calc_value(A,n_pre);
A_post = calc_value(A,n_post);

% Amax
[A_hor_max n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);
[A_max n_Amax] = calc_max_value(A,n_pre,n_post);

accel_angle_hor_vel_Ahormax = calc_value(accel_angle_hor_vel,n_Ahormax);
accel_angle_hor_vel_Amax = calc_value(accel_angle_hor_vel,n_Amax);

accel_angle_hor_body_Ahormax = calc_value(accel_angle_hor_body,n_Ahormax);
accel_angle_hor_body_Amax = calc_value(accel_angle_hor_body,n_Amax);

stim_angle_vel_Ahormax = calc_value(stim_angle_vel,n_Ahormax);
stim_angle_vel_Amax = calc_value(stim_angle_vel,n_Amax);

stim_angle_accel_Ahormax = calc_value(stim_angle_accel,n_Ahormax);
stim_angle_accel_Amax = calc_value(stim_angle_accel,n_Amax);

stim_angle_yaw_Ahormax = calc_value(stim_angle_yaw,n_Ahormax);
stim_angle_yaw_Amax = calc_value(stim_angle_yaw,n_Amax);

An_hor_Ahormax = calc_value(An_hor,n_Ahormax);
An_hor_Amax = calc_value(An_hor,n_Amax);

At_hor_Ahormax = calc_value(At_hor,n_Ahormax);
At_hor_Amax = calc_value(At_hor,n_Amax);

slip_Ahormax = calc_value(slip,n_Ahormax);
slip_Amax = calc_value(slip,n_Amax);

yaw_Ahormax = calc_value(yaw,n_Ahormax);
yaw_Amax = calc_value(yaw,n_Amax);

pitch_Ahormax = calc_value(pitch,n_Ahormax);
pitch_Amax = calc_value(pitch,n_Amax);

roll_Ahormax = calc_value(roll,n_Ahormax);
roll_Amax = calc_value(roll,n_Amax);

slip_aero_Ahormax = calc_value(slip_aero,n_Ahormax);
slip_aero_Amax = calc_value(slip_aero,n_Amax);

pitch_aero_Ahormax = calc_value(pitch_aero,n_Ahormax);
pitch_aero_Amax = calc_value(pitch_aero,n_Amax);

%% angles
% stim_angle_vel
heading_pre = calc_value(stim_angle_vel,n_pre);
heading_post = calc_value(stim_angle_vel,n_post);
stim_angle_vel_pre = calc_value(stim_angle_vel,n_pre);
stim_angle_vel_post = calc_value(stim_angle_vel,n_post);

% stim_angle_yaw
stim_angle_yaw_pre = calc_value(stim_angle_yaw,n_pre);
stim_angle_yaw_post = calc_value(stim_angle_yaw,n_post);

% stim F
stim_angle_F_pre = calc_value(stim_angle_F,n_pre);
stim_angle_F_post = calc_value(stim_angle_F,n_post);

% stim spn
stim_angle_spn_pre = calc_value(stim_angle_spn,n_pre);
stim_angle_spn_post = calc_value(stim_angle_spn,n_post);

% Adir_stimulus
Adir_pre = calc_value(stim_angle_accel,n_pre);
Adir_post = calc_value(stim_angle_accel,n_post);
[Adir_mean Adir_ul Adir_ll] = circ_calc_mean_value(stim_angle_accel,n_pre,n_post);

% Adir_body
Adir_body_pre = calc_value(accel_angle_hor_body,n_pre);
Adir_body_post = calc_value(accel_angle_hor_body,n_post);
[Adir_body_mean Adir_body_ul Adir_body_ll] = circ_calc_mean_value(accel_angle_hor_body,n_pre,n_post);

% Adir_vel
Adir_vel_pre = calc_value(accel_angle_hor_vel,n_pre);
Adir_vel_post = calc_value(accel_angle_hor_vel,n_post);
[Adir_vel_mean Adir_vel_ul Adir_vel_ll] = circ_calc_mean_value(accel_angle_hor_vel,n_pre,n_post);
A_mean = calc_mean_value(A,n_pre,n_post);

%% body angles
% slip
slip_pre = calc_value(slip,n_pre);
slip_post = calc_value(slip,n_post);
[slip_mean slip_ul slip_ll] = circ_calc_mean_value(slip,n_pre,n_post);

% yaw
yaw_pre = calc_value(yaw,n_pre);
yaw_post = calc_value(yaw,n_post);
[yaw_mean yaw_ul yaw_ll] = circ_calc_mean_value(yaw,n_pre,n_post);
yaw_max = calc_max_abs_value(yaw,n_pre,n_post);

% pitch
pitch_pre = calc_value(pitch,n_pre);
pitch_post = calc_value(pitch,n_post);
[pitch_mean pitch_ul pitch_ll] = circ_calc_mean_value(pitch,n_pre,n_post);
[pitch_steady_mean pitch_steady_ul pitch_steady_ll] = circ_calc_mean_value(pitch,n_first,n_pre);
pitch_max = calc_max_value(pitch,n_pre,n_post);

% roll
roll_pre = calc_value(roll,n_pre);
roll_post = calc_value(roll,n_post);
[roll_mean roll_ul roll_ll] = circ_calc_mean_value(roll,n_pre,n_post);
roll_max = calc_max_abs_value(roll,n_pre,n_post);

%% global body angles
roll_global_first = calc_value(roll_global,n_first);
slip_global_first = calc_value(slip_global,n_first);
pitch_global_first = calc_value(pitch_global,n_first);

roll_global_pre = calc_value(roll_global,n_pre);
slip_global_pre = calc_value(slip_global,n_pre);
pitch_global_pre = calc_value(pitch_global,n_pre);

[roll_global_steady_mean roll_global_steady_ul roll_global_steady_ll] = circ_calc_mean_value(roll_global,n_first,n_pre);
[slip_global_steady_mean slip_global_steady_ul slip_global_steady_ll] = circ_calc_mean_value(slip_global,n_first,n_pre);
[pitch_global_steady_mean pitch_global_steady_ul pitch_global_steady_ll] = circ_calc_mean_value(pitch_global,n_first,n_pre);

%% body angles DOT
yaw_dot_max = calc_max_value(abs(yaw_dot),n_pre,n_post);
pitch_dot_max = calc_max_value(pitch_dot,n_pre,n_post);
% roll_dot_max = calc_max_value(roll_dot,n_pre,n_post);

yaw_dot_mean = calc_mean_value(yaw_dot,n_pre,n_post);
pitch_dot_mean = calc_mean_value(pitch_dot,n_pre,n_post);
% roll_dot_mean = calc_mean_value(roll_dot,n_pre,n_post);


%% AERO body angles
% slip_aero
slip_aero_pre = calc_value(slip_aero,n_pre);
slip_aero_post = calc_value(slip_aero,n_post);
[slip_aero_mean slip_aero_ul slip_aero_ll] = circ_calc_mean_value(slip_aero,n_pre,n_post);

% pitch_aero
pitch_aero_pre = calc_value(pitch_aero,n_pre);
pitch_aero_post = calc_value(pitch_aero,n_post);
[pitch_aero_mean pitch_aero_ul pitch_aero_ll] = circ_calc_mean_value(pitch_aero,n_pre,n_post);
[pitch_aero_steady_mean pitch_aero_steady_ul pitch_aero_steady_ll] = circ_calc_mean_value(pitch_aero,n_first,n_pre);

%% forces
% F_hor norm
F_hor_pre = calc_value(F_hor,n_pre);
F_hor_post = calc_value(F_hor,n_post);
F_hor_mean = calc_mean_value(F_hor,n_pre,n_post);

% F norm
F_pre = calc_value(F,n_pre);
F_post = calc_value(F,n_post);
F_mean = calc_mean_value(F,n_pre,n_post);

% Fsp_pitch
Fsp_pitch_pre = calc_value(Fsp_pitch,n_pre);
Fsp_pitch_post = calc_value(Fsp_pitch,n_post);

% Fsp_roll
Fsp_roll_pre = calc_value(Fsp_roll,n_pre);
Fsp_roll_post = calc_value(Fsp_roll,n_post);

% Fb_pitch
Fb_pitch_pre = calc_value(Fb_pitch,n_pre);
Fb_pitch_post = calc_value(Fb_pitch,n_post);

% Fb_roll
Fb_roll_pre = calc_value(Fb_roll,n_pre);
Fb_roll_post = calc_value(Fb_roll,n_post);

% F
F_pre = calc_value(F,n_pre);
F_post = calc_value(F,n_post);
F_mean = calc_mean_value(F,n_pre,n_post);

% Fmax & Fhormax
[F_hor_max n_Fhormax] = calc_max_value(F_hor,n_pre,n_post);
[F_max n_Fmax] = calc_max_value(F,n_pre,n_post);

F_angle_hor_vel_Fhormax = calc_value(F_angle_hor_vel,n_Fhormax);
F_angle_hor_vel_Fmax = calc_value(F_angle_hor_vel,n_Fmax);

F_angle_hor_body_Fhormax = calc_value(F_angle_hor_body,n_Fhormax);
F_angle_hor_body_Fmax = calc_value(F_angle_hor_body,n_Fmax);

accel_angle_hor_vel_Fhormax = calc_value(accel_angle_hor_vel,n_Fhormax);
accel_angle_hor_vel_Fmax = calc_value(accel_angle_hor_vel,n_Fmax);

accel_angle_hor_body_Fhormax = calc_value(accel_angle_hor_body,n_Fhormax);
accel_angle_hor_body_Fmax = calc_value(accel_angle_hor_body,n_Fmax);

stim_angle_vel_Fhormax = calc_value(stim_angle_vel,n_Fhormax);
stim_angle_vel_Fmax = calc_value(stim_angle_vel,n_Fmax);

stim_angle_accel_Fhormax = calc_value(stim_angle_accel,n_Fhormax);
stim_angle_accel_Fmax = calc_value(stim_angle_accel,n_Fmax);

stim_angle_yaw_Fhormax = calc_value(stim_angle_yaw,n_Fhormax);
stim_angle_yaw_Fmax = calc_value(stim_angle_yaw,n_Fmax);

stim_angle_F_Fhormax = calc_value(stim_angle_F,n_Fhormax);
stim_angle_F_Fmax = calc_value(stim_angle_F,n_Fmax);

stim_angle_spn_Fhormax = calc_value(stim_angle_spn,n_Fhormax);
stim_angle_spn_Fmax = calc_value(stim_angle_spn,n_Fmax);

An_hor_Fhormax = calc_value(An_hor,n_Fhormax);
An_hor_Fmax = calc_value(An_hor,n_Fmax);

At_hor_Fhormax = calc_value(At_hor,n_Fhormax);
At_hor_Fmax = calc_value(At_hor,n_Fmax);

Fn_hor_Fhormax = calc_value(Fn_hor,n_Fhormax);
Fn_hor_Fmax = calc_value(Fn_hor,n_Fmax);

Ft_hor_Fhormax = calc_value(Ft_hor,n_Fhormax);
Ft_hor_Fmax = calc_value(Ft_hor,n_Fmax);

slip_Fhormax = calc_value(slip,n_Fhormax);
slip_Fmax = calc_value(slip,n_Fmax);

yaw_Fhormax = calc_value(yaw,n_Fhormax);
yaw_Fmax = calc_value(yaw,n_Fmax);

pitch_Fhormax = calc_value(pitch,n_Fhormax);
pitch_Fmax = calc_value(pitch,n_Fmax);

roll_Fhormax = calc_value(roll,n_Fhormax);
roll_Fmax = calc_value(roll,n_Fmax);

Fsp_pitch_Fhormax = calc_value(Fsp_pitch,n_Fhormax);
Fsp_pitch_Fmax = calc_value(Fsp_pitch,n_Fmax);

Fsp_roll_Fhormax = calc_value(Fsp_roll,n_Fhormax);
Fsp_roll_Fmax = calc_value(Fsp_roll,n_Fmax);

Fb_pitch_Fhormax = calc_value(Fb_pitch,n_Fhormax);
Fb_pitch_Fmax = calc_value(Fb_pitch,n_Fmax);

Fb_roll_Fhormax = calc_value(Fb_roll,n_Fhormax);
Fb_roll_Fmax = calc_value(Fb_roll,n_Fmax);


%% pre response n_first
F_hor_first = calc_value(F_hor,n_first);
F_first = calc_value(F,n_first);

F_angle_hor_vel_first = calc_value(F_angle_hor_vel,n_first);
F_angle_hor_body_first = calc_value(F_angle_hor_body,n_first);

accel_angle_hor_vel_first = calc_value(accel_angle_hor_vel,n_first);
accel_angle_hor_body_first = calc_value(accel_angle_hor_body,n_first);

stim_angle_vel_first = calc_value(stim_angle_vel,n_first);
stim_angle_accel_first = calc_value(stim_angle_accel,n_first);
stim_angle_yaw_first = calc_value(stim_angle_yaw,n_first);
stim_angle_F_first = calc_value(stim_angle_F,n_first);
stim_angle_spn_first = calc_value(stim_angle_spn,n_first);

An_hor_first = calc_value(An_hor,n_first);
At_hor_first = calc_value(At_hor,n_first);

Fn_hor_first = calc_value(Fn_hor,n_first);
Ft_hor_first = calc_value(Ft_hor,n_first);

slip_first = calc_value(slip,n_first);
yaw_first = calc_value(yaw,n_first);
pitch_first = calc_value(pitch,n_first);
roll_first = calc_value(roll,n_first);

Fsp_pitch_first = calc_value(Fsp_pitch,n_first);
Fsp_roll_first = calc_value(Fsp_roll,n_first);
Fb_pitch_first = calc_value(Fb_pitch,n_first);
Fb_roll_first = calc_value(Fb_roll,n_first);


