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

pitch_Ahormax = calc_value(pitch,n_Ahormax);
pitch_Amax = calc_value(pitch,n_Amax);

roll_Ahormax = calc_value(roll,n_Ahormax);
roll_Amax = calc_value(roll,n_Amax);



%% angles
% heading
heading_pre = calc_value(stim_angle_vel,n_pre);
heading_post = calc_value(stim_angle_vel,n_post);

% yaw
yaw_pre = calc_value(stim_angle_yaw,n_pre);
yaw_post = calc_value(stim_angle_yaw,n_post);

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
    
% pitch
pitch_pre = calc_value(pitch,n_pre);
pitch_post = calc_value(pitch,n_post);

% roll
roll_pre = calc_value(roll,n_pre);
roll_post = calc_value(roll,n_post);

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

% F norm
V_pre = calc_value(V,n_pre);
V_post = calc_value(V,n_post);
V_mean = calc_mean_value(V,n_pre,n_post);

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
pitch_first = calc_value(pitch,n_first);
roll_first = calc_value(roll,n_first);

Fsp_pitch_first = calc_value(Fsp_pitch,n_first);
Fsp_roll_first = calc_value(Fsp_roll,n_first);
Fb_pitch_first = calc_value(Fb_pitch,n_first);
Fb_roll_first = calc_value(Fb_roll,n_first);


