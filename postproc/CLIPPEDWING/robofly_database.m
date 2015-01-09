

%% ROLL
clear
clc
close all

load('MOD_norm_data_buttercut7.mat')
load(rollDB_all)
load('MOD_norm_data_buttercut7.mat')

skip = 10;

t = t_allNOfreq(1:skip:end,1);
t = t-min(t);
t = t/f_robo2fly;

%% kinematics
stroke_L = stroke_L_allNOfreq(1:skip:end,:);
stroke_R = stroke_R_allNOfreq(1:skip:end,:);

deviation_L = dev_L_allNOfreq(1:skip:end,:);
deviation_R = dev_R_allNOfreq(1:skip:end,:);

rotation_L = pitch_L_allNOfreq(1:skip:end,:);
rotation_R = pitch_R_allNOfreq(1:skip:end,:);

%% forces & torques

% Fx
Fx_all(:,:) = mean(Fx_allNOfreq,2);
Fx_all = Fx_all(1:skip:end,1:end-1);

Fx_str(:,:) = mean(Fx_stroke,2);
Fx_str = Fx_str(1:skip:end,1:end-1);

Fx_rotation(:,:) = mean(Fx_pitch,2);
Fx_rotation = Fx_rotation(1:skip:end,1:end-1);

Fx_deviation(:,:) = mean(Fx_dev,2);
Fx_deviation = Fx_deviation(1:skip:end,1:end-1);

% Fy
Fy_all(:,:) = mean(Fy_allNOfreq,2);
Fy_all = Fy_all(1:skip:end,1:end-1);

Fy_str(:,:) = mean(Fy_stroke,2);
Fy_str = Fy_str(1:skip:end,1:end-1);

Fy_rotation(:,:) = mean(Fy_pitch,2);
Fy_rotation = Fy_rotation(1:skip:end,1:end-1);

Fy_deviation(:,:) = mean(Fy_dev,2);
Fy_deviation = Fy_deviation(1:skip:end,1:end-1);

% Fz
Fz_all(:,:) = mean(Fz_allNOfreq,2);
Fz_all = Fz_all(1:skip:end,1:end-1);

Fz_str(:,:) = mean(Fz_stroke,2);
Fz_str = Fz_str(1:skip:end,1:end-1);

Fz_rotation(:,:) = mean(Fz_pitch,2);
Fz_rotation = Fz_rotation(1:skip:end,1:end-1);

Fz_deviation(:,:) = mean(Fz_dev,2);
Fz_deviation = Fz_deviation(1:skip:end,1:end-1);

% Mx
Mx_all(:,:) = mean(Mx_allNOfreq,2);
Mx_all = Mx_all(1:skip:end,1:end-1);

Mx_str(:,:) = mean(Mx_stroke,2);
Mx_str = Mx_str(1:skip:end,1:end-1);

Mx_rotation(:,:) = mean(Mx_pitch,2);
Mx_rotation = Mx_rotation(1:skip:end,1:end-1);

Mx_deviation(:,:) = mean(Mx_dev,2);
Mx_deviation = Mx_deviation(1:skip:end,1:end-1);

% My
My_all(:,:) = mean(My_allNOfreq,2);
My_all = My_all(1:skip:end,1:end-1);

My_str(:,:) = mean(My_stroke,2);
My_str = My_str(1:skip:end,1:end-1);

My_rotation(:,:) = mean(My_pitch,2);
My_rotation = My_rotation(1:skip:end,1:end-1);

My_deviation(:,:) = mean(My_dev,2);
My_deviation = My_deviation(1:skip:end,1:end-1);

% Mz
Mz_all(:,:) = mean(Mz_allNOfreq,2);
Mz_all = Mz_all(1:skip:end,1:end-1);

Mz_str(:,:) = mean(Mz_stroke,2);
Mz_str = Mz_str(1:skip:end,1:end-1);

Mz_rotation(:,:) = mean(Mz_pitch,2);
Mz_rotation = Mz_rotation(1:skip:end,1:end-1);

Mz_deviation(:,:) = mean(Mz_dev,2);
Mz_deviation = Mz_deviation(1:skip:end,1:end-1);

%% store data

roll.t = t(isnan(t)==0,:);
roll.roll_accel_norm = mod_value_allNOfreq*rollaccel_norm;

roll.stroke_L = stroke_L(isnan(t)==0,:);
roll.stroke_R = stroke_R(isnan(t)==0,:);

roll.deviation_L = deviation_L(isnan(t)==0,:);
roll.deviation_R = deviation_R(isnan(t)==0,:);

roll.rotation_L = rotation_L(isnan(t)==0,:);
roll.rotation_R = rotation_R(isnan(t)==0,:);

roll.Fx_norm_all = Fx_all(isnan(t)==0,:)/Mg_fly;
roll.Fx_norm_stroke = Fx_str(isnan(t)==0,:)/Mg_fly;
roll.Fx_norm_rotation = Fx_rotation(isnan(t)==0,:)/Mg_fly;
roll.Fx_norm_deviation = Fx_deviation(isnan(t)==0,:)/Mg_fly;

roll.Fy_norm_all = Fy_all(isnan(t)==0,:)/Mg_fly;
roll.Fy_norm_stroke = Fy_str(isnan(t)==0,:)/Mg_fly;
roll.Fy_norm_rotation = Fy_rotation(isnan(t)==0,:)/Mg_fly;
roll.Fy_norm_deviation = Fy_deviation(isnan(t)==0,:)/Mg_fly;

roll.Fz_norm_all = Fz_all(isnan(t)==0,:)/Mg_fly;
roll.Fz_norm_stroke = Fz_str(isnan(t)==0,:)/Mg_fly;
roll.Fz_norm_rotation = Fz_rotation(isnan(t)==0,:)/Mg_fly;
roll.Fz_norm_deviation = Fz_deviation(isnan(t)==0,:)/Mg_fly;

roll.Mx_norm_all = Mx_all(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mx_norm_stroke = Mx_str(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mx_norm_rotation = Mx_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mx_norm_deviation = Mx_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

roll.My_norm_all = My_all(isnan(t)==0,:)/Mg_fly/Lwing;
roll.My_norm_stroke = My_str(isnan(t)==0,:)/Mg_fly/Lwing;
roll.My_norm_rotation = My_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
roll.My_norm_deviation = My_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

roll.Mz_norm_all = Mz_all(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mz_norm_stroke = Mz_str(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mz_norm_rotation = Mz_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
roll.Mz_norm_deviation = Mz_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

%% save
save('DBkinematics_forces_torques.mat','roll','skip')

%% PITCH
clear
clc
close all

load('MOD_norm_data_buttercut7.mat')
load(pitchDB_all)
load('MOD_norm_data_buttercut7.mat')
load('DBkinematics_forces_torques.mat')

t = t_allNOfreq(1:skip:end,1);
t = t-min(t);
t = t/f_robo2fly;

%% kinematics
stroke = stroke_L_allNOfreq(1:skip:end,:);
deviation = dev_L_allNOfreq(1:skip:end,:);
rotation = pitch_L_allNOfreq(1:skip:end,:);

%% forces & torques

% Fx
Fx_all(:,:) = mean(Fx_allNOfreq,2);
Fx_all = Fx_all(1:skip:end,1:end-1);

Fx_str(:,:) = mean(Fx_stroke,2);
Fx_str = Fx_str(1:skip:end,1:end-1);

Fx_rotation(:,:) = mean(Fx_pitch,2);
Fx_rotation = Fx_rotation(1:skip:end,1:end-1);

Fx_deviation(:,:) = mean(Fx_dev,2);
Fx_deviation = Fx_deviation(1:skip:end,1:end-1);

% Fy
Fy_all(:,:) = mean(Fy_allNOfreq,2);
Fy_all = Fy_all(1:skip:end,1:end-1);

Fy_str(:,:) = mean(Fy_stroke,2);
Fy_str = Fy_str(1:skip:end,1:end-1);

Fy_rotation(:,:) = mean(Fy_pitch,2);
Fy_rotation = Fy_rotation(1:skip:end,1:end-1);

Fy_deviation(:,:) = mean(Fy_dev,2);
Fy_deviation = Fy_deviation(1:skip:end,1:end-1);

% Fz
Fz_all(:,:) = mean(Fz_allNOfreq,2);
Fz_all = Fz_all(1:skip:end,1:end-1);

Fz_str(:,:) = mean(Fz_stroke,2);
Fz_str = Fz_str(1:skip:end,1:end-1);

Fz_rotation(:,:) = mean(Fz_pitch,2);
Fz_rotation = Fz_rotation(1:skip:end,1:end-1);

Fz_deviation(:,:) = mean(Fz_dev,2);
Fz_deviation = Fz_deviation(1:skip:end,1:end-1);

% Mx
Mx_all(:,:) = mean(Mx_allNOfreq,2);
Mx_all = Mx_all(1:skip:end,1:end-1);

Mx_str(:,:) = mean(Mx_stroke,2);
Mx_str = Mx_str(1:skip:end,1:end-1);

Mx_rotation(:,:) = mean(Mx_pitch,2);
Mx_rotation = Mx_rotation(1:skip:end,1:end-1);

Mx_deviation(:,:) = mean(Mx_dev,2);
Mx_deviation = Mx_deviation(1:skip:end,1:end-1);

% My
My_all(:,:) = mean(My_allNOfreq,2);
My_all = My_all(1:skip:end,1:end-1);

My_str(:,:) = mean(My_stroke,2);
My_str = My_str(1:skip:end,1:end-1);

My_rotation(:,:) = mean(My_pitch,2);
My_rotation = My_rotation(1:skip:end,1:end-1);

My_deviation(:,:) = mean(My_dev,2);
My_deviation = My_deviation(1:skip:end,1:end-1);

% Mz
Mz_all(:,:) = mean(Mz_allNOfreq,2);
Mz_all = Mz_all(1:skip:end,1:end-1);

Mz_str(:,:) = mean(Mz_stroke,2);
Mz_str = Mz_str(1:skip:end,1:end-1);

Mz_rotation(:,:) = mean(Mz_pitch,2);
Mz_rotation = Mz_rotation(1:skip:end,1:end-1);

Mz_deviation(:,:) = mean(Mz_dev,2);
Mz_deviation = Mz_deviation(1:skip:end,1:end-1);

%% store data

pitch.t = t(isnan(t)==0,:);
pitch.pitch_accel_norm = mod_value_allNOfreq*pitchaccel_norm;

pitch.stroke = stroke(isnan(t)==0,:);
pitch.deviation = deviation(isnan(t)==0,:);
pitch.rotation = rotation(isnan(t)==0,:);

pitch.Fx_norm_all = Fx_all(isnan(t)==0,:)/Mg_fly;
pitch.Fx_norm_stroke = Fx_str(isnan(t)==0,:)/Mg_fly;
pitch.Fx_norm_rotation = Fx_rotation(isnan(t)==0,:)/Mg_fly;
pitch.Fx_norm_deviation = Fx_deviation(isnan(t)==0,:)/Mg_fly;

pitch.Fy_norm_all = Fy_all(isnan(t)==0,:)/Mg_fly;
pitch.Fy_norm_stroke = Fy_str(isnan(t)==0,:)/Mg_fly;
pitch.Fy_norm_rotation = Fy_rotation(isnan(t)==0,:)/Mg_fly;
pitch.Fy_norm_deviation = Fy_deviation(isnan(t)==0,:)/Mg_fly;

pitch.Fz_norm_all = Fz_all(isnan(t)==0,:)/Mg_fly;
pitch.Fz_norm_stroke = Fz_str(isnan(t)==0,:)/Mg_fly;
pitch.Fz_norm_rotation = Fz_rotation(isnan(t)==0,:)/Mg_fly;
pitch.Fz_norm_deviation = Fz_deviation(isnan(t)==0,:)/Mg_fly;

pitch.Mx_norm_all = Mx_all(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mx_norm_stroke = Mx_str(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mx_norm_rotation = Mx_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mx_norm_deviation = Mx_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

pitch.My_norm_all = My_all(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.My_norm_stroke = My_str(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.My_norm_rotation = My_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.My_norm_deviation = My_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

pitch.Mz_norm_all = Mz_all(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mz_norm_stroke = Mz_str(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mz_norm_rotation = Mz_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
pitch.Mz_norm_deviation = Mz_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

%% save
save('DBkinematics_forces_torques.mat','roll','pitch','skip')

%% FORCE
clear
clc
close all

load('MOD_norm_data_buttercut7.mat')
load(FDB_all)
load('MOD_norm_data_buttercut7.mat')
load('DBkinematics_forces_torques.mat')

t_INCfreq = t_freq(1:skip:end,:);
for i=1:size(t_freq,2)
    t_INCfreq(:,i) = t_INCfreq(:,i) - min(t_INCfreq(:,i));
end
t_INCfreq = t_INCfreq/f_robo2fly;

t_NOfreq = t_allNOfreq(1:skip:end,1);
t_NOfreq = t_NOfreq-min(t_NOfreq);
t_NOfreq = t_NOfreq/f_robo2fly;

%% kinematics
stroke = stroke_L_allNOfreq(1:skip:end,:);
deviation = dev_L_allNOfreq(1:skip:end,:);
rotation = pitch_L_allNOfreq(1:skip:end,:);

%% forces & torques

% Fx
Fx_allvars(:,:) = mean(Fx_all,2);
Fx_allvars = Fx_allvars(1:skip:end,1:end-1);
Fx_allvars(Fx_allvars==0)=nan;

Fx_f(:,:) = mean(Fx_freq,2);
Fx_f = Fx_f(1:skip:end,1:end-1);
Fx_f(Fx_f==0)=nan;

Fx_str(:,:) = mean(Fx_stroke,2);
Fx_str = Fx_str(1:skip:end,1:end-1);
Fx_str(Fx_str==0)=nan;

Fx_rotation(:,:) = mean(Fx_pitch,2);
Fx_rotation = Fx_rotation(1:skip:end,1:end-1);
Fx_rotation(Fx_rotation==0)=nan;

Fx_deviation(:,:) = mean(Fx_dev,2);
Fx_deviation = Fx_deviation(1:skip:end,1:end-1);
Fx_deviation(Fx_deviation==0)=nan;

% Fy
Fy_allvars(:,:) = mean(Fy_all,2);
Fy_allvars = Fy_allvars(1:skip:end,1:end-1);
Fy_allvars(Fy_allvars==0)=nan;

Fy_f(:,:) = mean(Fy_freq,2);
Fy_f = Fy_f(1:skip:end,1:end-1);
Fy_f(Fy_f==0)=nan;

Fy_str(:,:) = mean(Fy_stroke,2);
Fy_str = Fy_str(1:skip:end,1:end-1);
Fy_str(Fy_str==0)=nan;

Fy_rotation(:,:) = mean(Fy_pitch,2);
Fy_rotation = Fy_rotation(1:skip:end,1:end-1);
Fy_rotation(Fy_rotation==0)=nan;

Fy_deviation(:,:) = mean(Fy_dev,2);
Fy_deviation = Fy_deviation(1:skip:end,1:end-1);
Fy_deviation(Fy_deviation==0)=nan;

% Fz
Fz_allvars(:,:) = mean(Fz_all,2);
Fz_allvars = Fz_allvars(1:skip:end,1:end-1);
Fz_allvars(Fz_allvars==0)=nan;

Fz_f(:,:) = mean(Fz_freq,2);
Fz_f = Fz_f(1:skip:end,1:end-1);
Fz_f(Fz_f==0)=nan;

Fz_str(:,:) = mean(Fz_stroke,2);
Fz_str = Fz_str(1:skip:end,1:end-1);
Fz_str(Fz_str==0)=nan;

Fz_rotation(:,:) = mean(Fz_pitch,2);
Fz_rotation = Fz_rotation(1:skip:end,1:end-1);
Fz_rotation(Fz_rotation==0)=nan;

Fz_deviation(:,:) = mean(Fz_dev,2);
Fz_deviation = Fz_deviation(1:skip:end,1:end-1);
Fz_deviation(Fz_deviation==0)=nan;

% Mx
Mx_allvars(:,:) = mean(Mx_all,2);
Mx_allvars = Mx_allvars(1:skip:end,1:end-1);
Mx_allvars(Mx_allvars==0)=nan;

Mx_f(:,:) = mean(Mx_freq,2);
Mx_f = Mx_f(1:skip:end,1:end-1);
Mx_f(Mx_f==0)=nan;

Mx_str(:,:) = mean(Mx_stroke,2);
Mx_str = Mx_str(1:skip:end,1:end-1);
Mx_str(Mx_str==0)=nan;

Mx_rotation(:,:) = mean(Mx_pitch,2);
Mx_rotation = Mx_rotation(1:skip:end,1:end-1);
Mx_rotation(Mx_rotation==0)=nan;

Mx_deviation(:,:) = mean(Mx_dev,2);
Mx_deviation = Mx_deviation(1:skip:end,1:end-1);
Mx_deviation(Mx_deviation==0)=nan;

% My
My_allvars(:,:) = mean(My_all,2);
My_allvars = My_allvars(1:skip:end,1:end-1);
My_allvars(My_allvars==0)=nan;

My_f(:,:) = mean(My_freq,2);
My_f = My_f(1:skip:end,1:end-1);
My_f(My_f==0)=nan;

My_str(:,:) = mean(My_stroke,2);
My_str = My_str(1:skip:end,1:end-1);
My_str(My_str==0)=nan;

My_rotation(:,:) = mean(My_pitch,2);
My_rotation = My_rotation(1:skip:end,1:end-1);
My_rotation(My_rotation==0)=nan;

My_deviation(:,:) = mean(My_dev,2);
My_deviation = My_deviation(1:skip:end,1:end-1);
My_deviation(My_deviation==0)=nan;

% Mz
Mz_allvars(:,:) = mean(Mz_all,2);
Mz_allvars = Mz_allvars(1:skip:end,1:end-1);
Mz_allvars(Mz_allvars==0)=nan;

Mz_f(:,:) = mean(Mz_freq,2);
Mz_f = Mz_f(1:skip:end,1:end-1);
Mz_f(Mz_f==0)=nan;

Mz_str(:,:) = mean(Mz_stroke,2);
Mz_str = Mz_str(1:skip:end,1:end-1);
Mz_str(Mz_str==0)=nan;

Mz_rotation(:,:) = mean(Mz_pitch,2);
Mz_rotation = Mz_rotation(1:skip:end,1:end-1);
Mz_rotation(Mz_rotation==0)=nan;

Mz_deviation(:,:) = mean(Mz_dev,2);
Mz_deviation = Mz_deviation(1:skip:end,1:end-1);
Mz_deviation(Mz_deviation==0)=nan;

%% store data
t = t_NOfreq;

force.t_NOfreq = t_NOfreq(isnan(t_NOfreq)==0,:);
force.t_INCfreq = t_INCfreq;
force.force_norm = mod_value_all*Fenhance_norm+1;

force.stroke = stroke(isnan(t)==0,:);
force.deviation = deviation(isnan(t)==0,:);
force.rotation = rotation(isnan(t)==0,:);

force.Fx_norm_all = Fx_allvars/Mg_fly;
force.Fx_norm_freq = Fx_f/Mg_fly;
force.Fx_norm_stroke = Fx_str(isnan(t)==0,:)/Mg_fly;
force.Fx_norm_rotation = Fx_rotation(isnan(t)==0,:)/Mg_fly;
force.Fx_norm_deviation = Fx_deviation(isnan(t)==0,:)/Mg_fly;

force.Fy_norm_all = Fy_allvars/Mg_fly;
force.Fy_norm_freq = Fy_f/Mg_fly;
force.Fy_norm_stroke = Fy_str(isnan(t)==0,:)/Mg_fly;
force.Fy_norm_rotation = Fy_rotation(isnan(t)==0,:)/Mg_fly;
force.Fy_norm_deviation = Fy_deviation(isnan(t)==0,:)/Mg_fly;

force.Fz_norm_all = Fz_allvars/Mg_fly;
force.Fz_norm_freq = Fz_f/Mg_fly;
force.Fz_norm_stroke = Fz_str(isnan(t)==0,:)/Mg_fly;
force.Fz_norm_rotation = Fz_rotation(isnan(t)==0,:)/Mg_fly;
force.Fz_norm_deviation = Fz_deviation(isnan(t)==0,:)/Mg_fly;

force.Mx_norm_all = Mx_allvars/Mg_fly/Lwing;
force.Mx_norm_freq = Mx_f/Mg_fly/Lwing;
force.Mx_norm_stroke = Mx_str(isnan(t)==0,:)/Mg_fly/Lwing;
force.Mx_norm_rotation = Mx_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
force.Mx_norm_deviation = Mx_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

force.My_norm_all = My_allvars/Mg_fly/Lwing;
force.My_norm_freq = My_f/Mg_fly/Lwing;
force.My_norm_stroke = My_str(isnan(t)==0,:)/Mg_fly/Lwing;
force.My_norm_rotation = My_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
force.My_norm_deviation = My_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

force.Mz_norm_all = Mz_allvars/Mg_fly/Lwing;
force.Mz_norm_freq = Mz_f/Mg_fly/Lwing;
force.Mz_norm_stroke = Mz_str(isnan(t)==0,:)/Mg_fly/Lwing;
force.Mz_norm_rotation = Mz_rotation(isnan(t)==0,:)/Mg_fly/Lwing;
force.Mz_norm_deviation = Mz_deviation(isnan(t)==0,:)/Mg_fly/Lwing;

%% save
save('DBkinematics_forces_torques.mat','roll','pitch','force')

