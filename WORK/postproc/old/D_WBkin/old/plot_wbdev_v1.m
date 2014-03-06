clc
clear
close all

% load body kin
loadname = 'kinflightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3.mat';
load(loadname)

ddevLminR = kinDB.devL - kinDB.devR;
roll_dot_sp = pathDB.roll_dot_sp;

%% MIRROR ATTITUDES based on roll_max
roll_dot_pre_mirror = calc_value(roll_dot,n_pre);

slip_global_mirror = slip_global;
roll_global_mirror = roll_global;
yaw_mirror = yaw;
roll_mirror = roll;
dyaw_mirror = dyaw;
droll_mirror = droll;
yaw_dot_mirror = yaw_dot;
roll_dot_mirror = roll_dot;
Fsp_roll_mirror = Fsp_roll;

for i=1:length(roll_dot_pre_mirror)
    if roll_dot_pre_mirror(i) < 0
        slip_global_mirror(:,i) = -slip_global_mirror(:,i);
        roll_global_mirror(:,i) = -roll_global_mirror(:,i);
        yaw_mirror(:,i) = -yaw_mirror(:,i);
        roll_mirror(:,i) = -roll_mirror(:,i);
        dyaw_mirror(:,i) = -dyaw_mirror(:,i);
        droll_mirror(:,i) = -droll_mirror(:,i);
        yaw_dot_mirror(:,i) = -yaw_dot_mirror(:,i);
        roll_dot_mirror(:,i) = -roll_dot_mirror(:,i);
        Fsp_roll_mirror(:,i) = -Fsp_roll_mirror(:,i);
    end
end


for 