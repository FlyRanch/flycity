function [FM_trans,FM_trans_norm,FM_transNrot,FM_trans_normNrot,...
    FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_trans_normNrot,...
    FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_trans_normNrot] =...
    calc_QSforcesNmoments_steadyNmodNbodyKin(const,wing_model,body_model,freq,...
    stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
    stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
    stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R)

%% constants
Mfly = const.Mfly;
Mg_fly = const.Mg_fly;
stroke_angle = const.stroke_angle;

R_strk = body_model.R_strk;
u_strk = body_model.u_strk;
w_strk = body_model.w_strk;

wing_length = wing_model.length;
nr_sect = wing_model.nr_sect;
nr_timepoints = wing_model.nr_timepoints;

rot_lift_on = 1;
rot_lift_off = 0;

%% calc steady WB kin
[stroke_L_steady,stroke_R_steady,stroke_dot_L_steady,stroke_dot_R_steady,stroke_dot_dot_L_steady,stroke_dot_dot_R_steady,...
    dev_L_steady,dev_R_steady,dev_dot_L_steady,dev_dot_R_steady,dev_dot_dot_L_steady,dev_dot_dot_R_steady,...
    rot_L_steady,rot_R_steady,rot_dot_L_steady,rot_dot_R_steady,rot_dot_dot_L_steady,rot_dot_dot_R_steady,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs_symm(nr_timepoints,freq,stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady);

%% calc WBmod kin
[stroke_L_MOD,stroke_R_MOD,stroke_dot_L_MOD,stroke_dot_R_MOD,stroke_dot_dot_L_MOD,stroke_dot_dot_R_MOD,...
    dev_L_MOD,dev_R_MOD,dev_dot_L_MOD,dev_dot_R_MOD,dev_dot_dot_L_MOD,dev_dot_dot_R_MOD,...
    rot_L_MOD,rot_R_MOD,rot_dot_L_MOD,rot_dot_R_MOD,rot_dot_dot_L_MOD,rot_dot_dot_R_MOD,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs(nr_timepoints,freq,stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);

%% sum steady & MOD

stroke_L = stroke_L_steady + stroke_L_MOD;
stroke_dot_L = stroke_dot_L_steady + stroke_dot_L_MOD;
stroke_dot_dot_L = stroke_dot_dot_L_steady + stroke_dot_dot_L_MOD;

stroke_R = stroke_R_steady + stroke_R_MOD;
stroke_dot_R = stroke_dot_R_steady + stroke_dot_R_MOD;
stroke_dot_dot_R = stroke_dot_dot_R_steady + stroke_dot_dot_R_MOD;

dev_L = dev_L_steady + dev_L_MOD;
dev_dot_L = dev_dot_L_steady + dev_dot_L_MOD;
dev_dot_dot_L = dev_dot_dot_L_steady + dev_dot_dot_L_MOD;

dev_R = dev_R_steady + dev_R_MOD;
dev_dot_R = dev_dot_R_steady + dev_dot_R_MOD;
dev_dot_dot_R = dev_dot_dot_R_steady + dev_dot_dot_R_MOD;

rot_L = rot_L_steady + rot_L_MOD;
rot_dot_L = rot_dot_L_steady + rot_dot_L_MOD;
rot_dot_dot_L = rot_dot_dot_L_steady + rot_dot_dot_L_MOD;

rot_R = rot_R_steady + rot_R_MOD;
rot_dot_R = rot_dot_R_steady + rot_dot_R_MOD;
rot_dot_dot_R = rot_dot_dot_R_steady + rot_dot_dot_R_MOD;

%% calc WB kin velocities
wb_loc(1) = 1;
wb_loc(2) = length(stroke_L);
down_loc(1) = 1;
down_loc(2) = find(stroke_L == max(stroke_L));
up_loc(1) = find(stroke_L == max(stroke_L))+1;
up_loc(2) = length(stroke_L);

[kine] = calc_angular_velocities_INCdots(stroke_L,dev_L,rot_L,stroke_R,dev_R,rot_R,stroke_dot_L,dev_dot_L,rot_dot_L,stroke_dot_R,dev_dot_R,rot_dot_R,stroke_dot_dot_L,dev_dot_dot_L,rot_dot_dot_L,stroke_dot_dot_R,dev_dot_dot_R,rot_dot_dot_R,freq,R_strk);

kine.R_strk = R_strk;
kine.u_strk = u_strk;
kine.w_strk = w_strk;

wb.wb_loc          = wb_loc;
wb.down_loc        = down_loc;
wb.up_loc          = up_loc;
wb.dt = dt;

%% NO rotational lift
[ FM_trans, FM_L_trans, FM_R_trans ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
    Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_off);

FM_trans_norm(1:3,:) = FM_trans(1:3,:)./Mg_fly;
FM_trans_norm(4:6,:) = FM_trans(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

FM_L_trans_norm(1:3,:) = FM_L_trans(1:3,:)./Mg_fly;
FM_L_trans_norm(4:6,:) = FM_L_trans(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

FM_R_transNrot_norm(1:3,:) = FM_R_transNrot(1:3,:)./Mg_fly;
FM_R_transNrot_norm(4:6,:) = FM_R_transNrot(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

% Fx_trans = FM_trans(1,:)./Mg_fly;
% Fy_trans = FM_trans(2,:)./Mg_fly;
% Fz_trans = FM_trans(3,:)./Mg_fly;
% 
% Mx_trans = FM_trans(4,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% My_trans = FM_trans(5,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% Mz_trans = FM_trans(6,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% 
% F_trans = sqrt(Fx_trans.^2 + Fy_trans.^2 + Fz_trans.^2);


%% INC rotational lift
[ FM_transNrot, FM_L_transNrot, FM_R_transNrot ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_on);

FM_transNrot_norm(1:3,:) = FM_transNrot(1:3,:)./Mg_fly;
FM_transNrot_norm(4:6,:) = FM_transNrot(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

FM_L_transNrot_norm(1:3,:) = FM_L_transNrot(1:3,:)./Mg_fly;
FM_L_transNrot_norm(4:6,:) = FM_L_transNrot(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

FM_R_transNrot_norm(1:3,:) = FM_R_transNrot(1:3,:)./Mg_fly;
FM_R_transNrot_norm(4:6,:) = FM_R_transNrot(4:6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

% Fx_transNrot = FM_transNrot(1,:)./Mg_fly;
% Fy_transNrot = FM_transNrot(2,:)./Mg_fly;
% Fz_transNrot = FM_transNrot(3,:)./Mg_fly;
% 
% Mx_transNrot = FM_transNrot(4,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% My_transNrot = FM_transNrot(5,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% Mz_transNrot = FM_transNrot(6,:)./freq^2./(1e-3*wing_length)^5*1e-3;
% 
% F_transNrot = sqrt(Fx_transNrot.^2 + Fy_transNrot.^2 + Fz_transNrot.^2);

%% plot FnM
% figure
% subplot(3,2,1)
% hold on
% plot(t_norm,Fx_trans)
% 
% subplot(3,2,3)
% hold on
% plot(t_norm,Fy_trans)
% 
% subplot(3,2,5)
% hold on
% plot(t_norm,Fz_trans)
% 
% subplot(3,2,2)
% hold on
% plot(t_norm,Mx_trans)
% 
% subplot(3,2,4)
% hold on
% plot(t_norm,My_trans)
% 
% subplot(3,2,6)
% hold on
% plot(t_norm,Mz_trans)
% 
% subplot(3,2,1)
% hold on
% plot(t_norm,Fx_transNrot,'--r')
% grid on
% 
% subplot(3,2,3)
% hold on
% plot(t_norm,Fy_transNrot,'--r')
% grid on
% 
% subplot(3,2,5)
% hold on
% plot(t_norm,Fz_transNrot,'--r')
% grid on
% 
% subplot(3,2,2)
% hold on
% plot(t_norm,Mx_transNrot,'--r')
% grid on
% 
% subplot(3,2,4)
% hold on
% plot(t_norm,My_transNrot,'--r')
% grid on
% 
% subplot(3,2,6)
% hold on
% plot(t_norm,Mz_transNrot,'--r')
% grid on



