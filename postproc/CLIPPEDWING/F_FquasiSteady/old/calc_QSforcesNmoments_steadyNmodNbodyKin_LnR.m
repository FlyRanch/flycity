function [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
    FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
    FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
    calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(MOD_val,const,wing_model,body_model,freq,...
    stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
    stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
    stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R,plot_on)

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

%% calc WBmod kin LEFT
[stroke,stroke_dot,stroke_dot_dot,dev,dev_dot,dev_dot_dot,rot,rot_dot,rot_dot_dot,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs(nr_timepoints,freq,stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L);

WBkinMOD2MOD_L

%% calc WBmod kin LEFT
[stroke,stroke_dot,stroke_dot_dot,dev,dev_dot,dev_dot_dot,rot,rot_dot,rot_dot_dot,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs(nr_timepoints,freq,stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);

WBkinMOD2MOD_R

%% sum steady & MOD

stroke_L = stroke_L_steady + MOD_val*stroke_L_MOD;
stroke_dot_L = stroke_dot_L_steady + MOD_val*stroke_dot_L_MOD;
stroke_dot_dot_L = stroke_dot_dot_L_steady + MOD_val*stroke_dot_dot_L_MOD;

stroke_R = stroke_R_steady + MOD_val*stroke_R_MOD;
stroke_dot_R = stroke_dot_R_steady + MOD_val*stroke_dot_R_MOD;
stroke_dot_dot_R = stroke_dot_dot_R_steady + MOD_val*stroke_dot_dot_R_MOD;

dev_L = dev_L_steady + MOD_val*dev_L_MOD;
dev_dot_L = dev_dot_L_steady + MOD_val*dev_dot_L_MOD;
dev_dot_dot_L = dev_dot_dot_L_steady + MOD_val*dev_dot_dot_L_MOD;

dev_R = dev_R_steady + MOD_val*dev_R_MOD;
dev_dot_R = dev_dot_R_steady + MOD_val*dev_dot_R_MOD;
dev_dot_dot_R = dev_dot_dot_R_steady + MOD_val*dev_dot_dot_R_MOD;

rot_L = rot_L_steady + MOD_val*rot_L_MOD;
rot_dot_L = rot_dot_L_steady + MOD_val*rot_dot_L_MOD;
rot_dot_dot_L = rot_dot_dot_L_steady + MOD_val*rot_dot_dot_L_MOD;

rot_R = rot_R_steady + MOD_val*rot_R_MOD;
rot_dot_R = rot_dot_R_steady + MOD_val*rot_dot_R_MOD;
rot_dot_dot_R = rot_dot_dot_R_steady + MOD_val*rot_dot_dot_R_MOD;

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
% [ FM_trans, FM_L_trans, FM_R_trans ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
%     Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_off);

[ FM_trans, FM_L_trans, FM_R_trans ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
    Aerodynamic_forces( kine, body_model, wing_model, rot_lift_off);

% moments Nmm -> Nm
FM_trans(4:6,:) = FM_trans(4:6,:)*1e-3;
FM_L_trans(4:6,:) = FM_L_trans(4:6,:)*1e-3;
FM_R_trans(4:6,:) = FM_R_trans(4:6,:)*1e-3;

% norm
FM_trans_norm(1:3,:) = FM_trans(1:3,:)./Mg_fly;
FM_trans_norm(4:6,:) = FM_trans(4:6,:)./Mg_fly./(1e-3*wing_length);

FM_L_trans_norm(1:3,:) = FM_L_trans(1:3,:)./Mg_fly;
FM_L_trans_norm(4:6,:) = FM_L_trans(4:6,:)./Mg_fly./(1e-3*wing_length);

FM_R_trans_norm(1:3,:) = FM_R_trans(1:3,:)./Mg_fly;
FM_R_trans_norm(4:6,:) = FM_R_trans(4:6,:)./Mg_fly./(1e-3*wing_length);

%% INC rotational lift
% [ FM_transNrot, FM_L_transNrot, FM_R_transNrot ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
%     Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_on);

[ FM_transNrot, FM_L_transNrot, FM_R_transNrot ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
    Aerodynamic_forces( kine, body_model, wing_model, rot_lift_on);

% moments Nmm -> Nm
FM_transNrot(4:6,:) = FM_transNrot(4:6,:)*1e-3;
FM_L_transNrot(4:6,:) = FM_L_transNrot(4:6,:)*1e-3;
FM_R_transNrot(4:6,:) = FM_R_transNrot(4:6,:)*1e-3;

% norm
FM_transNrot_norm(1:3,:) = FM_transNrot(1:3,:)./Mg_fly;
FM_transNrot_norm(4:6,:) = FM_transNrot(4:6,:)./Mg_fly./(1e-3*wing_length);

FM_L_transNrot_norm(1:3,:) = FM_L_transNrot(1:3,:)./Mg_fly;
FM_L_transNrot_norm(4:6,:) = FM_L_transNrot(4:6,:)./Mg_fly./(1e-3*wing_length);

FM_R_transNrot_norm(1:3,:) = FM_R_transNrot(1:3,:)./Mg_fly;
FM_R_transNrot_norm(4:6,:) = FM_R_transNrot(4:6,:)./Mg_fly./(1e-3*wing_length);

%% plot
if plot_on ==1
    %% plot kinematics
    figure(1)
    subplot(3,3,1)
    hold off
    plot(t_norm,stroke_L)
    hold on
    plot(t_norm,stroke_R,'c')
    % plot(t_norm,stroke_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke','fontsize',10)

    subplot(3,3,2)
    hold off
    plot(t_norm,stroke_dot_L)
    hold on
    plot(t_norm,stroke_dot_R,'c')
    % plot(t_norm,stroke_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke dot','fontsize',10)

    subplot(3,3,3)
    hold off
    plot(t_norm,stroke_dot_dot_L)
    hold on
    plot(t_norm,stroke_dot_dot_R,'c')
    % plot(t_norm,stroke_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke ddot','fontsize',10)

    subplot(3,3,4)
    hold off
    plot(t_norm,dev_L)
    hold on
    plot(t_norm,dev_R,'c')
    % plot(t_norm,dev_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev','fontsize',10)

    subplot(3,3,5)
    hold off
    plot(t_norm,dev_dot_L)
    hold on
    plot(t_norm,dev_dot_R,'c')
    % plot(t_norm,dev_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev dot','fontsize',10)

    subplot(3,3,6)
    hold off
    plot(t_norm,dev_dot_dot_L)
    hold on
    plot(t_norm,dev_dot_dot_R,'c')
    % plot(t_norm,dev_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev ddot','fontsize',10)

    subplot(3,3,7)
    hold off
    plot(t_norm,rot_L)
    hold on
    plot(t_norm,rot_R,'c')
    % plot(t_norm,rot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot','fontsize',10)

    subplot(3,3,8)
    hold off
    plot(t_norm,rot_dot_L)
    hold on
    plot(t_norm,rot_dot_R,'c')
    % plot(t_norm,rot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot dot','fontsize',10)

    subplot(3,3,9)
    hold off
    plot(t_norm,rot_dot_dot_L)
    hold on
    plot(t_norm,rot_dot_dot_R,'c')
    % plot(t_norm,rot_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot ddot','fontsize',10)

    saveas(gca,['WBkin_quasisteadyModel_steadyWBnMOD',sprintf('%2.1f', MOD_val),'.png'])

    %% plot FnM
    figure(2)
    subplot(3,2,1)
    hold off
    plot(t_norm,FM_trans_norm(1,:))

    subplot(3,2,3)
    hold off
    plot(t_norm,FM_trans_norm(2,:))

    subplot(3,2,5)
    hold off
    plot(t_norm,FM_trans_norm(3,:))

    subplot(3,2,2)
    hold off
    plot(t_norm,FM_trans_norm(4,:))

    subplot(3,2,4)
    hold off
    plot(t_norm,FM_trans_norm(5,:))

    subplot(3,2,6)
    hold off
    plot(t_norm,FM_trans_norm(6,:))

    subplot(3,2,1)
    hold on
    plot(t_norm,FM_transNrot_norm(1,:),'--r')
    grid on

    subplot(3,2,3)
    hold on
    plot(t_norm,FM_transNrot_norm(2,:),'--r')
    grid on

    subplot(3,2,5)
    hold on
    plot(t_norm,FM_transNrot_norm(3,:),'--r')
    grid on

    subplot(3,2,2)
    hold on
    plot(t_norm,FM_transNrot_norm(4,:),'--r')
    grid on

    subplot(3,2,4)
    hold on
    plot(t_norm,FM_transNrot_norm(5,:),'--r')
    grid on

    subplot(3,2,6)
    hold on
    plot(t_norm,FM_transNrot_norm(6,:),'--r')
    grid on

    saveas(gca,['FnM_quasisteadyModel_steadyWBnMOD',sprintf('%2.1f', MOD_val),'.png'])
end

