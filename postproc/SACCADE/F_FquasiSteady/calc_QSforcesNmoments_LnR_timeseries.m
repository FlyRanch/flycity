%% WBmod kin LEFT
stroke = stroke_wb_L_seq_bins_mean_all;
dev = dev_wb_L_seq_bins_mean_all;
rot = pitch_wb_L_seq_bins_mean_all;
t = t_wb_seq_bins_mean_all;

[stroke,stroke_dot,stroke_dot_dot,dev,dev_dot,dev_dot_dot,rot,rot_dot,rot_dot_dot,dt] = ...
    calc_WBkin_timeseries(t,stroke,dev,rot);

WBkin2kinL

%% calc WBmod kin RIGHT
stroke = stroke_wb_R_seq_bins_mean_all;
dev = dev_wb_R_seq_bins_mean_all;
rot = pitch_wb_R_seq_bins_mean_all;
t = t_wb_seq_bins_mean_all;

[stroke,stroke_dot,stroke_dot_dot,dev,dev_dot,dev_dot_dot,rot,rot_dot,rot_dot_dot,dt] = ...
    calc_WBkin_timeseries(t,stroke,dev,rot);

WBkin2kinR


%% calc WB kin velocities
[kine] = calc_angular_velocities_INCdots_timeseries(stroke_L,dev_L,rot_L,stroke_R,dev_R,rot_R,stroke_dot_L,dev_dot_L,rot_dot_L,stroke_dot_R,dev_dot_R,rot_dot_R,stroke_dot_dot_L,dev_dot_dot_L,rot_dot_dot_L,stroke_dot_dot_R,dev_dot_dot_R,rot_dot_dot_R,t,R_strk);

kine.R_strk = R_strk;
kine.u_strk = u_strk;
kine.w_strk = w_strk;
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

%% wingbeat average forces & torques
FM_trans_WBmean = nan(length(t_wb_seq_mean_all),size(FM_trans_norm,1));
FM_transNrot_WBmean = nan(length(t_wb_seq_mean_all),size(FM_trans_norm,1));

FM_trans_norm_WBmean = nan(length(t_wb_seq_mean_all),size(FM_trans_norm,1));
FM_transNrot_norm_WBmean = nan(length(t_wb_seq_mean_all),size(FM_trans_norm,1));

for i = 1:length(t_wb_seq_mean_all)-1
    if isnan(t_wb_seq_mean_all(i))==0 && isnan(t_wb_seq_mean_all(i+1))==0
    
        subset = find(t_wb_seq_bins_mean_all >= t_wb_seq_mean_all(i) & t_wb_seq_bins_mean_all < t_wb_seq_mean_all(i+1));
        
        FM_trans_WBmean(i,:)        = nanmean(FM_trans(:,subset)');
        FM_transNrot_WBmean(i,:)    = nanmean(FM_transNrot(:,subset)');
        
        FM_trans_norm_WBmean(i,:)        = nanmean(FM_trans_norm(:,subset)');
        FM_transNrot_norm_WBmean(i,:)    = nanmean(FM_transNrot_norm(:,subset)');
    end
end

%% plot
if plot_on ==1
    %% plot kinematics
    figure(1)
    subplot(3,3,1)
    hold off
    plot(t,stroke_L)
    hold on
    plot(t,stroke_R,'c')
    % plot(t,stroke_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke','fontsize',10)

    subplot(3,3,2)
    hold off
    plot(t,stroke_dot_L)
    hold on
    plot(t,stroke_dot_R,'c')
    % plot(t,stroke_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke dot','fontsize',10)

    subplot(3,3,3)
    hold off
    plot(t,stroke_dot_dot_L)
    hold on
    plot(t,stroke_dot_dot_R,'c')
    % plot(t,stroke_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('stroke ddot','fontsize',10)

    subplot(3,3,4)
    hold off
    plot(t,dev_L)
    hold on
    plot(t,dev_R,'c')
    % plot(t,dev_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev','fontsize',10)

    subplot(3,3,5)
    hold off
    plot(t,dev_dot_L)
    hold on
    plot(t,dev_dot_R,'c')
    % plot(t,dev_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev dot','fontsize',10)

    subplot(3,3,6)
    hold off
    plot(t,dev_dot_dot_L)
    hold on
    plot(t,dev_dot_dot_R,'c')
    % plot(t,dev_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('dev ddot','fontsize',10)

    subplot(3,3,7)
    hold off
    plot(t,rot_L)
    hold on
    plot(t,rot_R,'c')
    % plot(t,rot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot','fontsize',10)

    subplot(3,3,8)
    hold off
    plot(t,rot_dot_L)
    hold on
    plot(t,rot_dot_R,'c')
    % plot(t,rot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot dot','fontsize',10)

    subplot(3,3,9)
    hold off
    plot(t,rot_dot_dot_L)
    hold on
    plot(t,rot_dot_dot_R,'c')
    % plot(t,rot_dot_dot_dps,'--r')
    xlabel('t*','fontsize',10) 
    ylabel('rot ddot','fontsize',10)

    saveas(gca,['WBkin_quasisteadyModel_timeseries.png'])

    %% plot FnM
    figure(2)
    subplot(3,2,1)
    hold off
    plot(t,FM_trans_norm(1,:))

    subplot(3,2,3)
    hold off
    plot(t,FM_trans_norm(2,:))

    subplot(3,2,5)
    hold off
    plot(t,FM_trans_norm(3,:))

    subplot(3,2,2)
    hold off
    plot(t,FM_trans_norm(4,:))

    subplot(3,2,4)
    hold off
    plot(t,FM_trans_norm(5,:))

    subplot(3,2,6)
    hold off
    plot(t,FM_trans_norm(6,:))

    subplot(3,2,1)
    hold on
    plot(t,FM_transNrot_norm(1,:),'--r')
    grid on

    subplot(3,2,3)
    hold on
    plot(t,FM_transNrot_norm(2,:),'--r')
    grid on

    subplot(3,2,5)
    hold on
    plot(t,FM_transNrot_norm(3,:),'--r')
    grid on

    subplot(3,2,2)
    hold on
    plot(t,FM_transNrot_norm(4,:),'--r')
    grid on

    subplot(3,2,4)
    hold on
    plot(t,FM_transNrot_norm(5,:),'--r')
    grid on

    subplot(3,2,6)
    hold on
    plot(t,FM_transNrot_norm(6,:),'--r')
    grid on

    saveas(gca,['FnM_quasisteadyModel_timeseries.png'])
    
    %% plot FnM WBmean
    figure(3)
    subplot(3,2,1)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,1),'-b.')

    subplot(3,2,3)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,2),'-b.')

    subplot(3,2,5)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,3),'-b.')

    subplot(3,2,2)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,4),'-b.')

    subplot(3,2,4)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,5),'-b.')

    subplot(3,2,6)
    hold off
    plot(t_wb_seq_mean_all,FM_trans_norm_WBmean(:,6),'-b.')

    subplot(3,2,1)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,1),'-ro')
    grid on

    subplot(3,2,3)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,2),'-ro')
    grid on

    subplot(3,2,5)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,3),'-ro')
    grid on

    subplot(3,2,2)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,4),'-ro')
    grid on

    subplot(3,2,4)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,5),'-ro')
    grid on

    subplot(3,2,6)
    hold on
    plot(t_wb_seq_mean_all,FM_transNrot_norm_WBmean(:,6),'-ro')
    grid on

    saveas(gca,['FnM_quasisteadyModel_timeseries_WBmean.png'])
end

