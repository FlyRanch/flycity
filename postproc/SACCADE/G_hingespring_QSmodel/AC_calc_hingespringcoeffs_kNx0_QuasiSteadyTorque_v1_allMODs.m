% clear;
clc;
close all
warning off

loc_set = 2
nr_timepoints = 5000;

load(['hingespring_Arot_vs_dTspanwise_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.mat'])

figure
hold on
plot(t_norm,rot_fwd_all(:,1)-90,'-k','linewidth',2)
plot(t_norm,rot_fwd_all(:,end)-90,'-r','linewidth',2)
plot(t_norm,rot_rwd_all(:,end)-90,'-b','linewidth',2)

plot(t_norm(n_rot_max_fwd_all(1)),rot_fwd_all(n_rot_max_fwd_all(1),1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_min_fwd_all(1)),rot_fwd_all(n_rot_min_fwd_all(1),1)-90,'ok','linewidth',2)

plot(t_norm(n_rot_max_fwd_all(end)),rot_fwd_all(n_rot_max_fwd_all(end),end)-90,'or','linewidth',2)
plot(t_norm(n_rot_min_fwd_all(end)),rot_fwd_all(n_rot_min_fwd_all(end),end)-90,'or','linewidth',2)

plot(t_norm(n_rot_max_rwd_all(end)),rot_rwd_all(n_rot_max_rwd_all(end),end)-90,'ob','linewidth',2)
plot(t_norm(n_rot_min_rwd_all(end)),rot_rwd_all(n_rot_min_rwd_all(end),end)-90,'ob','linewidth',2)

legend('steady','fwd','rwd')
axis tight
axis square
axis([0,1,-90,90])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-90 0 90])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')

saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% calc rot0 as function of Tyaw
rot0_max_fwd_all_rad = rot_max_fwd_all + Mwing_total_fwd_atRotmax/k_ArotdMtotal_fwd_norm_rad;
rot0_min_fwd_all_rad = rot_min_fwd_all + Mwing_total_fwd_atRotmin/k_ArotdMtotal_fwd_norm_rad;

rot0_max_rwd_all_rad = rot_max_rwd_all + Mwing_total_rwd_atRotmax/k_ArotdMtotal_rwd_norm_rad;
rot0_min_rwd_all_rad = rot_min_rwd_all + Mwing_total_rwd_atRotmin/k_ArotdMtotal_rwd_norm_rad;

rot_max_fwd_all_deg = rad2deg(rot_max_fwd_all);
rot_min_fwd_all_deg = rad2deg(rot_min_fwd_all);

rot_max_rwd_all_deg = rad2deg(rot_max_rwd_all);
rot_min_rwd_all_deg = rad2deg(rot_min_rwd_all);

rot0_max_fwd_all_deg = rot_max_fwd_all_deg + Mwing_total_fwd_atRotmax/k_ArotdMtotal_fwd_norm_deg;
rot0_min_fwd_all_deg = rot_min_fwd_all_deg + Mwing_total_fwd_atRotmin/k_ArotdMtotal_fwd_norm_deg;

rot0_max_rwd_all_deg = rot_max_rwd_all_deg + Mwing_total_rwd_atRotmax/k_ArotdMtotal_rwd_norm_deg;
rot0_min_rwd_all_deg = rot_min_rwd_all_deg + Mwing_total_rwd_atRotmin/k_ArotdMtotal_rwd_norm_deg;

% trendlines
rot0Tyaw_coeff_fwd_rotmax_deg = polyfit([YawTorques_all],[rot0_max_fwd_all_deg],1);
rot0Tyaw_coeff_rwd_rotmax_deg = polyfit([YawTorques_all],[rot0_max_rwd_all_deg],1);

rot0Tyaw_coeff_fwd_rotmin_deg = polyfit([YawTorques_all],[rot0_min_fwd_all_deg],1);
rot0Tyaw_coeff_rwd_rotmin_deg = polyfit([YawTorques_all],[rot0_min_rwd_all_deg],1);

    %% save data & plots
    save(['hingespring_Arot_vs_dTspanwise_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.mat'])
    
    mkdir('qsModel_FnM_n_Arot_VS_dTspanwise')
    cd('qsModel_FnM_n_Arot_VS_dTspanwise')

%% plot results
figure

% rot0 @ rot_max
rot0max_min = 40;
rot0max_max = 50;

subplot(2,1,1)
hold on
plot(YawTorques_all,rot0_max_fwd_all_deg,'dk','markerfacecolor','w','markersize',10)
plot(YawTorques_all,rot0_max_rwd_all_deg,'ok','markerfacecolor','w','markersize',10)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_fwd_rotmax_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_rwd_rotmax_deg,[min(YawTorques_all) max(YawTorques_all)]),'color',[.5 .5 .5],'linewidth',2)

legend('fwd@ds','rwd@ds')
axis tight
% axis square
axis([YawTorque_min,YawTorque_max,rot0max_min,rot0max_max])
set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
set(gca,'ytick',rot0max_min:(rot0max_max-rot0max_min):rot0max_max)
xlabel('normalized yaw torque')
ylabel('rot0 [deg]')
title('rot0 @ mid downstroke')

% rot0 @ rot_min
rot0min_min = 140;
rot0min_max = 150;

subplot(2,1,2)
hold on
plot(YawTorques_all,rot0_min_fwd_all_deg,'dk','markerfacecolor',[.5 .5 .5],'markersize',10)
plot(YawTorques_all,rot0_min_rwd_all_deg,'ok','markerfacecolor',[.5 .5 .5],'markersize',10)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_rwd_rotmin_deg,[min(YawTorques_all) max(YawTorques_all)]),'color',[.5 .5 .5],'linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_fwd_rotmin_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)

legend('fwd@us','rwd@us')
axis tight
% axis square
axis([YawTorque_min,YawTorque_max,rot0min_min,rot0min_max])
set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
set(gca,'ytick',rot0min_min:(rot0min_max-rot0min_min):rot0min_max)
xlabel('normalized yaw torque')
ylabel('rot0 [deg]')
title('rot0 @ mid upstroke')

saveas(gca,['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])
