%% NORMALIZED WINGBEATS SUBSETS NOMIRROR
if plot_WB_NOmirror == 1

    % calc wb values NOmirror
    calc_wb_data_NOmirror

    % ROLL ACCEL SUBSETS
    subset = roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_pos = 1.5*subset_std;

    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetRollAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetRollAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetRollAccel_NOmirror.svg')

    % YAW ACCEL SUBSETS
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .025*subset_std;
    subset_pos = 1.5*subset_std;

    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetYawAccel_NOmirror.svg')
    
    % PITCH ACCEL NO ROLL SUBSETS
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_pos = 1.5*subset_std;
    
    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetPitchAccel_NOmirror.svg')
    
    
    % FORCE NO ROLL SUBSETS
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_pos = 2.75*subset_std;
    
    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetF_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetF_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetF_NOmirror.svg')
end

%% NORMALIZED WINGBEATS ALL
if plot_WB_ALL == 1

    cmap_jet = jet(200);

    % ROLL ACCEL
    color_code = roll_dot_dot_mean_wb;
    color_code_mean = nanmean(color_code(:));
    color_code_std = nanstd(color_code(:));
    color_code_limit = abs(color_code_mean) + 2*color_code_std;
    color_code = 100*color_code/color_code_limit + 100;
    color_code(color_code<1) = 1;
    color_code(color_code>200) = 200;
    color_code = round(color_code);

    plot_wbNORM_wingatt_jetcolor
    saveas(gca,'WBnorm_wingkin_colorRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_colorRollAccel.png')
    plot2svg('WBnorm_wingkin_colorRollAccel.svg')

    % YAW ACCEL
    color_code = yaw_dot_dot_mean_wb;
    color_code_mean = nanmean(color_code(:));
    color_code_std = nanstd(color_code(:));
    color_code_limit = abs(color_code_mean) + 2*color_code_std;
    color_code = 100*color_code/color_code_limit + 100;
    color_code(color_code<1) = 1;
    color_code(color_code>200) = 200;
    color_code = round(color_code);

    plot_wbNORM_wingatt_jetcolor
    saveas(gca,'WBnorm_wingkin_colorYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_colorYawAccel.png')
    plot2svg('WBnorm_wingkin_colorYawAccel.svg')

    % PITCH ACCEL & FLIGHT FORCE
    % pitch
    color_code = pitch_dot_dot_mean_wb;
    color_code_mean = nanmean(color_code(:));
    color_code_std = nanstd(color_code(:));
    color_code_limit = abs(color_code_mean) + 2*color_code_std;
    color_code = 100*color_code/color_code_limit + 100;
    color_code(color_code<1) = 1;
    color_code(color_code>200) = 200;
    color_code = round(color_code);

    plot_wbNORM_wingatt_jetcolor
    saveas(gca,'WBnorm_wingkin_colorPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_colorPitchAccel.png')
    plot2svg('WBnorm_wingkin_colorPitchAccel.svg')

    % Force
    color_code = F_mean_wb-1;
    color_code_mean = nanmean(color_code(:));
    color_code_std = nanstd(color_code(:));
    color_code_limit = abs(color_code_mean) + 3*color_code_std;
    color_code = 100*color_code/color_code_limit + 100;
    color_code(color_code<1) = 1;
    color_code(color_code>200) = 200;
    color_code2 = round(color_code);

    plot_wbNORM_wingatt_jetcolor
    saveas(gca,'WBnorm_wingkin_colorF.fig')
    saveas(gca,'WBnorm_wingkin_colorF.png')
    plot2svg('WBnorm_wingkin_colorF.svg')
end

%% calc wb values INCmirror
calc_wb_data
% calc_wb_data_NOmirror

% ONLY STEADY
if plot_WB_STEADY == 1

    % STEADY LIMIT
    limit = .25;
    roll_limit = limit*nanstd(roll_dot_dot_mean_wb(:));
    pitch_limit = limit*nanstd(pitch_dot_dot_mean_wb(:));
    yaw_limit = limit*nanstd(yaw_dot_dot_mean_wb(:));
    F_limit = limit*nanstd(F_mean_wb(:));

%     color_code_now = 'b';
    color_code_now = [.5 .5 .5];
    color_mid = [.25 .25 .25];

    plot_wbNORM_wingatt_steady
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD.png'])
    plot2svg(['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD.svg'])
    
    plot_WBfunc_heatmap_steady
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.png'])
    plot2svg(['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.svg'])
    
    plot_WBfreqNdt_hist_steady_LnR
    saveas(gca,['WBfreqNdt_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.png'])
    plot2svg(['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.svg'])
    
    plot_wbNORM_wingatt_steady_minMEAN
%     plot_wbNORM_wingatt_steady_devNORM
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_minMEAN.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_minMEAN.png'])
    plot2svg(['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_minMEAN.svg'])
    
    plot_WBfunc_heatmap_steady_minMEAN
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap_minMEAN.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap_minMEAN.png'])
    plot2svg(['WBnorm_wingkin_',num2str(mid),'steadyWBs_limit',num2str(limit),'xSTD_heatmap_minMEAN.svg'])
    
    
    
end


if plot_WB_SUBSETS == 1
    
    color_mid = [.25 .25 .25];
    color_max = [0.5 0 0];
    color_min = [0 0 0.5];
    
    %% ROLL ACCEL SUBSETS
    subset = roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0.05*subset_std;
    subset_pos_max = 200000;
    subset_pos_min = 50000;
    subset_neg_min = -200000;
    subset_neg_max = -50000;
    
    subset_OFF = F_mean_wb-1;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;
    subset_OFF_mid = 0.2;

    % HISTOGRAM
    h_max = max(hist(subset(:),100));
    figure
    plot([subset_mid subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    hold on
    plot([-subset_mid -subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    plot([subset_pos_max subset_pos_max],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_pos_min subset_pos_min],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_neg_max subset_neg_max],[0 h_max],'color',color_min,'linewidth',5)
    plot([subset_neg_min subset_neg_min],[0 h_max],'color',color_min,'linewidth',5)
    hist(subset(:),100)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    saveas(gca,'WBnorm_hist_RollAccel.fig')
    saveas(gca,'WBnorm_hist_RollAccel.png')
    plot2svg('WBnorm_hist_RollAccel.svg')
    
    % POS NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_INCfreq
    figure(1)
    saveas(gca,'WBnorm_wingkin_subsetRollAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetRollAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetRollAccelnSTEADY.svg')
    figure(2)
    saveas(gca,'WBfreq_wingkin_subsetRollAccelnSTEADY.fig')
    saveas(gca,'WBfreq_wingkin_subsetRollAccelnSTEADY.png')
    plot2svg('WBfreq_wingkin_subsetRollAccelnSTEADY.svg')
    
    % POS & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccelnSTEADY.svg')
    
    % NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccelnSTEADY.svg')
    
    % POS
    subset_mid_ON = 0;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccel.svg')
   
    % NEG
    subset_mid_ON = 0;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccel.svg')
    
    
    % POS STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelNoF.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelNoF.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccelNoF.svg')
    figure(2)
    saveas(gca,'WBfreq_wingkin_subsetPosRollAccelNoF.fig')
    saveas(gca,'WBfreq_wingkin_subsetPosRollAccelNoF.png')
    plot2svg('WBfreq_wingkin_subsetPosRollAccelNoF.svg')
    
    % NEG STEADY & PITCH OFF
    subset_OFF = pitch_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;
    subset_OFF_mid = 25000;
    subset_neg_max = -25000;

    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelNOpitch.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelNOpitch.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccelNOpitch.svg')
    figure(2)
    saveas(gca,'WBfreq_wingkin_subsetNegRollAccelNOpitch.fig')
    saveas(gca,'WBfreq_wingkin_subsetNegRollAccelNOpitch.png')
    plot2svg('WBfreq_wingkin_subsetNegRollAccelNOpitch.svg')
    
    %% YAW ACCEL SUBSETS
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0.05*subset_std;
    subset_pos_max = 200000;
    subset_pos_min = 25000;
    subset_neg_min = -200000;
    subset_neg_max = -25000;

    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;
    subset_OFF_mid = 25000;

    % HISTOGRAM
    h_max = max(hist(subset(:),100));
    figure
    plot([subset_mid subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    hold on
    plot([-subset_mid -subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    plot([subset_pos_max subset_pos_max],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_pos_min subset_pos_min],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_neg_max subset_neg_max],[0 h_max],'color',color_min,'linewidth',5)
    plot([subset_neg_min subset_neg_min],[0 h_max],'color',color_min,'linewidth',5)
    hist(subset(:),100)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    saveas(gca,'WBnorm_hist_YawAccel.fig')
    saveas(gca,'WBnorm_hist_YawAccel.png')
    plot2svg('WBnorm_hist_YawAccel.svg')
    
    % POS NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_INCfreq
    figure(1)
   saveas(gca,'WBnorm_wingkin_subsetYawAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetYawAccelnSTEADY.svg')
     figure(2)
   saveas(gca,'WBfreq_wingkin_subsetYawAccelnSTEADY.fig')
    saveas(gca,'WBfreq_wingkin_subsetYawAccelnSTEADY.png')
    plot2svg('WBfreq_wingkin_subsetYawAccelnSTEADY.svg')
    
    % POS & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccelnSTEADY.svg')
    
    % NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccelnSTEADY.svg')
    
    % POS
    subset_mid_ON = 0;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccel.svg')
   
    % NEG
    subset_mid_ON = 0;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccel.svg')
    
    % POS STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccelNOroll.svg')
    
    % NEG STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccelNOroll.svg')

    %% PITCH ACCEL SUBSETS
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0.05*subset_std;
    subset_pos_max = 200000;
    subset_pos_min = 25000;
    subset_neg_min = -200000;
    subset_neg_max = -25000;

    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;
    subset_OFF_mid = 25000;

    % HISTOGRAM
    h_max = max(hist(subset(:),100));
    figure
    plot([subset_mid subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    hold on
    plot([-subset_mid -subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    plot([subset_pos_max subset_pos_max],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_pos_min subset_pos_min],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_neg_max subset_neg_max],[0 h_max],'color',color_min,'linewidth',5)
    plot([subset_neg_min subset_neg_min],[0 h_max],'color',color_min,'linewidth',5)
    hist(subset(:),100)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    saveas(gca,'WBnorm_hist_PitchAccel.fig')
    saveas(gca,'WBnorm_hist_PitchAccel.png')
    plot2svg('WBnorm_hist_PitchAccel.svg')
    
    % POS NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_INCfreq
    figure(1)
    saveas(gca,'WBnorm_wingkin_subsetPitchAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPitchAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPitchAccelnSTEADY.svg')
    figure(2)
    saveas(gca,'WBfreq_wingkin_subsetPitchAccelnSTEADY.fig')
    saveas(gca,'WBfreq_wingkin_subsetPitchAccelnSTEADY.png')
    plot2svg('WBfreq_wingkin_subsetPitchAccelnSTEADY.svg')
    
    % POS & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccelnSTEADY.svg')
    
    % NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccelnSTEADY.svg')
    
    % POS
    subset_mid_ON = 0;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccel.svg')
   
    % NEG
    subset_mid_ON = 0;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccel.svg')
        
    % POS STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccelNOroll.svg')
    
    % NEG STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccelNOroll.svg')
    
 
    %% FORCE SUBSETS
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0.05*subset_std;
    subset_pos_max = 1;
    subset_pos_min = .2;
    subset_neg_min = -.2;
    subset_neg_max = -1;

    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;
    subset_OFF_mid = 25000;

    % HISTOGRAM
    h_max = max(hist(subset(:),100));
    figure
    plot([subset_mid subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    hold on
    plot([-subset_mid -subset_mid],[0 h_max],'color',color_mid,'linewidth',5)
    plot([subset_pos_max subset_pos_max],[0 h_max],'color',color_max,'linewidth',5)
    plot([subset_pos_min subset_pos_min],[0 h_max],'color',color_max,'linewidth',5)
%     plot([subset_neg_max subset_neg_max],[0 h_max],'color',color_min,'linewidth',5)
%     plot([subset_neg_min subset_neg_min],[0 h_max],'color',color_min,'linewidth',5)
    hist(subset(:),100)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    saveas(gca,'WBnorm_hist_Force.fig')
    saveas(gca,'WBnorm_hist_Force.png')
    plot2svg('WBnorm_hist_Force.svg')
    
    % POS NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_INCfreq
    figure(1)
    saveas(gca,'WBnorm_wingkin_subsetForcenSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetForcenSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetForcenSTEADY.svg')
    figure(2)
    saveas(gca,'WBfreq_wingkin_subsetForcenSTEADY.fig')
    saveas(gca,'WBfreq_wingkin_subsetForcenSTEADY.png')
    plot2svg('WBfreq_wingkin_subsetForcenSTEADY.svg')
    
    % POS & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosForcenSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosForcenSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosForcenSTEADY.svg')
    
    % POS
    subset_mid_ON = 0;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetPosForce.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosForce.png')
    plot2svg('WBnorm_wingkin_subsetPosForce.svg')
   
    % POS STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 1;
    subset_neg_ON = 0;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF
 
   saveas(gca,'WBnorm_wingkin_subsetPosForceNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosForceNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPosForceNOroll.svg')
    
    % NEG & STEADY
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegForcenSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegForcenSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegForcenSTEADY.svg')
    
    % NEG
    subset_mid_ON = 0;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 0;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegForce.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegForce.png')
    plot2svg('WBnorm_wingkin_subsetNegForce.svg')
   
    % NEG STEADY & OFF
    subset_mid_ON = 1;
    subset_pos_ON = 0;
    subset_neg_ON = 1;
    subset_OFF_ON = 1;
    plot_wbNORM_wingatt_3subsetsINCoff_ONOFF

    saveas(gca,'WBnorm_wingkin_subsetNegForceNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegForceNOroll.png')
    plot2svg('WBnorm_wingkin_subsetNegForceNOroll.svg')
end
