%% NORMALIZED WINGBEATS SUBSETS NOMIRROR
if plot_WB_NOmirror == 1

    % calc wb values NOmirror
    calc_wb_data_NOmirror

    % ROLL ACCEL SUBSETS
    subset = roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.5*subset_std;

    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetRollAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetRollAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetRollAccel_NOmirror.svg')

    % YAW ACCEL SUBSETS
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .025*subset_std;
    subset_max = 1.5*subset_std;

    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetYawAccel_NOmirror.svg')
    
    % PITCH ACCEL NO ROLL SUBSETS
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.5*subset_std;
    
    plot_wbNORM_wingatt_3subsets_NOmirror
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel_NOmirror.fig')
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel_NOmirror.png')
    plot2svg('WBnorm_wingkin_subsetPitchAccel_NOmirror.svg')
    
    
    % FORCE NO ROLL SUBSETS
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 2.75*subset_std;
    
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
    
    % POS ROLL ACCEL SUBSETS
    color_max = [0.5 0 0];
    subset = roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 1.5*subset_std;
    subset_min = 5*subset_std;

    % NO OFF
    subset_OFF = F_mean_wb-1;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccel.svg')

    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccel_heatmap.svg')
    
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccelnSTEADY.svg')

    % INC OFF
    subset_max = 1*subset_std;
    subset_OFF_mid = .5*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
%     plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_neg_trendline
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelNoF.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosRollAccelNoF.png')
    plot2svg('WBnorm_wingkin_subsetPosRollAccelNoF.svg')

    % NEG ROLL ACCEL SUBSETS
    color_max = [0 0 0.5];
    subset = -roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 2*subset_std;

    % NO OFF
    subset_OFF = pitch_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccel.svg')

    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccel_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccelnSTEADY.svg')
    
    % INC OFF
    subset_max = 1.5*subset_std;
    subset_OFF_mid = .5*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
%     plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_neg_trendline
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelNOpitch.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegRollAccelNOpitch.png')
    plot2svg('WBnorm_wingkin_subsetNegRollAccelNOpitch.svg')


    % POS YAW ACCEL SUBSETS
    color_max = [0.5 0 0];
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 2*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccel.svg')

    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccel_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccelnSTEADY.svg')

    % INC OFF
    subset_max = 1*subset_std;
    subset_OFF_mid = .5*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosYawAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPosYawAccelNOroll.svg')

    % NEG YAW ACCEL SUBSETS
    color_max = [0 0 0.5];
    subset = -yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 1*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccel.svg')

    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetYawAccel_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccelnSTEADY.svg')
    
    % INC OFF
    subset_OFF_mid = 1*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegYawAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetNegYawAccelNOroll.svg')

    % POS PITCH ACCEL SUBSETS
    color_max = [0.5 0 0];
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 2*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccel.svg')
    
    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccel_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccelnSTEADY.svg')

    % INC OFF
    subset_max = 1*subset_std;
    subset_OFF_mid = .5*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPosPitchAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPosPitchAccelNOroll.svg')

    % NEG PITCH ACCEL SUBSETS
    color_max = [0 0 0.5];
    subset = -pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 1.5*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccel.svg')
    
    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccel_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccel_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccelnSTEADY.svg')

    % INC OFF
    subset_max = 1*subset_std;
    subset_OFF_mid = 1*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetNegPitchAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetNegPitchAccelNOroll.svg')

    % FORCE SUBSETS
    color_max = [0.5 0 0];
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = 0*subset_std;
    subset_max = 2.75*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetForce.fig')
    saveas(gca,'WBnorm_wingkin_subsetForce.png')
    plot2svg('WBnorm_wingkin_subsetForce.svg')
    
    plot_WBfunc_heatmap_hi
    saveas(gca,'WBnorm_wingkin_subsetForce_heatmap.fig')
    saveas(gca,'WBnorm_wingkin_subsetForce_heatmap.png')
    plot2svg('WBnorm_wingkin_subsetForce_heatmap.svg')
        
    % INC steady
    subset_mid = 0.05*subset_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetFnSTEADY.fig')
    saveas(gca,'WBnorm_wingkin_subsetFnSTEADY.png')
    plot2svg('WBnorm_wingkin_subsetFnSTEADY.svg')
    
    % INC OFF
    subset_max = 1.25*subset_std;
    subset_OFF_mid = .5*subset_OFF_std;
    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetForceNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetForceNOroll.png')
    plot2svg('WBnorm_wingkin_subsetForceNOroll.svg')

end