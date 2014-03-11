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

if plot_WB_SUBSETS == 1

    % ROLL ACCEL SUBSETS
    subset = -roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 2*subset_std;

    % NO OFF
    subset_OFF = pitch_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetRollAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetRollAccel.png')
    plot2svg('WBnorm_wingkin_subsetRollAccel.svg')

    % YAW ACCEL SUBSETS
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 2*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetYawAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccel.png')
    plot2svg('WBnorm_wingkin_subsetYawAccel.svg')

    % PITCH ACCEL NO ROLL SUBSETS
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.5*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel.fig')
    saveas(gca,'WBnorm_wingkin_subsetPitchAccel.png')
    plot2svg('WBnorm_wingkin_subsetPitchAccel.svg')
    
    
    % FORCE NO ROLL SUBSETS
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 2.75*subset_std;

    % NO OFF
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetF.fig')
    saveas(gca,'WBnorm_wingkin_subsetF.png')
    plot2svg('WBnorm_wingkin_subsetF.svg')
end
    
    %% INC filters
if plot_WB_SUBSETS_FILTER == 1
    
    % ROLL ACCEL NO PITCH SUBSETS
    subset = -roll_dot_dot_mean_wb;
%     subset = roll_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.5*subset_std;
    
    subset_OFF = pitch_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 1*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
%     plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_neg_trendline
    saveas(gca,'WBnorm_wingkin_subsetRollAccelNOpitch.fig')
    saveas(gca,'WBnorm_wingkin_subsetRollAccelNOpitch.png')
    plot2svg('WBnorm_wingkin_subsetRollAccelNOpitch.svg')

    % YAW ACCEL NO ROLL SUBSETS
    subset = yaw_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.25*subset_std;
    
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = 1*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetYawAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetYawAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetYawAccelNOroll.svg')

    
    % PITCH ACCEL NO ROLL SUBSETS
    subset = pitch_dot_dot_mean_wb;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1*subset_std;
    
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetPitchAccelNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetPitchAccelNOroll.png')
    plot2svg('WBnorm_wingkin_subsetPitchAccelNOroll.svg')
    
    
    % FORCE NO ROLL SUBSETS
    subset = F_mean_wb-1;
    subset_mean = nanmean(subset(:));
    subset_std = nanstd(subset(:));
    subset_mid = .05*subset_std;
    subset_max = 1.75*subset_std;
    
    subset_OFF = roll_dot_dot_mean_wb;
    subset_OFF_std = nanstd(subset_OFF(:));
    subset_OFF_mid = .5*subset_OFF_std;

    plot_wbNORM_wingatt_2subsetsINCoff_INCmirror_pos
    saveas(gca,'WBnorm_wingkin_subsetForceNOroll.fig')
    saveas(gca,'WBnorm_wingkin_subsetForceNOroll.png')
    plot2svg('WBnorm_wingkin_subsetForceNOroll.svg')
end
