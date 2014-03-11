    
    plot_wbNORM_wingatt_steady
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.png'])
    plot2svg(['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.svg'])
    
    plot_WBfunc_heatmap_steady
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.png'])
    plot2svg(['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD_heatmap.svg'])
    
    plot_WBfreqNpitchV_hist_steady_LnR
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.fig'])
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.png'])
    plot2svg(['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.svg'])
