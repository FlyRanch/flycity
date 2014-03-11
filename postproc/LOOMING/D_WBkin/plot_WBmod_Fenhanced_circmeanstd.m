binx = t_wb_Fenhance_bins(:,1);
t_bins = [t_wb_Fenhance_bins t_wb_Fenhance_bins];           


            subplot(3,1,1)
            var = [strokeMOD_wb_L_Fenhance_bins strokeMOD_wb_R_Fenhance_bins];
            plot_WBmeanstd
            
            subplot(3,1,2)
            var = [pitchMOD_wb_L_Fenhance_bins pitchMOD_wb_R_Fenhance_bins];
            plot_WBmeanstd

            subplot(3,1,3)
            var = [devMOD_wb_L_Fenhance_bins devMOD_wb_R_Fenhance_bins];
            plot_WBmeanstd

