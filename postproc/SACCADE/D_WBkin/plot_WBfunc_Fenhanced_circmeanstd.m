binx = t_wb_Fenhance_bins(:,1);
t_bins = [t_wb_Fenhance_bins t_wb_Fenhance_bins];           


            subplot(3,1,1)
            var = [stroke_wb_L_Fenhance_bins stroke_wb_R_Fenhance_bins];
            plot_WBmeanstd
            
            subplot(3,1,2)
            var = [pitch_wb_L_Fenhance_bins pitch_wb_R_Fenhance_bins];
            plot_WBmeanstd

            subplot(3,1,3)
            var = [dev_wb_L_Fenhance_bins dev_wb_R_Fenhance_bins];
            plot_WBmeanstd

