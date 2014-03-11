binx = t_wb_steady_bins(:,1);
t_bins = [t_wb_steady_bins t_wb_steady_bins];           


            subplot(3,1,1)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
            
            subplot(3,1,2)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd

            subplot(3,1,3)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd

