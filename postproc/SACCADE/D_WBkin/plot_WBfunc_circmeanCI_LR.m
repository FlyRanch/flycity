binx = t_steady_bins(:,1);
t_bins = t_steady_bins;           


            subplot(3,5,1)
            var = stroke_wb_L_steady_bins;
            plot_WBmeanCI
            stroke_wb_L_steady_bins_meanCI = var_meanCI;
            subplot(3,5,6)
            var = stroke_wb_R_steady_bins;
            plot_WBmeanCI
            stroke_wb_R_steady_bins_meanCI = var_meanCI;
            subplot(3,5,11)
            var = Dstroke_wb_steady_bins;
            plot_WBmeanCI
            Dstroke_wb_steady_bins_meanCI = var_meanCI;
            
            subplot(3,5,2)
            var = pitch_wb_L_steady_bins;
            plot_WBmeanCI
            pitch_wb_L_steady_bins_meanCI = var_meanCI;
            subplot(3,5,7)
            var = pitch_wb_R_steady_bins;
            plot_WBmeanCI
            pitch_wb_R_steady_bins_meanCI = var_meanCI;
            subplot(3,5,12)
            var = Dpitch_wb_steady_bins;
            plot_WBmeanCI
            Dpitch_wb_steady_bins_meanCI = var_meanCI;

            subplot(3,5,3)
            var = dev_wb_L_steady_bins;
            plot_WBmeanCI
            dev_wb_L_steady_bins_meanCI = var_meanCI;
            subplot(3,5,8)
            var = dev_wb_R_steady_bins;
            dev_wb_R_steady_bins_meanCI = var_meanCI;
            plot_WBmeanCI
            subplot(3,5,13)
            var = Ddev_wb_steady_bins;
            plot_WBmeanCI
            Ddev_wb_steady_bins_meanCI = var_meanCI;

             subplot(3,5,4)
            var = aoa_wb_L_steady_bins;
            plot_WBmeanCI
            subplot(3,5,9)
            var = aoa_wb_R_steady_bins;
            plot_WBmeanCI
            subplot(3,5,14)
            var = Daoa_wb_steady_bins;
            plot_WBmeanCI

            subplot(3,5,5)
            var = U_wb_L_steady_bins;
            plot_WBmeanCI
            subplot(3,5,10)
            var = U_wb_R_steady_bins;
            plot_WBmeanCI
            subplot(3,5,15)
            var = DU_wb_steady_bins;
            plot_WBmeanCI
