binx = t_steady_bins(:,1);
t_bins = [t_steady_bins t_steady_bins];           


            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_steady_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_steady_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_steady_bins_meanCIstd = var_meanCIstd;

