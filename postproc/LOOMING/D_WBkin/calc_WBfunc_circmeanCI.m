binx = t_steady_bins(:,1);
t_bins = [t_steady_bins t_steady_bins];           


            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                var_meanCI(i,:) = [mu ul ll];

            end
            stroke_wb_steady_bins_meanCI = var_meanCI;
            
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                var_meanCI(i,:) = [mu ul ll];

            end
            pitch_wb_steady_bins_meanCI = var_meanCI;

            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                var_meanCI(i,:) = [mu ul ll];

            end
            dev_wb_steady_bins_meanCI = var_meanCI;

