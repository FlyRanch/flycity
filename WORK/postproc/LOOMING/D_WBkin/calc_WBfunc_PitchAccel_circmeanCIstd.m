% wb mean
            t_bins = [t_wb_PitchAccel_bins t_wb_PitchAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_wb_L_PitchAccel_bins stroke_wb_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_PitchAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_L_PitchAccel_bins pitch_wb_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_PitchAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_L_PitchAccel_bins dev_wb_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_PitchAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_PitchAccel_bins t_ds_PitchAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_ds_L_PitchAccel_bins stroke_ds_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_ds_PitchAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_ds_L_PitchAccel_bins pitch_ds_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_ds_PitchAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_ds_L_PitchAccel_bins dev_ds_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_ds_PitchAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_PitchAccel_bins t_us_PitchAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_us_L_PitchAccel_bins stroke_us_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_us_PitchAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_us_L_PitchAccel_bins pitch_us_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_us_PitchAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_us_L_PitchAccel_bins dev_us_R_PitchAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_us_PitchAccel_bins_meanCIstd = var_meanCIstd;

