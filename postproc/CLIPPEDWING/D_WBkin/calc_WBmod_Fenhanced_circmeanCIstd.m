% wb mean
            t_bins = [t_wb_Fenhance_bins t_wb_Fenhance_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_wb_L_Fenhance_bins strokeMOD_wb_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_wb_Fenhance_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_wb_L_Fenhance_bins pitchMOD_wb_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_wb_Fenhance_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_wb_L_Fenhance_bins devMOD_wb_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_wb_Fenhance_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_Fenhance_bins t_ds_Fenhance_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_ds_L_Fenhance_bins strokeMOD_ds_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_ds_Fenhance_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_ds_L_Fenhance_bins pitchMOD_ds_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_ds_Fenhance_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_ds_L_Fenhance_bins devMOD_ds_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_ds_Fenhance_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_Fenhance_bins t_us_Fenhance_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_us_L_Fenhance_bins strokeMOD_us_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_us_Fenhance_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_us_L_Fenhance_bins pitchMOD_us_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_us_Fenhance_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_us_L_Fenhance_bins devMOD_us_R_Fenhance_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_us_Fenhance_bins_meanCIstd = var_meanCIstd;

