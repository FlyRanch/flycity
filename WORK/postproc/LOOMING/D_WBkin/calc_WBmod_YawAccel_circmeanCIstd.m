%% forwards
% wb mean
            t_bins = [t_wb_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_wb_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_wb_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_wb_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_wb_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_wb_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_wb_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_ds_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_ds_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_ds_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_ds_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_ds_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_ds_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_us_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_us_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_us_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_us_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_us_fwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_us_fwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% rearwards
% wb mean
            t_bins = [t_wb_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_wb_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_wb_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_wb_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_wb_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_wb_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_wb_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_ds_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_ds_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_ds_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_ds_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_ds_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_ds_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_us_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_us_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_us_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_us_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_us_rwd_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_us_rwd_YawAccel_bins_meanCIstd = var_meanCIstd;

%% for - rearwards
% wb mean
            t_bins = [t_wb_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_wb_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_wb_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_wb_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_wb_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_wb_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_wb_YawAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_ds_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_ds_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_ds_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_ds_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_ds_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_ds_YawAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_YawAccel_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_us_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_us_YawAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_us_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_us_YawAccel_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_us_YawAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_us_YawAccel_bins_meanCIstd = var_meanCIstd;

