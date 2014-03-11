%% upwards

% wb mean
            t_bins = [t_wb_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_wb_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_up_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_up_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_up_RollAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_ds_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_ds_up_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_ds_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_ds_up_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_ds_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_ds_up_RollAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_us_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_us_up_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_us_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_us_up_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_us_up_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_us_up_RollAccel_bins_meanCIstd = var_meanCIstd;

%% downwards

% wb mean
            t_bins = [t_wb_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_wb_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_down_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_down_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_down_RollAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_ds_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_ds_down_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_ds_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_ds_down_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_ds_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_ds_down_RollAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_us_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_us_down_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_us_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_us_down_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [dev_us_down_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_us_down_RollAccel_bins_meanCIstd = var_meanCIstd;

%% up - downwards

% wb mean
            t_bins = [t_wb_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_wb_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_wb_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_wb_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_wb_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_wb_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_wb_RollAccel_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_ds_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_ds_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_ds_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_ds_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_ds_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_ds_RollAccel_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_RollAccel_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_us_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_us_RollAccel_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_us_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_us_RollAccel_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_us_RollAccel_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                [mu ul ll] = circ_mean_deg_nonan(var_now);
                [std std0] = circ_std_deg_nonan(var_now);
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_us_RollAccel_bins_meanCIstd = var_meanCIstd;

