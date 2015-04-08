%% Left

% wb mean
            t_bins = [t_wb_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_wb_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_ds_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_ds_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_ds_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_ds_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_ds_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_ds_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_us_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_us_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_us_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_us_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_us_L_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_us_L_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% Right

% wb mean
            t_bins = [t_wb_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_wb_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_wb_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_wb_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_wb_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_wb_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_wb_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_ds_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_ds_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_ds_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_ds_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_ds_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_ds_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [stroke_us_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            stroke_us_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [pitch_us_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitch_us_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [dev_us_R_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;

                var_meanCIstd(i,:) = [mu ul ll std];
            end
            dev_us_R_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% L-R

% wb mean
            t_bins = [t_wb_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_wb_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_wb_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_wb_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_wb_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_wb_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_wb_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_ds_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_ds_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_ds_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_ds_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_ds_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_ds_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3funcClipped_bins];           
            binx = t_bins(:,1);
            
            var = [Dstroke_us_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dstroke_us_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;
            
            var = [Dpitch_us_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Dpitch_us_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

            var = [Ddev_us_S2S3funcClipped_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            Ddev_us_S2S3funcClipped_bins_meanCIstd = var_meanCIstd;

