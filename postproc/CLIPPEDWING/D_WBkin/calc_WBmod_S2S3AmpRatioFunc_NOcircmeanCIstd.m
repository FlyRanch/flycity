%% Left
% wb mean
            t_bins = [t_wb_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_wb_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_wb_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_wb_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;

                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_ds_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_ds_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_ds_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_ds_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_ds_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_ds_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_us_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_us_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_us_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_us_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_us_L_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_us_L_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% Right
% wb mean
            t_bins = [t_wb_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_wb_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_wb_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_wb_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_ds_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_ds_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_ds_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_ds_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_ds_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_ds_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [strokeMOD_us_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            strokeMOD_us_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [pitchMOD_us_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            pitchMOD_us_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [devMOD_us_R_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            devMOD_us_R_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% L-R
% wb mean
            t_bins = [t_wb_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_wb_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_wb_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_wb_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_wb_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_wb_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_wb_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% downstroke mean
            t_bins = [t_ds_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_ds_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_ds_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_ds_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_ds_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_ds_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_ds_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

%% upstroke mean
            t_bins = [t_us_S2S3AmpRatioFunc_bins];           
            binx = t_bins(:,1);
            
            var = [DstrokeMOD_us_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DstrokeMOD_us_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;
            
            var = [DpitchMOD_us_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DpitchMOD_us_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

            var = [DdevMOD_us_S2S3AmpRatioFunc_bins];
            for i = 1:length(binx)

                var_now = var(t_bins == binx(i));
                mu = nanmean(var_now);
                std = nanstd(var_now);
                CI = 1.96 * std / sqrt(length(var_now));
                ul = mu + CI;
                ll = mu - CI;
                
                var_meanCIstd(i,:) = [mu ul ll std];
            end
            DdevMOD_us_S2S3AmpRatioFunc_bins_meanCIstd = var_meanCIstd;

