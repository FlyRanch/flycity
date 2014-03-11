

%% max values
dh=gradient(stim_angle_vel_plot')';
dh(dh>45)=nan;
dh(dh<-45)=nan;
dhdt = dh/dt;
dhdms = dhdt/1000;

dhdt_max = max(abs(dhdt));
dhdWB_max = dhdt_max/freq_steady;
dhdms_max = dhdt_max/1000;

n_seqs = sum(isnan(dhdt_max)==0);

dhdtMaxMean = nanmean(dhdt_max)
dhdtMaxMedian = nanmedian(dhdt_max)
dhdtMaxMode = mode(dhdt_max)
dhdtMaxSTE = nanstd(dhdt_max)/sqrt(n_seqs)
dhdtMaxQuartiles = prctile(dhdt_max,[25 75])

dhdWBMaxMean = dhdtMaxMean/freq_steady
dhdWBMaxMedian = dhdtMaxMedian/freq_steady
dhdWBMaxMode = dhdtMaxMode/freq_steady
dhdWBMaxSTE = dhdtMaxSTE/freq_steady
dhdWBMaxQuartiles = dhdtMaxQuartiles/freq_steady

dhdmsMaxMean = dhdtMaxMean/1000
dhdmsMaxMedian = dhdtMaxMedian/1000
dhdmsMaxMode = dhdtMaxMode/1000
dhdmsMaxSTE = dhdtMaxSTE/1000
dhdmsMaxQuartiles = dhdtMaxQuartiles/1000

%% turn mean
% dt_turn = (n_post-n_pre)*dt;
% dheadingdt = turn_angle_vel_mirror./dt_turn
% dheadingdms = dheadingdt/1000;
% dheadingdWB = dheadingdt/freq_steady;
% 
% dheadingdt_mean = nanmean(abs(dheadingdt))
% dheadingdms_mean = nanmean(abs(dheadingdms))
% dheadingdWB_mean = nanmean(abs(dheadingdWB))
% 
% dheadingdt_ste = nanstd(abs(dheadingdt))/sqrt(n_seqs)
% dheadingdms_ste = nanstd(abs(dheadingdms))/sqrt(n_seqs)
% dheadingdWB_ste = nanstd(abs(dheadingdWB))/sqrt(n_seqs)

%% plot
% 
% figure,boxplot(dhdms_max)
% figure,hist(log10(dhdt_max))
% 
% figure, plot(dhdt_max)
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),dhdms(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -50 300])

