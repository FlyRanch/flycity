dyawdWB=gradient(yaw_dot_mirror')';
dyawdt = dyawdWB*freq_steady;
dyawdms = dyawdt/1000;

dyawdt_max = max(abs(dyawdt));
n_seqs = sum(isnan(dyawdt_max)==0);

dyawdtMaxMean = nanmean(dyawdt_max)
dyawdtMaxMedian = nanmedian(dyawdt_max)
dyawdtMaxMode = mode(dyawdt_max)
dyawdtMaxSTE = nanstd(dyawdt_max)/sqrt(n_seqs)

dyawdWBMaxMean = dyawdtMaxMean/freq_steady
dyawdWBMaxMedian = dyawdtMaxMedian/freq_steady
dyawdWBMaxMode = dyawdtMaxMode/freq_steady
dyawdWBMaxSTE = dyawdtMaxSTE/freq_steady

dyawdmsMaxMean = dyawdtMaxMean/1000
dyawdmsMaxMedian = dyawdtMaxMedian/1000
dyawdmsMaxMode = dyawdtMaxMode/1000
dyawdmsMaxSTE = dyawdtMaxSTE/1000

%% plot
% figure,hist(dyawdt_max)
% 
% figure, plot(dyawdt_max)
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),dyawdWB(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -15 15])
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),dyawdms(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -5 5])
