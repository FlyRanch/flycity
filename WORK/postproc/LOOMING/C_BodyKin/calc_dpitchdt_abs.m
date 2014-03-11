dpitchdWB=gradient(pitch_dot')';
dpitchdt = dpitchdWB*freq_steady;
dpitchdms = dpitchdt/1000;

dpitchdt_max = max(abs(dpitchdt));
n_seqs = sum(isnan(dpitchdt_max)==0);

dpitchdtMaxMean = nanmean(dpitchdt_max)
dpitchdtMaxMedian = nanmedian(dpitchdt_max)
dpitchdtMaxMode = mode(dpitchdt_max)
dpitchdtMaxSTE = nanstd(dpitchdt_max)/sqrt(n_seqs)

dpitchdWBMaxMean = dpitchdtMaxMean/freq_steady
dpitchdWBMaxMedian = dpitchdtMaxMedian/freq_steady
dpitchdWBMaxMode = dpitchdtMaxMode/freq_steady
dpitchdWBMaxSTE = dpitchdtMaxSTE/freq_steady

dpitchdmsMaxMean = dpitchdtMaxMean/1000
dpitchdmsMaxMedian = dpitchdtMaxMedian/1000
dpitchdmsMaxMode = dpitchdtMaxMode/1000
dpitchdmsMaxSTE = dpitchdtMaxSTE/1000

%% plot
% figure,hist(dpitchdt_max)
% 
% figure, plot(dpitchdt_max)
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),dpitchdWB(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -15 15])
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),dpitchdms(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -5 5])
