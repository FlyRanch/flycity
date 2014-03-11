drolldWB=gradient(roll_dot_mirror')';
drolldt = drolldWB*freq_steady;
drolldms = drolldt/1000;

drolldt_max = max(abs(drolldt));
n_seqs = sum(isnan(drolldt_max)==0);

drolldtMaxMean = nanmean(drolldt_max)
drolldtMaxMedian = nanmedian(drolldt_max)
drolldtMaxMode = mode(drolldt_max)
drolldtMaxSTE = nanstd(drolldt_max)/sqrt(n_seqs)

drolldWBMaxMean = drolldtMaxMean/freq_steady
drolldWBMaxMedian = drolldtMaxMedian/freq_steady
drolldWBMaxMode = drolldtMaxMode/freq_steady
drolldWBMaxSTE = drolldtMaxSTE/freq_steady

drolldmsMaxMean = drolldtMaxMean/1000
drolldmsMaxMedian = drolldtMaxMedian/1000
drolldmsMaxMode = drolldtMaxMode/1000
drolldmsMaxSTE = drolldtMaxSTE/1000

%% plot
% figure,hist(drolldt_max)
% 
% figure, plot(drolldt_max)
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),drolldWB(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -15 15])
% 
% figure,hold on
% for i=1:length(color_var)
%     if isnan(color_var(i)) == 0 && color_var(i)~=0
%         plot(t-t_shift(i),drolldms(:,i),'-','color',grey_color,'linewidth',.25)
%     end
% end
% axis([t_start t_stop -5 5])
