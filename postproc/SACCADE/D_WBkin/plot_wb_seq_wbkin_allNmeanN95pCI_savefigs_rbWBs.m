% wb kin

%% freq
figure
subplot(2,1,1)
hold on
% plot(t_wb_seq_pre_mean_all,f_wb_seq_pre,'.k')
% plot(t_wb_seq_post_mean_all,f_wb_seq_post,'.k')

plot(t_wb_seq_pre_mean_all,f_wb_seq_pre,'ok','MarkerFaceColor','k','MarkerSize',2)
plot(t_wb_seq_post_mean_all,f_wb_seq_post,'ok','MarkerFaceColor','k','MarkerSize',2)

% plot(t_wb_seq_pre_mean_all,f_wb_seq_pre,'.','color',[0 0 1])
% plot(t_wb_seq_post_mean_all,f_wb_seq_post,'.','color',[1 0 0])
% 
% plot(t_wb_seq_pre_mean_all,f_wb_seq_pre_mean_all,'color',[0 0 1])
% plot(t_wb_seq_post_mean_all,f_wb_seq_post_mean_all,'color',[1 0 0])

ciplot(f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all,'-k','linewidth',1)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('freq [1/sec]')
ylim([160 240])

saveas(gca,['WBvsTimeSeqs_freq_95pCI.fig'])
saveas(gca,['WBvsTimeSeqs_freq_95pCI.png'])
plot2svg(['WBvsTimeSeqs_freq_95pCI.svg'])

%% angles ALL
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])
% 
% plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all,'color',[0 0 1],'linewidth',2)
% plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all,'color',[1 0 0],'linewidth',2)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

% plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all,'color',[0 0 1],'linewidth',2)
% plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all,'color',[1 0 0],'linewidth',2)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all,'color',[0 0 1],'linewidth',2)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all,'color',[1 0 0],'linewidth',2)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_all.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_all.png'])
plot2svg(['WBvsTimeSeqs_WBangles_all.svg'])


%% angles Mean & 95% CI
figure
subplot(3,1,1)
hold on
% plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all,'--','color',[0 0 1],'linewidth',2)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all,'--','color',[1 0 0],'linewidth',2)

ciplot(stroke_wb_L_seq_bins_mean_all-1.96*stroke_wb_L_seq_bins_ste_all,...
    stroke_wb_L_seq_bins_mean_all+1.96*stroke_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
ciplot(stroke_wb_R_seq_bins_mean_all-1.96*stroke_wb_R_seq_bins_ste_all,...
    stroke_wb_R_seq_bins_mean_all+1.96*stroke_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
% plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all,'--','color',[0 0 1],'linewidth',2)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all,'--','color',[1 0 0],'linewidth',2)

ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
    dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
    dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

% xlabel('time [sec]')
xlim([-.05 .06])
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all,'--','color',[0 0 1],'linewidth',2)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all,'--','color',[1 0 0],'linewidth',2)

ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
    pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
    pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_mean_allN95pCI.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_mean_allN95pCI.png'])
plot2svg(['WBvsTimeSeqs_WBangles_mean_allN95pCI.svg'])


%% angles Only Mean
figure
subplot(3,1,1)
hold on
% plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all,'-','color',[0 0 1],'linewidth',1)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all,'-','color',[1 0 0],'linewidth',1)

% ciplot(stroke_wb_L_seq_bins_mean_all-1.96*stroke_wb_L_seq_bins_ste_all,...
%     stroke_wb_L_seq_bins_mean_all+1.96*stroke_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
% ciplot(stroke_wb_R_seq_bins_mean_all-1.96*stroke_wb_R_seq_bins_ste_all,...
%     stroke_wb_R_seq_bins_mean_all+1.96*stroke_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
% plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all,'-','color',[0 0 1],'linewidth',1)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all,'-','color',[1 0 0],'linewidth',1)

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

% xlabel('time [sec]')
xlim([-.05 .06])
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'linewidth',.25,'color',[0 0 1])
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'linewidth',.25,'color',[1 0 0])

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all,'-','color',[0 0 1],'linewidth',1)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all,'-','color',[1 0 0],'linewidth',1)

% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 0 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_mean_all.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_mean_all.png'])
plot2svg(['WBvsTimeSeqs_WBangles_mean_all.svg'])
