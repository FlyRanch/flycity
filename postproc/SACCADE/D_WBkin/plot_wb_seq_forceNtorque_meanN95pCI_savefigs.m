
figure

% force
subplot(4,1,1)
hold on

n_nonan = find(isnan(F_mean_wb_seq_mean_all)==0);
ciplot(F_mean_wb_seq_mean_all(n_nonan)-1.96*F_mean_wb_seq_ste_all(n_nonan),...
    F_mean_wb_seq_mean_all(n_nonan)+1.96*F_mean_wb_seq_ste_all(n_nonan),...
    t_wb_seq_mean_all(n_nonan),[.5 .5 .5])
plot(t_wb_seq_mean_all,F_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('F/mg')
xlim([-.05 .06]) 
ylim([.75 1.25])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[0:.25:2]) 

% M_R
subplot(4,1,2)
hold on

n_nonan = find(isnan(M_R_mean_wb_seq_mean_all)==0);
ciplot(M_R_mean_wb_seq_mean_all(n_nonan)-1.96*M_R_mean_wb_seq_ste_all(n_nonan),...
    M_R_mean_wb_seq_mean_all(n_nonan)+1.96*M_R_mean_wb_seq_ste_all(n_nonan),...
    t_wb_seq_mean_all(n_nonan),[.5 .5 .5])
plot(t_wb_seq_mean_all,M_R_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('M_R')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.03:.03:.03]) 

% M_L
subplot(4,1,3)
hold on

n_nonan = find(isnan(M_L_mean_wb_seq_mean_all)==0);
ciplot(M_L_mean_wb_seq_mean_all(n_nonan)-1.96*M_L_mean_wb_seq_ste_all(n_nonan),...
    M_L_mean_wb_seq_mean_all(n_nonan)+1.96*M_L_mean_wb_seq_ste_all(n_nonan),...
    t_wb_seq_mean_all(n_nonan),[.5 .5 .5])
plot(t_wb_seq_mean_all,M_L_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('M_L')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.03:.03:.03]) 

% Myaw
subplot(4,1,4)
hold on

n_nonan = find(isnan(Myaw_mean_wb_seq_mean_all)==0);
ciplot(Myaw_mean_wb_seq_mean_all(n_nonan)-1.96*Myaw_mean_wb_seq_ste_all(n_nonan),...
    Myaw_mean_wb_seq_mean_all(n_nonan)+1.96*Myaw_mean_wb_seq_ste_all(n_nonan),...
    t_wb_seq_mean_all(n_nonan),[.5 .5 .5])
plot(t_wb_seq_mean_all,Myaw_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('Myaw')
xlim([-.05 .06]) 
ylim([0 .06])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-.03:.03:.06]) 

saveas(gca,['TorqueVsTime_rotLnRnYaw_meanN95pConfInt.fig'])
saveas(gca,['TorqueVsTime_rotLnRnYaw_meanN95pConfInt.png'])
plot2svg(['TorqueVsTime_rotLnRnYaw_meanN95pConfInt.svg'])
