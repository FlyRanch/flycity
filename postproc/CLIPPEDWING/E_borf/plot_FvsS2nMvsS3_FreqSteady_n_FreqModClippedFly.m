% F-S2
figure
subplot(1,2,1)
hold on

% non-cut
plot(S2_ratio(cut_type==0),Fx_norm(cut_type==0),'sk','markerfacecolor','b')
plot(S2_ratio(cut_type==0),Fy_norm(cut_type==0),'sk','markerfacecolor','r')
plot(S2_ratio(cut_type==0),Fz_norm(cut_type==0),'sk','markerfacecolor','g')

plot(S2_ratio(cut_type==0),Fx_norm_freqMod(cut_type==0),'sk','markerfacecolor','c')
plot(S2_ratio(cut_type==0),Fy_norm_freqMod(cut_type==0),'sk','markerfacecolor','m')
plot(S2_ratio(cut_type==0),Fz_norm_freqMod(cut_type==0),'sk','markerfacecolor','y')

% steady freq
plot(S2_ratio(cut_type==1),Fx_norm(cut_type==1),'ok','markerfacecolor','b')
plot(S2_ratio(cut_type==2),Fx_norm(cut_type==2),'dk','markerfacecolor','b')

plot(S2_ratio(cut_type==1),Fy_norm(cut_type==1),'ok','markerfacecolor','r')
plot(S2_ratio(cut_type==2),Fy_norm(cut_type==2),'dk','markerfacecolor','r')

plot(S2_ratio(cut_type==1),Fz_norm(cut_type==1),'ok','markerfacecolor','g')
plot(S2_ratio(cut_type==2),Fz_norm(cut_type==2),'dk','markerfacecolor','g')

% clipped fly freq
plot(S2_ratio(cut_type==1),Fx_norm_freqMod(cut_type==1),'ok','markerfacecolor','c')
plot(S2_ratio(cut_type==2),Fx_norm_freqMod(cut_type==2),'dk','markerfacecolor','c')

plot(S2_ratio(cut_type==1),Fy_norm_freqMod(cut_type==1),'ok','markerfacecolor','m')
plot(S2_ratio(cut_type==2),Fy_norm_freqMod(cut_type==2),'dk','markerfacecolor','m')

plot(S2_ratio(cut_type==1),Fz_norm_freqMod(cut_type==1),'ok','markerfacecolor','y')
plot(S2_ratio(cut_type==2),Fz_norm_freqMod(cut_type==2),'dk','markerfacecolor','y')

% linear fits
plot([min(S2_ratio) max(S2_ratio)],polyval(Fx_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fy_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fz_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')

plot([min(S2_ratio) max(S2_ratio)],polyval(Fx_S2_fit_freqMod,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fy_S2_fit_freqMod,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fz_S2_fit_freqMod,[min(S2_ratio) max(S2_ratio)]),'k')

legend('Fx steady','Fy steady','Fz steady','Fx freqMod','Fy freqMod','Fz freqMod','location','SW')
xlabel('S2 ratio')
ylabel('F/mg')
axis([0 1 -1.5 .25])

% M-S3
subplot(1,2,2)
hold on

% steady freq
plot(S3_ratio(cut_type==0),Mx_norm(cut_type==0)-Mx_norm(cut_type==0),'sk','markerfacecolor','b')
plot(S3_ratio(cut_type==1),Mx_norm(cut_type==1)-Mx_norm(cut_type==0),'ok','markerfacecolor','b')
plot(S3_ratio(cut_type==2),Mx_norm(cut_type==2)-Mx_norm(cut_type==0),'dk','markerfacecolor','b')

plot(S3_ratio(cut_type==0),My_norm_CoM(cut_type==0),'sk','markerfacecolor','r')
plot(S3_ratio(cut_type==1),My_norm_CoM(cut_type==1),'ok','markerfacecolor','r')
plot(S3_ratio(cut_type==2),My_norm_CoM(cut_type==2),'dk','markerfacecolor','r')

plot(S3_ratio(cut_type==0),Mz_norm(cut_type==0)-Mz_norm(cut_type==0),'sk','markerfacecolor','g')
plot(S3_ratio(cut_type==1),Mz_norm(cut_type==1)-Mz_norm(cut_type==0),'ok','markerfacecolor','g')
plot(S3_ratio(cut_type==2),Mz_norm(cut_type==2)-Mz_norm(cut_type==0),'dk','markerfacecolor','g')

plot([min(S3_ratio) max(S3_ratio)],polyval(MxMinSteady_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(My_CoM_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(MzMinSteady_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')

% clipped fly freq
plot(S3_ratio(cut_type==0),Mx_norm_freqMod(cut_type==0)-Mx_norm_freqMod(cut_type==0),'sk','markerfacecolor','c')
plot(S3_ratio(cut_type==1),Mx_norm_freqMod(cut_type==1)-Mx_norm_freqMod(cut_type==0),'ok','markerfacecolor','c')
plot(S3_ratio(cut_type==2),Mx_norm_freqMod(cut_type==2)-Mx_norm_freqMod(cut_type==0),'dk','markerfacecolor','c')

plot(S3_ratio(cut_type==0),My_norm_freqMod_CoM(cut_type==0),'sk','markerfacecolor','m')
plot(S3_ratio(cut_type==1),My_norm_freqMod_CoM(cut_type==1),'ok','markerfacecolor','m')
plot(S3_ratio(cut_type==2),My_norm_freqMod_CoM(cut_type==2),'dk','markerfacecolor','m')

plot(S3_ratio(cut_type==0),Mz_norm_freqMod(cut_type==0)-Mz_norm_freqMod(cut_type==0),'sk','markerfacecolor','y')
plot(S3_ratio(cut_type==1),Mz_norm_freqMod(cut_type==1)-Mz_norm_freqMod(cut_type==0),'ok','markerfacecolor','y')
plot(S3_ratio(cut_type==2),Mz_norm_freqMod(cut_type==2)-Mz_norm_freqMod(cut_type==0),'dk','markerfacecolor','y')

plot([min(S3_ratio) max(S3_ratio)],polyval(MxMinSteady_S3_fit_freqMod,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(My_CoM_S3_fit_freqMod,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(MzMinSteady_S3_fit_freqMod,[min(S3_ratio) max(S3_ratio)]),'k')

legend('Intact wing','Tip cut','Trailing Edge')
xlabel('S3 ratio')
ylabel('T/mgl')
axis([0 1 -.1 .5])
