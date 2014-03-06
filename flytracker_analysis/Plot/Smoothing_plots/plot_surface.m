
addpath('/home/matt/Dropbox/flytracker_analysis/');
addpath('/home/matt/Dropbox/flytracker_analysis/Plot');

i = 104;
pathDB.qL1(1:200,i) = nan;
pathDB.qL2(1:200,i) = nan;
pathDB.qL3(1:200,i) =nan;
pathDB.qL4(1:200,i) = nan;
pathDB.qR1(1:200,i) = nan;
pathDB.qR2(1:200,i) = nan;
pathDB.qR3(1:200,i) = nan;
pathDB.qR4(1:200,i) =nan;
pathDB.qL1_filt2(1:200,i) = nan;
pathDB.qL2_filt2(1:200,i) = nan;
pathDB.qL3_filt2(1:200,i) = nan;
pathDB.qL4_filt2(1:200,i) = nan;
pathDB.qR1_filt2(1:200,i) = nan;
pathDB.qR2_filt2(1:200,i) = nan;
pathDB.qR3_filt2(1:200,i) = nan;
pathDB.qR4_filt2(1:200,i) = nan;
% pathDB.wing_l(1:80,i) = 0;



for i = 104

settings = settings;
seq_nr = i;
qL1_r = pathDB.qL1(:,i);
qL2_r = pathDB.qL2(:,i);
qL3_r = pathDB.qL3(:,i);
qL4_r = pathDB.qL4(:,i);
qR1_r = pathDB.qR1(:,i);
qR2_r = pathDB.qR2(:,i);
qR3_r = pathDB.qR3(:,i);
qR4_r = pathDB.qR4(:,i);
qL1_f2 = pathDB.qL1_filt2(:,i);
qL2_f2 = pathDB.qL2_filt2(:,i);
qL3_f2 = pathDB.qL3_filt2(:,i);
qL4_f2 = pathDB.qL4_filt2(:,i);
qR1_f2 = pathDB.qR1_filt2(:,i);
qR2_f2 = pathDB.qR2_filt2(:,i);
qR3_f2 = pathDB.qR3_filt2(:,i);
qR4_f2 = pathDB.qR4_filt2(:,i);
wing_l = pathDB.wing_l(:,i);
fig_nr1 = i;
save_on_off = 0;


Sphere_plot_raw_filt2(settings,seq_nr,qL1_r,qL2_r,qL3_r,qL4_r,qR1_r,qR2_r,qR3_r,qR4_r,qL1_f2,qL2_f2,qL3_f2,qL4_f2,qR1_f2,qR2_f2,qR3_f2,qR4_f2,wing_l,fig_nr1,save_on_off)

end
% 
% % figure(1)
% % hold on
% % plot(pathDB.theta_R(:,1))
% % plot(pathDB.theta_L(:,1),'-r')
% % hold off
% % 
% figure(103)
% hold on
% plot(phi_R(:,110),theta_R(:,110))
% plot(phi_L(:,110),theta_L(:,110),'-r')
% %plot(pathDB.phi_L(:,101),'-r')
% hold off
% axis equal
% 
% plot(pathDB.phi_R(:,1))


figure(1)
grid off; box on; hold on
title('steady flight')
ylabel('stroke deviation (deg)')
xlabel('stroke angle (deg)')
% xlim([-80 120]); ylim([-30 50]);
axis equal
% set(gca,'XTick',[-20 60],'XTickLabel',{'down stroke';'up stroke'},'TickLength',[0 0],'LineWidth', 1.2);
% x = [1:200];
for i = 104
    for j = [3:4,14:length(stroke_wb_R_MATT_bins(1,:,i))]
        
plot(stroke_wb_R_MATT_bins(:,j,i),dev_wb_R_MATT_bins(:,j,i))
plot(stroke_wb_L_MATT_bins(:,j,i),dev_wb_L_MATT_bins(:,j,i),'-r')
% pause
% display(j)
%plot(pathDB.phi_L(:,101),'-r')

% x = [1:200]' + max(x);
    
    end
end

hold off

% print -dpng -r600 turning_stroke_v_dev.png

% figure(3)
% plot(pathDB.theta_R(:,101)-pathDB.theta_R(:,101))

clear wing wing_l x i j
clear qL1_f2 qL1_r qL2_f2 qL2_r qL3_f2 qL3_r qL4_f2 qL4_r qR1_f2 qR1_r qR2_f2 qR2_r qR3_f2 qR3_r qR4_f2 qR4_r