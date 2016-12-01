clear
clc
close all

mkdir('MODsNstats_bodyNfreqNvel_indiv')

Eqname=dir('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq*')
Eqname=Eqname.name;
load(Eqname)

loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

%% S2S3 Amplitude Ratio func settings
for i = 1: length(SecondMomentRatio)
%     counter = length(SecondMomentRatio)-i
    
    sol = subs(solAdAiRatio,S2,SecondMomentRatio(i));
    sol = subs(sol,S3,ThirdMomentRatio(i));
    S2S3AmpRatioFunc(i,1) = eval(sol);
end

%% intact flies
f_WBsteady_intact = f_wb_steady;
pitch_WBsteady_intact = pitch_global_steady;
V_WBsteady_intact = V_steady;

Aplus_WBsteady_intact = ones(size(f_WBsteady_intact));

%% steady wb's only
clip_type_WBsteady = clip_type(steady_nr_mean_wb==1);

Aplus_WBsteady = S2S3AmpRatioFunc(steady_nr_mean_wb==1);
RS2_WBsteady = SecondMomentRatio(steady_nr_mean_wb==1);
RS3_WBsteady = ThirdMomentRatio(steady_nr_mean_wb==1);

f_wb = mean([f_wb_L,f_wb_R]')';
f_WBsteady = f_wb(steady_nr_mean_wb==1);
pitch_WBsteady = pitch_mean_wb(steady_nr_mean_wb==1);
roll_WBsteady = roll_mean_wb(steady_nr_mean_wb==1);
slip_WBsteady = slip_mean_wb(steady_nr_mean_wb==1);
V_WBsteady = V_mean_wb;

%% individual data
% intact flies
seq_intact_unique = unique(seq_nr_steady);
for i = 1:length(seq_intact_unique)
    n_now = find(seq_nr_steady==seq_intact_unique(i));
    
    f_indiv_intact(i,1) = mean(f_WBsteady_intact(n_now));
    pitch_indiv_intact(i,1) = mean(pitch_WBsteady_intact(n_now));
    V_indiv_intact(i,1) = mean(V_WBsteady_intact(n_now));
end
Aplus_indiv_intact = ones(size(f_indiv_intact));
RS2_indiv_intact = ones(size(f_indiv_intact));
RS3_indiv_intact = ones(size(f_indiv_intact));
clip_type_indiv_intact = zeros(size(f_indiv_intact));
roll_indiv_intact = zeros(size(f_indiv_intact));
slip_indiv_intact = zeros(size(f_indiv_intact));

% damaged flies
RS2_WBsteady_unique = unique(RS2_WBsteady);
for i = 1:length(RS2_WBsteady_unique)
    n_now = find(RS2_WBsteady==RS2_WBsteady_unique(i));
    n_indiv(i,1) = length(n_now);
    
    clip_type_indiv(i,1) = mean(clip_type_WBsteady(n_now));
    
    Aplus_indiv(i,1) = mean(Aplus_WBsteady(n_now));
    RS2_indiv(i,1) = mean(RS2_WBsteady(n_now));
    RS3_indiv(i,1) = mean(RS3_WBsteady(n_now));
    
    f_indiv(i,1) = mean(f_WBsteady(n_now));
    pitch_indiv(i,1) = mean(pitch_WBsteady(n_now));
    roll_indiv(i,1) = mean(roll_WBsteady(n_now));
    slip_indiv(i,1) = mean(slip_WBsteady(n_now));
    V_indiv(i,1) = mean(V_WBsteady(n_now));
end

%% wb based statistics
% V_all = [V_WBsteady_intact;V_WBsteady];
% V_groups(1:length(V_WBsteady_intact),1)=1;
% V_groups(length(V_WBsteady_intact)+1:length(V_all),1)=2;
% [p_V_all,tbl_V_all,stats_V_all] = kruskalwallis(V_all,V_groups);
% 
% f_all = [f_WBsteady_intact;f_WBsteady];
% f_groups(1:length(f_WBsteady_intact),1)=1;
% f_groups(length(f_WBsteady_intact)+1:length(f_all),1)=2;
% [p_freq_all,tbl_freq_all,stats_freq_all] = kruskalwallis(f_all,f_groups);
% 
% [pval_pitch_all, med_pitch_all, P_pitch_all]=circ_cmtest(deg2rad(pitch_WBsteady_intact),deg2rad(pitch_WBsteady));
% % figure
% % rose(deg2rad(pitch_WBsteady_intact),360)
% % hold on
% % h=rose(deg2rad(pitch_WBsteady),360)
% % set(h,'color','r')

%% individual based statistics 
Aplus_indiv_all = [Aplus_indiv_intact;Aplus_indiv];
f_indiv_all = [f_indiv_intact;f_indiv];
pitch_indiv_all = [pitch_indiv_intact;pitch_indiv];
roll_indiv_all = [roll_indiv_intact;roll_indiv];
slip_indiv_all = [slip_indiv_intact;slip_indiv];
V_indiv_all = [V_indiv_intact;V_indiv];

%% V stats
V_indiv_groups(1:length(V_indiv_intact),1)=1;
V_indiv_groups(length(V_indiv_intact)+1:length(V_indiv_all),1)=2;
[kruskalwallis_test.p_V_indiv,kruskalwallis_test.tbl_V_indiv,kruskalwallis_test.stats_V_indiv] = kruskalwallis(V_indiv_all,V_indiv_groups);
[ranksum_test.p_V_indiv,ranksum_test.h_V_indiv,ranksum_test.stats_V_indiv] = ranksum(V_indiv_intact,V_indiv);

%% freq stats
f_indiv_groups(1:length(f_indiv_intact),1)=1;
f_indiv_groups(length(f_indiv_intact)+1:length(f_indiv_all),1)=2;
[kruskalwallis_test.p_freq_indiv,kruskalwallis_test.tbl_freq_indiv,kruskalwallis_test.stats_freq_indiv] = kruskalwallis(f_indiv_all,f_indiv_groups);
[ranksum_test.p_freq_indiv,ranksum_test.h_freq_indiv,ranksum_test.stats_freq_indiv] = ranksum(f_indiv_intact,f_indiv);

[freq_S2S3AmpRatioFunc_indiv_asympFit1, freq_S2S3AmpRatioFunc_indiv_asympFit1_gof]      = freq_S2S3AmpRatio_asymp_fit1_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit2, freq_S2S3AmpRatioFunc_indiv_asympFit2_gof]      = freq_S2S3AmpRatio_asymp_fit2_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit3, freq_S2S3AmpRatioFunc_indiv_asympFit3_gof]      = freq_S2S3AmpRatio_asymp_fit3_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit4, freq_S2S3AmpRatioFunc_indiv_asympFit4_gof]      = freq_S2S3AmpRatio_asymp_fit4_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit5, freq_S2S3AmpRatioFunc_indiv_asympFit5_gof]      = freq_S2S3AmpRatio_asymp_fit5_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit9, freq_S2S3AmpRatioFunc_indiv_asympFit9_gof]      = freq_S2S3AmpRatio_asymp_fit9_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit10, freq_S2S3AmpRatioFunc_indiv_asympFit10_gof]    = freq_S2S3AmpRatio_asymp_fit10_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit11, freq_S2S3AmpRatioFunc_indiv_asympFit11_gof]    = freq_S2S3AmpRatio_asymp_fit11_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit15, freq_S2S3AmpRatioFunc_indiv_asympFit15_gof]    = freq_S2S3AmpRatio_asymp_fit15_fsteady(Aplus_indiv_all , f_indiv_all, 0);

[freq_S2S3AmpRatioFunc_indiv_linear, freq_S2S3AmpRatioFunc_indiv_linear_gof]      = fit(Aplus_indiv_all , f_indiv_all, 'poly1');
[freq_S2S3AmpRatioFunc_indiv_power2, freq_S2S3AmpRatioFunc_indiv_power2_gof]      = fit(Aplus_indiv_all , f_indiv_all, 'power2');
[freq_S2S3AmpRatioFunc_indiv_smooth999, freq_S2S3AmpRatioFunc_indiv_smooth999_gof]= fit(Aplus_indiv_all , f_indiv_all, 'smoothingspline','SmoothingParam',.999);

%% pitch stats
[pitch_indiv_intact_mean pitch_indiv_intact_ul pitch_indiv_intact_ll]=circ_mean_deg_nonan(pitch_indiv_intact);
[pitch_indiv_mean pitch_indiv_ul pitch_indiv_ll]=circ_mean_deg_nonan(pitch_indiv);

[pitch_indiv_intact_med]=rad2deg(circ_median(deg2rad(pitch_indiv_intact)));
[pitch_indiv_med]=rad2deg(circ_median(deg2rad(pitch_indiv)));

[pval_pitch_indiv med_pitch_indiv P_pitch_indiv]=circ_cmtest(deg2rad(pitch_indiv_intact),deg2rad(pitch_indiv));

%% roll stats
[roll_indiv_mean roll_indiv_ul roll_indiv_ll]=circ_mean_deg_nonan(roll_indiv);
[roll_indiv_med]=rad2deg(circ_median(deg2rad(roll_indiv)));

[h_roll_indiv mu_roll_indiv ul_roll_indiv ll_roll_indiv] = circ_mtest(deg2rad(roll_indiv),0);
mu_roll_indiv = rad2deg(mu_roll_indiv)
ul_roll_indiv = rad2deg(ul_roll_indiv)
ll_roll_indiv = rad2deg(ll_roll_indiv)

[roll_S2S3AmpRatioFunc_indiv_asympFit1, roll_S2S3AmpRatioFunc_indiv_asympFit1_gof]      = roll_S2S3AmpRatio_asymp_fit1_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit2, roll_S2S3AmpRatioFunc_indiv_asympFit2_gof]      = roll_S2S3AmpRatio_asymp_fit2_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit3, roll_S2S3AmpRatioFunc_indiv_asympFit3_gof]      = roll_S2S3AmpRatio_asymp_fit3_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit4, roll_S2S3AmpRatioFunc_indiv_asympFit4_gof]      = roll_S2S3AmpRatio_asymp_fit4_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit5, roll_S2S3AmpRatioFunc_indiv_asympFit5_gof]      = roll_S2S3AmpRatio_asymp_fit5_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit9, roll_S2S3AmpRatioFunc_indiv_asympFit9_gof]      = roll_S2S3AmpRatio_asymp_fit9_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit10, roll_S2S3AmpRatioFunc_indiv_asympFit10_gof]    = roll_S2S3AmpRatio_asymp_fit10_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit11, roll_S2S3AmpRatioFunc_indiv_asympFit11_gof]    = roll_S2S3AmpRatio_asymp_fit11_fsteady(Aplus_indiv_all , roll_indiv_all, 0);
[roll_S2S3AmpRatioFunc_indiv_asympFit15, roll_S2S3AmpRatioFunc_indiv_asympFit15_gof]    = roll_S2S3AmpRatio_asymp_fit15_fsteady(Aplus_indiv_all , roll_indiv_all, 0);

[roll_S2S3AmpRatioFunc_indiv_linear, roll_S2S3AmpRatioFunc_indiv_linear_gof]      = fit(Aplus_indiv_all , roll_indiv_all, 'poly1');
[roll_S2S3AmpRatioFunc_indiv_power2, roll_S2S3AmpRatioFunc_indiv_power2_gof]      = fit(Aplus_indiv_all , roll_indiv_all, 'power2');
[roll_S2S3AmpRatioFunc_indiv_smooth999, roll_S2S3AmpRatioFunc_indiv_smooth999_gof]= fit(Aplus_indiv_all , roll_indiv_all, 'smoothingspline','SmoothingParam',.999);

%% slip stats
[slip_indiv_mean slip_indiv_ul slip_indiv_ll]=circ_mean_deg_nonan(slip_indiv);
[slip_indiv_med]=rad2deg(circ_median(deg2rad(slip_indiv)));

[h_slip_indiv mu_slip_indiv ul_slip_indiv ll_slip_indiv] = circ_mtest(deg2rad(slip_indiv),0);
mu_slip_indiv = rad2deg(mu_slip_indiv)
ul_slip_indiv = rad2deg(ul_slip_indiv)
ll_slip_indiv = rad2deg(ll_slip_indiv)

%% save data
save_freqNvelNbodyKin_statsNmods

%% figures
cd('MODsNstats_bodyNfreqNvel_indiv')

%% freq stats fig
figure(2)
saveas(gca,['freq_indiv_stats.fig'])
saveas(gca,['freq_indiv_stats.png'])
plot2svg(['freq_indiv_stats.svg'])

figure(4)
saveas(gca,['vel_indiv_stats.fig'])
saveas(gca,['vel_indiv_stats.png'])
plot2svg(['vel_indiv_stats.svg'])

%% freq fits fig
figure
subplot(2,2,1)
plot(freq_S2S3AmpRatioFunc_indiv_asympFit10,Aplus_indiv_all, f_indiv_all)
legend off
xlabel('A+')
ylabel('freq')
title('asymp10')

subplot(2,2,2)
plot(freq_S2S3AmpRatioFunc_indiv_linear,Aplus_indiv_all, f_indiv_all)
legend off
xlabel('A+')
ylabel('freq')
title('linear')

subplot(2,2,3)
plot(freq_S2S3AmpRatioFunc_indiv_power2,Aplus_indiv_all, f_indiv_all)
legend off
xlabel('A+')
ylabel('freq')
title('power2')

subplot(2,2,4)
plot(freq_S2S3AmpRatioFunc_indiv_smooth999,Aplus_indiv_all, f_indiv_all)
legend off
xlabel('A+')
ylabel('freq')
title('smooth999')

saveas(gca,['freq_indiv_fits.fig'])
saveas(gca,['freq_indiv_fits.png'])
plot2svg(['freq_indiv_fits.svg'])

%% MSfig freq fit
xmin = 1;
xmax = 1.3;
y_min = 150;
y_max = 250;
N_min = 1;
N_max = 50;

cmap_Nwb=cbrewer('seq','Blues',N_max);

figure(10)
subplot(2,2,3)
hold on

plot([xmin:(xmax-xmin)/99:xmax],feval(freq_S2S3AmpRatioFunc_indiv_asympFit10,[xmin:(xmax-xmin)/99:xmax]),'-k','linewidth',3)
plot([xmin:(xmax-xmin)/99:xmax],feval(freq_S2S3AmpRatioFunc_indiv_smooth999,[xmin:(xmax-xmin)/99:xmax]),'--b','linewidth',3)

% errorbar(1, mean(f_indiv_intact),std(f_indiv_intact)/sqrt(length(f_indiv_intact)),'sk','markersize',7,'markerfacecolor','w')
errorbar(1, mean(f_indiv_intact),std(f_indiv_intact),'sk','markersize',7,'markerfacecolor','w')
for i = 1:length(f_indiv)
    N_now = n_indiv(i);
    color_nr = round(24/(N_max-N_min)*(N_now-N_min)+1);
    if color_nr > N_max
        color_nr = N_max;
    end

    if clip_type_indiv(i) == 1
        plot(Aplus_indiv(i), f_indiv(i),'ok','markersize',7,'markerfacecolor',cmap_Nwb(color_nr,:))
    else
        plot(Aplus_indiv(i), f_indiv(i),'dk','markersize',7,'markerfacecolor',cmap_Nwb(color_nr,:))
    end
end

%     legend('NONclipped steady WB','Clipped WBs','asymptotic fit','location','SE')
axis([xmin xmax y_min y_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',y_min:(y_max-y_min)/2:y_max) 
xlabel('A+','fontsize',10)
ylabel('frequency','fontsize',10)
axis square

subplot(2,2,1)
colormap(cmap_Nwb)
caxis([N_min N_max])
h = colorbar('location','northoutside'); 
title(h,'#wingbeats')
set(h,'xtick',N_min:(N_max-N_min)/2:N_max)

saveas(gca,['WBmod_freq_Aplus_indiv_spline999_asympFit10_MSfig.fig'])
saveas(gca,['WBmod_freq_Aplus_indiv_spline999_asympFit10_MSfig.png'])
plot2svg(['WBmod_freq_Aplus_indiv_spline999_asympFit10_MSfig.svg'])
    
%% roll stats fig
figure
subplot(2,2,1)
plot(roll_S2S3AmpRatioFunc_indiv_asympFit10,Aplus_indiv_all, roll_indiv_all)
legend off
xlabel('A+')
ylabel('roll')
title('asymp10')

subplot(2,2,2)
plot(roll_S2S3AmpRatioFunc_indiv_linear,Aplus_indiv_all, roll_indiv_all)
legend off
xlabel('A+')
ylabel('roll')
title('linear')

subplot(2,2,3)
plot(roll_S2S3AmpRatioFunc_indiv_power2,Aplus_indiv_all, roll_indiv_all)
legend off
xlabel('A+')
ylabel('roll')
title('power2')

subplot(2,2,4)
plot(roll_S2S3AmpRatioFunc_indiv_smooth999,Aplus_indiv_all, roll_indiv_all)
legend off
xlabel('A+')
ylabel('roll')
title('smooth999')

saveas(gca,['roll_indiv_fits.fig'])
saveas(gca,['roll_indiv_fits.png'])
plot2svg(['roll_indiv_fits.svg'])

%% MSfig freq fit
xmin = 1;
xmax = 1.3;
y_min = -15;
y_max = 30;
N_min = 1;
N_max = 50;

cmap_Nwb=cbrewer('seq','Blues',N_max);

figure(10)
subplot(2,2,4)
hold on
plot([xmin:(xmax-xmin)/99:xmax],feval(roll_S2S3AmpRatioFunc_indiv_asympFit10,[xmin:(xmax-xmin)/99:xmax]),'-k','linewidth',3)
plot([xmin:(xmax-xmin)/99:xmax],feval(roll_S2S3AmpRatioFunc_indiv_smooth999,[xmin:(xmax-xmin)/99:xmax]),'--b','linewidth',3)
plot(1, 0,'sk','markersize',7,'markerfacecolor','w')

for i = 1:length(roll_indiv)
    N_now = n_indiv(i);
    color_nr = round(24/(N_max-N_min)*(N_now-N_min)+1);
    if color_nr > N_max
        color_nr = N_max;
    end

    if clip_type_indiv(i) == 1
        plot(Aplus_indiv(i), roll_indiv(i),'ok','markersize',7,'markerfacecolor',cmap_Nwb(color_nr,:))
    else
        plot(Aplus_indiv(i), roll_indiv(i),'dk','markersize',7,'markerfacecolor',cmap_Nwb(color_nr,:))
    end
end

%     legend('NONclipped steady WB','Clipped WBs','asymptotic fit','location','SE')
axis([xmin xmax y_min y_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',y_min:(y_max-y_min)/3:y_max) 
xlabel('A+','fontsize',10)
ylabel('roll angle','fontsize',10)
axis square

subplot(2,2,2)
colormap(cmap_Nwb)
caxis([N_min N_max])
h = colorbar('location','northoutside'); 
title(h,'#wingbeats')
set(h,'xtick',N_min:(N_max-N_min)/2:N_max)

saveas(gca,['WBmod_freq_Aplus_indiv_asympFit_order10_MSfig.fig'])
saveas(gca,['WBmod_freq_Aplus_indiv_asympFit_order10_MSfig.png'])
plot2svg(['WBmod_freq_Aplus_indiv_asympFit_order10_MSfig.svg'])

%% rose plots of body orientation
% plot intact pitch
figure

angles = deg2rad(pitch_indiv_intact);
angle_mean = deg2rad(pitch_indiv_intact_mean);
angle_ll = deg2rad(pitch_indiv_intact_ll);
angle_ul = deg2rad(pitch_indiv_intact_ul);

c_now = [0 0 1];

subplot(2,2,1)
h=rose(angles,180);
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'b') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot damaged pitch
angles = deg2rad(pitch_indiv);
angle_mean = deg2rad(pitch_indiv_mean);
angle_ll = deg2rad(pitch_indiv_ll);
angle_ul = deg2rad(pitch_indiv_ul);

c_now = [1 0 0];

subplot(2,2,1)
h=rose(angles,180);
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot roll
angles = deg2rad(roll_indiv);
angle_mean = deg2rad(roll_indiv_mean);
angle_ll = deg2rad(roll_indiv_ll);
angle_ul = deg2rad(roll_indiv_ul);

c_now = [1 0 0];

subplot(2,2,2)
h=rose(angles,90);
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot slip
angles = deg2rad(slip_indiv);
angle_mean = deg2rad(slip_indiv_mean);
angle_ll = deg2rad(slip_indiv_ll);
angle_ul = deg2rad(slip_indiv_ul);

c_now = [1 0 0];

subplot(2,2,3)
h=rose(angles,90);
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% save rose plots
saveas(gca,['BodyKin_indiv_MSfig.fig'])
saveas(gca,['BodyKin_indiv_MSfig.png'])
plot2svg(['BodyKin_indiv_MSfig.svg'])

cd ..







