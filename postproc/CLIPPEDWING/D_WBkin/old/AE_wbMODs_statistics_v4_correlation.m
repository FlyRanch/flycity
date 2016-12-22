clear
clc

mkdir('MODsNstats_bodyNfreq_indiv')

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

%% individual data
% intact flies
seq_intact_unique = unique(seq_nr_steady);
for i = 1:length(seq_intact_unique)
    n_now = find(seq_nr_steady==seq_intact_unique(i));
    
    f_indiv_intact(i,1) = mean(f_WBsteady_intact(n_now));
    pitch_indiv_intact(i,1) = mean(pitch_WBsteady_intact(n_now));
end
Aplus_indiv_intact = ones(size(f_indiv_intact));
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
end

%% wb based statistics
f_all = [f_WBsteady_intact;f_WBsteady];
f_groups(1:length(f_WBsteady_intact),1)=1;
f_groups(length(f_WBsteady_intact)+1:length(f_all),1)=2;
[p_freq_all,tbl_freq_all,stats_freq_all] = kruskalwallis(f_all,f_groups)

[pval_pitch_all, med_pitch_all, P_pitch_all]=circ_cmtest(deg2rad(pitch_WBsteady_intact),deg2rad(pitch_WBsteady))
% figure
% rose(deg2rad(pitch_WBsteady_intact),360)
% hold on
% h=rose(deg2rad(pitch_WBsteady),360)
% set(h,'color','r')

%% individual based statistics 
Aplus_indiv_all = [Aplus_indiv_intact;Aplus_indiv];
f_indiv_all = [f_indiv_intact;f_indiv];
pitch_indiv_all = [pitch_indiv_intact;pitch_indiv];
roll_indiv_all = [roll_indiv_intact;roll_indiv];
slip_indiv_all = [slip_indiv_intact;slip_indiv];

% freq
f_indiv_groups(1:length(f_indiv_intact),1)=1;
f_indiv_groups(length(f_indiv_intact)+1:length(f_indiv_all),1)=2;
[p_freq_indiv,tbl_freq_indiv,stats_freq_indiv] = kruskalwallis(f_indiv_all,f_indiv_groups)

[freq_S2S3AmpRatioFunc_indiv_asympFit1, freq_S2S3AmpRatioFunc_indiv_asympFit1_gof]      = freq_S2S3AmpRatio_asymp_fit1_fsteady(Aplus_indiv_all , f_indiv_all, plot_on);
[freq_S2S3AmpRatioFunc_indiv_asympFit2, freq_S2S3AmpRatioFunc_indiv_asympFit2_gof]      = freq_S2S3AmpRatio_asymp_fit2_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit3, freq_S2S3AmpRatioFunc_indiv_asympFit3_gof]      = freq_S2S3AmpRatio_asymp_fit3_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit4, freq_S2S3AmpRatioFunc_indiv_asympFit4_gof]      = freq_S2S3AmpRatio_asymp_fit4_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit5, freq_S2S3AmpRatioFunc_indiv_asympFit5_gof]      = freq_S2S3AmpRatio_asymp_fit5_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit9, freq_S2S3AmpRatioFunc_indiv_asympFit9_gof]      = freq_S2S3AmpRatio_asymp_fit9_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit10, freq_S2S3AmpRatioFunc_indiv_asympFit10_gof]    = freq_S2S3AmpRatio_asymp_fit10_fsteady(Aplus_indiv_all , f_indiv_all, plot_on);
[freq_S2S3AmpRatioFunc_indiv_asympFit11, freq_S2S3AmpRatioFunc_indiv_asympFit11_gof]    = freq_S2S3AmpRatio_asymp_fit11_fsteady(Aplus_indiv_all , f_indiv_all, 0);
[freq_S2S3AmpRatioFunc_indiv_asympFit15, freq_S2S3AmpRatioFunc_indiv_asympFit15_gof]    = freq_S2S3AmpRatio_asymp_fit15_fsteady(Aplus_indiv_all , f_indiv_all, 0);

[freq_S2S3AmpRatioFunc_indiv_linear, freq_S2S3AmpRatioFunc_indiv_linear_gof]      = fit(Aplus_indiv_all , f_indiv_all, 'poly1');
[freq_S2S3AmpRatioFunc_indiv_power2, freq_S2S3AmpRatioFunc_indiv_power2_gof]      = fit(Aplus_indiv_all , f_indiv_all, 'power2');
[freq_S2S3AmpRatioFunc_indiv_smooth999, freq_S2S3AmpRatioFunc_indiv_smooth999_gof]      = fit(Aplus_indiv_all , f_indiv_all, 'smoothingspline','SmoothingParam',.999);

figure
subplot(2,2,1)
plot(freq_S2S3AmpRatioFunc_indiv_asympFit10,Aplus_indiv_all, f_indiv_all)
subplot(2,2,2)
plot(freq_S2S3AmpRatioFunc_indiv_linear,Aplus_indiv_all, f_indiv_all)
subplot(2,2,3)
plot(freq_S2S3AmpRatioFunc_indiv_power2,Aplus_indiv_all, f_indiv_all)
subplot(2,2,4)
plot(freq_S2S3AmpRatioFunc_indiv_smooth999,Aplus_indiv_all, f_indiv_all)

% pitch
[pitch_indiv_intact_mean pitch_intact_ul pitch_intact_ll]=circ_mean_deg_nonan(pitch_indiv_intact)
[pitch_indiv_mean pitch_indiv_ul pitch_indiv_ll]=circ_mean_deg_nonan(pitch_indiv)

[pitch_indiv_intact_med]=rad2deg(circ_median(deg2rad(pitch_indiv_intact)))
[pitch_indiv_med]=rad2deg(circ_median(deg2rad(pitch_indiv)))

[pval_pitch_indiv med_pitch_indiv P_pitch_indiv]=circ_cmtest(deg2rad(pitch_indiv_intact),deg2rad(pitch_indiv))

% roll
[roll_indiv_mean roll_indiv_ul roll_indiv_ll]=circ_mean_deg_nonan(roll_indiv)
pval_roll_indiv = circ_medtest(deg2rad(roll_indiv),0)
[h_roll_indiv mu_roll_indiv ul_roll_indiv ll_roll_indiv] = circ_mtest(deg2rad(roll_indiv),0)
figure
rose(deg2rad(roll_indiv),360)

% slip
[slip_indiv_mean slip_indiv_ul slip_indiv_ll]=circ_mean_deg_nonan(slip_indiv)
pval_slip_indiv=circ_medtest(deg2rad(slip_indiv),0)
[h_slip_indiv mu_slip_indiv ul_slip_indiv ll_slip_indiv] = circ_mtest(deg2rad(slip_indiv),0)

%% plot intact pitch
figure

angles = deg2rad(pitch_indiv_intact);
angle_mean = deg2rad(pitch_indiv_intact_mean);
angle_ll = deg2rad(pitch_intact_ll);
angle_ul = deg2rad(pitch_intact_ul);

c_now = [0 0 1];

subplot(2,2,1)
h=rose(angles,180)
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
h=rose(angles,180)
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
h=rose(angles,90)
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
h=rose(angles,90)
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
cd('MODsNstats_bodyNfreq_indiv')

saveas(gca,['BodyKin_indiv_MSfig.fig'])
saveas(gca,['BodyKin_indiv_MSfig.png'])
plot2svg(['BodyKin_indiv_MSfig.svg'])

cd ..







