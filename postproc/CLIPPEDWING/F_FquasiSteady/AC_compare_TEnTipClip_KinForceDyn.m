clear;
clc;
close all
warning off

%% plot all results
freq_asymFitNr = 10;
figure(2)
subplot(2,2,2)
hold on
subplot(2,2,4)
hold on

load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])
subplot(2,2,2)
plot(S2ratios,-Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','w')
subplot(2,2,4)
plot(S3ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','w')

load(['allMODs_TEClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])
subplot(2,2,2)
plot(S2ratios,-Fz_damaged_mean_all,'d-k','markersize',10,'markerfacecolor','w')
subplot(2,2,4)
plot(S3ratios,Mx_damaged_mean_all,'d-k','markersize',10,'markerfacecolor','w')

%% do analysis
clear;
load('bodyNwingModel_4qsModel.mat')

loadname=dir('WBdataset_steadyNclipMods_S2S3AmpRatioFunc.mat')
loadname = loadname.name;
load(loadname)

plot_on = 1
% plot_on = 0

rot_on=1;
% rot_on=0;

%% constants
freq_asymFitNr = 10;
fit_order = 2;

nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

Mg_fly = body_model.Mg_fly;
l_wing = wing_model.length/1000;
    
t_norm = 0:1/(nr_timepoints-1):1;

%% steady WB
f_steady = f_wb_steady_meanCIstd(1);

% fouriers coeffs
stroke_coeffs_steady = stroke_steady_fourier_coeffs_binmean;
dev_coeffs_steady = dev_steady_fourier_coeffs_binmean;
rot_coeffs_steady = pitch_steady_fourier_coeffs_binmean;

[stroke_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))';
[dev_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))';
[rot_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))';

%% damaged fly WB kin MODs
f_clipped_fly = f_wb_S2S3AmpRatioFunc_meanCIstd(1);
freqMOD = freqMOD_wb_S2S3AmpRatioFunc_meanCIstd(1);

% fouriers coeffs
strokeMOD_coeffs_intact = strokeMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;
devMOD_coeffs_intact = devMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;
rotMOD_coeffs_intact = pitchMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;

[strokeMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_intact,0))';
[devMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_intact,0))';
[rotMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_intact,0))';

strokeMOD_coeffs_damaged = strokeMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;
devMOD_coeffs_damaged = devMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;
rotMOD_coeffs_damaged = pitchMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;

[strokeMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_damaged,0))';
[devMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_damaged,0))';
[rotMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_damaged,0))';


%% TIP CLIP
% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_R';
c_max = max(wing_model.chords_R');

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% loop with different cuts
sect_min = 16;
sect_max = nr_sect;
for i = 1:(sect_max-sect_min+1)

    N_cut = sect_max -i +1
    
    span_ratio = N_cut/nr_sect;
    
    % S2 & S3 of damaged wing
    y_sect = wing_model.y_sect_R(1:N_cut,2);
    c_sect = wing_model.chords_R(1:N_cut)';
    
    S2_sect = c_sect .* h_sect .* y_sect.^2;
    S3_sect = c_sect .* h_sect .* y_sect.^3;

    S2_damaged = sum(S2_sect);
    S3_damaged = sum(S3_sect);
    
    S2ratio = S2_damaged/S2_intact;
    S3ratio = S3_damaged/S3_intact;
    
    % store data
    span_ratios_tip(i,1) = span_ratio;
    S2ratios_tip(i,1) = S2ratio;
    S3ratios_tip(i,1) = S3ratio;
end

% trendline
S2span_fit = polyfit(span_ratios_tip,S2ratios_tip,fit_order)
S3span_fit = polyfit(span_ratios_tip,S3ratios_tip,fit_order)
S2S3_fit_tip = polyfit(S2ratios_tip,S3ratios_tip,fit_order)

% plot
figure(1)
subplot(2,2,1)
hold on
plot(span_ratios_tip,S2ratios_tip,'ok','markersize',10)
plot(span_ratios_tip,polyval(S2span_fit,span_ratios_tip),'-k')

subplot(2,2,2)
hold on
plot(span_ratios_tip,S3ratios_tip,'ok','markersize',10)
plot(span_ratios_tip,polyval(S3span_fit,span_ratios_tip),'-k')

subplot(2,2,3)
hold on
plot(S2ratios_tip,S3ratios_tip,'ok','markersize',10)
plot(S2ratios_tip,polyval(S2S3_fit_tip,S2ratios_tip),'-k')

%% TE clip
% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_R';
c_max = max(wing_model.chords_R');

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% loop with different cuts
TEclip_ratio_min = 0.45;
TEclip_Dratio = -(1-TEclip_ratio_min)/4;
TEclip_ratios = [1:TEclip_Dratio:TEclip_ratio_min]';

for i = 1:length(TEclip_ratios)

    %% S2 & S3 of damaged wing
    TEclip_ratio_now = TEclip_ratios(i)
    c_cut = TEclip_ratio_now * c_max;
    
    c_sect_cut = c_sect;
    c_sect_cut(c_sect_cut>c_cut) = c_cut;
    
    S2_sect = c_sect_cut .* h_sect .* y_sect.^2;
    S3_sect = c_sect_cut .* h_sect .* y_sect.^3;

    S2_damaged = sum(S2_sect);
    S3_damaged = sum(S3_sect);
    
    S2ratio = S2_damaged/S2_intact;
    S3ratio = S3_damaged/S3_intact;
    
    % store data
    S2ratios_TE(i,1) = S2ratio;
    S3ratios_TE(i,1) = S3ratio;
end

% trendline
cord_ratios_TE = TEclip_ratios;
S2cord_fit = polyfit(cord_ratios_TE,S2ratios_TE,fit_order);
S3cord_fit = polyfit(cord_ratios_TE,S3ratios_TE,fit_order);
S2S3_fit_TE = polyfit(S2ratios_TE,S3ratios_TE,fit_order);

cordS2_fit = polyfit(S2ratios_TE,cord_ratios_TE,fit_order);
cordS3_fit = polyfit(S3ratios_TE,cord_ratios_TE,fit_order);

% plot
figure(1)
subplot(2,2,1)
hold on
plot(cord_ratios_TE,S2ratios_TE,'dk','markersize',10)
plot(cord_ratios_TE,polyval(S2cord_fit,cord_ratios_TE),'-k')

xlabel('cut ratio')
ylabel('S2 ratio')
axis equal
axis([0.45 1 0.5 1])
set(gca,'xtick',0:.5:1)
set(gca,'ytick',0:.5:1)

subplot(2,2,2)
hold on
plot(cord_ratios_TE,S3ratios_TE,'dk','markersize',10)
plot(cord_ratios_TE,polyval(S3cord_fit,cord_ratios_TE),'-k')

xlabel('cut ratio')
ylabel('S3 ratio')
axis equal
axis([0.45 1 0.5 1])
set(gca,'xtick',0:.5:1)
set(gca,'ytick',0:.5:1)

subplot(2,2,3)
hold on
plot(S2ratios_TE,S3ratios_TE,'dk','markersize',10)
plot(S2ratios_TE,polyval(S2S3_fit_TE,S2ratios_TE),'-k')

xlabel('S2 ratio')
ylabel('S3 ratio')
axis equal
axis([0.5 1 0.5 1])
set(gca,'xtick',0:.5:1)
set(gca,'ytick',0:.5:1)

%% select cut percentages
N_tip = length(span_ratios_tip)-1;
span_ratio_tip = span_ratios_tip(N_tip);
S2_tip = S2ratios_tip(N_tip);
S3_tip = S3ratios_tip(N_tip);

% cord_ratio_4S2tip = polyval(cordS2_fit,S2_tip)
r = S2cord_fit;
r(end) = r(end)-S2_tip;
rootr = roots(r)
cord_ratio_4S2tip = rootr(end)
S2TE_4S2tip = polyval(S2cord_fit,cord_ratio_4S2tip)
S3TE_4S2tip = polyval(S3cord_fit,cord_ratio_4S2tip)

% cord_ratio_4S3tip = polyval(cordS3_fit,S3_tip)
r = S3cord_fit;
r(end) = r(end)-S3_tip;
rootr = roots(r)
cord_ratio_4S3tip = rootr(end)
S2TE_4S3tip = polyval(S2cord_fit,cord_ratio_4S3tip)
S3TE_4S3tip = polyval(S3cord_fit,cord_ratio_4S3tip)

subplot(2,2,1)
plot([span_ratio_tip cord_ratio_4S2tip],[S2_tip S2TE_4S2tip],'+-r','markersize',10)
    
subplot(2,2,2)
plot([span_ratio_tip cord_ratio_4S3tip],[S3_tip S3TE_4S3tip],'x-r','markersize',10)
    
subplot(2,2,3)
plot([S2_tip S2TE_4S2tip],[S3_tip S3TE_4S2tip],'+-r','markersize',10)
plot([S2_tip S2TE_4S3tip],[S3_tip S3TE_4S3tip],'x-r','markersize',10)

%% TIP CLIP
% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_R';
c_max = max(wing_model.chords_R');

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% tip cut wing
N_cut = nr_sect -N_tip +1;
span_ratio = N_cut/nr_sect;

% S2 & S3 of damaged wing
y_sect = wing_model.y_sect_R(1:N_cut,2);
c_sect = wing_model.chords_R(1:N_cut)';

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_damaged = sum(S2_sect);
S3_damaged = sum(S3_sect);

S2ratio = S2_damaged/S2_intact;
S3ratio = S3_damaged/S3_intact;

S2ratio_tip = S2ratio;
S3ratio_tip = S3ratio;
subplot(2,2,3)
plot(S2ratio_tip,S3ratio_tip,'or','markersize',10)

%% WB kin MOD based on stroke amplitude RATIO
sol = subs(solAdAiRatio,S2,S2ratio);
sol = subs(sol,S3,S3ratio);
S2S3AmpRatioFunc_now = eval(sol);

% asymptotic freq fit
if freq_asymFitNr == 1
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit1,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 2
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit2,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 3
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit3,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 4
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit4,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 5
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit5,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 9
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit9,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 10
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit10,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 11
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit11,S2S3AmpRatioFunc_now);
elseif freq_asymFitNr == 15
    freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit15,S2S3AmpRatioFunc_now);
end        

stroke_intact_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_intact + stroke_steady;    
dev_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_intact    + dev_steady;    
rot_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_intact    + rot_steady;    

stroke_damaged_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_damaged + stroke_steady;    
dev_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_damaged    + dev_steady;    
rot_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_damaged    + rot_steady;    

    %% TIPCUT WB kin NO MODs
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_steady = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_steady = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_steady = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_steady = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_steady = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_steady = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_steady(N_cut+1:end,:) = nan;
    Fy_damaged_steady(N_cut+1:end,:) = nan;
    Fz_damaged_steady(N_cut+1:end,:) = nan;
    
    Mx_damaged_steady(N_cut+1:end,:) = nan;
    My_damaged_steady(N_cut+1:end,:) = nan;
    Mz_damaged_steady(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_steady = (nansum(Fx_damaged_steady,1))' / Mg_fly;
    Fy_damagedwing_steady = (nansum(Fy_damaged_steady,1))' / Mg_fly;
    Fz_damagedwing_steady = (nansum(Fz_damaged_steady,1))' / Mg_fly;
    
    Mx_damagedwing_steady = (nansum(Mx_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_steady = (nansum(My_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_steady = (nansum(Mz_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_steady = (nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_intactwing_steady = (nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_intactwing_steady = (nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_intactwing_steady = (nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_steady = (nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_steady = (nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_steady_tip = (nansum(Fx_damaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_damaged_total_steady_tip = (nansum(Fy_damaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_damaged_total_steady_tip = (nansum(Fz_damaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_damaged_total_steady_tip = (nansum(Mx_damaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_steady_tip = (nansum(My_damaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_steady_tip = (nansum(Mz_damaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_steady_tip = (nansum(Fx_NONdamaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_NONdamaged_total_steady_tip = (nansum(Fy_NONdamaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_NONdamaged_total_steady_tip = (nansum(Fz_NONdamaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_NONdamaged_total_steady_tip = (nansum(Mx_NONdamaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_steady_tip = (nansum(My_NONdamaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_steady_tip = (nansum(Mz_NONdamaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    % wingbeat average forces & torques: steady wb
    Fx_damaged_mean_steady_tip = nanmean(Fx_damaged_total_steady_tip);
    Fy_damaged_mean_steady_tip = nanmean(Fy_damaged_total_steady_tip);
    Fz_damaged_mean_steady_tip = nanmean(Fz_damaged_total_steady_tip);
    
    Mx_damaged_mean_steady_tip = nanmean(Mx_damaged_total_steady_tip);
    My_damaged_mean_steady_tip = nanmean(My_damaged_total_steady_tip);
    Mz_damaged_mean_steady_tip = nanmean(Mz_damaged_total_steady_tip);
    
    Fx_NONdamaged_mean_steady_tip = nanmean(Fx_NONdamaged_total_steady_tip);
    Fy_NONdamaged_mean_steady_tip = nanmean(Fy_NONdamaged_total_steady_tip);
    Fz_NONdamaged_mean_steady_tip = nanmean(Fz_NONdamaged_total_steady_tip);
    
    Mx_NONdamaged_mean_steady_tip = nanmean(Mx_NONdamaged_total_steady_tip);
    My_NONdamaged_mean_steady_tip = nanmean(My_NONdamaged_total_steady_tip);
    Mz_NONdamaged_mean_steady_tip = nanmean(Mz_NONdamaged_total_steady_tip);
    
    %% TIPCUT WB kin ALL MODs
    freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_all = (nansum(Fx_damaged_all,1))' / Mg_fly;
    Fy_damagedwing_all = (nansum(Fy_damaged_all,1))' / Mg_fly;
    Fz_damagedwing_all = (nansum(Fz_damaged_all,1))' / Mg_fly;
    
    Mx_damagedwing_all = (nansum(Mx_damaged_all,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_all = (nansum(My_damaged_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_all = (nansum(Mz_damaged_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_all = (nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_intactwing_all = (nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_intactwing_all = (nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_intactwing_all = (nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_all = (nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_all = (nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_all_tip = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all_tip = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all_tip = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all_tip = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_all_tip = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_all_tip = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_all_tip = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all_tip = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all_tip = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all_tip = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_all_tip = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_all_tip = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    % wingbeat average forces & torques: ALL MODs
    Fx_damaged_mean_all_tip = nanmean(Fx_damaged_total_all_tip);
    Fy_damaged_mean_all_tip = nanmean(Fy_damaged_total_all_tip);
    Fz_damaged_mean_all_tip = nanmean(Fz_damaged_total_all_tip);
    
    Mx_damaged_mean_all_tip = nanmean(Mx_damaged_total_all_tip);
    My_damaged_mean_all_tip = nanmean(My_damaged_total_all_tip);
    Mz_damaged_mean_all_tip = nanmean(Mz_damaged_total_all_tip);
    
    Fx_NONdamaged_mean_all_tip = nanmean(Fx_NONdamaged_total_all_tip);
    Fy_NONdamaged_mean_all_tip = nanmean(Fy_NONdamaged_total_all_tip);
    Fz_NONdamaged_mean_all_tip = nanmean(Fz_NONdamaged_total_all_tip);
    
    Mx_NONdamaged_mean_all_tip = nanmean(Mx_NONdamaged_total_all_tip);
    My_NONdamaged_mean_all_tip = nanmean(My_NONdamaged_total_all_tip);
    Mz_NONdamaged_mean_all_tip = nanmean(Mz_NONdamaged_total_all_tip);

%% TE CUT 4S2tip
% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_R';
c_max = max(wing_model.chords_R');

S2_sect_element = 1/12 .* c_sect .* h_sect.^3;
S3_sect_element = 1/36 .* c_sect .* h_sect.^4;

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% TE CUT morph
TEclip_ratio_now = cord_ratio_4S2tip;
c_cut = TEclip_ratio_now * c_max;

c_sect_cut = c_sect;
c_sect_cut(c_sect_cut>c_cut) = c_cut;

S2_sect = c_sect_cut .* h_sect .* y_sect.^2;
S3_sect = c_sect_cut .* h_sect .* y_sect.^3;

S2_damaged = sum(S2_sect);
S3_damaged = sum(S3_sect);

S2ratio = S2_damaged/S2_intact;
S3ratio = S3_damaged/S3_intact;

S2ratio_TE4S2tip = S2ratio;
S3ratio_TE4S2tip = S3ratio;
subplot(2,2,3)
plot(S2ratio_TE4S2tip,S3ratio_TE4S2tip,'dr','markersize',10)

% write damaged wing data to wingmodel
wing_model.chords_R = c_sect_cut';
    
    %% TE CUT 4S2tip wingshape & steady wingbeat
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_steady = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_steady = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_steady = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_steady = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_steady = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_steady = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % total force & torque damaged & intact wing
    Fx_damagedwing_steady = (nansum(Fx_damaged_steady,1))' / Mg_fly;
    Fy_damagedwing_steady = (nansum(Fy_damaged_steady,1))' / Mg_fly;
    Fz_damagedwing_steady = (nansum(Fz_damaged_steady,1))' / Mg_fly;
    
    Mx_damagedwing_steady = (nansum(Mx_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_steady = (nansum(My_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_steady = (nansum(Mz_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_steady = (nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_intactwing_steady = (nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_intactwing_steady = (nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_intactwing_steady = (nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_steady = (nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_steady = (nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_steady_TE4S2tip = (nansum(Fx_damaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_damaged_total_steady_TE4S2tip = (nansum(Fy_damaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_damaged_total_steady_TE4S2tip = (nansum(Fz_damaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_damaged_total_steady_TE4S2tip = (nansum(Mx_damaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_steady_TE4S2tip = (nansum(My_damaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_steady_TE4S2tip = (nansum(Mz_damaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_steady_TE4S2tip = (nansum(Fx_NONdamaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_NONdamaged_total_steady_TE4S2tip = (nansum(Fy_NONdamaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_NONdamaged_total_steady_TE4S2tip = (nansum(Fz_NONdamaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_NONdamaged_total_steady_TE4S2tip = (nansum(Mx_NONdamaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_steady_TE4S2tip = (nansum(My_NONdamaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_steady_TE4S2tip = (nansum(Mz_NONdamaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    % wingbeat average forces & torques: steady wb
    Fx_damaged_mean_steady_TE4S2tip = nanmean(Fx_damaged_total_steady_TE4S2tip);
    Fy_damaged_mean_steady_TE4S2tip = nanmean(Fy_damaged_total_steady_TE4S2tip);
    Fz_damaged_mean_steady_TE4S2tip = nanmean(Fz_damaged_total_steady_TE4S2tip);
    
    Mx_damaged_mean_steady_TE4S2tip = nanmean(Mx_damaged_total_steady_TE4S2tip);
    My_damaged_mean_steady_TE4S2tip = nanmean(My_damaged_total_steady_TE4S2tip);
    Mz_damaged_mean_steady_TE4S2tip = nanmean(Mz_damaged_total_steady_TE4S2tip);
    
    Fx_NONdamaged_mean_steady_TE4S2tip = nanmean(Fx_NONdamaged_total_steady_TE4S2tip);
    Fy_NONdamaged_mean_steady_TE4S2tip = nanmean(Fy_NONdamaged_total_steady_TE4S2tip);
    Fz_NONdamaged_mean_steady_TE4S2tip = nanmean(Fz_NONdamaged_total_steady_TE4S2tip);
    
    Mx_NONdamaged_mean_steady_TE4S2tip = nanmean(Mx_NONdamaged_total_steady_TE4S2tip);
    My_NONdamaged_mean_steady_TE4S2tip = nanmean(My_NONdamaged_total_steady_TE4S2tip);
    Mz_NONdamaged_mean_steady_TE4S2tip = nanmean(Mz_NONdamaged_total_steady_TE4S2tip);

    %% TE CUT 4S2tip wingshape & TIP CUT wingbeat
    freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % total force & torque damaged & intact wing
    Fx_damagedwing_all = (nansum(Fx_damaged_all,1))' / Mg_fly;
    Fy_damagedwing_all = (nansum(Fy_damaged_all,1))' / Mg_fly;
    Fz_damagedwing_all = (nansum(Fz_damaged_all,1))' / Mg_fly;
    
    Mx_damagedwing_all = (nansum(Mx_damaged_all,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_all = (nansum(My_damaged_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_all = (nansum(Mz_damaged_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_all = (nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_intactwing_all = (nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_intactwing_all = (nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_intactwing_all = (nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_all = (nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_all = (nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_all_TE4S2tip = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all_TE4S2tip = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all_TE4S2tip = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all_TE4S2tip = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_all_TE4S2tip = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_all_TE4S2tip = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_all_TE4S2tip = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all_TE4S2tip = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all_TE4S2tip = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all_TE4S2tip = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_all_TE4S2tip = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_all_TE4S2tip = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    % wingbeat average forces & torques: ALL MODs
    Fx_damaged_mean_all_TE4S2tip = nanmean(Fx_damaged_total_all_TE4S2tip);
    Fy_damaged_mean_all_TE4S2tip = nanmean(Fy_damaged_total_all_TE4S2tip);
    Fz_damaged_mean_all_TE4S2tip = nanmean(Fz_damaged_total_all_TE4S2tip);
    
    Mx_damaged_mean_all_TE4S2tip = nanmean(Mx_damaged_total_all_TE4S2tip);
    My_damaged_mean_all_TE4S2tip = nanmean(My_damaged_total_all_TE4S2tip);
    Mz_damaged_mean_all_TE4S2tip = nanmean(Mz_damaged_total_all_TE4S2tip);
    
    Fx_NONdamaged_mean_all_TE4S2tip = nanmean(Fx_NONdamaged_total_all_TE4S2tip);
    Fy_NONdamaged_mean_all_TE4S2tip = nanmean(Fy_NONdamaged_total_all_TE4S2tip);
    Fz_NONdamaged_mean_all_TE4S2tip = nanmean(Fz_NONdamaged_total_all_TE4S2tip);
    
    Mx_NONdamaged_mean_all_TE4S2tip = nanmean(Mx_NONdamaged_total_all_TE4S2tip);
    My_NONdamaged_mean_all_TE4S2tip = nanmean(My_NONdamaged_total_all_TE4S2tip);
    Mz_NONdamaged_mean_all_TE4S2tip = nanmean(Mz_NONdamaged_total_all_TE4S2tip);


%% TE CUT 4S3tip
% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_L';
c_max = max(wing_model.chords_L');

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% TE CUT morph
TEclip_ratio_now = cord_ratio_4S3tip;
c_cut = TEclip_ratio_now * c_max;

c_sect_cut = c_sect;
c_sect_cut(c_sect_cut>c_cut) = c_cut;

S2_sect = c_sect_cut .* h_sect .* y_sect.^2;
S3_sect = c_sect_cut .* h_sect .* y_sect.^3;

S2_damaged = sum(S2_sect);
S3_damaged = sum(S3_sect);

S2ratio = S2_damaged/S2_intact;
S3ratio = S3_damaged/S3_intact;

S2ratio_TE4S3tip = S2ratio;
S3ratio_TE4S3tip = S3ratio;
subplot(2,2,3)
plot(S2ratio_TE4S3tip,S3ratio_TE4S3tip,'dr','markersize',10)

% write damaged wing data to wingmodel
wing_model.chords_R = c_sect_cut';
    
    %% TE CUT 4S3tip wingshape & steady wingbeat
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_steady = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_steady = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_steady = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_steady = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_steady = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_steady = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % total force & torque damaged & intact wing
    Fx_damagedwing_steady = (nansum(Fx_damaged_steady,1))' / Mg_fly;
    Fy_damagedwing_steady = (nansum(Fy_damaged_steady,1))' / Mg_fly;
    Fz_damagedwing_steady = (nansum(Fz_damaged_steady,1))' / Mg_fly;
    
    Mx_damagedwing_steady = (nansum(Mx_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_steady = (nansum(My_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_steady = (nansum(Mz_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_steady = (nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_intactwing_steady = (nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_intactwing_steady = (nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_intactwing_steady = (nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_steady = (nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_steady = (nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_steady_TE4S3tip = (nansum(Fx_damaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_damaged_total_steady_TE4S3tip = (nansum(Fy_damaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_damaged_total_steady_TE4S3tip = (nansum(Fz_damaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_damaged_total_steady_TE4S3tip = (nansum(Mx_damaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_steady_TE4S3tip = (nansum(My_damaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_steady_TE4S3tip = (nansum(Mz_damaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_steady_TE4S3tip = (nansum(Fx_NONdamaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_NONdamaged_total_steady_TE4S3tip = (nansum(Fy_NONdamaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_NONdamaged_total_steady_TE4S3tip = (nansum(Fz_NONdamaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_NONdamaged_total_steady_TE4S3tip = (nansum(Mx_NONdamaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_steady_TE4S3tip = (nansum(My_NONdamaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_steady_TE4S3tip = (nansum(Mz_NONdamaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    % wingbeat average forces & torques: steady wb
    Fx_damaged_mean_steady_TE4S3tip = nanmean(Fx_damaged_total_steady_TE4S3tip);
    Fy_damaged_mean_steady_TE4S3tip = nanmean(Fy_damaged_total_steady_TE4S3tip);
    Fz_damaged_mean_steady_TE4S3tip = nanmean(Fz_damaged_total_steady_TE4S3tip);
    
    Mx_damaged_mean_steady_TE4S3tip = nanmean(Mx_damaged_total_steady_TE4S3tip);
    My_damaged_mean_steady_TE4S3tip = nanmean(My_damaged_total_steady_TE4S3tip);
    Mz_damaged_mean_steady_TE4S3tip = nanmean(Mz_damaged_total_steady_TE4S3tip);
    
    Fx_NONdamaged_mean_steady_TE4S3tip = nanmean(Fx_NONdamaged_total_steady_TE4S3tip);
    Fy_NONdamaged_mean_steady_TE4S3tip = nanmean(Fy_NONdamaged_total_steady_TE4S3tip);
    Fz_NONdamaged_mean_steady_TE4S3tip = nanmean(Fz_NONdamaged_total_steady_TE4S3tip);
    
    Mx_NONdamaged_mean_steady_TE4S3tip = nanmean(Mx_NONdamaged_total_steady_TE4S3tip);
    My_NONdamaged_mean_steady_TE4S3tip = nanmean(My_NONdamaged_total_steady_TE4S3tip);
    Mz_NONdamaged_mean_steady_TE4S3tip = nanmean(Mz_NONdamaged_total_steady_TE4S3tip);

    %% TE CUT 4S3tip wingshape & TIP CUT wingbeat
    freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % total force & torque damaged & intact wing
    Fx_damagedwing_all = (nansum(Fx_damaged_all,1))' / Mg_fly;
    Fy_damagedwing_all = (nansum(Fy_damaged_all,1))' / Mg_fly;
    Fz_damagedwing_all = (nansum(Fz_damaged_all,1))' / Mg_fly;
    
    Mx_damagedwing_all = (nansum(Mx_damaged_all,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_all = (nansum(My_damaged_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_all = (nansum(Mz_damaged_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_all = (nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_intactwing_all = (nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_intactwing_all = (nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_intactwing_all = (nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_all = (nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_all = (nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_all_TE4S3tip = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all_TE4S3tip = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all_TE4S3tip = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all_TE4S3tip = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_all_TE4S3tip = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_all_TE4S3tip = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_all_TE4S3tip = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all_TE4S3tip = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all_TE4S3tip = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all_TE4S3tip = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_all_TE4S3tip = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_all_TE4S3tip = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

    %% wingbeat average forces & torques: ALL MODs
    Fx_damaged_mean_all_TE4S3tip = nanmean(Fx_damaged_total_all_TE4S3tip);
    Fy_damaged_mean_all_TE4S3tip = nanmean(Fy_damaged_total_all_TE4S3tip);
    Fz_damaged_mean_all_TE4S3tip = nanmean(Fz_damaged_total_all_TE4S3tip);
    
    Mx_damaged_mean_all_TE4S3tip = nanmean(Mx_damaged_total_all_TE4S3tip);
    My_damaged_mean_all_TE4S3tip = nanmean(My_damaged_total_all_TE4S3tip);
    Mz_damaged_mean_all_TE4S3tip = nanmean(Mz_damaged_total_all_TE4S3tip);
    
    Fx_NONdamaged_mean_all_TE4S3tip = nanmean(Fx_NONdamaged_total_all_TE4S3tip);
    Fy_NONdamaged_mean_all_TE4S3tip = nanmean(Fy_NONdamaged_total_all_TE4S3tip);
    Fz_NONdamaged_mean_all_TE4S3tip = nanmean(Fz_NONdamaged_total_all_TE4S3tip);
    
    Mx_NONdamaged_mean_all_TE4S3tip = nanmean(Mx_NONdamaged_total_all_TE4S3tip);
    My_NONdamaged_mean_all_TE4S3tip = nanmean(My_NONdamaged_total_all_TE4S3tip);
    Mz_NONdamaged_mean_all_TE4S3tip = nanmean(Mz_NONdamaged_total_all_TE4S3tip);

    %% save data
    save(['allMODs_TE4TipClip_freqAsym',num2str(freq_asymFitNr),'_roton',num2str(rot_on),'.mat'])
    
    %% plot
    if plot_on == 1
    mkdir(['qsModel_FnM_TE4TipClip_roton',num2str(rot_on)])
    cd(['qsModel_FnM_TE4TipClip_roton',num2str(rot_on)])

    % with all syms
    figure(2)
    
    % S2-S3
    subplot(2,2,1)
    hold on
    plot(S2ratios_tip,polyval(S2S3_fit_tip,S2ratios_tip),'-k')
    plot(S2ratios_TE,polyval(S2S3_fit_TE,S2ratios_TE),'-k')

    plot(S2ratios_tip,S3ratios_tip,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios_TE,S3ratios_TE,'dk','markersize',10,'markerfacecolor','w')

    plot([S2ratio_tip S2ratio_TE4S2tip],[S3ratio_tip S3ratio_TE4S2tip],'-k')
    plot([S2ratio_tip S2ratio_TE4S3tip],[S3ratio_tip S3ratio_TE4S3tip],'-k')
    
    plot(S2ratio_tip,S3ratio_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratio_TE4S2tip,S3ratio_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratio_TE4S3tip,S3ratio_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    xlabel('S2ratio')
    ylabel('S3ratio')
    axis equal
    axis([.5 1 .5 1])
    set(gca,'xtick',.5:.5:1)
    set(gca,'ytick',.5:.5:1)
    
    % Fz
    subplot(2,2,2)
    hold on
    plot([0 5],[1 1],'--k')
    plot(S2ratio_tip,-Fz_damaged_mean_steady_tip,'ok','markersize',10,'markerfacecolor','y')
    plot(S2ratio_tip,-Fz_damaged_mean_all_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratio_TE4S2tip,-Fz_damaged_mean_all_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratio_TE4S3tip,-Fz_damaged_mean_all_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    xlabel('S2ratio')
    ylabel('Fz/mg')
    axis([.5 1 .5 1.5])
    set(gca,'xtick',.5:.5:1)
    set(gca,'ytick',.5:.5:1.5)    
    
    % Tx
    subplot(2,2,4)
    hold on
    plot([0 5],[0 0],'--k')
    plot(S3ratio_tip,Mx_damaged_mean_steady_tip,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratio_tip,Mx_damaged_mean_all_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(S3ratio_TE4S2tip,Mx_damaged_mean_all_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(S3ratio_TE4S3tip,Mx_damaged_mean_all_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    xlabel('S3ratio')
    ylabel('Tx/mgl')
    axis([.5 1 -.1 .1])
    set(gca,'xtick',.5:.5:1)
    set(gca,'ytick',-.1:.1:.1)   
    
    % legend
    subplot(2,2,3)
    hold on
    plot(1,1,'ok','markersize',10,'markerfacecolor','y')
    plot(1,1,'ok','markersize',10,'markerfacecolor','g')
    plot(1,1,'dk','markersize',10,'markerfacecolor','b')
    plot(1,1,'dk','markersize',10,'markerfacecolor','r')
    
    axis off
    legend('tip cut normal','tip cut WDR','TE cut same S2','TE cut same S3','location','E')
    axis([0 .5 0 .5])
    
    saveas(gca,['FzTx_vs_S2S3_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzTx_vs_S2S3_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzTx_vs_S2S3_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

    %% list
    figure(3)
    % S2-S3
    subplot(2,2,1)
    hold on
    plot(S2ratios_tip,polyval(S2S3_fit_tip,S2ratios_tip),'-k')
    plot(S2ratios_TE,polyval(S2S3_fit_TE,S2ratios_TE),'-k')

    plot(S2ratios_tip,S3ratios_tip,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios_TE,S3ratios_TE,'dk','markersize',10,'markerfacecolor','w')

    plot([S2ratio_tip S2ratio_TE4S2tip],[S3ratio_tip S3ratio_TE4S2tip],'-k')
    plot([S2ratio_tip S2ratio_TE4S3tip],[S3ratio_tip S3ratio_TE4S3tip],'-k')
    
    plot(S2ratio_tip,S3ratio_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratio_TE4S2tip,S3ratio_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratio_TE4S3tip,S3ratio_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    xlabel('S2ratio')
    ylabel('S3ratio')
    axis equal
    axis([.5 1 .5 1])
    set(gca,'xtick',.5:.5:1)
    set(gca,'ytick',.5:.5:1)
    
    % Fz
    subplot(2,2,2)
    hold on
    plot([0 5],[1 1],'--k')
    plot(1,-Fz_damaged_mean_steady_tip,'ok','markersize',10,'markerfacecolor','y')
    plot(2,-Fz_damaged_mean_all_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(3,-Fz_damaged_mean_all_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(4,-Fz_damaged_mean_all_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    ylabel('Fz/mg')
    axis([0 5 .5 1.5])
    set(gca,'xtick',[])
    set(gca,'ytick',.5:.5:1.5)    
    

    % Tx
    subplot(2,2,4)
    hold on
    plot([0 5],[0 0],'--k')
    plot(1,Mx_damaged_mean_steady_tip,'ok','markersize',10,'markerfacecolor','y')
    plot(2,Mx_damaged_mean_all_tip,'ok','markersize',10,'markerfacecolor','g')
    plot(3,Mx_damaged_mean_all_TE4S2tip,'dk','markersize',10,'markerfacecolor','b')
    plot(4,Mx_damaged_mean_all_TE4S3tip,'dk','markersize',10,'markerfacecolor','r')
    
    ylabel('Tx/mgl')
    axis([0 5 -.1 .1])
    set(gca,'xtick',[])
    set(gca,'ytick',-.1:.1:.1)   

    % legend
    subplot(2,2,3)
    hold on
    plot(1,1,'ok','markersize',10,'markerfacecolor','y')
    plot(1,1,'ok','markersize',10,'markerfacecolor','g')
    plot(1,1,'dk','markersize',10,'markerfacecolor','b')
    plot(1,1,'dk','markersize',10,'markerfacecolor','r')
    
    axis off
    legend('tip cut normal','tip cut WDR','TE cut same S2','TE cut same S3','location','E')
    axis([0 .5 0 .5])

    saveas(gca,['Fz_Tx_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Tx_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Tx_TE4TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])
    cd ..
end
