
clear
clc
close all

%% const
% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% robo data
Lrobo = .23;
crobo = .065;
ARrobo = Lrobo/crobo;

% fluid data
v_air = 1.568e-5; %kinematic viscocity in m^2/s
v_oil = 115.0*1e-6; %kinematic viscocity in m^2/s
% rho_air = 1.22521; %kg/m^3
rho_oil = 880.0; %kg/m^3

%% Elzinga scaling
% l_scale_up = crobo/c_fly;
c_fly = Lwing/ARwing_fly; %mean chord of actual fly in m
f_robo2fly = (v_air*crobo*Lrobo)/(c_fly*Lwing*v_oil);
F_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^4)*((f_robo2fly)^2);
M_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^5)*((f_robo2fly)^2);

%% cut wing geometry (FROM BW IMAGES OF SOLIDWORKS PROJECTIONS)
geom_file = dir('BorfMorphCutDatabase*');
load(geom_file.name)
BorfMorphCutData.cut_ratio = BorfMorphCutData.cut_perc/100;

%% run data
Nwb = 10;
wb_start = 2;
wb_stop = 10;

%% cali CF for all (others have ERROR: Mx is assymetric)
% cali_template = 'cali_matrix_interp_*.mat';
cali_file = 'cali_matrix_interp_L_Full_R_Full_CF.mat';
load(cali_file)

%% load & save intact
cut_type_now = 0;  % no cut
cut_perc_now = 100;

dir_now = dir('L_Full_R_Full_Acrylic*');

if dir_now.isdir == 1
    cd(dir_now.name)
    m=0;
    
    %% full amp reduce
    file_template = '*reduce*.mat';
    file_names = dir(file_template);
    Amp_type_now = 2;

    for i = 1:length(file_names)
    
        file_now = file_names(i).name
        load(file_now)
        
        m=m+1;
        
        Amp_type(m,1) = Amp_type_now;

        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);
        
        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now);
        
        cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
        cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
        cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);
        
        l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
        CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
        A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
        S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
        S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
        S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);
        
        CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
        A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
        S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
        S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
        S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);
        
        l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
        CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
        S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
        S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
        S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);

        extract_NONclipped_FnM_AmpReduce_MIRROR_INCcali

    end
end

cd ..

%% remove offset
n_Amp_increase = find(Amp_type==1);
n_Amp_reduce = find(Amp_type==2);

n_Amp_increase_steady = find(Amp_type==1&Amp_ratio==1);
n_Amp_reduce_steady = find(Amp_type==2&Amp_ratio==1);

% % fix Fz jump
% Fz_norm_fix(n_Amp_increase,1) = Fz_norm(n_Amp_increase) - Fz_norm(n_Amp_increase_steady) - 1;
% Fz_norm_fix(n_Amp_reduce,1) = Fz_norm(n_Amp_reduce) - Fz_norm(n_Amp_reduce_steady) - 1;

% Fy equilibrium
FyMinSteady_norm(n_Amp_increase,1) = Fy_norm(n_Amp_increase) - Fy_norm(n_Amp_increase_steady);
FyMinSteady_norm(n_Amp_reduce,1) = Fy_norm(n_Amp_reduce) - Fy_norm(n_Amp_reduce_steady);

% roll equilibrium
MxMinSteady_norm(n_Amp_increase,1) = Mx_norm(n_Amp_increase) - Mx_norm(n_Amp_increase_steady);
MxMinSteady_norm(n_Amp_reduce,1) = Mx_norm(n_Amp_reduce) - Mx_norm(n_Amp_reduce_steady);

% pitch equilibrium
MyMinSteady_norm(n_Amp_increase,1) = My_norm(n_Amp_increase) - My_norm(n_Amp_increase_steady);
MyMinSteady_norm(n_Amp_reduce,1) = My_norm(n_Amp_reduce) - My_norm(n_Amp_reduce_steady);

% yaw equilibrium
MzMinSteady_norm(n_Amp_increase,1) = Mz_norm(n_Amp_increase) - Mz_norm(n_Amp_increase_steady);
MzMinSteady_norm(n_Amp_reduce,1) = Mz_norm(n_Amp_reduce) - Mz_norm(n_Amp_reduce_steady);

%% calc linear fits
% F-Aratio
[Fx_Amp_fit, Fx_Amp_fit_error] = polyfit(Amp_ratio,Fx_norm,1);
[Fy_Amp_fit, Fy_Amp_fit_error] = polyfit(Amp_ratio,Fy_norm,1);
[Fz_Amp_fit, Fz_Amp_fit_error] = polyfit(Amp_ratio,Fz_norm,1);

[FyMinSteady_Amp_fit, FyMinSteady_Amp_fit_error] = polyfit(Amp_ratio,FyMinSteady_norm,1);
% [Fz_fix_Amp_fit, Fz_fix_Amp_fit_error] = polyfit(Amp_ratio,Fz_norm_fix,1);

% M-S
[Mx_Amp_fit, Mx_Amp_fit_error] = polyfit(Amp_ratio,Mx_norm,1);
[My_Amp_fit, My_Amp_fit_error] = polyfit(Amp_ratio,My_norm,1);
[Mz_Amp_fit, Mz_Amp_fit_error] = polyfit(Amp_ratio,Mz_norm,1);

[MxMinSteady_Amp_fit, MxMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MxMinSteady_norm,1);
[MyMinSteady_Amp_fit, MyMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MyMinSteady_norm,1);
[MzMinSteady_Amp_fit, MzMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MzMinSteady_norm,1);

[My_CoM_Amp_fit, My_CoM_Amp_fit_error] = polyfit(Amp_ratio,My_CoM_norm,1);

%% calc linear fits with Amplitude Squared
% F-Aratio
[Fx_Amp2_fit, Fx_Amp2_fit_error] = polyfit(Amp_ratio.^2,Fx_norm,1);
[Fy_Amp2_fit, Fy_Amp2_fit_error] = polyfit(Amp_ratio.^2,Fy_norm,1);
[Fz_Amp2_fit, Fz_Amp2_fit_error] = polyfit(Amp_ratio.^2,Fz_norm,1);

[FyMinSteady_Amp2_fit, FyMinSteady_Amp2_fit_error] = polyfit(Amp_ratio.^2,FyMinSteady_norm,1);
% [Fz_fix_Amp2_fit, Fz_fix_Amp2_fit_error] = polyfit(Amp_ratio.^2,Fz_norm_fix,1);

% M-S
[Mx_Amp2_fit, Mx_Amp2_fit_error] = polyfit(Amp_ratio.^2,Mx_norm,1);
[My_Amp2_fit, My_Amp2_fit_error] = polyfit(Amp_ratio.^2,My_norm,1);
[Mz_Amp2_fit, Mz_Amp2_fit_error] = polyfit(Amp_ratio.^2,Mz_norm,1);

[MxMinSteady_Amp2_fit, MxMinSteady_Amp2_fit_error] = polyfit(Amp_ratio.^2,MxMinSteady_norm,1);
[MyMinSteady_Amp2_fit, MyMinSteady_Amp2_fit_error] = polyfit(Amp_ratio.^2,MyMinSteady_norm,1);
[MzMinSteady_Amp2_fit, MzMinSteady_Amp2_fit_error] = polyfit(Amp_ratio.^2,MzMinSteady_norm,1);

[My_CoM_Amp2_fit, My_CoM_Amp2_fit_error] = polyfit(Amp_ratio.^2,My_CoM_norm,1);

%% calc parabola fits
% F-Aratio
[Fx_Amp_fit2, Fx_Amp_fit2_error] = polyfit(Amp_ratio,Fx_norm,2);
[Fy_Amp_fit2, Fy_Amp_fit2_error] = polyfit(Amp_ratio,Fy_norm,2);
[Fz_Amp_fit2, Fz_Amp_fit2_error] = polyfit(Amp_ratio,Fz_norm,2);

[FyMinSteady_Amp_fit2, FyMinSteady_Amp_fit2_error] = polyfit(Amp_ratio,FyMinSteady_norm,2);
% [Fz_fix_Amp_fit2, Fz_fix_Amp_fit2_error] = polyfit(Amp_ratio,Fz_norm_fix,2);

% M-S
[Mx_Amp_fit2, Mx_Amp_fit2_error] = polyfit(Amp_ratio,Mx_norm,2);
[My_Amp_fit2, My_Amp_fit2_error] = polyfit(Amp_ratio,My_norm,2);
[Mz_Amp_fit2, Mz_Amp_fit2_error] = polyfit(Amp_ratio,Mz_norm,2);

[MxMinSteady_Amp_fit2, MxMinSteady_Amp_fit2_error] = polyfit(Amp_ratio,MxMinSteady_norm,2);
[MyMinSteady_Amp_fit2, MyMinSteady_Amp_fit2_error] = polyfit(Amp_ratio,MyMinSteady_norm,2);
[MzMinSteady_Amp_fit2, MzMinSteady_Amp_fit2_error] = polyfit(Amp_ratio,MzMinSteady_norm,2);

[My_CoM_Amp_fit2, My_CoM_Amp_fit2_error] = polyfit(Amp_ratio,My_CoM_norm,2);

%% plot F-Amp & M-Amp (My@"CoM")
% F-Amp
figure
subplot(1,2,1)
hold on
plot(Amp_ratio,Fx_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio,Fy_norm,'ok','markerfacecolor','r')
plot(Amp_ratio,FyMinSteady_norm,'ok','markerfacecolor','r')
plot(Amp_ratio,Fz_norm,'ok','markerfacecolor','g')

% linear fits
plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fx_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fy_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(FyMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fz_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

% legend('Fx','Fy','Fz','location','SW')
legend('Fx','Fy','Fz','location','E')
xlabel('Amp ratio')
ylabel('F/mg')
axis([0.75 1 -1 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)

% M-Amp
subplot(1,2,2)
hold on

plot(Amp_ratio,MxMinSteady_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio,MyMinSteady_norm,'ok','markerfacecolor','g')
plot(Amp_ratio,My_CoM_norm,'dk','markerfacecolor','r')
plot(Amp_ratio,MzMinSteady_norm,'ok','markerfacecolor','g')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(MxMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(MyMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(MzMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(My_CoM_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

legend('Mx','My@CoM','Mz','location','NE')
xlabel('Amp ratio')
ylabel('T/mgl')
axis([0.75 1 -.05 .2])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1:.05:1)

mkdir('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')
cd('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.png'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.svg'])

cd ..

%% plot F-Amp^2 & M-Amp^2 (My@"CoM")
% F-Amp^2
figure
subplot(1,2,1)
hold on
plot(Amp_ratio.^2,Fx_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio.^2,Fy_norm,'ok','markerfacecolor','r')
plot(Amp_ratio.^2,FyMinSteady_norm,'ok','markerfacecolor','r')
plot(Amp_ratio.^2,Fz_norm,'ok','markerfacecolor','g')

% linear fits
plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(Fx_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')
% plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(Fy_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')
plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(FyMinSteady_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')
plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(Fz_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')

% legend('Fx','Fy','Fz','location','SW')
legend('Fx','Fy','Fz','location','E')
xlabel('Amp Squared Ratio')
ylabel('F/mg')
axis([0.5 1 -1 .25])
set(gca,'xtick',0.5:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)

% M-Amp
subplot(1,2,2)
hold on

plot(Amp_ratio.^2,MxMinSteady_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio.^2,MyMinSteady_norm,'ok','markerfacecolor','g')
plot(Amp_ratio.^2,My_CoM_norm,'dk','markerfacecolor','r')
plot(Amp_ratio.^2,MzMinSteady_norm,'ok','markerfacecolor','g')

plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(MxMinSteady_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')
% plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(MyMinSteady_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')
plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(MzMinSteady_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')

plot([min(Amp_ratio.^2) max(Amp_ratio.^2)],polyval(My_CoM_Amp2_fit,[min(Amp_ratio.^2) max(Amp_ratio.^2)]),'k')

legend('Mx','My@CoM','Mz','location','NE')
xlabel('Amp Squared Ratio')
ylabel('T/mgl')
axis([0.5 1 -.05 .2])
set(gca,'xtick',0.5:.25:1.25)
set(gca,'ytick',-1:.05:1)

mkdir('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')
cd('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_AmpSquareFit.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_AmpSquareFit.png'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_AmpSquareFit.svg'])

cd ..

%% plot F-Amp & M-Amp (My@"CoM") PARABOLA FIT
% F-Amp
figure
subplot(1,2,1)
hold on
plot(Amp_ratio,Fx_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio,Fy_norm,'ok','markerfacecolor','r')
plot(Amp_ratio,FyMinSteady_norm,'ok','markerfacecolor','r')
plot(Amp_ratio,Fz_norm,'ok','markerfacecolor','g')

% linear fits
plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(Fx_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')
% plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(Fy_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')
plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(FyMinSteady_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')
plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(Fz_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')

% legend('Fx','Fy','Fz','location','SW')
legend('Fx','Fy','Fz','location','E')
xlabel('Amp ratio')
ylabel('F/mg')
axis([0.75 1 -1 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)

% M-Amp
subplot(1,2,2)
hold on

plot(Amp_ratio,MxMinSteady_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio,MyMinSteady_norm,'ok','markerfacecolor','g')
plot(Amp_ratio,My_CoM_norm,'dk','markerfacecolor','r')
plot(Amp_ratio,MzMinSteady_norm,'ok','markerfacecolor','g')

plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(MxMinSteady_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')
% plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(MyMinSteady_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')
plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(MzMinSteady_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')

plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(My_CoM_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')

legend('Mx','My@CoM','Mz','location','NE')
xlabel('Amp ratio')
ylabel('T/mgl')
axis([0.75 1 -.05 .2])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1:.05:1)

mkdir('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')
cd('figures_NONcutWing_FnM_vs_ReducedStrokeAmplitude')

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_ParabFit.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_ParabFit.png'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_ParabFit.svg'])

cd ..
%% save data
save('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF',...
    'Amp_ratio',...
    'Fx_norm',...
    'Fy_norm',...
    'Fz_norm',...
    'Mx_norm',...
    'My_norm',...
    'Mz_norm',...
    'FyMinSteady_norm',...
    'My_CoM_norm',...
    'MxMinSteady_norm',...
    'MyMinSteady_norm',...
    'MzMinSteady_norm',...
    ...
    'Fx_Amp_fit',...
    'Fx_Amp_fit_error',...
    'Fy_Amp_fit',...
    'Fy_Amp_fit_error',...
    'Fz_Amp_fit',...
    'Fz_Amp_fit_error',...
    'FyMinSteady_Amp_fit',...
    'FyMinSteady_Amp_fit_error',...
    'Mx_Amp_fit',...
    'Mx_Amp_fit_error',...
    'My_Amp_fit',...
    'My_Amp_fit_error',...
    'Mz_Amp_fit',...
    'Mz_Amp_fit_error',...
    'MxMinSteady_Amp_fit',...
    'MxMinSteady_Amp_fit_error',...
    'MyMinSteady_Amp_fit',...
    'MyMinSteady_Amp_fit_error',...
    'MzMinSteady_Amp_fit',...
    'MzMinSteady_Amp_fit_error',...
    ...
    'Fx_Amp2_fit',...
    'Fx_Amp2_fit_error',...
    'Fy_Amp2_fit',...
    'Fy_Amp2_fit_error',...
    'Fz_Amp2_fit',...
    'Fz_Amp2_fit_error',...
    'FyMinSteady_Amp2_fit',...
    'FyMinSteady_Amp2_fit_error',...
    'Mx_Amp2_fit',...
    'Mx_Amp2_fit_error',...
    'My_Amp2_fit',...
    'My_Amp2_fit_error',...
    'Mz_Amp2_fit',...
    'Mz_Amp2_fit_error',...
    'MxMinSteady_Amp2_fit',...
    'MxMinSteady_Amp2_fit_error',...
    'MyMinSteady_Amp2_fit',...
    'MyMinSteady_Amp2_fit_error',...
    'MzMinSteady_Amp2_fit',...
    'MzMinSteady_Amp2_fit_error',...
    ...
    'Fx_Amp_fit2',...
    'Fx_Amp_fit2_error',...
    'Fy_Amp_fit2',...
    'Fy_Amp_fit2_error',...
    'Fz_Amp_fit2',...
    'Fz_Amp_fit2_error',...
    'FyMinSteady_Amp_fit2',...
    'FyMinSteady_Amp_fit2_error',...
    'Mx_Amp_fit2',...
    'Mx_Amp_fit2_error',...
    'My_Amp_fit2',...
    'My_Amp_fit2_error',...
    'Mz_Amp_fit2',...
    'Mz_Amp_fit2_error',...
    'MxMinSteady_Amp_fit2',...
    'MxMinSteady_Amp_fit2_error',...
    'MyMinSteady_Amp_fit2',...
    'MyMinSteady_Amp_fit2_error',...
    'MzMinSteady_Amp_fit2',...
    'MzMinSteady_Amp_fit2_error');
