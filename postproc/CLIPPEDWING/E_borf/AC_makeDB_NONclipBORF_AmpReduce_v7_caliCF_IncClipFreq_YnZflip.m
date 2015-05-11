
clear
clc
close all
warning off

%% const
% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% Fenhance_by_freqMod_4SteadyWB
load('Fenhance_by_freqMod_4SteadyWB.mat')

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
    
        file_now = file_names(i).name;
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
%         extract_NONclipped_FnM_AmpReduce_MIRROR_NOcali

    end
end

cd ..

%% remove offset
FyMinSteady_norm = Fy_norm - Fy_norm(Amp_ratio==1);
MxMinSteady_norm = Mx_norm - Mx_norm(Amp_ratio==1);
MyMinSteady_norm = My_norm - My_norm(Amp_ratio==1);
MzMinSteady_norm = Mz_norm - Mz_norm(Amp_ratio==1);

%% calc Forces & Torques at wingbeat frequency of cut wing flies
Fx_norm_freqMod = Fx_norm*Fnorm_clip_steady;
Fy_norm_freqMod = Fy_norm*Fnorm_clip_steady;
Fz_norm_freqMod = Fz_norm*Fnorm_clip_steady;

Mx_norm_freqMod = Mx_norm*Fnorm_clip_steady;
My_norm_freqMod = My_norm*Fnorm_clip_steady;
Mz_norm_freqMod = Mz_norm*Fnorm_clip_steady;

FyMinSteady_norm_freqMod = FyMinSteady_norm*Fnorm_clip_steady;
MxMinSteady_norm_freqMod = MxMinSteady_norm*Fnorm_clip_steady;
MyMinSteady_norm_freqMod = MyMinSteady_norm*Fnorm_clip_steady;
MzMinSteady_norm_freqMod = MzMinSteady_norm*Fnorm_clip_steady;

My_CoM_norm_freqMod = My_CoM_norm*Fnorm_clip_steady;

%% calc linear fits for steady wb frequency
% F-Aratio
[Fx_Amp_fit, Fx_Amp_fit_error] = polyfit(Amp_ratio,Fx_norm,1);
[Fy_Amp_fit, Fy_Amp_fit_error] = polyfit(Amp_ratio,Fy_norm,1);
[Fz_Amp_fit, Fz_Amp_fit_error] = polyfit(Amp_ratio,Fz_norm,1);

[FyMinSteady_Amp_fit, FyMinSteady_Amp_fit_error] = polyfit(Amp_ratio,FyMinSteady_norm,1);

% M-Aratio
[Mx_Amp_fit, Mx_Amp_fit_error] = polyfit(Amp_ratio,Mx_norm,1);
[My_Amp_fit, My_Amp_fit_error] = polyfit(Amp_ratio,My_norm,1);
[Mz_Amp_fit, Mz_Amp_fit_error] = polyfit(Amp_ratio,Mz_norm,1);

[MxMinSteady_Amp_fit, MxMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MxMinSteady_norm,1);
[MyMinSteady_Amp_fit, MyMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MyMinSteady_norm,1);
[MzMinSteady_Amp_fit, MzMinSteady_Amp_fit_error] = polyfit(Amp_ratio,MzMinSteady_norm,1);

[My_CoM_Amp_fit, My_CoM_Amp_fit_error] = polyfit(Amp_ratio,My_CoM_norm,1);

%% calc linear fits for clipped fly wb frequency
% F-Aratio
[Fx_Amp_fit_freqMod, Fx_Amp_fit_error_freqMod] = polyfit(Amp_ratio,Fx_norm_freqMod,1);
[Fy_Amp_fit_freqMod, Fy_Amp_fit_error_freqMod] = polyfit(Amp_ratio,Fy_norm_freqMod,1);
[Fz_Amp_fit_freqMod, Fz_Amp_fit_error_freqMod] = polyfit(Amp_ratio,Fz_norm_freqMod,1);

[FyMinSteady_Amp_fit_freqMod, FyMinSteady_Amp_fit_error_freqMod] = polyfit(Amp_ratio,FyMinSteady_norm_freqMod,1);

% M-Aratio
[Mx_Amp_fit_freqMod, Mx_Amp_fit_error_freqMod] = polyfit(Amp_ratio,Mx_norm_freqMod,1);
[My_Amp_fit_freqMod, My_Amp_fit_error_freqMod] = polyfit(Amp_ratio,My_norm_freqMod,1);
[Mz_Amp_fit_freqMod, Mz_Amp_fit_error_freqMod] = polyfit(Amp_ratio,Mz_norm_freqMod,1);

[MxMinSteady_Amp_fit_freqMod, MxMinSteady_Amp_fit_error_freqMod] = polyfit(Amp_ratio,MxMinSteady_norm_freqMod,1);
[MyMinSteady_Amp_fit_freqMod, MyMinSteady_Amp_fit_error_freqMod] = polyfit(Amp_ratio,MyMinSteady_norm_freqMod,1);
[MzMinSteady_Amp_fit_freqMod, MzMinSteady_Amp_fit_error_freqMod] = polyfit(Amp_ratio,MzMinSteady_norm_freqMod,1);

[My_CoM_Amp_fit_freqMod, My_CoM_Amp_fit_error_freqMod] = polyfit(Amp_ratio,My_CoM_norm_freqMod,1);

%% plot F-Amp & M-Amp (My@"CoM")

% F-Amp
figure
subplot(1,2,1)
hold on
plot(Amp_ratio,Fx_norm,'sk','markersize',7,'markerfacecolor','b')
% plot(Amp_ratio,-Fy_norm,'sk','markersize',7,'markerfacecolor','r')
plot(Amp_ratio,-FyMinSteady_norm,'sk','markersize',7,'markerfacecolor','r')
plot(Amp_ratio,-Fz_norm,'sk','markersize',7,'markerfacecolor','g')

plot(Amp_ratio,Fx_norm_freqMod,'sk','markersize',7,'markerfacecolor','c')
% plot(Amp_ratio,-Fy_norm_freqMod,'sk','markersize',7,'markerfacecolor',[1 .5 0])
plot(Amp_ratio,-FyMinSteady_norm_freqMod,'sk','markersize',7,'markerfacecolor',[1 .5 0])
plot(Amp_ratio,-Fz_norm_freqMod,'sk','markersize',7,'markerfacecolor',[0 .5 0])

% linear fits
plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fx_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(-Fy_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-FyMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-Fz_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fx_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(-Fy_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-FyMinSteady_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-Fz_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')

legend('Fx steady','-Fy steady','-Fz steady','Fx freqMod','-Fy freqMod','-Fz freqMod','location','E')
xlabel('Amp ratio')
ylabel('F/mg')
axis([0.75 1 0 1.5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)

% M-Amp
subplot(1,2,2)
hold on

plot(Amp_ratio,MxMinSteady_norm,'sk','markersize',7,'markerfacecolor','b')
% plot(Amp_ratio,-MyMinSteady_norm,'sk','markersize',7,'markerfacecolor','r')
plot(Amp_ratio,-My_CoM_norm,'sk','markersize',7,'markerfacecolor','r')
plot(Amp_ratio,-MzMinSteady_norm,'sk','markersize',7,'markerfacecolor','g')

plot(Amp_ratio,MxMinSteady_norm_freqMod,'sk','markersize',7,'markerfacecolor','c')
% plot(Amp_ratio,-MyMinSteady_norm_freqMod,'sk','markersize',7,'markerfacecolor',[1 .5 0])
plot(Amp_ratio,-My_CoM_norm_freqMod,'sk','markersize',7,'markerfacecolor',[1 .5 0])
plot(Amp_ratio,-MzMinSteady_norm_freqMod,'sk','markersize',7,'markerfacecolor',[0 .5 0])

plot([min(Amp_ratio) max(Amp_ratio)],polyval(MxMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(-MyMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-My_CoM_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-MzMinSteady_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(MxMinSteady_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(-MyMinSteady_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-My_CoM_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(-MzMinSteady_Amp_fit_freqMod,[min(Amp_ratio) max(Amp_ratio)]),'k')

legend('Mx','-My@CoM','-Mz','location','NE')
xlabel('Amp ratio')
ylabel('T/mgl')
axis([0.75 1 -.1 .4])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1:.1:1)

mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit_YnZflip.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit_YnZflip.png'])
% saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit_YnZflip.svg'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit_YnZflip.svg'])

cd ..

%% save data
save('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF_INCclippedWingFreq',...
    'Amp_ratio',...
    ...
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
    'Fx_norm_freqMod',...
    'Fy_norm_freqMod',...
    'Fz_norm_freqMod',...
    'Mx_norm_freqMod',...
    'My_norm_freqMod',...
    'Mz_norm_freqMod',...
    'FyMinSteady_norm_freqMod',...
    'My_CoM_norm_freqMod',...
    'MxMinSteady_norm_freqMod',...
    'MyMinSteady_norm_freqMod',...
    'MzMinSteady_norm_freqMod',...
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
    'My_CoM_Amp_fit',...
    'My_CoM_Amp_fit_error',...
    ...
    'Fx_Amp_fit_freqMod',...
    'Fx_Amp_fit_error_freqMod',...
    'Fy_Amp_fit_freqMod',...
    'Fy_Amp_fit_error_freqMod',...
    'Fz_Amp_fit_freqMod',...
    'Fz_Amp_fit_error_freqMod',...
    'FyMinSteady_Amp_fit_freqMod',...
    'FyMinSteady_Amp_fit_error_freqMod',...
    'Mx_Amp_fit_freqMod',...
    'Mx_Amp_fit_error_freqMod',...
    'My_Amp_fit_freqMod',...
    'My_Amp_fit_error_freqMod',...
    'Mz_Amp_fit_freqMod',...
    'Mz_Amp_fit_error_freqMod',...
    'MxMinSteady_Amp_fit_freqMod',...
    'MxMinSteady_Amp_fit_error_freqMod',...
    'MyMinSteady_Amp_fit_freqMod',...
    'MyMinSteady_Amp_fit_error_freqMod',...
    'MzMinSteady_Amp_fit_freqMod',...
    'MzMinSteady_Amp_fit_error_freqMod',...
    'My_CoM_Amp_fit_freqMod',...
    'My_CoM_Amp_fit_error_freqMod');

