
clear
clc
close all

%% const
% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steadyWB F&M&CoM data
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

m=0;

%% cali CF for all (others have ERROR: Mx is assymetric)
% cali_template = 'cali_matrix_interp_*.mat';
cali_file = 'cali_matrix_interp_L_Full_R_Full_CF.mat';
load(cali_file)

%% load & save intact
cut_type_now = 0;  % no cut

% dirs = dir('L_Full_R_Full*')
% dirs = dir('L_Full_R_Full_CF*')
dirs = dir('L_Full_R_Full_Acrylic*');

% file_template = 'F_*_0.mat';
% file_template = 'F_*increase_0.mat';
file_template = 'F_*reduce_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name(1).name)
        
        m=m+1;
        cut_perc_now = 100;
        
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

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save distal data
cut_type_now = 1;  % TIP

dirs = dir('L_Full_R_Distal*');

file_template = 'F_*increase_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        m=m+1;
        cut_perc_now = str2num(dir_now.name(end-1:end));
        
        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now & BorfMorphCutData.cut_ratio == cut_perc(m,1)/100);
        
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

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save trailing edge data
cut_type_now = 2;  % TE

dirs = dir('L_Full_R_Trailing*');

file_template = 'F_*increase_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        m=m+1;
        cut_perc_now = str2num(dir_now.name(end-1:end));
        
        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now & BorfMorphCutData.cut_ratio == cut_perc(m,1)/100);
        
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
        
        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% calc My at CoG;
My_norm_CoG = My_norm - d_norm_steady*Fz_norm;

%% calc linear fits
% F-S
[Fx_S1_fit, Fx_S1_fit_error] = polyfit(S1_ratio,Fx_norm,1);
[Fy_S1_fit, Fy_S1_fit_error] = polyfit(S1_ratio,Fy_norm,1);
[Fz_S1_fit, Fz_S1_fit_error] = polyfit(S1_ratio,Fz_norm,1);

[Fx_S2_fit, Fx_S2_fit_error] = polyfit(S2_ratio,Fx_norm,1);
[Fy_S2_fit, Fy_S2_fit_error] = polyfit(S2_ratio,Fy_norm,1);
[Fz_S2_fit, Fz_S2_fit_error] = polyfit(S2_ratio,Fz_norm,1);

[Fx_S3_fit, Fx_S3_fit_error] = polyfit(S3_ratio,Fx_norm,1);
[Fy_S3_fit, Fy_S3_fit_error] = polyfit(S3_ratio,Fy_norm,1);
[Fz_S3_fit, Fz_S3_fit_error] = polyfit(S3_ratio,Fz_norm,1);

[FzPlus1_S1_fit, FzPlus1_S1_fit_error] = polyfit(S1_ratio,Fz_norm+1,1);
[FzPlus1_S2_fit, FzPlus1_S2_fit_error] = polyfit(S2_ratio,Fz_norm+1,1);
[FzPlus1_S3_fit, FzPlus1_S3_fit_error] = polyfit(S3_ratio,Fz_norm+1,1);

% M-S
[Mx_S1_fit, Mx_S1_fit_error] = polyfit(S1_ratio,Mx_norm,1);
[My_S1_fit, My_S1_fit_error] = polyfit(S1_ratio,My_norm,1);
[Mz_S1_fit, Mz_S1_fit_error] = polyfit(S1_ratio,Mz_norm,1);

[Mx_S2_fit, Mx_S2_fit_error] = polyfit(S2_ratio,Mx_norm,1);
[My_S2_fit, My_S2_fit_error] = polyfit(S2_ratio,My_norm,1);
[Mz_S2_fit, Mz_S2_fit_error] = polyfit(S2_ratio,Mz_norm,1);

[Mx_S3_fit, Mx_S3_fit_error] = polyfit(S3_ratio,Mx_norm,1);
[My_S3_fit, My_S3_fit_error] = polyfit(S3_ratio,My_norm,1);
[Mz_S3_fit, Mz_S3_fit_error] = polyfit(S3_ratio,Mz_norm,1);

[My_CoG_S1_fit, My_CoG_S1_fit_error] = polyfit(S1_ratio,My_norm_CoG,1);
[My_CoG_S2_fit, My_CoG_S2_fit_error] = polyfit(S2_ratio,My_norm_CoG,1);
[My_CoG_S3_fit, My_CoG_S3_fit_error] = polyfit(S3_ratio,My_norm_CoG,1);

[MxMinSteadyMx_S1_fit, MxMinSteadyMx_S1_fit_error] = polyfit(S1_ratio,Mx_norm-Mx_norm(cut_type==0),1);
[MyMinSteadyMy_S1_fit, MyMinSteadyMy_S1_fit_error] = polyfit(S1_ratio,My_norm-My_norm(cut_type==0),1);
[MzMinSteadyMz_S1_fit, MzMinSteadyMz_S1_fit_error] = polyfit(S1_ratio,Mz_norm-Mz_norm(cut_type==0),1);

[MxMinSteadyMx_S2_fit, MxMinSteadyMx_S2_fit_error] = polyfit(S2_ratio,Mx_norm-Mx_norm(cut_type==0),1);
[MyMinSteadyMy_S2_fit, MyMinSteadyMy_S2_fit_error] = polyfit(S2_ratio,My_norm-My_norm(cut_type==0),1);
[MzMinSteadyMz_S2_fit, MzMinSteadyMz_S2_fit_error] = polyfit(S2_ratio,Mz_norm-Mz_norm(cut_type==0),1);

[MxMinSteadyMx_S3_fit, MxMinSteadyMx_S3_fit_error] = polyfit(S3_ratio,Mx_norm-Mx_norm(cut_type==0),1);
[MyMinSteadyMy_S3_fit, MyMinSteadyMy_S3_fit_error] = polyfit(S3_ratio,My_norm-My_norm(cut_type==0),1);
[MzMinSteadyMz_S3_fit, MzMinSteadyMz_S3_fit_error] = polyfit(S3_ratio,Mz_norm-Mz_norm(cut_type==0),1);

[MyMinSteadyMy_CoG_S1_fit, MyMinSteadyMy_CoG_S1_fit_error] = polyfit(S1_ratio,My_norm_CoG-My_norm_CoG(cut_type==0),1);
[MyMinSteadyMy_CoG_S2_fit, MyMinSteadyMy_CoG_S2_fit_error] = polyfit(S2_ratio,My_norm_CoG-My_norm_CoG(cut_type==0),1);
[MyMinSteadyMy_CoG_S3_fit, MyMinSteadyMy_CoG_S3_fit_error] = polyfit(S3_ratio,My_norm_CoG-My_norm_CoG(cut_type==0),1);

%% plot data
mkdir('figures_CutWing_FnM_vs_S2nS3_steadyWB')
cd('figures_CutWing_FnM_vs_S2nS3_steadyWB')

%% plot F-S2 & M-S3 (My@"CoG")
% F-S2
figure
subplot(1,2,1)
hold on
plot(S2_ratio(cut_type==0),Fx_norm(cut_type==0),'sk','markerfacecolor','b')
plot(S2_ratio(cut_type==0),Fy_norm(cut_type==0),'sk','markerfacecolor','r')
plot(S2_ratio(cut_type==0),Fz_norm(cut_type==0),'sk','markerfacecolor','g')

plot(S2_ratio(cut_type==1),Fx_norm(cut_type==1),'ok','markerfacecolor','b')
plot(S2_ratio(cut_type==2),Fx_norm(cut_type==2),'dk','markerfacecolor','b')

plot(S2_ratio(cut_type==1),Fy_norm(cut_type==1),'ok','markerfacecolor','r')
plot(S2_ratio(cut_type==2),Fy_norm(cut_type==2),'dk','markerfacecolor','r')

plot(S2_ratio(cut_type==1),Fz_norm(cut_type==1),'ok','markerfacecolor','g')
plot(S2_ratio(cut_type==2),Fz_norm(cut_type==2),'dk','markerfacecolor','g')

% linear fits
plot([min(S2_ratio) max(S2_ratio)],polyval(Fx_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fy_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fz_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')

legend('x-axis','y-axis','z-axis','location','E')
xlabel('S2 ratio')
ylabel('F/mg')
axis([0 1 -1 .2])

% M-S3
subplot(1,2,2)
hold on

plot(S3_ratio(cut_type==0),Mx_norm(cut_type==0)-Mx_norm(cut_type==0),'sk','markerfacecolor','b')
plot(S3_ratio(cut_type==1),Mx_norm(cut_type==1)-Mx_norm(cut_type==0),'ok','markerfacecolor','b')
plot(S3_ratio(cut_type==2),Mx_norm(cut_type==2)-Mx_norm(cut_type==0),'dk','markerfacecolor','b')

plot(S3_ratio(cut_type==0),My_norm_CoG(cut_type==0),'sk','markerfacecolor','r')
plot(S3_ratio(cut_type==1),My_norm_CoG(cut_type==1),'ok','markerfacecolor','r')
plot(S3_ratio(cut_type==2),My_norm_CoG(cut_type==2),'dk','markerfacecolor','r')

plot(S3_ratio(cut_type==0),Mz_norm(cut_type==0)-Mz_norm(cut_type==0),'sk','markerfacecolor','g')
plot(S3_ratio(cut_type==1),Mz_norm(cut_type==1)-Mz_norm(cut_type==0),'ok','markerfacecolor','g')
plot(S3_ratio(cut_type==2),Mz_norm(cut_type==2)-Mz_norm(cut_type==0),'dk','markerfacecolor','g')

plot([min(S3_ratio) max(S3_ratio)],polyval(MxMinSteadyMx_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(My_CoG_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(MzMinSteadyMz_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')

legend('Intact wing','Tip cut','Trailing Edge')
xlabel('S3 ratio')
ylabel('T/mgl')
axis([0 1 -.05 .3])

saveas(gcf,['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoG.fig'])
saveas(gcf,['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoG.png'])
plot2svg(['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoG.svg'])

%% plot&calc linear fit for F-S & T-S
% F-S1
figure
subplot(2,3,1)
hold on
plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[Fx_norm(cut_type==1);Fx_norm(cut_type==0)],'ok','markerfacecolor','b')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[Fx_norm(cut_type==2);Fx_norm(cut_type==0)],'dk','markerfacecolor','b')

plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[Fy_norm(cut_type==1);Fy_norm(cut_type==0)],'ok','markerfacecolor','r')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[Fy_norm(cut_type==2);Fy_norm(cut_type==0)],'dk','markerfacecolor','r')

plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[Fz_norm(cut_type==1);Fz_norm(cut_type==0)]+1,'ok','markerfacecolor','g')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[Fz_norm(cut_type==2);Fz_norm(cut_type==0)]+1,'dk','markerfacecolor','g')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz+1 Tip','Fz+1 TE')
% legend('x','y','z')
legend('Tip','TE')

plot([min(S1_ratio) max(S1_ratio)],polyval(Fx_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')
plot([min(S1_ratio) max(S1_ratio)],polyval(Fy_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')
plot([min(S1_ratio) max(S1_ratio)],polyval(FzPlus1_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')

% xlabel('S1 ratio')
ylabel('F/mg')
axis([0 1 -.1 .6])

% F-S2
subplot(2,3,2)
hold on
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fx_norm(cut_type==1);Fx_norm(cut_type==0)],'ok','markerfacecolor','b')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fx_norm(cut_type==2);Fx_norm(cut_type==0)],'dk','markerfacecolor','b')

plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fy_norm(cut_type==1);Fy_norm(cut_type==0)],'ok','markerfacecolor','r')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fy_norm(cut_type==2);Fy_norm(cut_type==0)],'dk','markerfacecolor','r')

plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fz_norm(cut_type==1);Fz_norm(cut_type==0)]+1,'ok','markerfacecolor','g')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fz_norm(cut_type==2);Fz_norm(cut_type==0)]+1,'dk','markerfacecolor','g')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')

% linear fits
plot([min(S2_ratio) max(S2_ratio)],polyval(Fx_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(Fy_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(FzPlus1_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')

% xlabel('S2 ratio')
% ylabel('F/mg')
axis([0 1 -.1 .6])

% F-S3
subplot(2,3,3)
hold on
plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Fx_norm(cut_type==1);Fx_norm(cut_type==0)],'ok','markerfacecolor','b')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Fx_norm(cut_type==2);Fx_norm(cut_type==0)],'dk','markerfacecolor','b')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Fy_norm(cut_type==1);Fy_norm(cut_type==0)],'ok','markerfacecolor','r')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Fy_norm(cut_type==2);Fy_norm(cut_type==0)],'dk','markerfacecolor','r')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Fz_norm(cut_type==1);Fz_norm(cut_type==0)]+1,'ok','markerfacecolor','g')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Fz_norm(cut_type==2);Fz_norm(cut_type==0)]+1,'dk','markerfacecolor','g')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')

% linear fits
plot([min(S3_ratio) max(S3_ratio)],polyval(Fx_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(Fy_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(FzPlus1_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')

% xlabel('S3 ratio')
% ylabel('F/mg')
axis([0 1 -.1 .6])

%M-S1
subplot(2,3,4)
hold on
plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[Mx_norm(cut_type==1);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'ok','markerfacecolor','b')
plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[My_norm(cut_type==1);My_norm(cut_type==0)]-My_norm(cut_type==0),'ok','markerfacecolor','r')
plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[Mz_norm(cut_type==1);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'ok','markerfacecolor','g')
plot([S1_ratio(cut_type==1);S1_ratio(cut_type==0)],[My_norm_CoG(cut_type==1);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'ok','markerfacecolor','c')

plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[Mx_norm(cut_type==2);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'dk','markerfacecolor','b')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[My_norm(cut_type==2);My_norm(cut_type==0)]-My_norm(cut_type==0),'dk','markerfacecolor','r')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[Mz_norm(cut_type==2);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'dk','markerfacecolor','g')
plot([S1_ratio(cut_type==2);S1_ratio(cut_type==0)],[My_norm_CoG(cut_type==2);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'dk','markerfacecolor','c')

% legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
% legend('x','y','z','y_c_g')

% linear fits
plot([min(S1_ratio) max(S1_ratio)],polyval(MxMinSteadyMx_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')
plot([min(S1_ratio) max(S1_ratio)],polyval(MyMinSteadyMy_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')
plot([min(S1_ratio) max(S1_ratio)],polyval(MzMinSteadyMz_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')

plot([min(S1_ratio) max(S1_ratio)],polyval(MyMinSteadyMy_CoG_S1_fit,[min(S1_ratio) max(S1_ratio)]),'k')

xlabel('S1 ratio')
ylabel('T/mgl')
axis([0 1 -.05 .3])

%M-S2
subplot(2,3,5)
hold on
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Mx_norm(cut_type==1);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'ok','markerfacecolor','b')
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[My_norm(cut_type==1);My_norm(cut_type==0)]-My_norm(cut_type==0),'ok','markerfacecolor','r')
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Mz_norm(cut_type==1);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'ok','markerfacecolor','g')
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[My_norm_CoG(cut_type==1);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'ok','markerfacecolor','c')

plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Mx_norm(cut_type==2);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'dk','markerfacecolor','b')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[My_norm(cut_type==2);My_norm(cut_type==0)]-My_norm(cut_type==0),'dk','markerfacecolor','r')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Mz_norm(cut_type==2);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'dk','markerfacecolor','g')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[My_norm_CoG(cut_type==2);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'dk','markerfacecolor','c')

% legend('Tip','TE')
legend('x','y','z','y_c_g')

% linear fits
plot([min(S2_ratio) max(S2_ratio)],polyval(MxMinSteadyMx_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(MyMinSteadyMy_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')
plot([min(S2_ratio) max(S2_ratio)],polyval(MzMinSteadyMz_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')

plot([min(S2_ratio) max(S2_ratio)],polyval(MyMinSteadyMy_CoG_S2_fit,[min(S2_ratio) max(S2_ratio)]),'k')

xlabel('S2 ratio')
% ylabel('T/mgl')
axis([0 1 -.05 .3])

%M-S3
subplot(2,3,6)
hold on
plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mx_norm(cut_type==1);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'ok','markerfacecolor','b')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mx_norm(cut_type==2);Mx_norm(cut_type==0)]-Mx_norm(cut_type==0),'dk','markerfacecolor','b')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[My_norm(cut_type==1);My_norm(cut_type==0)]-My_norm(cut_type==0),'ok','markerfacecolor','r')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[My_norm(cut_type==2);My_norm(cut_type==0)]-My_norm(cut_type==0),'dk','markerfacecolor','r')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mz_norm(cut_type==1);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'ok','markerfacecolor','g')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mz_norm(cut_type==2);Mz_norm(cut_type==0)]-Mz_norm(cut_type==0),'dk','markerfacecolor','g')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[My_norm_CoG(cut_type==1);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'ok','markerfacecolor','c')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[My_norm_CoG(cut_type==2);My_norm_CoG(cut_type==0)]-My_norm_CoG(cut_type==0),'dk','markerfacecolor','c')

% legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')


plot([min(S3_ratio) max(S3_ratio)],polyval(MxMinSteadyMx_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(MyMinSteadyMy_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')
plot([min(S3_ratio) max(S3_ratio)],polyval(MzMinSteadyMz_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')

plot([min(S3_ratio) max(S3_ratio)],polyval(MyMinSteadyMy_CoG_S3_fit,[min(S3_ratio) max(S3_ratio)]),'k')

xlabel('S3 ratio')
% ylabel('T/mgl')
axis([0 1 -.05 .3])

saveas(gcf,['FnMvsS1S2S3_robofly_CutWing_steadyWB_INCcali.fig'])
saveas(gcf,['FnMvsS1S2S3_robofly_CutWing_steadyWB_INCcali.png'])
plot2svg(['FnMvsS1S2S3_robofly_CutWing_steadyWB_INCcali.svg'])

cd ..

%% save data
save('roboflyDB_CutWing_steadyWB_INCcaliCF',...
    'time',...
    'time_norm',...
    'fx_norm',...
    'fy_norm',...
    'fz_norm',...
    'mx_norm',...
    'my_norm',...
    'mz_norm',...
    'Fx_norm',...
    'Fy_norm',...
    'Fz_norm',...
    'Mx_norm',...
    'My_norm',...
    'Mz_norm',...
    'cut_ratio',...
    'cut_perc',...
    'cut_type',...
    'cut_perc_geom',...
    'cut_ratio_geom',...
    'cut_type_geom',...
    'l_ratio',...
    'CoA_ratio',...
    'A_ratio',...
    'S1_ratio',...
    'S2_ratio',...
    'S3_ratio',...
    'CoA_normL',...
    'A_normL',...
    'S1_normL',...
    'S2_normL',...
    'S3_normL',...
    'CoA_normA',...
    'l_normA',...
    'S1_normA',...
    'S2_normA',...
    'S3_normA',...
    'Fx_S1_fit',...
    'Fx_S1_fit_error',...
    'Fy_S1_fit',...
    'Fy_S1_fit_error',...
    'Fz_S1_fit',...
    'Fz_S1_fit_error',...
    'FzPlus1_S1_fit',...
    'FzPlus1_S1_fit_error',...
    'Fx_S2_fit',...
    'Fx_S2_fit_error',...
    'Fy_S2_fit',...
    'Fy_S2_fit_error',...
    'Fz_S2_fit',...
    'Fz_S2_fit_error',...
    'FzPlus1_S2_fit',...
    'FzPlus1_S2_fit_error',...
    'Fx_S3_fit',...
    'Fx_S3_fit_error',...
    'Fy_S3_fit',...
    'Fy_S3_fit_error',...
    'Fz_S3_fit',...
    'Fz_S3_fit_error',...
    'FzPlus1_S3_fit',...
    'FzPlus1_S3_fit_error',...
    'Mx_S1_fit',...
    'Mx_S1_fit_error',...
    'My_S1_fit',...
    'My_S1_fit_error',...
    'Mz_S1_fit',...
    'Mz_S1_fit_error',...
    'MxMinSteadyMx_S1_fit',...
    'MxMinSteadyMx_S1_fit_error',...
    'MyMinSteadyMy_S1_fit',...
    'MyMinSteadyMy_S1_fit_error',...
    'MzMinSteadyMz_S1_fit',...
    'MzMinSteadyMz_S1_fit_error',...
    'Mx_S2_fit',...
    'Mx_S2_fit_error',...
    'My_S2_fit',...
    'My_S2_fit_error',...
    'Mz_S2_fit',...
    'Mz_S2_fit_error',...
    'MxMinSteadyMx_S2_fit',...
    'MxMinSteadyMx_S2_fit_error',...
    'MyMinSteadyMy_S2_fit',...
    'MyMinSteadyMy_S2_fit_error',...
    'MzMinSteadyMz_S2_fit',...
    'MzMinSteadyMz_S2_fit_error',...
    'Mx_S3_fit',...
    'Mx_S3_fit_error',...
    'My_S3_fit',...
    'My_S3_fit_error',...
    'Mz_S3_fit',...
    'Mz_S3_fit_error',...
    'MxMinSteadyMx_S3_fit',...
    'MxMinSteadyMx_S3_fit_error',...
    'MyMinSteadyMy_S3_fit',...
    'MyMinSteadyMy_S3_fit_error',...
    'MzMinSteadyMz_S3_fit',...
    'MzMinSteadyMz_S3_fit_error');

