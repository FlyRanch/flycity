
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

% dir_now = dir('L_Full_R_Full_Acrylic*');
dir_now = dir('L_Full_R_Full_CF*');

if dir_now.isdir == 1
    cd(dir_now.name)
    m=0;
    
    %% clipped amp increase
%     file_template = '*increase*.mat';
%     file_names = dir(file_template);
%     Amp_type_now = 1;
%     
%     for i = 1:length(file_names)
%     
%         file_now = file_names(i).name
%         load(file_now)
%         
%         m=m+1;
%         
%         Amp_type(m,1) = Amp_type_now;
%         cut_ratio(m,1) = cut_perc_now/100;
%         cut_perc(m,1) = cut_perc_now;
%         cut_type(m,1) = cut_type_now;
% 
%         n_total = length(t);
%         n_wb = n_total/Nwb;
%         n_start = round(wb_start*n_wb+1);
%         n_stop = round(wb_stop*n_wb);
%         
%         % borf geometry
%         n_geom = find(BorfMorphCutData.cut_type == cut_type_now);
%         
%         cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
%         cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
%         cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);
%         
%         l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
%         CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
%         A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
%         S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
%         S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
%         S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);
%         
%         CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
%         A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
%         S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
%         S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
%         S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);
%         
%         l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
%         CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
%         S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
%         S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
%         S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);
% 
%         extract_NONclipped_FnM_AmpIncrease_INCcali
% 
%     end
    
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

% fix Fz jump
Fz_norm_fix(n_Amp_increase,1) = Fz_norm(n_Amp_increase) - Fz_norm(n_Amp_increase_steady) - 1;
Fz_norm_fix(n_Amp_reduce,1) = Fz_norm(n_Amp_reduce) - Fz_norm(n_Amp_reduce_steady) - 1;

% Fy equilibrium
FyMinSteadyFy_norm(n_Amp_increase,1) = Fy_norm(n_Amp_increase) - Fy_norm(n_Amp_increase_steady);
FyMinSteadyFy_norm(n_Amp_reduce,1) = Fy_norm(n_Amp_reduce) - Fy_norm(n_Amp_reduce_steady);

% roll equilibrium
MxMinSteadyMx_norm(n_Amp_increase,1) = Mx_norm(n_Amp_increase) - Mx_norm(n_Amp_increase_steady);
MxMinSteadyMx_norm(n_Amp_reduce,1) = Mx_norm(n_Amp_reduce) - Mx_norm(n_Amp_reduce_steady);

% pitch equilibrium
MyMinSteadyMy_norm(n_Amp_increase,1) = My_norm(n_Amp_increase) - My_norm(n_Amp_increase_steady);
MyMinSteadyMy_norm(n_Amp_reduce,1) = My_norm(n_Amp_reduce) - My_norm(n_Amp_reduce_steady);

% yaw equilibrium
MzMinSteadyMz_norm(n_Amp_increase,1) = Mz_norm(n_Amp_increase) - Mz_norm(n_Amp_increase_steady);
MzMinSteadyMz_norm(n_Amp_reduce,1) = Mz_norm(n_Amp_reduce) - Mz_norm(n_Amp_reduce_steady);

%% calc linear fits
% F-Aratio
[Fx_Amp_fit, Fx_Amp_fit_error] = polyfit(Amp_ratio,Fx_norm,1);
[Fy_Amp_fit, Fy_Amp_fit_error] = polyfit(Amp_ratio,Fy_norm,1);
[Fz_Amp_fit, Fz_Amp_fit_error] = polyfit(Amp_ratio,Fz_norm,1);

[FyMinSteadyFy_Amp_fit, FyMinSteadyFy_Amp_fit_error] = polyfit(Amp_ratio,FyMinSteadyFy_norm,1);
[Fz_fix_Amp_fit, Fz_fix_Amp_fit_error] = polyfit(Amp_ratio,Fz_norm_fix,1);

% M-S
[Mx_Amp_fit, Mx_Amp_fit_error] = polyfit(Amp_ratio,Mx_norm,1);
[My_Amp_fit, My_Amp_fit_error] = polyfit(Amp_ratio,My_norm,1);
[Mz_Amp_fit, Mz_Amp_fit_error] = polyfit(Amp_ratio,Mz_norm,1);

[MxMinSteadyMx_Amp_fit, MxMinSteadyMx_Amp_fit_error] = polyfit(Amp_ratio,MxMinSteadyMx_norm,1);
[MyMinSteadyMy_Amp_fit, MyMinSteadyMy_Amp_fit_error] = polyfit(Amp_ratio,MyMinSteadyMy_norm,1);
[MzMinSteadyMz_Amp_fit, MzMinSteadyMz_Amp_fit_error] = polyfit(Amp_ratio,MzMinSteadyMz_norm,1);

% My second order fit
[My_Amp_fit2, My_Amp_fit2_error] = polyfit(Amp_ratio,My_norm,2);
[My_CoG_Amp_fit2, My_CoG_Amp_fit2_error] = polyfit(Amp_ratio,My_norm_CoG,2);
[MyMinSteadyMy_Amp_fit2, MyMinSteadyMy_Amp_fit2_error] = polyfit(Amp_ratio,MyMinSteadyMy_norm,2);

%% plot F-Amp & M-Amp (My@"CoG")
% F-Amp
figure
subplot(1,2,1)
hold on
plot(Amp_ratio,Fx_norm,'ok','markerfacecolor','b')
% plot(Amp_ratio,Fy_norm,'ok','markerfacecolor','r')
% plot(Amp_ratio,Fz_norm,'dk','markerfacecolor','g')

plot(Amp_ratio,FyMinSteadyFy_norm,'ok','markerfacecolor','r')
plot(Amp_ratio,Fz_norm_fix,'ok','markerfacecolor','g')

% linear fits
plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fx_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fy_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fz_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(FyMinSteadyFy_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(Fz_fix_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

% legend('Fx','Fy','Fz','location','SW')
legend('Fx','Fy','Fz fix','location','E')
xlabel('Amp ratio')
ylabel('F/mg')
axis([0.75 1 -1 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)

% M-Amp
subplot(1,2,2)
hold on

plot(Amp_ratio,MxMinSteadyMx_norm,'ok','markerfacecolor','b')

% plot(Amp_ratio,MyMinSteadyMy_norm,'ok','markerfacecolor','g')
plot(Amp_ratio,My_norm_CoG,'dk','markerfacecolor','r')

plot(Amp_ratio,MzMinSteadyMz_norm,'ok','markerfacecolor','g')

plot([min(Amp_ratio) max(Amp_ratio)],polyval(MxMinSteadyMx_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
% plot([min(Amp_ratio) max(Amp_ratio)],polyval(My_CoG_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')
plot([min(Amp_ratio) max(Amp_ratio)],polyval(MzMinSteadyMz_Amp_fit,[min(Amp_ratio) max(Amp_ratio)]),'k')

plot([min(Amp_ratio):.01:max(Amp_ratio)],polyval(My_CoG_Amp_fit2,[min(Amp_ratio):.01:max(Amp_ratio)]),'k')

legend('Mx','My@CoM','Mz','location','NE')
xlabel('Amp ratio')
ylabel('T/mgl')
axis([0.75 1 -.05 .2])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1:.05:1)

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoG_robofly_NONcutWing_INCcali_REDUCE_CF.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoG_robofly_NONcutWing_INCcali_REDUCE_CF.png'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoG_robofly_NONcutWing_INCcali_REDUCE_CF.svg'])

%% save data
save('roboflyDB_AmpStrokeRatio_NONcutWing_INCcali_REDUCE_CF',...
    'Amp_ratio',...
    'Fx_norm',...
    'Fy_norm',...
    'Fz_norm',...
    'Mx_norm',...
    'My_norm',...
    'Mz_norm',...
    'Fz_norm_fix',...
    'FyMinSteadyFy_norm',...
    'My_norm_CoG',...
    'MxMinSteadyMx_norm',...
    'MyMinSteadyMy_norm',...
    'MzMinSteadyMz_norm',...
    'Fx_Amp_fit',...
    'Fx_Amp_fit_error',...
    'Fy_Amp_fit',...
    'Fy_Amp_fit_error',...
    'Fz_Amp_fit',...
    'Fz_Amp_fit_error',...
    'FyMinSteadyFy_Amp_fit',...
    'FyMinSteadyFy_Amp_fit_error',...
    'Fz_fix_Amp_fit',...
    'Fz_fix_Amp_fit_error',...
    'Mx_Amp_fit',...
    'Mx_Amp_fit_error',...
    'My_Amp_fit',...
    'My_Amp_fit_error',...
    'Mz_Amp_fit',...
    'Mz_Amp_fit_error',...
    'MxMinSteadyMx_Amp_fit',...
    'MxMinSteadyMx_Amp_fit_error',...
    'MyMinSteadyMy_Amp_fit',...
    'MyMinSteadyMy_Amp_fit_error',...
    'MzMinSteadyMz_Amp_fit',...
    'MzMinSteadyMz_Amp_fit_error',...
    'My_Amp_fit2',...
    'My_Amp_fit2_error',...
    'MyMinSteadyMy_Amp_fit2',...
    'MyMinSteadyMy_Amp_fit2_error',...
    'My_CoG_Amp_fit2',...
    'My_CoG_Amp_fit2_error');

