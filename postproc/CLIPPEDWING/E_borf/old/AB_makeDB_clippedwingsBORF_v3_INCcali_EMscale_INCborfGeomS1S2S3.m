
clear
clc
% close all

%% const
% fly data
var_file = dir('flyVar*');
load(var_file.name)

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
cali_template = 'cali_matrix_interp_*.mat';

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
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
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
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
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
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
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

save('roboflyDB_CutWing_INCcali',...
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
    'S3_normA');

%% F-S2 & T-S3
figure
subplot(1,2,1)
hold on
plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fx_norm(cut_type==1);Fx_norm(cut_type==0)],'o-')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fx_norm(cut_type==2);Fx_norm(cut_type==0)],'*-')

plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fy_norm(cut_type==1);Fy_norm(cut_type==0)],'or-')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fy_norm(cut_type==2);Fy_norm(cut_type==0)],'*r-')

plot([S2_ratio(cut_type==1);S2_ratio(cut_type==0)],[Fz_norm(cut_type==1);Fz_norm(cut_type==0)]+1,'og-')
plot([S2_ratio(cut_type==2);S2_ratio(cut_type==0)],[Fz_norm(cut_type==2);Fz_norm(cut_type==0)]+1,'*g-')

legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
xlabel('S2 ratio')
ylabel('F/mg')

subplot(1,2,2)
hold on

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mx_norm(cut_type==1);Mx_norm(cut_type==0)],'o-')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mx_norm(cut_type==2);Mx_norm(cut_type==0)],'*-')
% plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mx_norm(cut_type==1);Mx_norm(cut_type==0)]-mean(Mx_norm(cut_type==0)),'o-')
% plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mx_norm(cut_type==2);Mx_norm(cut_type==0)]-mean(Mx_norm(cut_type==0)),'*-')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[My_norm(cut_type==1);My_norm(cut_type==0)]-mean(My_norm(cut_type==0)),'or-')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[My_norm(cut_type==2);My_norm(cut_type==0)]-mean(My_norm(cut_type==0)),'*r-')

plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mz_norm(cut_type==1);Mz_norm(cut_type==0)],'og-')
plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mz_norm(cut_type==2);Mz_norm(cut_type==0)],'*g-')
% plot([S3_ratio(cut_type==1);S3_ratio(cut_type==0)],[Mz_norm(cut_type==1);Mz_norm(cut_type==0)]-mean(Mz_norm(cut_type==0)),'og-')
% plot([S3_ratio(cut_type==2);S3_ratio(cut_type==0)],[Mz_norm(cut_type==2);Mz_norm(cut_type==0)]-mean(Mz_norm(cut_type==0)),'*g-')

legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
xlabel('S3 ratio')
ylabel('T/mgl')

