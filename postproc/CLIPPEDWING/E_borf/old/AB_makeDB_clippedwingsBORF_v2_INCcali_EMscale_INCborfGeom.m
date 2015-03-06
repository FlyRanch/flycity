
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

%% clipped wing geometry
geom_file = dir('borf_clipped_wing_geometry*');
borf_geom = load(geom_file.name)

borf_geom.CoA_norm = borf_geom.CoA_mm/(Lrobo*1e3);
borf_geom.A_norm = borf_geom.A_mm2/(Lrobo*1e3)^2;
borf_geom.S1_norm = borf_geom.S1_mm3/(Lrobo*1e3)^3;
borf_geom.S2_norm = borf_geom.S2_mm4/(Lrobo*1e3)^4;

%% run data
Nwb = 10;
wb_start = 2;
wb_stop = 9;

m=0;
file_template = 'F_*_0.mat';
file_template = 'F_*_increase_1.mat';
cali_template = 'cali_matrix_interp_*.mat';

%% load & save intact
clip_type_now = 0;  % TE
% dirs = dir('L_Full_R_Full*')
dirs = dir('L_Full_R_Full_Acrylic*')

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name(1).name)
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
        m=m+1;
        clip_perc_now = 100;
        
        clip_perc(m,1) = clip_perc_now;
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);
        
        % borf geometry
        n_geom = find(borf_geom.clip_type == clip_type_now);
        
        clip_ratio_geom(m,1) = borf_geom.clip_ratio(n_geom);
        clip_type_geom(m,1) = borf_geom.clip_type(n_geom);
        
        CoA_ratio(m,1) = borf_geom.CoA_ratio(n_geom);
        A_ratio(m,1) = borf_geom.A_ratio(n_geom);
        S1_ratio(m,1) = borf_geom.S1_ratio(n_geom);
        S2_ratio(m,1) = borf_geom.S2_ratio(n_geom);
        
        CoA_norm(m,1) = borf_geom.CoA_norm(n_geom);
        A_norm(m,1) = borf_geom.A_norm(n_geom);
        S1_norm(m,1) = borf_geom.S1_norm(n_geom);
        S2_norm(m,1) = borf_geom.S2_norm(n_geom);

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save distal data
clip_type_now = 1;  % TIP
dirs = dir('L_Full_R_Distal*');

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
        m=m+1;
        clip_perc_now = str2num(dir_now.name(end-1:end));
        
        clip_perc(m,1) = clip_perc_now;
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(borf_geom.clip_type == clip_type_now & borf_geom.clip_ratio == clip_perc(m,1)/100);
        
        clip_ratio_geom(m,1) = borf_geom.clip_ratio(n_geom);
        clip_type_geom(m,1) = borf_geom.clip_type(n_geom);
        
        CoA_ratio(m,1) = borf_geom.CoA_ratio(n_geom);
        A_ratio(m,1) = borf_geom.A_ratio(n_geom);
        S1_ratio(m,1) = borf_geom.S1_ratio(n_geom);
        S2_ratio(m,1) = borf_geom.S2_ratio(n_geom);
        
        CoA_norm(m,1) = borf_geom.CoA_norm(n_geom);
        A_norm(m,1) = borf_geom.A_norm(n_geom);
        S1_norm(m,1) = borf_geom.S1_norm(n_geom);
        S2_norm(m,1) = borf_geom.S2_norm(n_geom);

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save trailing edge data
clip_type_now = 2;  % TE
dirs = dir('L_Full_R_Trailing*')

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        cali_name = dir(cali_template);
        load(cali_name(1).name)
        
        m=m+1;
        clip_perc_now = str2num(dir_now.name(end-1:end));
        
        clip_perc(m,1) = clip_perc_now;
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(borf_geom.clip_type == clip_type_now & borf_geom.clip_ratio == clip_perc(m,1)/100);
        
        clip_ratio_geom(m,1) = borf_geom.clip_ratio(n_geom);
        clip_type_geom(m,1) = borf_geom.clip_type(n_geom);
        
        CoA_ratio(m,1) = borf_geom.CoA_ratio(n_geom);
        A_ratio(m,1) = borf_geom.A_ratio(n_geom);
        S1_ratio(m,1) = borf_geom.S1_ratio(n_geom);
        S2_ratio(m,1) = borf_geom.S2_ratio(n_geom);
        
        CoA_norm(m,1) = borf_geom.CoA_norm(n_geom);
        A_norm(m,1) = borf_geom.A_norm(n_geom);
        S1_norm(m,1) = borf_geom.S1_norm(n_geom);
        S2_norm(m,1) = borf_geom.S2_norm(n_geom);
        
        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

save('roboflyDB_ClippedWing_INCcali',...
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
    'clip_perc',...
    'clip_type',...
    'clip_ratio_geom',...
    'clip_type_geom',...
    'CoA_ratio',...
    'A_ratio',...
    'S1_ratio',...
    'S2_ratio',...
    'CoA_norm',...
    'A_norm',...
    'S1_norm',...
    'S2_norm');

%% figures

% %% CoA
% figure
% hold on
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')
% 
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')
% 
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)],'og-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)],'*g-')
% 
% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
% xlabel('CoA ratio')
% ylabel('F/mg')
% 
% figure
% hold on
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'o-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'*-')
% 
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')
% 
% plot([CoA_ratio(clip_type==1);CoA_ratio(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'og-')
% plot([CoA_ratio(clip_type==2);CoA_ratio(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'*g-')
% 
% legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
% xlabel('nanmean ratio')
% ylabel('T/mgl')
% 
% %% A
% figure
% hold on
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')
% 
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')
% 
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)],'og-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)],'*g-')
% 
% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
% xlabel('A ratio')
% ylabel('F/mg')
% 
% figure
% hold on
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'o-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'*-')
% 
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')
% 
% plot([A_ratio(clip_type==1);A_ratio(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'og-')
% plot([A_ratio(clip_type==2);A_ratio(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'*g-')
% 
% legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
% xlabel('A ratio')
% ylabel('T/mgl')

%% S1
figure
subplot(2,2,1)
hold on
plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)]+1,'og-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)]+1,'*g-')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
xlabel('S1 ratio')
ylabel('F/mg')

subplot(2,2,2)
hold on
plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'o-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'*-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'og-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'*g-')

% legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
xlabel('S1 ratio')
ylabel('T/mgl')

%% S2
subplot(2,2,3)
hold on
plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')

plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')

plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)]+1,'og-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)]+1,'*g-')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
xlabel('S2 ratio')
ylabel('F/mg')

subplot(2,2,4)
hold on
plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'o-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'*-')

plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')

plot([S2_ratio(clip_type==1);S2_ratio(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'og-')
plot([S2_ratio(clip_type==2);S2_ratio(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'*g-')

% legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
legend('F&Mx Tip','F&Mx TE','F&My Tip','F&My TE','F&Mz Tip','F&Mz TE')
xlabel('S2 ratio')
ylabel('T/mgl')

