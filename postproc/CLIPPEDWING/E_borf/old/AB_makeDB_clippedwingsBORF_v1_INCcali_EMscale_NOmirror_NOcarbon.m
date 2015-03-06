
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

%% run data
Nwb = 10;
wb_start = 2;
wb_stop = 9;

m=0;
file_template = 'F_*_0.mat';
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
        
        clip_perc(m,1) = 100;
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save distal data
clip_type_now = 1;  % TIP
dirs = dir('L_Full_R_Distal*')

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        m=m+1;
        
        clip_perc(m,1) = str2num(dir_now.name(end-1:end));
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);


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
        
        m=m+1;
        
        clip_perc(m,1) = str2num(dir_now.name(end-1:end));
        clip_type(m,1) = clip_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);


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
    'clip_type');

figure
hold on
plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)],'og-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)],'*g-')

legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
xlabel('percentage clipped')
ylabel('F/mg')

figure
hold on
plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)],'o-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)],'*-')

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)]-mean(My_norm(clip_type==0)),[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)]-mean(My_norm(clip_type==0)),[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)],'og-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)],'*g-')

legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
xlabel('percentage clipped')
ylabel('T/mgl')

% 
% %% plot
% figure
% hold on
% plot(time,stroke_intact,'b')
% plot(time,stroke_clipped,'r')
% plot(time,rot_intact,'b')
% plot(time,rot_clipped,'r')
% plot(time,dev_intact,'b')
% plot(time,dev_clipped,'r')


