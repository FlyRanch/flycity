
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

%         time = t(n_start:n_stop);
% 
%         stroke_intact = kine(n_start:n_stop,5);
%         stroke_clipped = kine(n_start:n_stop,6);
%         rot_intact = kine(n_start:n_stop,1);
%         rot_clipped = kine(n_start:n_stop,4);
%         dev_intact = kine(n_start:n_stop,2);
%         dev_clipped = kine(n_start:n_stop,3);

        %% time        
        time_robo = t(1:n_wb);
        time(:,m) = time_robo/f_robo2fly;
        time_norm(:,m) = [0:1/(n_wb-1):1]';
        
        %% temporal dynamics
        fx_series =  ft(:,1);  % thrust, fwd pos
        fy_series =  ft(:,3);  % thrust, fwd pos
        fz_series =  -ft(:,2);  % thrust, fwd pos

        mx_series =  ft(:,4);  % thrust, fwd pos
        my_series =  ft(:,6);  % thrust, fwd pos
        mz_series =  -ft(:,5);  % thrust, fwd pos
        
        for wb = wb_start:wb_stop-1
            fx_wb(:,wb-wb_start+1) = fx_series((wb-1)*n_wb+1:wb*n_wb);
            fy_wb(:,wb-wb_start+1) = fy_series((wb-1)*n_wb+1:wb*n_wb);
            fz_wb(:,wb-wb_start+1) = fz_series((wb-1)*n_wb+1:wb*n_wb);
            
            mx_wb(:,wb-wb_start+1) = mx_series((wb-1)*n_wb+1:wb*n_wb);
            my_wb(:,wb-wb_start+1) = my_series((wb-1)*n_wb+1:wb*n_wb);
            mz_wb(:,wb-wb_start+1) = mz_series((wb-1)*n_wb+1:wb*n_wb);
        end
        fx_robo = nanmean(fx_wb,2);
        fy_robo = nanmean(fy_wb,2);
        fz_robo = nanmean(fz_wb,2);
        
        mx_robo = nanmean(mx_wb,2);
        my_robo = nanmean(my_wb,2);
        mz_robo = nanmean(mz_wb,2);

        fx_fly = F_robo2fly*fx_robo;
        fy_fly = F_robo2fly*fy_robo;
        fz_fly = F_robo2fly*fz_robo;

        mx_fly = M_robo2fly*mx_robo;
        my_fly = M_robo2fly*my_robo;
        mz_fly = M_robo2fly*mz_robo;

        fx_norm(:,m) = fx_fly / Mg_fly;
        fy_norm(:,m) = fy_fly / Mg_fly;
        fz_norm(:,m) = fz_fly / Mg_fly;

        mx_norm(:,m) = mx_fly / Mg_fly / Lwing;
        my_norm(:,m) = my_fly / Mg_fly / Lwing;
        mz_norm(:,m) = mz_fly / Mg_fly / Lwing;
        
        %% mean F&M
        Fx_robo =  mean(ft(n_start:n_stop,1));  % thrust, fwd pos
        Fy_robo =  mean(ft(n_start:n_stop,3));  % sideways, right pos
        Fz_robo = -mean(ft(n_start:n_stop,2));  % vertical, down pos

        Mx_robo =  mean(ft(n_start:n_stop,4));  % Roll, right pos
        My_robo =  mean(ft(n_start:n_stop,6));  % pitch, nose up pos
        Mz_robo = -mean(ft(n_start:n_stop,5));  % yaw, right pos

        Fx_fly = F_robo2fly*Fx_robo;
        Fy_fly = F_robo2fly*Fy_robo;
        Fz_fly = F_robo2fly*Fz_robo;

        Mx_fly = M_robo2fly*Mx_robo;
        My_fly = M_robo2fly*My_robo;
        Mz_fly = M_robo2fly*Mz_robo;

        Fx_norm(m,1) = Fx_fly / Mg_fly;
        Fy_norm(m,1) = Fy_fly / Mg_fly;
        Fz_norm(m,1) = Fz_fly / Mg_fly;

        Mx_norm(m,1) = Mx_fly / Mg_fly / Lwing;
        My_norm(m,1) = My_fly / Mg_fly / Lwing;
        Mz_norm(m,1) = Mz_fly / Mg_fly / Lwing;
        
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

%         time = t(n_start:n_stop);
% 
%         stroke_intact = kine(n_start:n_stop,5);
%         stroke_clipped = kine(n_start:n_stop,6);
%         rot_intact = kine(n_start:n_stop,1);
%         rot_clipped = kine(n_start:n_stop,4);
%         dev_intact = kine(n_start:n_stop,2);
%         dev_clipped = kine(n_start:n_stop,3);

        %% time        
        time_robo = t(1:n_wb);
        time(:,m) = time_robo/f_robo2fly;
        time_norm(:,m) = [0:1/(n_wb-1):1]';
        
        %% temporal dynamics
        fx_series =  ft(:,1);  % thrust, fwd pos
        fy_series =  ft(:,3);  % thrust, fwd pos
        fz_series =  -ft(:,2);  % thrust, fwd pos

        mx_series =  ft(:,4);  % thrust, fwd pos
        my_series =  ft(:,6);  % thrust, fwd pos
        mz_series =  -ft(:,5);  % thrust, fwd pos
        
        for wb = wb_start:wb_stop-1
            fx_wb(:,wb-wb_start+1) = fx_series((wb-1)*n_wb+1:wb*n_wb);
            fy_wb(:,wb-wb_start+1) = fy_series((wb-1)*n_wb+1:wb*n_wb);
            fz_wb(:,wb-wb_start+1) = fz_series((wb-1)*n_wb+1:wb*n_wb);
            
            mx_wb(:,wb-wb_start+1) = mx_series((wb-1)*n_wb+1:wb*n_wb);
            my_wb(:,wb-wb_start+1) = my_series((wb-1)*n_wb+1:wb*n_wb);
            mz_wb(:,wb-wb_start+1) = mz_series((wb-1)*n_wb+1:wb*n_wb);
        end
        fx_robo = nanmean(fx_wb,2);
        fy_robo = nanmean(fy_wb,2);
        fz_robo = nanmean(fz_wb,2);
        
        mx_robo = nanmean(mx_wb,2);
        my_robo = nanmean(my_wb,2);
        mz_robo = nanmean(mz_wb,2);

        fx_fly = F_robo2fly*fx_robo;
        fy_fly = F_robo2fly*fy_robo;
        fz_fly = F_robo2fly*fz_robo;

        mx_fly = M_robo2fly*mx_robo;
        my_fly = M_robo2fly*my_robo;
        mz_fly = M_robo2fly*mz_robo;

        fx_norm(:,m) = fx_fly / Mg_fly;
        fy_norm(:,m) = fy_fly / Mg_fly;
        fz_norm(:,m) = fz_fly / Mg_fly;

        mx_norm(:,m) = mx_fly / Mg_fly / Lwing;
        my_norm(:,m) = my_fly / Mg_fly / Lwing;
        mz_norm(:,m) = mz_fly / Mg_fly / Lwing;
        
        %% mean F&M
        Fx_robo =  mean(ft(n_start:n_stop,1));  % thrust, fwd pos
        Fy_robo =  mean(ft(n_start:n_stop,3));  % sideways, right pos
        Fz_robo = -mean(ft(n_start:n_stop,2));  % vertical, down pos

        Mx_robo =  mean(ft(n_start:n_stop,4));  % Roll, right pos
        My_robo =  mean(ft(n_start:n_stop,6));  % pitch, nose up pos
        Mz_robo = -mean(ft(n_start:n_stop,5));  % yaw, right pos

        Fx_fly = F_robo2fly*Fx_robo;
        Fy_fly = F_robo2fly*Fy_robo;
        Fz_fly = F_robo2fly*Fz_robo;

        Mx_fly = M_robo2fly*Mx_robo;
        My_fly = M_robo2fly*My_robo;
        Mz_fly = M_robo2fly*Mz_robo;

        Fx_norm(m,1) = Fx_fly / Mg_fly;
        Fy_norm(m,1) = Fy_fly / Mg_fly;
        Fz_norm(m,1) = Fz_fly / Mg_fly;

        Mx_norm(m,1) = Mx_fly / Mg_fly / Lwing;
        My_norm(m,1) = My_fly / Mg_fly / Lwing;
        Mz_norm(m,1) = Mz_fly / Mg_fly / Lwing;
        
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
        
%         time = t(n_start:n_stop);
% 
%         stroke_intact = kine(n_start:n_stop,5);
%         stroke_clipped = kine(n_start:n_stop,6);
%         rot_intact = kine(n_start:n_stop,1);
%         rot_clipped = kine(n_start:n_stop,4);
%         dev_intact = kine(n_start:n_stop,2);
%         dev_clipped = kine(n_start:n_stop,3);

        %% time        
        time_robo = t(1:n_wb);
        time(:,m) = time_robo/f_robo2fly;
        time_norm(:,m) = [0:1/(n_wb-1):1]';
        
        %% temporal dynamics
        fx_series =  ft(:,1);  % thrust, fwd pos
        fy_series =  ft(:,3);  % thrust, fwd pos
        fz_series =  -ft(:,2);  % thrust, fwd pos

        mx_series =  ft(:,4);  % thrust, fwd pos
        my_series =  ft(:,6);  % thrust, fwd pos
        mz_series =  -ft(:,5);  % thrust, fwd pos
        
        for wb = wb_start:wb_stop-1
            fx_wb(:,wb-wb_start+1) = fx_series((wb-1)*n_wb+1:wb*n_wb);
            fy_wb(:,wb-wb_start+1) = fy_series((wb-1)*n_wb+1:wb*n_wb);
            fz_wb(:,wb-wb_start+1) = fz_series((wb-1)*n_wb+1:wb*n_wb);
            
            mx_wb(:,wb-wb_start+1) = mx_series((wb-1)*n_wb+1:wb*n_wb);
            my_wb(:,wb-wb_start+1) = my_series((wb-1)*n_wb+1:wb*n_wb);
            mz_wb(:,wb-wb_start+1) = mz_series((wb-1)*n_wb+1:wb*n_wb);
        end
        fx_robo = nanmean(fx_wb,2);
        fy_robo = nanmean(fy_wb,2);
        fz_robo = nanmean(fz_wb,2);
        
        mx_robo = nanmean(mx_wb,2);
        my_robo = nanmean(my_wb,2);
        mz_robo = nanmean(mz_wb,2);

        fx_fly = F_robo2fly*fx_robo;
        fy_fly = F_robo2fly*fy_robo;
        fz_fly = F_robo2fly*fz_robo;

        mx_fly = M_robo2fly*mx_robo;
        my_fly = M_robo2fly*my_robo;
        mz_fly = M_robo2fly*mz_robo;

        fx_norm(:,m) = fx_fly / Mg_fly;
        fy_norm(:,m) = fy_fly / Mg_fly;
        fz_norm(:,m) = fz_fly / Mg_fly;

        mx_norm(:,m) = mx_fly / Mg_fly / Lwing;
        my_norm(:,m) = my_fly / Mg_fly / Lwing;
        mz_norm(:,m) = mz_fly / Mg_fly / Lwing;
        
        %% mean F&M
        Fx_robo =  mean(ft(n_start:n_stop,1));  % thrust, fwd pos
        Fy_robo =  mean(ft(n_start:n_stop,3));  % sideways, right pos
        Fz_robo = -mean(ft(n_start:n_stop,2));  % vertical, down pos

        Mx_robo =  mean(ft(n_start:n_stop,4));  % Roll, right pos
        My_robo =  mean(ft(n_start:n_stop,6));  % pitch, nose up pos
        Mz_robo = -mean(ft(n_start:n_stop,5));  % yaw, right pos

        Fx_fly = F_robo2fly*Fx_robo;
        Fy_fly = F_robo2fly*Fy_robo;
        Fz_fly = F_robo2fly*Fz_robo;

        Mx_fly = M_robo2fly*Mx_robo;
        My_fly = M_robo2fly*My_robo;
        Mz_fly = M_robo2fly*Mz_robo;

        Fx_norm(m,1) = Fx_fly / Mg_fly;
        Fy_norm(m,1) = Fy_fly / Mg_fly;
        Fz_norm(m,1) = Fz_fly / Mg_fly;

        Mx_norm(m,1) = Mx_fly / Mg_fly / Lwing;
        My_norm(m,1) = My_fly / Mg_fly / Lwing;
        Mz_norm(m,1) = Mz_fly / Mg_fly / Lwing;
        
        cd ..

    end
end

save('roboflyDB_ClippedWing_NOcali',...
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

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')

plot([clip_perc(clip_type==1);clip_perc(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)],'og-')
plot([clip_perc(clip_type==2);clip_perc(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)],'*g-')

legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
xlabel('percentage clipped')
ylabel('T/mgl')

figure
hold on
plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fx_norm(clip_type==1);Fx_norm(clip_type==0)],'o-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fx_norm(clip_type==2);Fx_norm(clip_type==0)],'*-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fy_norm(clip_type==1);Fy_norm(clip_type==0)],'or-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fy_norm(clip_type==2);Fy_norm(clip_type==0)],'*r-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Fz_norm(clip_type==1);Fz_norm(clip_type==0)],'og-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Fz_norm(clip_type==2);Fz_norm(clip_type==0)],'*g-')

legend('Fx Tip','Fx TE','Fy Tip','Fy TE','Fz Tip','Fz TE')
xlabel('S1 ratio')
ylabel('F/mg')

figure
hold on
plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Mx_norm(clip_type==1);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'o-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Mx_norm(clip_type==2);Mx_norm(clip_type==0)]-mean(Mx_norm(clip_type==0)),'*-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[My_norm(clip_type==1);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'or-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[My_norm(clip_type==2);My_norm(clip_type==0)]-mean(My_norm(clip_type==0)),'*r-')

plot([S1_ratio(clip_type==1);S1_ratio(clip_type==0)],[Mz_norm(clip_type==1);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'og-')
plot([S1_ratio(clip_type==2);S1_ratio(clip_type==0)],[Mz_norm(clip_type==2);Mz_norm(clip_type==0)]-mean(Mz_norm(clip_type==0)),'*g-')

legend('Mx Tip','Mx TE','My Tip','My TE','Mz Tip','Mz TE')
xlabel('S1 ratio')
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


