% make cali file
clear
clc
close all

%% settings
start = 6000;
stop = 2000;
step_interp = 1;

settings.start=start;
settings.stop=stop;
settings.step_interp=step_interp;

dirs = dir('L_Full_R*');

for d = 1:length(dirs)
    if dirs(d).isdir==1
        cd(dirs(d).name)
        name=['cali_matrix_interp',dirs(d).name(9:end)];

%% left wing: fly ref NOT robo ref
left_name = dir('g_cal_left*.mat');
load(left_name.name)

% left wing stroke
str_run = radtodeg(kine(:,6));

% convert to double
ft = double(ft);

cali_L.stroke=radtodeg(stroke)';

cali_L.rotation=radtodeg(rotation);
cali_L.deviation=radtodeg(deviation);

cali_L.kin_list(:,1)=cali_L.stroke;
cali_L.kin_list(:,2)=cali_L.rotation;
cali_L.kin_list(:,3)=cali_L.deviation;

cali_L.Fx_run = ft(:,1);
cali_L.Fy_run = ft(:,3);
cali_L.Fz_run = -ft(:,2);

cali_L.Mx_run = ft(:,4);
cali_L.My_run = ft(:,6);
cali_L.Mz_run = -ft(:,5);

for i=1:length(cali_L.stroke)
    
    str_now = cali_L.stroke(i);
    
    str_loc=find(str_run==str_now);
    str_loc_start = str_loc(start);
    str_loc_stop = str_loc(end-stop);
    
    if str_run(str_loc_start:str_loc_stop) == str_now
    else
        ERROR=1
        return
    end
    
    cali_L.Fx(i,1) = mean(cali_L.Fx_run(str_loc_start:str_loc_stop));
    cali_L.Fy(i,1) = mean(cali_L.Fy_run(str_loc_start:str_loc_stop));
    cali_L.Fz(i,1) = mean(cali_L.Fz_run(str_loc_start:str_loc_stop));
    
    cali_L.Mx(i,1) = mean(cali_L.Mx_run(str_loc_start:str_loc_stop));
    cali_L.My(i,1) = mean(cali_L.My_run(str_loc_start:str_loc_stop));
    cali_L.Mz(i,1) = mean(cali_L.Mz_run(str_loc_start:str_loc_stop));
end
%% interp
cali_L.stroke_interp = [min(cali_L.stroke):step_interp:max(cali_L.stroke)]';

cali_L.Fx_interp = interp1(cali_L.stroke,cali_L.Fx,cali_L.stroke_interp,'spline');
cali_L.Fy_interp = interp1(cali_L.stroke,cali_L.Fy,cali_L.stroke_interp,'spline');
cali_L.Fz_interp = interp1(cali_L.stroke,cali_L.Fz,cali_L.stroke_interp,'spline');

cali_L.Mx_interp = interp1(cali_L.stroke,cali_L.Mx,cali_L.stroke_interp,'spline');
cali_L.My_interp = interp1(cali_L.stroke,cali_L.My,cali_L.stroke_interp,'spline');
cali_L.Mz_interp = interp1(cali_L.stroke,cali_L.Mz,cali_L.stroke_interp,'spline');

figure
hold on
plot(cali_L.stroke,cali_L.Fx,'ob')
plot(cali_L.stroke,cali_L.Fy,'or')
plot(cali_L.stroke,cali_L.Fz,'og')
plot(cali_L.stroke,cali_L.Mx,'oc')
plot(cali_L.stroke,cali_L.My,'om')
plot(cali_L.stroke,cali_L.Mz,'ok')

legend Fx Fy Fz Mx My Mz
title([name,'_L'])

plot(cali_L.stroke_interp,cali_L.Fx_interp,'-b')
plot(cali_L.stroke_interp,cali_L.Fy_interp,'-r')
plot(cali_L.stroke_interp,cali_L.Fz_interp,'-g')
plot(cali_L.stroke_interp,cali_L.Mx_interp,'-c')
plot(cali_L.stroke_interp,cali_L.My_interp,'-m')
plot(cali_L.stroke_interp,cali_L.Mz_interp,'-k')


clear loc i stroke rotation deviation kine ft t str_run str_now str_loc str_loc_start str_loc_stop

%% right wing: fly ref NOT robo ref
right_name = dir('g_cal_right*.mat');
load(right_name.name)

% right wing stroke
str_run = radtodeg(kine(:,5));

% convert to double
ft = double(ft);

cali_R.stroke=radtodeg(stroke)';
cali_R.rotation=radtodeg(rotation)';
cali_R.deviation=radtodeg(deviation)';

cali_R.kin_list(:,1)=cali_R.stroke;
cali_R.kin_list(:,2)=cali_R.rotation;
cali_R.kin_list(:,3)=cali_R.deviation;

cali_R.Fx_run = ft(:,1);
cali_R.Fy_run = ft(:,3);
cali_R.Fz_run = -ft(:,2);

cali_R.Mx_run = ft(:,4);
cali_R.My_run = ft(:,6);
cali_R.Mz_run = -ft(:,5);

for i=1:length(cali_R.stroke)
    
    str_now = cali_R.stroke(i);
    
    str_loc=find(str_run==str_now);
    str_loc_start = str_loc(start);
    str_loc_stop = str_loc(end-stop);
    
    if str_run(str_loc_start:str_loc_stop) == str_now
    else
        ERROR=1
        return
    end
    
    cali_R.Fx(i,1) = mean(cali_R.Fx_run(str_loc_start:str_loc_stop));
    cali_R.Fy(i,1) = mean(cali_R.Fy_run(str_loc_start:str_loc_stop));
    cali_R.Fz(i,1) = mean(cali_R.Fz_run(str_loc_start:str_loc_stop));
    
    cali_R.Mx(i,1) = mean(cali_R.Mx_run(str_loc_start:str_loc_stop));
    cali_R.My(i,1) = mean(cali_R.My_run(str_loc_start:str_loc_stop));
    cali_R.Mz(i,1) = mean(cali_R.Mz_run(str_loc_start:str_loc_stop));
end

%% interp
cali_R.stroke_interp = [min(cali_R.stroke):step_interp:max(cali_R.stroke)]';

cali_R.Fx_interp = interp1(cali_R.stroke,cali_R.Fx,cali_R.stroke_interp,'spline');
cali_R.Fy_interp = interp1(cali_R.stroke,cali_R.Fy,cali_R.stroke_interp,'spline');
cali_R.Fz_interp = interp1(cali_R.stroke,cali_R.Fz,cali_R.stroke_interp,'spline');

cali_R.Mx_interp = interp1(cali_R.stroke,cali_R.Mx,cali_R.stroke_interp,'spline');
cali_R.My_interp = interp1(cali_R.stroke,cali_R.My,cali_R.stroke_interp,'spline');
cali_R.Mz_interp = interp1(cali_R.stroke,cali_R.Mz,cali_R.stroke_interp,'spline');


figure
hold on
plot(cali_R.stroke,cali_R.Fx,'ob')
plot(cali_R.stroke,cali_R.Fy,'or')
plot(cali_R.stroke,cali_R.Fz,'og')
plot(cali_R.stroke,cali_R.Mx,'oc')
plot(cali_R.stroke,cali_R.My,'om')
plot(cali_R.stroke,cali_R.Mz,'ok')

legend Fx Fy Fz Mx My Mz
title([name,'_R'])

plot(cali_R.stroke_interp,cali_R.Fx_interp,'-b')
plot(cali_R.stroke_interp,cali_R.Fy_interp,'-r')
plot(cali_R.stroke_interp,cali_R.Fz_interp,'-g')
plot(cali_R.stroke_interp,cali_R.Mx_interp,'-c')
plot(cali_R.stroke_interp,cali_R.My_interp,'-m')
plot(cali_R.stroke_interp,cali_R.Mz_interp,'-k')

clear loc i stroke rotation deviation kine ft t str_run str_now str_loc str_loc_start str_loc_stop

save(name,'cali_L','cali_R','settings')
cd ..
    end
end
