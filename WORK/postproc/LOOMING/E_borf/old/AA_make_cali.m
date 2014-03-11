% make cali file
clear
clc
close all

loc = cd;
name=['cali_matrix_',loc(end-7:end)]
start = 3000;

%% left wing: fly ref NOT robo ref
list = dir('g_cal_left*.mat')
for i=1:length(list)
    load(list(i).name)
    
    cali_L.str_list(i,1)=rad2deg(stroke);
    cali_L.rot_list(i,1)=rad2deg(rotation);
    cali_L.dev_list(i,1)=rad2deg(deviation);
    
    cali_L.kin_list(i,1)=rad2deg(stroke);
    cali_L.kin_list(i,2)=rad2deg(rotation);
    cali_L.kin_list(i,3)=rad2deg(deviation);
    
    cali_L.Fx_list(i,1)=mean(ft(start:end,1));
    cali_L.Fx_list(i,2)=std(ft(start:end,1));
    
    cali_L.Fy_list(i,1)=mean(ft(start:end,3));
    cali_L.Fy_list(i,2)=std(ft(start:end,3));
    
    cali_L.Fz_list(i,1)=mean(-ft(start:end,2));
    cali_L.Fz_list(i,2)=std(-ft(start:end,2));
    
    cali_L.Mx_list(i,1)=mean(ft(start:end,4));
    cali_L.Mx_list(i,2)=std(ft(start:end,4));
    
    cali_L.My_list(i,1)=mean(ft(start:end,6));
    cali_L.My_list(i,2)=std(ft(start:end,6));
    
    cali_L.Mz_list(i,1)=mean(-ft(start:end,5));
    cali_L.Mz_list(i,2)=std(-ft(start:end,5));
end    
    
cali_L.str_steps = sort(unique(cali_L.str_list));
cali_L.rot_steps = sort(unique(cali_L.rot_list));
cali_L.dev_steps = sort(unique(cali_L.dev_list));

[cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh] = meshgrid(cali_L.str_steps,cali_L.rot_steps,cali_L.dev_steps);

for i=1:length(cali_L.str_list)
    
    str_now = cali_L.str_list(i);
    rot_now = cali_L.rot_list(i);
    dev_now = cali_L.dev_list(i);
    
    str_loc=find(cali_L.str_steps==str_now);
    rot_loc=find(cali_L.rot_steps==rot_now);
    dev_loc=find(cali_L.dev_steps==dev_now);
    
%     cali_L.str_array(rot_loc,str_loc,dev_loc) = cali_L.str(i,1);
%     cali_L.rot_array(rot_loc,str_loc,dev_loc) = cali_L.rot(i,1);
%     cali_L.dev_array(rot_loc,str_loc,dev_loc) = cali_L.dev(i,1);
    
    cali_L.Fx_array(rot_loc,str_loc,dev_loc) = cali_L.Fx_list(i,1);
    cali_L.Fy_array(rot_loc,str_loc,dev_loc) = cali_L.Fy_list(i,1);
    cali_L.Fz_array(rot_loc,str_loc,dev_loc) = cali_L.Fz_list(i,1);
    
    cali_L.Mx_array(rot_loc,str_loc,dev_loc) = cali_L.Mx_list(i,1);
    cali_L.My_array(rot_loc,str_loc,dev_loc) = cali_L.My_list(i,1);
    cali_L.Mz_array(rot_loc,str_loc,dev_loc) = cali_L.Mz_list(i,1);
    
end

%% right wing: fly ref NOT robo ref
list = dir('g_cal_right*.mat')
for i=1:length(list)
    load(list(i).name)
    
    cali_R.str_list(i,1)=rad2deg(stroke);
    cali_R.rot_list(i,1)=rad2deg(rotation);
    cali_R.dev_list(i,1)=rad2deg(deviation);
    
    cali_R.kin_list(i,1)=rad2deg(stroke);
    cali_R.kin_list(i,2)=rad2deg(rotation);
    cali_R.kin_list(i,3)=rad2deg(deviation);
    
    cali_R.Fx_list(i,1)=mean(ft(start:end,1));
    cali_R.Fx_list(i,2)=std(ft(start:end,1));
    
    cali_R.Fy_list(i,1)=mean(ft(start:end,2));
    cali_R.Fy_list(i,2)=std(ft(start:end,2));
    
    cali_R.Fz_list(i,1)=mean(ft(start:end,3));
    cali_R.Fz_list(i,2)=std(ft(start:end,3));
    
    cali_R.Mx_list(i,1)=mean(ft(start:end,4));
    cali_R.Mx_list(i,2)=std(ft(start:end,4));
    
    cali_R.My_list(i,1)=mean(ft(start:end,5));
    cali_R.My_list(i,2)=std(ft(start:end,5));
    
    cali_R.Mz_list(i,1)=mean(ft(start:end,6));
    cali_R.Mz_list(i,2)=std(ft(start:end,6));
end    
    
cali_R.str_steps = sort(unique(cali_R.str_list));
cali_R.rot_steps = sort(unique(cali_R.rot_list));
cali_R.dev_steps = sort(unique(cali_R.dev_list));

[cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh] = meshgrid(cali_R.str_steps,cali_R.rot_steps,cali_R.dev_steps);

for i=1:length(cali_R.str_list)
    
    str_now = cali_R.str_list(i);
    rot_now = cali_R.rot_list(i);
    dev_now = cali_R.dev_list(i);
    
    str_loc=find(cali_R.str_steps==str_now);
    rot_loc=find(cali_R.rot_steps==rot_now);
    dev_loc=find(cali_R.dev_steps==dev_now);
    
%     cali_R.str_array(rot_loc,str_loc,dev_loc) = cali_R.str(i,1);
%     cali_R.rot_array(rot_loc,str_loc,dev_loc) = cali_R.rot(i,1);
%     cali_R.dev_array(rot_loc,str_loc,dev_loc) = cali_R.dev(i,1);
    
    cali_R.Fx_array(rot_loc,str_loc,dev_loc) = cali_R.Fx_list(i,1);
    cali_R.Fy_array(rot_loc,str_loc,dev_loc) = cali_R.Fy_list(i,1);
    cali_R.Fz_array(rot_loc,str_loc,dev_loc) = cali_R.Fz_list(i,1);
    
    cali_R.Mx_array(rot_loc,str_loc,dev_loc) = cali_R.Mx_list(i,1);
    cali_R.My_array(rot_loc,str_loc,dev_loc) = cali_R.My_list(i,1);
    cali_R.Mz_array(rot_loc,str_loc,dev_loc) = cali_R.Mz_list(i,1);
    
end

clear list loc i stroke rotation deviation ft t ans str_loc rot_loc dev_loc str_now rot_now dev_now
cd ..
save(name)

% plot cali
figure(1)
subplot(3,2,1)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fx_array,[0],[0],[0])
title('Fx_L')
subplot(3,2,3)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fy_array,[0],[0],[0])
title('Fy_L')
subplot(3,2,5)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fz_array,[0],[0],[0])
title('Fz_L')

subplot(3,2,2)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fx_array,[0],[0],[0])
title('Fx_R')
subplot(3,2,4)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fy_array,[0],[0],[0])
title('Fy_R')
subplot(3,2,6)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fz_array,[0],[0],[0])
title('Fz_R')
saveas(gca,[name,'_F.png'])

figure(2)
subplot(3,2,1)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mx_array,[0],[0],[0])
title('Mx_L')
subplot(3,2,3)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.My_array,[0],[0],[0])
title('My_L')
subplot(3,2,5)
slice(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mz_array,[0],[0],[0])
title('Mz_L')

subplot(3,2,2)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mx_array,[0],[0],[0])
title('Mx_R')
subplot(3,2,4)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.My_array,[0],[0],[0])
title('My_R')
subplot(3,2,6)
slice(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mz_array,[0],[0],[0])
title('Mz_R')
saveas(gca,[name,'_M.png'])



% Fx = interp3(str_mesh,rot_mesh,dev_mesh,Fx_array,0,0,0)
  