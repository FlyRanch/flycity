% make cali file
clear
clc
close all

% load_name='cali_matrix_08132013.mat'
% load_name='cali_matrix_09292013.mat'

load_name=dir('cali_matrix*')
load_name = load_name.name;

load(load_name)
name = ['interp_',load_name]

%% cali interp
step_interp = 1;

str_steps_interp = [round(min(cali_R.str_steps)):step_interp:round(max(cali_R.str_steps))]';
rot_steps_interp = [round(min(cali_R.rot_steps)):step_interp:round(max(cali_R.rot_steps))]';
dev_steps_interp = [round(min(cali_R.dev_steps)):step_interp:round(max(cali_R.dev_steps))]';

str_steps_interp = [-100:step_interp:100]';
rot_steps_interp = [-80:step_interp:80]';
dev_steps_interp = [-35:step_interp:35]';

[rot_interp,str_interp,dev_interp] = meshgrid(rot_steps_interp,str_steps_interp,dev_steps_interp);
% [cali_L.rot_interp,cali_L.str_interp,cali_L.dev_interp] = meshgrid(rot_steps_interp,str_steps_interp,dev_steps_interp);
% [cali_R.rot_interp,cali_R.str_interp,cali_R.dev_interp] = meshgrid(rot_steps_interp,str_steps_interp,dev_steps_interp);

% cali_L.Fx_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fx_array,str_interp,rot_interp,dev_interp,'spline');
% cali_L.Fy_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fy_array,str_interp,rot_interp,dev_interp,'spline');
% cali_L.Fz_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fz_array,str_interp,rot_interp,dev_interp,'spline');
% cali_L.Mx_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mx_array,str_interp,rot_interp,dev_interp,'spline');
% cali_L.My_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.My_array,str_interp,rot_interp,dev_interp,'spline');
% cali_L.Mz_interp = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mz_array,str_interp,rot_interp,dev_interp,'spline');
% 
% cali_R.Fx_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fx_array,str_interp,rot_interp,dev_interp,'spline');
% cali_R.Fy_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fy_array,str_interp,rot_interp,dev_interp,'spline');
% cali_R.Fz_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fz_array,str_interp,rot_interp,dev_interp,'spline');
% cali_R.Mx_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mx_array,str_interp,rot_interp,dev_interp,'spline');
% cali_R.My_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.My_array,str_interp,rot_interp,dev_interp,'spline');
% cali_R.Mz_interp = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mz_array,str_interp,rot_interp,dev_interp,'spline');

cali_L.Fx_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fx_array,rot_interp,str_interp,dev_interp,'spline');
cali_L.Fy_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fy_array,rot_interp,str_interp,dev_interp,'spline');
cali_L.Fz_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fz_array,rot_interp,str_interp,dev_interp,'spline');
cali_L.Mx_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Mx_array,rot_interp,str_interp,dev_interp,'spline');
cali_L.My_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.My_array,rot_interp,str_interp,dev_interp,'spline');
cali_L.Mz_interp = interp3(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Mz_array,rot_interp,str_interp,dev_interp,'spline');

cali_R.Fx_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fx_array,rot_interp,str_interp,dev_interp,'spline');
cali_R.Fy_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fy_array,rot_interp,str_interp,dev_interp,'spline');
cali_R.Fz_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fz_array,rot_interp,str_interp,dev_interp,'spline');
cali_R.Mx_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Mx_array,rot_interp,str_interp,dev_interp,'spline');
cali_R.My_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.My_array,rot_interp,str_interp,dev_interp,'spline');
cali_R.Mz_interp = interp3(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Mz_array,rot_interp,str_interp,dev_interp,'spline');

%% save data
save(name)

%% plot cali
mkdir('cali_figs')
cd('cali_figs')

figure(1)
subplot(3,2,1)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fx_array,[0],[0],[0])
title('Fx_L')
subplot(3,2,3)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fy_array,[0],[0],[0])
title('Fy_L')
subplot(3,2,5)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Fz_array,[0],[0],[0])
title('Fz_L')

subplot(3,2,2)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fx_array,[0],[0],[0])
title('Fx_R')
subplot(3,2,4)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fy_array,[0],[0],[0])
title('Fy_R')
subplot(3,2,6)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Fz_array,[0],[0],[0])
title('Fz_R')
saveas(gca,[name,'_F.png'])

figure(2)
subplot(3,2,1)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Mx_array,[0],[0],[0])
title('Mx_L')
subplot(3,2,3)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.My_array,[0],[0],[0])
title('My_L')
subplot(3,2,5)
slice(cali_L.rot_mesh,cali_L.str_mesh,cali_L.dev_mesh,cali_L.Mz_array,[0],[0],[0])
title('Mz_L')

subplot(3,2,2)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Mx_array,[0],[0],[0])
title('Mx_R')
subplot(3,2,4)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.My_array,[0],[0],[0])
title('My_R')
subplot(3,2,6)
slice(cali_R.rot_mesh,cali_R.str_mesh,cali_R.dev_mesh,cali_R.Mz_array,[0],[0],[0])
title('Mz_R')
saveas(gca,[name,'_M.png'])

% plot cali interp
figure(3)
subplot(3,2,1)
h=slice(rot_interp,str_interp,dev_interp,cali_L.Fx_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_L.Fx_array(:)) max(cali_L.Fx_array(:))])
title('Fx_L')
subplot(3,2,3)
h=slice(rot_interp,str_interp,dev_interp,cali_L.Fy_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_L.Fy_array(:)) max(cali_L.Fy_array(:))])
title('Fy_L')
subplot(3,2,5)
h=slice(rot_interp,str_interp,dev_interp,cali_L.Fz_interp,[0],[0],[0]);
caxis([min(cali_L.Fz_array(:)) max(cali_L.Fz_array(:))])
set(h,'EdgeColor','none')
title('Fz_L')

subplot(3,2,2)
h=slice(rot_interp,str_interp,dev_interp,cali_R.Fx_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
title('Fx_R')
caxis([min(cali_R.Fx_array(:)) max(cali_R.Fx_array(:))])
subplot(3,2,4)
h=slice(rot_interp,str_interp,dev_interp,cali_R.Fy_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_R.Fy_array(:)) max(cali_R.Fy_array(:))])
title('Fy_R')
subplot(3,2,6)
h=slice(rot_interp,str_interp,dev_interp,cali_R.Fz_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_R.Fz_array(:)) max(cali_R.Fz_array(:))])
title('Fz_R')
saveas(gca,[name,'_Finterp.png'])

figure(4)
subplot(3,2,1)
h=slice(rot_interp,str_interp,dev_interp,cali_L.Mx_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_L.Mx_array(:)) max(cali_L.Mx_array(:))])
title('Mx_L')
subplot(3,2,3)
h=slice(rot_interp,str_interp,dev_interp,cali_L.My_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_L.My_array(:)) max(cali_L.My_array(:))])
title('My_L')
subplot(3,2,5)
h=slice(rot_interp,str_interp,dev_interp,cali_L.Mz_interp,[0],[0],[0]);
caxis([min(cali_L.Mz_array(:)) max(cali_L.Mz_array(:))])
set(h,'EdgeColor','none')
title('Mz_L')

subplot(3,2,2)
h=slice(rot_interp,str_interp,dev_interp,cali_R.Mx_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
title('Mx_R')
caxis([min(cali_R.Mx_array(:)) max(cali_R.Mx_array(:))])
subplot(3,2,4)
h=slice(rot_interp,str_interp,dev_interp,cali_R.My_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_R.My_array(:)) max(cali_R.My_array(:))])
title('My_R')
subplot(3,2,6)
h=slice(rot_interp,str_interp,dev_interp,cali_R.Mz_interp,[0],[0],[0]);
set(h,'EdgeColor','none')
caxis([min(cali_R.Mz_array(:)) max(cali_R.Mz_array(:))])
title('Mz_R')
saveas(gca,[name,'_Minterp.png'])

cd ..

  