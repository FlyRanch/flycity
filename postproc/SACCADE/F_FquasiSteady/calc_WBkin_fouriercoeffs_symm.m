function [stroke_L,stroke_R,stroke_dot_L,stroke_dot_R,stroke_dot_dot_L,stroke_dot_dot_R,...
    dev_L,dev_R,dev_dot_L,dev_dot_R,dev_dot_dot_L,dev_dot_dot_R,...
    rot_L,rot_R,rot_dot_L,rot_dot_R,rot_dot_dot_L,rot_dot_dot_R,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs_symm(nr_timepoints,freq,stroke_coeffs,dev_coeffs,rot_coeffs)

%  calc time vars
dt = (1/freq)/(nr_timepoints-1);
t = 0:dt:(dt*(nr_timepoints-1)); 
t_norm = [0:1/(nr_timepoints-1):1];
dt_norm = t_norm(2)-t_norm(1);

if isempty(stroke_coeffs) == 0
    [stroke] = calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs,0);
    [stroke_dot] = calc_val_dot_fourier_series_4thN8th_order(t_norm,stroke_coeffs,0);
    [stroke_dot_dot] = calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,stroke_coeffs,0);
else
    [stroke] = zeros(1,nr_timepoints);
    [stroke_dot] = zeros(1,nr_timepoints);
    [stroke_dot_dot] = zeros(1,nr_timepoints);
end
    
if isempty(dev_coeffs) == 0
    [dev] = calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs,0);
    [dev_dot] = calc_val_dot_fourier_series_4thN8th_order(t_norm,dev_coeffs,0);
    [dev_dot_dot] = calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,dev_coeffs,0);
else
    [dev] = zeros(1,nr_timepoints);
    [dev_dot] = zeros(1,nr_timepoints);
    [dev_dot_dot] = zeros(1,nr_timepoints);
end
    
if isempty(rot_coeffs) == 0
    [rot] = calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs,0);
    [rot_dot] = calc_val_dot_fourier_series_4thN8th_order(t_norm,rot_coeffs,0);
    [rot_dot_dot] = calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,rot_coeffs,0);
else
    [rot] = zeros(1,nr_timepoints);
    [rot_dot] = zeros(1,nr_timepoints);
    [rot_dot_dot] = zeros(1,nr_timepoints);
end

%% derivatives
stroke_dot     = stroke_dot./dt*dt_norm;
dev_dot     = dev_dot./dt*dt_norm;
rot_dot     = rot_dot./dt*dt_norm;

stroke_dot_dot     = stroke_dot_dot./dt./dt*dt_norm*dt_norm;
dev_dot_dot     = dev_dot_dot./dt./dt*dt_norm*dt_norm;
rot_dot_dot     = rot_dot_dot./dt./dt*dt_norm*dt_norm;

%% deg2rad
stroke = deg2rad(stroke);
stroke_dot = deg2rad(stroke_dot);
stroke_dot_dot = deg2rad(stroke_dot_dot);

dev = deg2rad(dev);
dev_dot = deg2rad(dev_dot);
dev_dot_dot = deg2rad(dev_dot_dot);

rot = deg2rad(rot);
rot_dot = deg2rad(rot_dot);
rot_dot_dot = deg2rad(rot_dot_dot);

%% L&R
WBkin2kinL
WBkin2kinR

%% plot kinematics
% figure
% subplot(3,3,1)
% plot(t_norm,stroke_L)
% hold on
% plot(t_norm,stroke_R,'c')
% % plot(t_norm,stroke_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('stroke','fontsize',10)
% 
% subplot(3,3,2)
% plot(t_norm,stroke_dot_L)
% hold on
% plot(t_norm,stroke_dot_R,'c')
% % plot(t_norm,stroke_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('stroke dot','fontsize',10)
% 
% subplot(3,3,3)
% plot(t_norm,stroke_dot_dot_L)
% hold on
% plot(t_norm,stroke_dot_dot_R,'c')
% % plot(t_norm,stroke_dot_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('stroke ddot','fontsize',10)
% 
% subplot(3,3,4)
% plot(t_norm,dev_L)
% hold on
% plot(t_norm,dev_R,'c')
% % plot(t_norm,dev_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('dev','fontsize',10)
% 
% subplot(3,3,5)
% plot(t_norm,dev_dot_L)
% hold on
% plot(t_norm,dev_dot_R,'c')
% % plot(t_norm,dev_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('dev dot','fontsize',10)
% 
% subplot(3,3,6)
% plot(t_norm,dev_dot_dot_L)
% hold on
% plot(t_norm,dev_dot_dot_R,'c')
% % plot(t_norm,dev_dot_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('dev ddot','fontsize',10)
% 
% subplot(3,3,7)
% plot(t_norm,rot_L)
% hold on
% plot(t_norm,rot_R,'c')
% % plot(t_norm,rot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('rot','fontsize',10)
% 
% subplot(3,3,8)
% plot(t_norm,rot_dot_L)
% hold on
% plot(t_norm,rot_dot_R,'c')
% % plot(t_norm,rot_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('rot dot','fontsize',10)
% 
% subplot(3,3,9)
% plot(t_norm,rot_dot_dot_L)
% hold on
% plot(t_norm,rot_dot_dot_R,'c')
% % plot(t_norm,rot_dot_dot_dps,'--r')
% xlabel('t*','fontsize',10) 
% ylabel('rot ddot','fontsize',10)
% 
% saveas(gca,'WBkin_quasisteadyModel_steadyWB.png')
