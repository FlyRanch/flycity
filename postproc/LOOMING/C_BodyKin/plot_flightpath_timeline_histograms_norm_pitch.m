%% plot flightpath timeline histograms ON

figure

% heading
var = stim_angle_vel_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,1)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -180 180])
set(gca,'XTick',-.05:.025:.5)   
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12)
xlabel('time','fontsize',18) 
ylabel('heading','fontsize',18) 

% yaw
var = stim_angle_yaw_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,4)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -180 180])
set(gca,'XTick',-.05:.025:.5)   
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12)
xlabel('time','fontsize',18) 
ylabel('yaw','fontsize',18) 


% Adir
var = stim_angle_accel_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,4)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -180 180])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('A direction','fontsize',18) 


% roll
var = roll_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,2)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -180 180])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('roll','fontsize',18) 

% slip
var = slip_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,5)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -90 90])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-90:90:90,'YTicklabel',90:-90:-90,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('slip','fontsize',18) 

% pitch
var = pitch_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,8)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -90 90])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-90:90:90,'YTicklabel',90:-90:-90,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('pitch','fontsize',18) 


% speed DIFFERENCE
% var = V_plot;
var = dV_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,3)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -.5 .25])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-.5:.25:.25],'YTicklabel',[.5:-.25:-.25],'fontsize',12) 
xlabel('time','fontsize',18)
ylabel('dV','fontsize',18) 


% An
var = An_hor_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,6)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -20 5])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-20;0;5],'YTicklabel',[20;0;-5],'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('An','fontsize',18) 


% At
var = At_hor_plot;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(3,3,9)
imagesc(x_hist,-y_hist,var_hist)
axis([-.025 .05 -20 15])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-20;0;15],'YTicklabel',[20;0;-15],'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('At','fontsize',18) 



