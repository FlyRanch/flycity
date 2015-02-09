%% plot flightpath timeline histograms ON

% t_hist = [];
% for i=1:size(stim_angle_vel,2)    
%     t_hist = [t_hist;t-t_shift(i)];
% end

% heading
var = stim_angle_vel;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,1)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -180 180])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12)
xlabel('time','fontsize',18) 
ylabel('heading','fontsize',18) 

% yaw
var = stim_angle_yaw;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,2)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -180 180])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('yaw','fontsize',18) 


% slip
var = slip;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,3)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -90 90])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',-180:180:180,'YTicklabel',90:-90:-90,'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('slip','fontsize',18) 


% speed
var = V;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,4)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -.8 0])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',[-.8;-.4;0],'YTicklabel',[.8;.4;0],'fontsize',12) 
xlabel('time','fontsize',18)
ylabel('V','fontsize',18) 


% An
var = abs(An_hor);

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,5)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -15 0])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',[-15:5:0],'YTicklabel',[15:-5:0],'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('An','fontsize',18) 


% At
var = At_hor;

var = var(:);
[N_histx,x_hist] = hist(t_hist,nx);
[N_histy, y_hist] = hist(var,ny);
var_hist = hist3([var,t_hist], [ny,nx]);
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i));
end

subplot(6,1,6)
imagesc(x_hist,-y_hist,var_hist)
axis([-.3 .15 -15 5])
set(gca,'XTick',-.45:.15:.15) 
set(gca,'YTick',[-15;0;5],'YTicklabel',[15;0;-5],'fontsize',12) 
xlabel('time','fontsize',18) 
ylabel('At','fontsize',18) 



