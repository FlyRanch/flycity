%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% heading
biny_min = -180;
biny_max = 180;
var = stim_angle_vel_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',8)
% xlabel('time','fontsize',10) 
ylabel('heading','fontsize',10) 

% speed DIFFERENCE
biny_min = -.25;
biny_max = .5;
% var = V_plot;
var = dV_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -.5 .25])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-.5:.25:.25],'YTicklabel',[.5:-.25:-.25],'fontsize',8) 
% xlabel('time','fontsize',10)
ylabel('dV','fontsize',10) 


% Adir
biny_min = -180;
biny_max = 180;
var = stim_angle_accel_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -180 180])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('A direction','fontsize',10) 



% A
biny_min = 0;
biny_max = 30;
var = A_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -20 0])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-20;0],'YTicklabel',[20;0],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('A','fontsize',10) 


% An
biny_min = -15;
biny_max = 15;
var = An_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -20 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-20;0;20],'YTicklabel',[20;0;-20],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('An','fontsize',10) 


% At
biny_min = -15;
biny_max = 15;
var = At_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -20 20])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-20;0;20],'YTicklabel',[20;0;-20],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('At','fontsize',10) 



