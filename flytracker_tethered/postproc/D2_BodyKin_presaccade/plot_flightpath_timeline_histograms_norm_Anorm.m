%% plot flightpath timeline histograms ON

figure

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

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -180 180])
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

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -.5 .25])
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

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -180 180])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('A direction','fontsize',10) 


% roll
biny_min = -180;
biny_max = 180;
var = roll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',-180:180:180,'YTicklabel',180:-180:-180,'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll','fontsize',10) 

% slip
biny_min = -90;
biny_max = 90;
var = slip_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -90 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',-90:90:90,'YTicklabel',90:-90:-90,'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('slip','fontsize',10) 

% pitch
biny_min = -90;
biny_max = 90;
var = pitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -90 90])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-90;0;90],'YTicklabel',[90;0;-90],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('pitch','fontsize',10) 


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

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -20 0])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-20;0],'YTicklabel',[20;0],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('A','fontsize',10) 


% An
biny_min = -30;
biny_max = 30;
var = An_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -20 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-20;0;20],'YTicklabel',[20;0;-20],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('An','fontsize',10) 


% At
biny_min = -30;
biny_max = 30;
var = At_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -20 10])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-20;0;10],'YTicklabel',[20;0;-10],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('At','fontsize',10) 



