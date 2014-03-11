%% plot flightpath timeline histograms ON

figure

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% heading
biny_min = -180;
biny_max = 180;
var = stim_angle_spn_plot;

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
ylabel('strokeplane normal','fontsize',10) 

% Fsp_pitch
biny_min = -30;
biny_max = 30;
var = Fsp_pitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -30 30])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',-180:30:180,'YTicklabel',180:-30:-180,'fontsize',8) 
%     xlabel('time','fontsize',10)
    ylabel('Fsp pitch','fontsize',10)


% Fsp_roll
biny_min = -30;
biny_max = 30;
var = Fsp_roll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -30 30])
    set(gca,'XTick',-0.5:.025:.5) 
set(gca,'YTick',-180:30:180,'YTicklabel',180:-30:-180,'fontsize',8) 
xlabel('time','fontsize',10) 
    ylabel('Fsp roll','fontsize',10)


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


% F
biny_min = 0;
biny_max = 2;
var = F_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -2 0])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-2;0],'YTicklabel',[2;0],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('F','fontsize',10) 


% Fn
biny_min = -2;
biny_max = 2;
var = Fn_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -2 2])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-2;0;2],'YTicklabel',[2;0;-2],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Fn','fontsize',10) 


% Ft
biny_min = -2;
biny_max = 2;
var = Ft_hor_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([-.025 .05 -2 1])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-2;0;1],'YTicklabel',[2;0;-1],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('Ft','fontsize',10) 



