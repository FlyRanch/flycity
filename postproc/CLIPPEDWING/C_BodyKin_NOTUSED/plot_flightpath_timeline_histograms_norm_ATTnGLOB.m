%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

%% strokeplane angles
% roll
biny_min = -90;
biny_max = 90;
var = roll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll','fontsize',10) 
title('strokeplane angles','fontsize',10) 

% pitch
biny_min = -90;
biny_max = 90;
var = pitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
% set(gca,'XTick',-.05:-t_start:t_stop)  
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('pitch','fontsize',10) 

% yaw
biny_min = -90;
biny_max = 90;
var = yaw_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('yaw','fontsize',10) 

%% euler angles
% roll_global
biny_min = -90;
biny_max = 90;
var = roll_global_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',[],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('euler roll','fontsize',10) 
title('body euler angles','fontsize',10) 

% pitch_global
biny_min = -90;
biny_max = 90;
var = pitch_global_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
% set(gca,'XTick',-.05:-t_start:t_stop)  
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',[],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('euler pitch','fontsize',10) 

% yaw_global
biny_min = -90;
biny_max = 90;
var = slip_global_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,2,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',[],'fontsize',8) 
xlabel('time','fontsize',10) 
% ylabel('euler slip','fontsize',10) 

