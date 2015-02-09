%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

%% attitudes
% roll
biny_min = -90;
biny_max = 90;
var = droll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll','fontsize',10) 

% pitch
biny_min = -90;
biny_max = 90;
var = dpitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,4)
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
var = dyaw_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('yaw','fontsize',10) 


%% angular vel
% roll_dot
biny_min = -3000;
biny_max = 3000;
var = roll_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll rate','fontsize',10) 

% pitch_dot
biny_min = -3000;
biny_max = 3000;
var = pitch_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
% set(gca,'XTick',-.05:-t_start:t_stop)  
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('pitch rate','fontsize',10) 


% yaw_dot
biny_min = -3000;
biny_max = 3000;
var = yaw_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('yaw rate','fontsize',10) 


%% anguklar accel
% roll_dot_dot
biny_min = -2e5;
biny_max = 2e5;
var = roll_dot_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll accel','fontsize',10) 

% pitch_dot_dot
biny_min = -2e5;
biny_max = 2e5;
var = pitch_dot_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
% set(gca,'XTick',-.05:-t_start:t_stop)  
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('pitch accel','fontsize',10) 


% yaw_dot_dot
biny_min = -2e5;
biny_max = 2e5;
var = yaw_dot_dot_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
end
var_hist(var_hist>1)=1;

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('yaw accel','fontsize',10) 

