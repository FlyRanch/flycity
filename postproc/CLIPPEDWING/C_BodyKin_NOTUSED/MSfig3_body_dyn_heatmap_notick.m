%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

%% attitudes
% roll
biny_min = -45;
biny_max = 90;
var = roll_plot;

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
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
ylabel('roll [deg]','fontsize',10) 

% pitch
biny_min = -45;
biny_max = 90;
var = pitch_plot;

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
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
ylabel('pitch [deg]','fontsize',10) 

% yaw
biny_min = -45;
biny_max = 90;
var = yaw_plot;

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
set(gca,'XTick',[]) 
set(gca,'YTick',[])
xlabel('time [s]','fontsize',10) 
ylabel('yaw [deg]','fontsize',10) 

