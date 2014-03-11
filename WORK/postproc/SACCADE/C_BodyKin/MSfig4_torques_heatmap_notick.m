%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% Mroll
biny_min = -.005;
biny_max = .005;
var = Mroll;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -.005 .005])
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
ylabel('Mroll [Nm]','fontsize',10) 

% Mpitch
biny_min = -.005;
biny_max = .005;
var = Mpitch;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -.005 .005])
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10)
ylabel('Mpitch [Nm]','fontsize',10) 

% Myaw
biny_min = -.005;
biny_max = .005;
var = Myaw;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -.005 .005])
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10)
ylabel('Myaw [Nm]','fontsize',10) 

% 
% % F
% biny_min = .5;
% biny_max = 2;
% var = F_plot;
% 
% var = var(:);
% biny = biny_min: (biny_max - biny_min)/ny :biny_max;
% var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(2:end-1,i)); % normalize per time bin
% end
% var_hist(var_hist>1)=1;
% 
% subplot(3,3,7)
% imagesc(binx,-biny,var_hist)
% axis([t_start t_stop -biny_max -biny_min])
% %     set(gca,'XTick',t_start:-t_start:t_stop) 
% set(gca,'XTick',[]) 
% set(gca,'YTick',[])
% % xlabel('time','fontsize',10) 
% ylabel('F/Mg [-]','fontsize',10) 
% 

