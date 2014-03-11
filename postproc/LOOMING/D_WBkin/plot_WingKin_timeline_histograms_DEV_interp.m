%% plot flightpath timeline histograms ON
figure

% bins
nx=1000;
ny=100;

[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% left ds max
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_max_udsPREV_L);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('left dev','fontsize',10) 
title('downstroke max')

% left ds min
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_min_ds_L);

% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left dev','fontsize',10) 
title('downstroke min')

% left up max
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_max_dus_L);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left dev','fontsize',10) 
title('upstroke max')

% left min
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_min_us_L);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left dev','fontsize',10) 
title('upstroke min')





% right ds max
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_max_udsPREV_R);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('right dev','fontsize',10) 
% title('downstroke max')

% right ds min
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_min_ds_R);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right dev','fontsize',10) 
% title('downstroke min')

% right up max
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_max_dus_R);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right dev','fontsize',10) 
% title('upstroke max')

% right min
biny_min = -45;
biny_max = 45;
var = rad2deg(dev_min_us_R);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,8)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right dev','fontsize',10) 
% title('upstroke min')








% L-R ds max
biny_min = -45;
biny_max = 45;
var = rad2deg(ddev_max_udsPREV);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,9)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('dev diff L-R','fontsize',10) 
% title('downstroke max')

% L-R ds min
biny_min = -45;
biny_max = 45;
var = rad2deg(ddev_min_ds);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,10)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R dev','fontsize',10) 
% title('downstroke min')

% L-R up max
biny_min = -45;
biny_max = 45;
var = rad2deg(ddev_max_dus);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,11)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R dev','fontsize',10) 
% title('upstroke max')

% L-R min
biny_min = -45;
biny_max = 45;
var = rad2deg(ddev_min_us);
% var = var(:);

var_interp = nan(size(t_hist_fr));
for i=1:size(var,2)
    x = t_hist_wb(:,i);
    y = var(:,i);
    xi = t_hist_fr(:,i);
    
    x = x(isnan(x)==0);
    y = y(isnan(y)==0);
    
    if isempty(x)==0
        yi = interp1(x,y,xi);
        var_interp(:,i) = yi;
    end
end
var = var_interp(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,4,12)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R dev','fontsize',10) 
% title('upstroke min')

