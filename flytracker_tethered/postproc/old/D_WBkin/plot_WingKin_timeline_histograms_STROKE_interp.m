%% plot flightpath timeline histograms ON
figure

% bins
nx=1000;
ny=100;

[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% left  max
biny_min = 45;
biny_max = 135;
var = rad2deg(stroke_max_wb_L);

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

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('left wingstroke','fontsize',10) 
title('wingstroke max')

% left  min
biny_min = -90;
biny_max = 0;
var = rad2deg(stroke_min_wb_L);

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

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left wingstroke','fontsize',10) 
title('wingstroke min')

% left A
biny_min = 90;
biny_max = 180;
var = rad2deg(Astroke_wb_L);
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

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left wingstroke','fontsize',10) 
title('A wingstroke')







% right  max
biny_min = 45;
biny_max = 135;
var = rad2deg(stroke_max_wb_R);
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

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('right wingstroke','fontsize',10) 
% title('wingstroke max')

% right  min
biny_min = -90;
biny_max = 0;
var = rad2deg(stroke_min_wb_R);
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

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right wingstroke','fontsize',10) 
% title('wingstroke min')

% right A
biny_min = 90;
biny_max = 180;
var = rad2deg(Astroke_wb_R);
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

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right wingstroke','fontsize',10) 
% title('A wingstroke')





% L-R  max
biny_min = -30;
biny_max = 30;
var = rad2deg(dstroke_max_wb);
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

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('stroke diff L-R','fontsize',10) 
% title('wingstroke max')

% L-R  min
biny_min = -30;
biny_max = 30;
var = rad2deg(dstroke_min_wb);
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

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R','fontsize',10) 
% title('wingstroke min')

% L-R up max
biny_min = -30;
biny_max = 30;
var = rad2deg(dAstroke_wb);
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

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R','fontsize',10) 
% title('upstroke max')
