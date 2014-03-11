%% plot flightpath timeline histograms ON
figure

% bins
nx=1000;
ny=100;

[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% downstroke dt
biny_min = .001;
biny_max = .004;
var = dt_ds_L;

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
ylabel('left wingpitch','fontsize',10) 
title('downstroke dt')

% upstroke dt
biny_min = .001;
biny_max = .004;
var = dt_us_L;

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
% ylabel('left wingpitch','fontsize',10) 
title('upstroke dt')

% wb freq
biny_min = 150;
biny_max = 250;
var = f_wb_L;

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
% ylabel('left wingpitch','fontsize',10) 
title('wingbeat frequency')







% downstroke dt
biny_min = .001;
biny_max = .004;
var = dt_ds_R;

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
ylabel('right wingpitch','fontsize',10) 
% title('wingpitch max')

% upstroke dt
biny_min = .001;
biny_max = .004;
var = dt_us_R;

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
% ylabel('right wingpitch','fontsize',10) 
% title('wingpitch min')

% wb freq
biny_min = 150;
biny_max = 250;
var = f_wb_R;

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
% ylabel('right wingpitch','fontsize',10) 
% title('A wingpitch')





% L-R  ds dt
biny_min = -.0015;
biny_max = .0015;
var = ddt_ds;

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
ylabel('pitch diff L-R','fontsize',10) 
% title('wingpitch max')

% L-R  us dt
biny_min = -.0015;
biny_max = .0015;
var = ddt_us;

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
% title('wingpitch min')

% L-R wb freq
biny_min = -50;
biny_max = 50;
var = df_wb;

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
% title('uppitch max')
