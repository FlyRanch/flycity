% plot wingbeat modifiactions for roll accel

t_hist = t_roll(:);

biny_min = -180;
biny_max = 180;

y_min = -45;
y_max = 45;

% stroke_L
var = Dstroke_rollaccel_L(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% pitch_L
var = Dpitch_rollaccel_L(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% dev_L
var = Ddev_rollaccel_L(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    line_now = var_hist(:,i);
    [meanx,meany] = ait_centroid(line_now)
    center_hist(i,1) = meany;
    
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

% plot(center_hist,'*')

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% stroke_R
var = Dstroke_rollaccel_R(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod right')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% pitch_R
var = Dpitch_rollaccel_R(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% dev_R
var = Ddev_rollaccel_R(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')

hold on
func = csaps(t_hist,-var)
fnplt(func)


% dstroke
var = Ddstroke_rollaccel(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod right')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% dpitch
var = Ddpitch_rollaccel(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

% ddev
var = Dddev_rollaccel(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')

hold on
func = csaps(t_hist,-var)
fnplt(func)

