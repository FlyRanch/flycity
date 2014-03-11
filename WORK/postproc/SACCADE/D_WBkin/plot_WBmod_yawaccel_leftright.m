% plot wingbeat modifiactions for yaw accel

t_hist = t_yaw(:);

biny_min = -180;
biny_max = 180;

y_min = -45;
y_max = 45;

% stroke_right
var = -Dstroke_yawaccel_right(:);

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
title('yaw accel mod left')

% pitch_right
var = -Dpitch_yawaccel_right(:);

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
% title('yaw accel mod left')

% dev_right
var = -Ddev_yawaccel_right(:);

biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('yaw accel left')

% stroke_left
var = -Dstroke_yawaccel_left(:);

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
title('yaw accel mod right')

% pitch_left
var = -Dpitch_yawaccel_left(:);

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
% title('yaw accel mod left')

% dev_left
var = -Ddev_yawaccel_left(:);

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
% title('yaw accel left')



% dstroke
var = -Ddstroke_yawaccel(:);

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
title('yaw accel mod right')

% dpitch
var = -Ddpitch_yawaccel(:);

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
% title('yaw accel mod left')

% ddev
var = -Dddev_yawaccel(:);

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
% title('yaw accel left')

