%% plot flightpath timeline histograms ON
figure

% bins
nx=100;
ny=25;

[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% left  max
biny_min = 90;
biny_max = 225;
var = rad2deg(pitch_max_wb_L(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('left wingpitch','fontsize',10) 
title('wingpitch max')

% left  min
biny_min = -45;
biny_max = 90;
var = rad2deg(pitch_min_wb_L(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left wingpitch','fontsize',10) 
title('wingpitch min')

% left A
biny_min = 90;
biny_max = 225;
var = rad2deg(Apitch_wb_L(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('left wingpitch','fontsize',10) 
title('A wingpitch')







% right  max
biny_min = 90;
biny_max = 225;
var = rad2deg(pitch_max_wb_R(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('right wingpitch','fontsize',10) 
% title('wingpitch max')

% right  min
biny_min = -45;
biny_max = 90;
var = rad2deg(pitch_min_wb_R(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right wingpitch','fontsize',10) 
% title('wingpitch min')

% right A
biny_min = 90;
biny_max = 225;
var = rad2deg(Apitch_wb_R(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('right wingpitch','fontsize',10) 
% title('A wingpitch')





% L-R  max
biny_min = -90;
biny_max = 90;
var = rad2deg(dpitch_max_wb(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('pitch diff L-R','fontsize',10) 
% title('wingpitch max')

% L-R  min
biny_min = -90;
biny_max = 90;
var = rad2deg(dpitch_min_wb(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R','fontsize',10) 
% title('wingpitch min')

% L-R up max
biny_min = -90;
biny_max = 90;
var = rad2deg(dApitch_wb(:));

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
% for i = 1:nx
%     var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
% end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',t_start:-t_start:t_stop) 
%     set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('L-R','fontsize',10) 
% title('uppitch max')
