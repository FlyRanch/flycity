%% plot flightpath timeline histograms ON

figure

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% heading
biny_min = -180;
biny_max = 180;
var = stim_angle_spn_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('strokeplane normal','fontsize',10) 

% Fsp_pitch
biny_min = -30;
biny_max = 30;
var = Fsp_pitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,3)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
%     xlabel('time','fontsize',10)
    ylabel('Fsp pitch','fontsize',10)


% Fsp_roll
biny_min = -30;
biny_max = 30;
var = Fsp_roll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,5)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
    set(gca,'XTick',-0.5:.025:.5) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
    ylabel('Fsp roll','fontsize',10)


% roll
biny_min = -90;
biny_max = 90;
var = roll_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,2)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('roll','fontsize',10) 

% yaw
biny_min = -90;
biny_max = 90;
var = yaw_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('yaw','fontsize',10) 

% pitch
biny_min = -45;
biny_max = 135;
var = pitch_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,2,6)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -biny_max -biny_min])
set(gca,'XTick',-.05:.025:.5)  
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
xlabel('time','fontsize',10) 
ylabel('pitch','fontsize',10) 


