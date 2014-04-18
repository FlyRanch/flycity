%% plot flightpath timeline histograms ON

figure
colormap(cmap_bw)

% time bins
[N_histx,x_hist] = hist(t_hist,nx);
binx = x_hist;

% heading
biny_min = -180;
biny_max = 180;
var = stim_angle_vel_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
ylabel('heading [deg]','fontsize',10) 

% speed DIFFERENCE
biny_min = -.25;
biny_max = .5;
% var = V_plot;
var = dV_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -.5 .25])
%     set(gca,'XTick',-0.5:.025:.5) 
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10)
ylabel('dV [m/s]','fontsize',10) 


% Fhor dir
biny_min = -180;
biny_max = 180;
var = stim_angle_accel_plot;

var = var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([t_start t_stop -180 180])
set(gca,'XTick',[]) 
set(gca,'YTick',[])
xlabel('time [s]','fontsize',10) 
ylabel('dir F_h_o_r [deg]','fontsize',10) 

% F
biny_min = .5;
biny_max = 2;
var = F_plot;

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
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
ylabel('F/Mg [-]','fontsize',10) 


% Fsp_roll
biny_min = -30;
biny_max = 30;
var = Fsp_roll_plot;

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
%     set(gca,'XTick',t_start:-t_start:t_stop) 
set(gca,'XTick',[]) 
set(gca,'YTick',[])
% xlabel('time','fontsize',10) 
    ylabel('F/Mg_r_o_l_l [deg]','fontsize',10)


% Fsp_pitch
biny_min = -30;
biny_max = 30;
var = Fsp_pitch_plot;

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
set(gca,'XTick',[]) 
set(gca,'YTick',[])
    xlabel('time [s]','fontsize',10)
    ylabel('F/Mg_p_i_t_c_h [deg]','fontsize',10)



