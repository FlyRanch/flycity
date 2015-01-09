
figure
subplot(3,3,1)
% title('downstroke dt')
hold on
subplot(3,3,2)
% title('upstroke dt')
hold on
subplot(3,3,3)
% title('flap freq')
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

mid=0;
hi=0;

dt_ds_L_steady = [];
dt_ds_R_steady = [];
Ddt_ds_steady = [];
dt_us_L_steady = [];
dt_us_R_steady = [];
Ddt_us_steady = [];
f_wb_L_steady = [];
f_wb_R_steady = [];
Df_wb_steady = [];

dt_ds_L_hi = [];
dt_ds_R_hi = [];
Ddt_ds_hi = [];
dt_us_L_hi = [];
dt_us_R_hi = [];
Ddt_us_hi = [];
f_wb_L_hi = [];
f_wb_R_hi = [];
Df_wb_hi = [];

dt_ds_L_low = [];
dt_ds_R_low = [];
Ddt_ds_low = [];
dt_us_L_low = [];
dt_us_R_low = [];
Ddt_us_low = [];
f_wb_L_low = [];
f_wb_R_low = [];
Df_wb_low = [];

for seq = 1:size(n_down_start_L,2)
    counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        if abs(roll_dot_dot_mean_wb(wb,seq)) < roll_limit &&...
                abs(pitch_dot_dot_mean_wb(wb,seq)) < pitch_limit &&...
                abs(yaw_dot_dot_mean_wb(wb,seq)) < yaw_limit &&...
                abs(F_mean_wb(wb,seq)-1) < F_limit
            
                mid=mid+1;

%         % only within steady&maneuver
%         if n_down_start_L_now > n_first(seq) && n_up_stop_L_now < n_post(seq) && ...
%                 n_down_start_R_now > n_first(seq) && n_up_stop_R_now < n_post(seq)
            
            dt_ds_L_now = dt_ds_L(wb,seq);
            dt_ds_R_now = dt_ds_R(wb,seq);
            dt_us_L_now = dt_us_L(wb,seq);
            dt_us_R_now = dt_us_R(wb,seq);
            f_wb_L_now = f_wb_L(wb,seq);
            f_wb_R_now = f_wb_R(wb,seq);
            
            Ddt_ds_now = dt_ds_L_now - dt_ds_R_now;
            Ddt_us_now = dt_us_L_now - dt_us_R_now;
            Df_wb_now = f_wb_L_now - f_wb_R_now;
            
            % store data
            dt_ds_L_steady(end+1,1) = dt_ds_L_now;
            dt_ds_R_steady(end+1,1) = dt_ds_R_now;
            Ddt_ds_steady(end+1,1) = Ddt_ds_now;
            dt_us_L_steady(end+1,1) = dt_us_L_now;
            dt_us_R_steady(end+1,1) = dt_us_R_now;
            Ddt_us_steady(end+1,1) = Ddt_us_now;
            f_wb_L_steady(end+1,1) = f_wb_L_now;
            f_wb_R_steady(end+1,1) = f_wb_R_now;
            Df_wb_steady(end+1,1) = Df_wb_now;

%             pause(.001)
        end
    end
end


subplot(3,3,1)
h=hist(dt_ds_L_steady,15);
hist(dt_ds_L_steady,15)
axis([.002 .004 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
%     grid on

subplot(3,3,4)
hist(dt_ds_R_steady,15)
axis([.002 .004 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Right wing','fontsize',10)
%     grid on

subplot(3,3,7)
hist(Ddt_ds_steady,15)
axis([-.001 .001 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
    xlabel('downstroke dt','fontsize',10)
    ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,3,2)
h=hist(dt_us_L_steady,15);
hist(dt_us_L_steady,15)
axis([.002 .004 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
%     grid on

subplot(3,3,5)
hist(dt_us_R_steady,15)
axis([.002 .004 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
%     grid on

subplot(3,3,8)
hist(Ddt_us_steady,15)
axis([-.001 .001 0 200])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:200,'fontsize',8)
    xlabel('downstroke dt','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,3,3)
hist(f_wb_L_steady,15)
axis([150 250 0 200])
    set(gca,'XTick',0:50:300) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
%     grid on

subplot(3,3,6)
hist(f_wb_R_steady,15)
axis([150 250 0 200])
    set(gca,'XTick',0:50:300) 
    set(gca,'YTick',0:100:200,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
%     grid on

subplot(3,3,9)
hist(Df_wb_steady,15)
axis([-10 10 0 200])
    set(gca,'XTick',-10:5:300) 
    set(gca,'YTick',0:100:200,'fontsize',8)
    xlabel('flap freq','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
            
mid
hi
