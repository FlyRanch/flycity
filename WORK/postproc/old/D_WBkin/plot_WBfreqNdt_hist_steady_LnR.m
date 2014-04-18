
figure
subplot(2,2,1)
% title('downstroke dt')
hold on
subplot(2,2,2)
% title('upstroke dt')
hold on
subplot(2,2,3)
% title('flap freq')
hold on
subplot(2,2,4)
hold on

mid=0;
hi=0;

dt_ds_steady = [];
dt_us_steady = [];
f_wb_steady = [];
Rds_steady = [];

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
            
            dt_ds_now = ([dt_ds_L(wb,seq),dt_ds_R(wb,seq)]);
            dt_us_now = ([dt_us_L(wb,seq),dt_us_R(wb,seq)]);
            f_wb_now = ([f_wb_L(wb,seq),f_wb_R(wb,seq)]);
            Rds_now = dt_ds_now ./ (dt_ds_now + dt_us_now);
            
            % store data
            dt_ds_steady(end+1:end+2,1) = dt_ds_now;
            dt_us_steady(end+1:end+2,1) = dt_us_now;
            f_wb_steady(end+1:end+2,1) = f_wb_now;
            Rds_steady(end+1:end+2,1) = Rds_now;

%             pause(.001)
        end
    end
end

subplot(2,2,1)
h=hist(dt_ds_steady,unique(dt_ds_steady));
hist(dt_ds_steady,unique(dt_ds_steady))
axis([.001 .004 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('downstroke dt','fontsize',10)
%     grid on

subplot(2,2,2)
h=hist(dt_us_steady,unique(dt_us_steady));
hist(dt_us_steady,unique(dt_us_steady))
axis([.001 .004 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.001:1) 
    set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('upstroke dt','fontsize',10)
%     grid on

subplot(2,2,3)
h=hist(Rds_steady,10);
hist(Rds_steady,10)
axis([0.4 .7 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.1:1) 
    set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('downstroke ratio','fontsize',10)
%     grid on
            
            
subplot(2,2,4)
h=hist(f_wb_steady,unique(f_wb_steady));
hist(f_wb_steady,unique(f_wb_steady))
axis([100 300 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:100:300) 
    set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('flap freq','fontsize',10)
%     grid on
            
            
mid
hi
