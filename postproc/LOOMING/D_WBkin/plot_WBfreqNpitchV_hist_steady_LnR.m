
figure
subplot(2,3,1)
% title('downstroke dt')
hold on
subplot(2,3,2)
title('steady wingbeat','fontsize',12)
hold on
subplot(2,3,3)
% title('flap freq')
hold on
subplot(2,3,4)
hold on
subplot(2,3,5)
% title('flap freq')
hold on
subplot(2,3,6)
hold on

subplot(2,3,1)
h=hist(dt_ds_steady,unique(dt_ds_steady));
hist(dt_ds_steady,unique(dt_ds_steady))
axis([.001 .004 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.001:1) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('downstroke dt','fontsize',10)
%     grid on

subplot(2,3,2)
h=hist(dt_us_steady,unique(dt_us_steady));
hist(dt_us_steady,unique(dt_us_steady))
axis([.001 .004 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.001:1) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('upstroke dt','fontsize',10)
%     grid on

subplot(2,3,4)
h=hist(Rds_steady,10);
hist(Rds_steady,10)
axis([0.4 .7 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.1:1) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('downstroke ratio','fontsize',10)
%     grid on
            
            
subplot(2,3,5)
h=hist(f_wb_steady,unique(f_wb_steady));
hist(f_wb_steady,unique(f_wb_steady))
axis([100 300 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:100:300) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('flap freq','fontsize',10)
%     grid on
            
subplot(2,3,3)
h=hist(V_steady,10);
hist(V_steady,10)
axis([-.1 .5 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:.25:1) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('flight speed','fontsize',10)
%     grid on
            
subplot(2,3,6)
h=hist(pitch_global_steady,10);
hist(pitch_global_steady,10)
axis([0 90 0 100*ceil(max(h)/100)])
    set(gca,'XTick',0:45:90) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('body pitch','fontsize',10)
%     grid on
            
