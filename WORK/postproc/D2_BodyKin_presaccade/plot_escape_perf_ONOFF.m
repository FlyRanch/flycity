heading_post = responseDB.stim_angle_vel_post_turn;
yaw_post = responseDB.stim_angle_yaw_post_turn;
t_resp = responseDB.t_resp;

V_post = (max([responseDB.V_post_resp responseDB.V_post_accel]'))';
An_max = abs(responseDB.An_hor_max);
At_max = responseDB.At_hor_max;

% ON/OFF
k=0;
l=0;

clear H1 Y1 T1 V1 An1 At1 
clear H2 Y2 T2 V2 An2 At2 

figure
subplot(3,2,1)
hold on
subplot(3,2,2)
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

for i=1:length(heading_post)
    if  settings.expansion.OFF(i) == 0
        subplot(3,2,1)
        plot(1,heading_post(i),'ko','MarkerSize',5)
        subplot(3,2,3)
        plot(1,yaw_post(i),'ko','MarkerSize',5)
        subplot(3,2,5)
        plot(1,t_resp(i),'ko','MarkerSize',5)
        
        subplot(3,2,2)
        plot(1,V_post(i),'ko','MarkerSize',5)
        subplot(3,2,4)
        plot(1,An_max(i),'ko','MarkerSize',5)
        subplot(3,2,6)
        plot(1,At_max(i),'ko','MarkerSize',5)
        
        k = k+1;
        H1(k) = heading_post(i);
        Y1(k) = yaw_post(i);
        T1(k) = t_resp(i);
        
        V1(k) = V_post(i);
        An1(k) = An_max(i);
        At1(k) = At_max(i);
        
    elseif  settings.expansion.OFF(i) == 1
        subplot(3,2,1)
        plot(2,heading_post(i),'ko','MarkerSize',5)
        subplot(3,2,3)
        plot(2,yaw_post(i),'ko','MarkerSize',5)
        subplot(3,2,5)
        plot(2,t_resp(i),'ko','MarkerSize',5)
        
        subplot(3,2,2)
        plot(2,V_post(i),'ko','MarkerSize',5)
        subplot(3,2,4)
        plot(2,An_max(i),'ko','MarkerSize',5)
        subplot(3,2,6)
        plot(2,At_max(i),'ko','MarkerSize',5)
        
        l = l+1;
        H2(l) = heading_post(i);
        Y2(l) = yaw_post(i);
        T2(l) = t_resp(i);
        
        V2(l) = V_post(i);
        An2(l) = An_max(i);
        At2(l) = At_max(i);
    end
end

        subplot(3,2,1)
errorbar(1,nanmean(H1),nanstd(H1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(H2),nanstd(H2),'dg','linewidth',2,'MarkerSize',10)
ylabel('escape heading','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[-180:180:180],'fontsize',12)
axis([0 3 -180 180])
% grid on

        subplot(3,2,3)
errorbar(1,nanmean(Y1),nanstd(Y1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(Y2),nanstd(Y2),'dg','linewidth',2,'MarkerSize',10)
ylabel('escape yaw','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[-180:180:180],'fontsize',12)
axis([0 3 -180 180])
% grid on

        subplot(3,2,5)
errorbar(1,nanmean(T1),nanstd(T1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(T2),nanstd(T2),'dg','linewidth',2,'MarkerSize',10)
ylabel('response time','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[0:.05:.1],'fontsize',12)
axis([0 3 0 .1])

        subplot(3,2,2)
errorbar(1,nanmean(V1),nanstd(V1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(V2),nanstd(V2),'dg','linewidth',2,'MarkerSize',10)
ylabel('escape speed','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[0:.4:.8],'fontsize',12)
axis([0 3 0 .8])

        subplot(3,2,4)
errorbar(1,nanmean(An1),nanstd(An1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(An2),nanstd(An2),'dg','linewidth',2,'MarkerSize',10)
ylabel('An max','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[0:10:20],'fontsize',12)
axis([0 3 0 20])

        subplot(3,2,6)
errorbar(1,nanmean(At1),nanstd(At1),'dr','linewidth',2,'MarkerSize',10)
errorbar(2,nanmean(At2),nanstd(At2),'dg','linewidth',2,'MarkerSize',10)
ylabel('At max','fontsize',18) 
set(gca,'XTick',[1:2])
set(gca,'XTickLabel',{'ON','OFF'})
set(gca,'YTick',[0:10:20],'fontsize',12)
axis([0 3 0 20])