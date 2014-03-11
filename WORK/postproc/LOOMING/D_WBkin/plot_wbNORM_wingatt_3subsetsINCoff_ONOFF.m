
figure
subplot(3,5,1)
title('wing stroke')
hold on
subplot(3,5,2)
title('wing pitch')
hold on
subplot(3,5,3)
title('stroke deviation')
hold on
subplot(3,5,4)
title('angle of attack')
hold on
subplot(3,5,5)
title('wing speed')
hold on
subplot(3,5,6)
hold on
subplot(3,5,7)
hold on
subplot(3,5,8)
hold on
subplot(3,5,9)
hold on
subplot(3,5,10)
hold on
subplot(3,5,11)
hold on
subplot(3,5,12)
hold on
subplot(3,5,13)
hold on
subplot(3,5,14)
hold on
subplot(3,5,15)
hold on

mid=0;
hi=0;
low=0;

t_steady = [];
stroke_wb_L_steady = [];
stroke_wb_R_steady = [];
pitch_wb_L_steady = [];
pitch_wb_R_steady = [];
dev_wb_L_steady = [];
dev_wb_R_steady = [];
aoa_wb_L_steady = [];
aoa_wb_R_steady = [];
U_wb_L_steady = [];
U_wb_R_steady = [];
Dstroke_wb_steady = [];
Dpitch_wb_steady = [];
Ddev_wb_steady = [];
Daoa_wb_steady = [];
DU_wb_steady = [];

t_hi = [];
stroke_wb_L_hi = [];
stroke_wb_R_hi = [];
pitch_wb_L_hi = [];
pitch_wb_R_hi = [];
dev_wb_L_hi = [];
dev_wb_R_hi = [];
aoa_wb_L_hi = [];
aoa_wb_R_hi = [];
U_wb_L_hi = [];
U_wb_R_hi = [];
Dstroke_wb_hi = [];
Dpitch_wb_hi = [];
Ddev_wb_hi = [];
Daoa_wb_hi = [];
DU_wb_hi = [];

t_low = [];
stroke_wb_L_low = [];
stroke_wb_R_low = [];
pitch_wb_L_low = [];
pitch_wb_R_low = [];
dev_wb_L_low = [];
dev_wb_R_low = [];
aoa_wb_L_low = [];
aoa_wb_R_low = [];
U_wb_L_low = [];
U_wb_R_low = [];
Dstroke_wb_low = [];
Dpitch_wb_low = [];
Ddev_wb_low = [];
Daoa_wb_low = [];
DU_wb_low = [];


for seq = 1:size(n_down_start_L,2)
    counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        subset_now = subset(wb,seq);
        subset_OFF_now = subset_OFF(wb,seq);
        
        color_code_now = 0;
        if subset_OFF_ON == 0 || abs(subset_OFF_now) < subset_OFF_mid
            if subset_mid_ON == 1 && abs(subset_now) < subset_mid
                color_code_now = color_mid + 0.25;
                mid=mid+1;
            elseif subset_pos_ON == 1 && subset_now > subset_pos_min && subset_now < subset_pos_max
                color_code_now = color_max + 0.5;
                hi=hi+1;
            elseif subset_neg_ON == 1 && subset_now > subset_neg_min && subset_now < subset_neg_max
                color_code_now = color_min + 0.5;
                low=low+1;
            end
        end
        
        if color_code_now ~= 0

            n_down_start_L_now = n_down_start_L(wb,seq);
            n_up_stop_L_now = n_down_start_L(wb+1,seq);
            n_down_start_R_now = n_down_start_R(wb,seq);
            n_up_stop_R_now = n_down_start_R(wb+1,seq);

%         % only within steady&maneuver
%         if n_down_start_L_now > n_first(seq) && n_up_stop_L_now < n_post(seq) && ...
%                 n_down_start_R_now > n_first(seq) && n_up_stop_R_now < n_post(seq)
            
            stroke_wb_L_now = stroke_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            stroke_wb_R_now = stroke_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            pitch_wb_L_now = unwrap(pitch_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq));
            pitch_wb_R_now = unwrap(pitch_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq));
            dev_wb_L_now = dev_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            dev_wb_R_now = dev_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            U_wb_L_now = U_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            U_wb_R_now = U_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            aoa_wb_L_now = aoa_L_mirror(n_down_start_L_now:n_up_stop_L_now,seq);
            aoa_wb_R_now = aoa_R_mirror(n_down_start_R_now:n_up_stop_R_now,seq);
            
            xi = [0:1/100:1]';
            yi = [1:(length(stroke_wb_L_now)-1)/100:length(stroke_wb_L_now)]';
%             xi_L = [1:length(stroke_wb_L_now)/100:length(stroke_wb_L_now)];
%             xi_R = [1:length(stroke_wb_R_now)/100:length(stroke_wb_R_now)];

            stroke_wb_L_interp = rad2deg(interp1(stroke_wb_L_now,yi));
            stroke_wb_R_interp = rad2deg(interp1(stroke_wb_R_now,yi));
            pitch_wb_L_interp = rad2deg(interp1(pitch_wb_L_now,yi));
            pitch_wb_R_interp = rad2deg(interp1(pitch_wb_R_now,yi));
            dev_wb_L_interp = rad2deg(interp1(dev_wb_L_now,yi));
            dev_wb_R_interp = rad2deg(interp1(dev_wb_R_now,yi));
            aoa_wb_L_interp = rad2deg(interp1(aoa_wb_L_now,yi));
            aoa_wb_R_interp = rad2deg(interp1(aoa_wb_R_now,yi));
            U_wb_L_interp = interp1(U_wb_L_now,yi)/1000;
            U_wb_R_interp = interp1(U_wb_R_now,yi)/1000;
            
            Dstroke_wb_interp = stroke_wb_L_interp - stroke_wb_R_interp;
            Dpitch_wb_interp = pitch_wb_L_interp - pitch_wb_R_interp;
            Ddev_wb_interp = dev_wb_L_interp - dev_wb_R_interp;
            Daoa_wb_interp = aoa_wb_L_interp - aoa_wb_R_interp;
            DU_wb_interp = U_wb_L_interp - U_wb_R_interp;

            % store data
            if abs(subset_now) < subset_mid
                t_steady(end+1:end+length(xi),1) = xi;
                stroke_wb_L_steady(end+1:end+length(xi),1) = stroke_wb_L_interp;
                stroke_wb_R_steady(end+1:end+length(xi),1) = stroke_wb_R_interp;
                pitch_wb_L_steady(end+1:end+length(xi),1) = pitch_wb_L_interp;
                pitch_wb_R_steady(end+1:end+length(xi),1) = pitch_wb_R_interp;
                dev_wb_L_steady(end+1:end+length(xi),1) = dev_wb_L_interp;
                dev_wb_R_steady(end+1:end+length(xi),1) = dev_wb_R_interp;
                aoa_wb_L_steady(end+1:end+length(xi),1) = aoa_wb_L_interp;
                aoa_wb_R_steady(end+1:end+length(xi),1) = aoa_wb_R_interp;
                U_wb_L_steady(end+1:end+length(xi),1) = U_wb_L_interp;
                U_wb_R_steady(end+1:end+length(xi),1) = U_wb_R_interp;
                
                Dstroke_wb_steady(end+1:end+length(xi),1) = Dstroke_wb_interp;
                Dpitch_wb_steady(end+1:end+length(xi),1) = Dpitch_wb_interp;
                Ddev_wb_steady(end+1:end+length(xi),1) = Ddev_wb_interp;
                Daoa_wb_steady(end+1:end+length(xi),1) = Daoa_wb_interp;
                DU_wb_steady(end+1:end+length(xi),1) = DU_wb_interp;
                
            elseif subset_now > subset_pos_min
                t_hi(end+1:end+length(xi),1) = xi;
                stroke_wb_L_hi(end+1:end+length(xi),1) = stroke_wb_L_interp;
                stroke_wb_R_hi(end+1:end+length(xi),1) = stroke_wb_R_interp;
                pitch_wb_L_hi(end+1:end+length(xi),1) = pitch_wb_L_interp;
                pitch_wb_R_hi(end+1:end+length(xi),1) = pitch_wb_R_interp;
                dev_wb_L_hi(end+1:end+length(xi),1) = dev_wb_L_interp;
                dev_wb_R_hi(end+1:end+length(xi),1) = dev_wb_R_interp;
                aoa_wb_L_hi(end+1:end+length(xi),1) = aoa_wb_L_interp;
                aoa_wb_R_hi(end+1:end+length(xi),1) = aoa_wb_R_interp;
                U_wb_L_hi(end+1:end+length(xi),1) = U_wb_L_interp;
                U_wb_R_hi(end+1:end+length(xi),1) = U_wb_R_interp;
                
                Dstroke_wb_hi(end+1:end+length(xi),1) = Dstroke_wb_interp;
                Dpitch_wb_hi(end+1:end+length(xi),1) = Dpitch_wb_interp;
                Ddev_wb_hi(end+1:end+length(xi),1) = Ddev_wb_interp;
                Daoa_wb_hi(end+1:end+length(xi),1) = Daoa_wb_interp;
                DU_wb_hi(end+1:end+length(xi),1) = DU_wb_interp;
                
            elseif subset_now < subset_neg_max
                t_low(end+1:end+length(xi),1) = xi;
                stroke_wb_L_low(end+1:end+length(xi),1) = stroke_wb_L_interp;
                stroke_wb_R_low(end+1:end+length(xi),1) = stroke_wb_R_interp;
                pitch_wb_L_low(end+1:end+length(xi),1) = pitch_wb_L_interp;
                pitch_wb_R_low(end+1:end+length(xi),1) = pitch_wb_R_interp;
                dev_wb_L_low(end+1:end+length(xi),1) = dev_wb_L_interp;
                dev_wb_R_low(end+1:end+length(xi),1) = dev_wb_R_interp;
                aoa_wb_L_low(end+1:end+length(xi),1) = aoa_wb_L_interp;
                aoa_wb_R_low(end+1:end+length(xi),1) = aoa_wb_R_interp;
                U_wb_L_low(end+1:end+length(xi),1) = U_wb_L_interp;
                U_wb_R_low(end+1:end+length(xi),1) = U_wb_R_interp;
                
                Dstroke_wb_low(end+1:end+length(xi),1) = Dstroke_wb_interp;
                Dpitch_wb_low(end+1:end+length(xi),1) = Dpitch_wb_interp;
                Ddev_wb_low(end+1:end+length(xi),1) = Ddev_wb_interp;
                Daoa_wb_low(end+1:end+length(xi),1) = Daoa_wb_interp;
                DU_wb_low(end+1:end+length(xi),1) = DU_wb_interp;
            end
            
            subplot(3,5,1)
            plot(xi,stroke_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,6)
            plot(xi,stroke_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,11)
            plot(xi,stroke_wb_L_interp-stroke_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,2)
            plot(xi,pitch_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,7)
            plot(xi,pitch_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,12)
            plot(xi,pitch_wb_L_interp-pitch_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,3)
            plot(xi,dev_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,8)
            plot(xi,dev_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,13)
            plot(xi,dev_wb_L_interp-dev_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)

             subplot(3,5,4)
            plot(xi,aoa_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,9)
            plot(xi,aoa_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,14)
            plot(xi,aoa_wb_L_interp-aoa_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,5)
            plot(xi,U_wb_L_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,10)
            plot(xi,U_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,15)
            plot(xi,U_wb_L_interp-U_wb_R_interp,'-','color',color_code_now,'linewidth',linewidth_timelines)
            
%             pause(.001)
        end
    end
end

calc_WBfunc_csaps    
plot_WBfunc_csaps_LRdiff
    
subplot(3,5,1)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,6)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Right wing','fontsize',10)
    grid on
subplot(3,5,11)
axis([0 1 -45 45])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
    ylabel('Left - Right','fontsize',10)
    grid on
            

subplot(3,5,2)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,7)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,12)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,5,3)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,8)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,13)
axis([0 1 -45 45])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
    

subplot(3,5,4)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,9)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
    grid on
subplot(3,5,14)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,5,5)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,10)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,15)
axis([0 1 -3 3])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-3:3:3,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
mid
hi
low
