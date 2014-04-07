
figure
subplot(2,5,1)
title('wing stroke')
hold on
subplot(2,5,2)
title('wing pitch')
hold on
subplot(2,5,3)
title('stroke deviation')
hold on
subplot(2,5,4)
title('angle of attack')
hold on
subplot(2,5,5)
title('wing speed')
hold on
subplot(2,5,6)
hold on
subplot(2,5,7)
hold on
subplot(2,5,8)
hold on
subplot(2,5,9)
hold on
subplot(2,5,10)
hold on


for seq = 1:size(n_down_start_L,2)
    for wb = 2:size(n_down_start_L,1)-1
        
        counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        n_down_start_L_now = n_down_start_L(wb,seq);
        n_up_stop_L_now = n_down_start_L(wb+1,seq);
        n_down_start_R_now = n_down_start_R(wb,seq);
        n_up_stop_R_now = n_down_start_R(wb+1,seq);

        % only within steady&maneuver
        if n_down_start_L_now > n_first(seq) && n_up_stop_L_now < n_post(seq) && ...
                n_down_start_R_now > n_first(seq) && n_up_stop_R_now < n_post(seq)
            
            color_code_now1 = color_code1(wb,seq);
            color_code_now2 = color_code2(wb,seq);

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

            subplot(2,5,1)
            plot(xi,stroke_wb_L_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
%             plot(xi,stroke_wb_R_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
            subplot(2,5,6)
%             plot(xi,stroke_wb_L_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)
            plot(xi,stroke_wb_R_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)

            subplot(2,5,2)
            plot(xi,pitch_wb_L_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
            subplot(2,5,7)
            plot(xi,pitch_wb_R_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)

            subplot(2,5,3)
            plot(xi,dev_wb_L_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
            subplot(2,5,8)
            plot(xi,dev_wb_R_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)

             subplot(2,5,4)
            plot(xi,aoa_wb_L_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
            subplot(2,5,9)
            plot(xi,aoa_wb_R_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)

            subplot(2,5,5)
            plot(xi,U_wb_L_interp,'-','color',cmap_jet(color_code_now1,:),'linewidth',linewidth_timelines)
            subplot(2,5,10)
            plot(xi,U_wb_R_interp,'-','color',cmap_jet(color_code_now2,:),'linewidth',linewidth_timelines)
        end
    end
end

subplot(2,5,1)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('body pitch','fontsize',10)
    grid on
subplot(2,5,6)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
    ylabel('flight force','fontsize',10)
    grid on
            
subplot(2,5,2)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(2,5,7)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
            
subplot(2,5,3)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(2,5,8)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
            
subplot(2,5,4)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(2,5,9)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
    grid on
            
subplot(2,5,5)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(2,5,10)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
