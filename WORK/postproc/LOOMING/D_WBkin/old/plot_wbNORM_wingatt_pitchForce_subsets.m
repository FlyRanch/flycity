
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

plot_pitch = 0;
mid_pitch=0;
hi_pitch=0;
low_pitch=0;

plot_F = 0;
mid_F=0;
hi_F=0;
low_F=0;

for seq = 1:size(n_down_start_L,2)
    counter = (size(n_down_start_L,2)-seq)
    for wb = 2:size(n_down_start_L,1)-1
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
        
        subset_pitch_now = subset_pitch(wb,seq);
        subset_F_now = subset_F(wb,seq);
        
        
        if         abs(subset_pitch_now) < subset_pitch_mid || abs(subset_pitch_now) > subset_pitch_max ||...
                   abs(subset_F_now) < subset_F_mid || abs(subset_F_now) > subset_F_max
            
            % plot pitch?
            if abs(subset_pitch_now) < subset_pitch_mid
                plot_pitch = 1;
                color_code_pitch_now = 'g';
                mid_pitch=mid_pitch+1;
            elseif subset_pitch_now > subset_pitch_max
                plot_pitch = 1;
                color_code_pitch_now = 'r';
                hi_pitch=hi_pitch+1;
            elseif subset_pitch_now < -subset_pitch_max
                plot_pitch = 1;
                color_code_pitch_now = 'b';
                low_pitch=low_pitch+1;
            else
                plot_pitch = 0;
            end
            
            % plot F?
            if abs(subset_F_now) < subset_F_mid
                plot_F = 1;
                color_code_F_now = 'g';
                mid_F=mid_F+1;
            elseif subset_F_now > subset_F_max
                plot_F = 1;
                color_code_F_now = 'r';
                hi_F=hi_F+1;
            elseif subset_F_now < -subset_F_max
                plot_F = 1;
                color_code_F_now = 'b';
                low_F=low_F+1;
            else
                plot_F = 0;
            end
            
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

            if plot_pitch ==1
                subplot(2,5,1)
                plot(xi,stroke_wb_L_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                plot(xi,stroke_wb_R_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                subplot(2,5,2)
                plot(xi,pitch_wb_L_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                plot(xi,pitch_wb_R_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                subplot(2,5,3)
                plot(xi,dev_wb_L_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                plot(xi,dev_wb_R_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                 subplot(2,5,4)
                plot(xi,aoa_wb_L_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                plot(xi,aoa_wb_R_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                subplot(2,5,5)
                plot(xi,U_wb_L_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
                plot(xi,U_wb_R_interp,'-','color',color_code_pitch_now,'linewidth',linewidth_timelines)
            end
            
            if plot_F ==1
                subplot(2,5,6)
                plot(xi,stroke_wb_L_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                plot(xi,stroke_wb_R_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                subplot(2,5,7)
                plot(xi,pitch_wb_L_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                plot(xi,pitch_wb_R_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                subplot(2,5,8)
                plot(xi,dev_wb_L_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                plot(xi,dev_wb_R_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                subplot(2,5,9)
                plot(xi,aoa_wb_L_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                plot(xi,aoa_wb_R_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                subplot(2,5,10)
                plot(xi,U_wb_L_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
                plot(xi,U_wb_R_interp,'-','color',color_code_F_now,'linewidth',linewidth_timelines)
            end
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

mid_pitch
hi_pitch
low_pitch
mid_F
hi_F
low_F
