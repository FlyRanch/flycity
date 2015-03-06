
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

mkdir('wb_seqs')
cd('wb_seqs')
for seq = 1:size(n_down_start_L,2)
    
    
    for wb = 2:size(n_down_start_L,1)-1
        
%         counter = 100*(size(n_down_start_L,2)-seq) +size(n_down_start_L,1) -wb
%         
%         color_code = round(200*wb/size(n_down_start_L,1));
        
        n_down_start_L_now = n_down_start_L(wb,seq);
        n_up_stop_L_now = n_down_start_L(wb+1,seq);
        n_down_start_R_now = n_down_start_R(wb,seq);
        n_up_stop_R_now = n_down_start_R(wb+1,seq);

%         % only within steady&maneuver
%         if n_down_start_L_now > n_first(seq) && n_up_stop_L_now < n_post(seq) && ...
%                 n_down_start_R_now > n_first(seq) && n_up_stop_R_now < n_post(seq)
            rollaccel_now = roll_dot_dot_mean_wb(wb,seq);

          if isnan(n_down_start_L_now)==0 && isnan(n_up_stop_L_now)==0 && isnan(rollaccel_now)==0
              
%             color_code_now = color_code(wb,seq);
%             color_code_now = round(200*wb/find(isnan(n_down_start_L(:,seq))==0, 1, 'last' ));
            color_code_now = round(100+100*rollaccel_now/max(abs(roll_dot_dot_mean_wb(:,seq))));
            if color_code_now>200
                color_code_now = 200
            end
            
            if color_code_now<1
                color_code_now = 1
            end
            
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

            subplot(3,5,1)
            plot(xi,stroke_wb_L_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,6)
            plot(xi,stroke_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,11)
            plot(xi,stroke_wb_L_interp-stroke_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on

            subplot(3,5,2)
            plot(xi,pitch_wb_L_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,7)
            plot(xi,pitch_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,12)
            plot(xi,pitch_wb_L_interp-pitch_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on

            subplot(3,5,3)
            plot(xi,dev_wb_L_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,8)
            plot(xi,dev_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,13)
            plot(xi,dev_wb_L_interp-dev_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on

             subplot(3,5,4)
            plot(xi,aoa_wb_L_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,9)
            plot(xi,aoa_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,14)
            plot(xi,aoa_wb_L_interp-aoa_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on

            subplot(3,5,5)
            plot(xi,U_wb_L_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,10)
            plot(xi,U_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
            subplot(3,5,15)
            plot(xi,U_wb_L_interp-U_wb_R_interp,'-','color',cmap_jet(color_code_now,:),'linewidth',linewidth_timelines)
            hold on
        end
    end
    
subplot(3,5,1)
hold off
subplot(3,5,2)
hold off
subplot(3,5,3)
hold off
subplot(3,5,4)
hold off
subplot(3,5,5)
hold off
subplot(3,5,6)
hold off
subplot(3,5,7)
hold off
subplot(3,5,8)
hold off
subplot(3,5,9)
hold off
subplot(3,5,10)
hold off
subplot(3,5,11)
hold off
subplot(3,5,12)
hold off
subplot(3,5,13)
hold off
subplot(3,5,14)
hold off
subplot(3,5,15)
hold off

seq_now = settings.seq(seq,:)
    saveas(gca,['wingkin_rollaccel_',num2str(seq_now(1)),'seq',num2str(seq_now(2)),'.fig'])
    saveas(gca,['wingkin_rollaccel_',num2str(seq_now(1)),'seq',num2str(seq_now(2)),'.png'])
    plot2svg(['wingkin_rollaccel_',num2str(seq_now(1)),'seq',num2str(seq_now(2)),'.svg'])
    
%     pause
end

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
            