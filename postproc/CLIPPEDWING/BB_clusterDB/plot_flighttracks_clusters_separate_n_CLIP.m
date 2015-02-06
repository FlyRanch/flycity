% plot timelines

    subplot(3,3,1)
    hold off
    plot(t,stim_angle_vel,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -180 180])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-180:90:180)
    xlabel('time')
    ylabel('heading angle')
    grid on
    
    subplot(3,3,4)
    hold off
    plot(t,stim_angle_accel,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -180 180])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-180:90:180)
    xlabel('time')
    ylabel('accel angle')
    grid on
    
    subplot(3,3,7)
    hold off
    plot(t,V,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) 0 .6])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',0:.2:6)
    xlabel('time')
    ylabel('V')
    grid on
    
    subplot(3,3,2)
    hold off
    plot(t,roll,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -90 90])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-180:45:180)
    xlabel('time')
    ylabel('roll')
    grid on
    
    subplot(3,3,5)
    hold off
    plot(t,slip,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -90 90])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-180:45:180)
    xlabel('time')
    ylabel('yaw')
    grid on
    
    subplot(3,3,8)
    hold off
    plot(t,pitch,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) 0 90])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-180:30:180)
    xlabel('time')
    ylabel('pitch')
    grid on
    
    subplot(3,3,3)
    hold off
    plot(t,A,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -6 6])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-10:2:10)
    xlabel('time')
    ylabel('A')
    grid on
    
    subplot(3,3,6)
    hold off
    plot(t,A_hor,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -6 6])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',-10:2:10)
    xlabel('time')
    ylabel('horizontal A')
    grid on
    
    subplot(3,3,9)
    hold off
    plot(t,A_ver,'-k')
    hold on
    axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -6 6])
    % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
    set(gca,'YTick',[-10:2:10])
    xlabel('time')
    ylabel('vertical A')
    grid on
%     
%     subplot(3,3,6)
%     hold off
%     plot(t,An_hor,'-k')
%     hold on
%     axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -6 6])
%     % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
%     set(gca,'YTick',-10:2:10)
%     xlabel('time')
%     ylabel('An')
%     grid on
%     
%     subplot(3,3,9)
%     hold off
%     plot(t,At_hor,'-k')
%     hold on
%     axis([min(t(isnan(V)==0)) max(t(isnan(V)==0)) -6 6])
%     % set(gca,'XTick',min(t(isnan(V)==0)):max(t(isnan(V)==0)))
%     set(gca,'YTick',[-10:2:10])
%     xlabel('time')
%     ylabel('At')
%     grid on
    
    for j = start_frame:skip:length(steady_frames)
        length(steady_frames)-j;
        if isnan(steady_frames(j))==0
            
            subplot(3,3,1)
            plot(t(j),stim_angle_vel(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',5)
            subplot(3,3,4)
            plot(t(j),stim_angle_accel(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',5)
            subplot(3,3,7)
            plot(t(j),(V(j)),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)

            subplot(3,3,2)
            plot(t(j),roll(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',5)
            subplot(3,3,5)
            plot(t(j),slip(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',5)
            subplot(3,3,8)
            plot(t(j),pitch(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',5)

            subplot(3,3,3)
            plot(t(j),A(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)
            subplot(3,3,6)
            plot(t(j),(A_hor(j)),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)
            subplot(3,3,9)
            plot(t(j),A_ver(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)
%             subplot(3,3,6)
%             plot(t(j),(An_hor(j)),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)
%             subplot(3,3,9)
%             plot(t(j),At_hor(j),'.','color',cmap_k(2-steady_frames(j),:),'MarkerSize',10)
        end
    end
