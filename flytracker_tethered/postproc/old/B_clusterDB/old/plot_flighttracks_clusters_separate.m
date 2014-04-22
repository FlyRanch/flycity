% plot timelines
    subplot(3,3,1)
    hold off
    plot(pathDB.t,patternDB.teta(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) 0 180])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',0:90:180)
    xlabel('time')
    ylabel('optical angle')
    grid on
    subplot(3,3,4)
    hold off
    plot(pathDB.t,pathDB.stim_angle_vel(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -180 180])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-180:180:180)
    xlabel('time')
    ylabel('heading angle')
    grid on
    subplot(3,3,7)
    hold off
    plot(pathDB.t,pathDB.stim_angle_accel(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -180 180])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-180:180:180)
    xlabel('time')
    ylabel('A dir')
    grid on
    
    
    
    subplot(3,3,2)
    hold off
    plot(pathDB.t,pathDB.roll(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -180 180])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-180:180:180)
    xlabel('time')
    ylabel('roll')
    grid on
    subplot(3,3,5)
    hold off
    plot(pathDB.t,pathDB.slip(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -90 90])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-180:90:180)
    xlabel('time')
    ylabel('yaw')
    grid on
    subplot(3,3,8)
    hold off
    plot(pathDB.t,pathDB.pitch(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -45 90])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-180:45:180)
    xlabel('time')
    ylabel('pitch')
    grid on
    

    subplot(3,3,3)
    hold off
    plot(pathDB.t,pathDB.V(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) 0 .8])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',0:.4:.8)
    xlabel('time')
    ylabel('V')
    grid on
    subplot(3,3,6)
    hold off
    plot(pathDB.t,pathDB.An_hor(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -20 20])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',-20:10:20)
    xlabel('time')
    ylabel('An')
    grid on
    subplot(3,3,9)
    hold off
    plot(pathDB.t,pathDB.At_hor(:,i),'-k')
    hold on
    axis([0 max(t(isnan(V)==0)) -15 20])
    set(gca,'XTick',0:.05:.5)
    set(gca,'YTick',[-15;0;20])
    xlabel('time')
    ylabel('At')
    grid on
    
    for j = 1:skip:length(IDX)
        length(IDX)-j;
        if isnan(IDX(j))==0
            
            subplot(3,3,1)
            plot(t(j),teta(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            subplot(3,3,4)
            plot(t(j),stim_angle_vel(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            subplot(3,3,7)
            plot(t(j),stim_angle_accel(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)

            
            subplot(3,3,2)
            plot(t(j),roll(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            subplot(3,3,5)
            plot(t(j),slip(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            subplot(3,3,8)
            plot(t(j),pitch(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)

            
            subplot(3,3,3)
            plot(t(j),V(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)
            subplot(3,3,6)
            plot(t(j),(An_hor(j)),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)
            subplot(3,3,9)
            plot(t(j),At_hor(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)
        end
    end

    % calc start
    if isnan(n_first) == 0
        subplot(3,3,1)
        plot(t_first,teta(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_first,stim_angle_vel(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_first,stim_angle_accel(n_first),'ok','linewidth',2,'MarkerSize',5)
        
        
        subplot(3,3,2)
        plot(t_first,roll(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_first,slip(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_first,pitch(n_first),'ok','linewidth',2,'MarkerSize',5)
        
        
        subplot(3,3,3)
        plot(t_first,V(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_first,(An_hor(n_first)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_first,At_hor(n_first),'ok','linewidth',2,'MarkerSize',5)
    end
    
    % resp start
    if isnan(n_resp) == 0
        
        subplot(3,3,1)
        plot(t_resp,teta(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_resp,stim_angle_vel(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_resp,stim_angle_accel(n_resp),'ok','linewidth',2,'MarkerSize',5)
        
        subplot(3,3,2)
        plot(t_resp,roll(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_resp,slip(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_resp,pitch(n_resp),'ok','linewidth',2,'MarkerSize',5)
        
        subplot(3,3,3)
        plot(t_resp,V(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_resp,(An_hor(n_resp)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_resp,At_hor(n_resp),'ok','linewidth',2,'MarkerSize',5)
    end



    % turn start
    if isnan(n_turn_start) == 0

        subplot(3,3,1)
        plot(t_turn_start,teta(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_turn_start,stim_angle_vel(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_turn_start,stim_angle_accel(n_turn_start),'xk','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_turn_start,roll(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_turn_start,slip(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_turn_start,pitch(n_turn_start),'xk','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_turn_start,V(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_turn_start,(An_hor(n_turn_start)),'xk','linewidth',2,'MarkerSize',5)
    end

    % turn stop
    if isnan(n_turn_stop) == 0

        subplot(3,3,1)
        plot(t_turn_stop,teta(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_turn_stop,stim_angle_vel(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_turn_stop,stim_angle_accel(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_turn_stop,roll(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_turn_stop,slip(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_turn_stop,pitch(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_turn_stop,V(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_turn_stop,(An_hor(n_turn_stop)),'xk','linewidth',2,'MarkerSize',5)
    end



    % turn max
    if isnan(n_turn_max) == 0
        subplot(3,3,1)
        plot(t_turn_max,teta(n_turn_max),'dk','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_turn_max,stim_angle_vel(n_turn_max),'dk','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_turn_max,stim_angle_accel(n_turn_max),'dk','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_turn_max,roll(n_turn_max),'dk','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_turn_max,slip(n_turn_max),'dk','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_turn_max,pitch(n_turn_max),'dk','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_turn_max,V(n_turn_max),'dk','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_turn_max,(An_hor_max),'dk','linewidth',2,'MarkerSize',5)
    end
                    
                    
                    
  % accel start
    if isnan(n_accel_start) == 0

        subplot(3,3,1)
        plot(t_accel_start,teta(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_accel_start,stim_angle_vel(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_accel_start,stim_angle_accel(n_accel_start),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_accel_start,roll(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_accel_start,slip(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_accel_start,pitch(n_accel_start),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_accel_start,V(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_accel_start,At_hor(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
    end

    % accel stop
    if isnan(n_accel_stop) == 0

        subplot(3,3,1)
        plot(t_accel_stop,teta(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
         subplot(3,3,4)
       plot(t_accel_stop,stim_angle_vel(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_accel_stop,stim_angle_accel(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_accel_stop,roll(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_accel_stop,slip(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_accel_stop,pitch(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_accel_stop,V(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_accel_stop,At_hor(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
    end
                    

    % accel max
    if isnan(n_accel_max) == 0

        subplot(3,3,1)
        plot(t_accel_max,teta(n_accel_max),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_accel_max,stim_angle_vel(n_accel_max),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_accel_max,stim_angle_accel(n_accel_max),'sk','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_accel_max,roll(n_accel_max),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_accel_max,slip(n_accel_max),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_accel_max,pitch(n_accel_max),'sk','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_accel_max,V(n_accel_max),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_accel_max,At_hor_max,'sk','linewidth',2,'MarkerSize',5)
    end

    
  % decel start
    if isnan(n_decel_start) == 0

        subplot(3,3,1)
        plot(t_decel_start,teta(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_decel_start,stim_angle_vel(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_decel_start,stim_angle_accel(n_decel_start),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_decel_start,roll(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_decel_start,slip(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_decel_start,pitch(n_decel_start),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_decel_start,V(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_decel_start,At_hor(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
    end

    % decel stop
    if isnan(n_decel_stop) == 0

        subplot(3,3,1)
        plot(t_decel_stop,teta(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
         subplot(3,3,4)
       plot(t_decel_stop,stim_angle_vel(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_decel_stop,stim_angle_accel(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_decel_stop,roll(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_decel_stop,slip(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_decel_stop,pitch(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_decel_stop,V(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_decel_stop,At_hor(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
    end
                    

    % decel max
    if isnan(n_decel_min) == 0

        subplot(3,3,1)
        plot(t_decel_min,teta(n_decel_min),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_decel_min,stim_angle_vel(n_decel_min),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_decel_min,stim_angle_accel(n_decel_min),'sk','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_decel_min,roll(n_decel_min),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_decel_min,slip(n_decel_min),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_decel_min,pitch(n_decel_min),'sk','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_decel_min,V(n_decel_min),'sk','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_decel_min,At_hor_min,'sk','linewidth',2,'MarkerSize',5)
    end

    % resp end
    if isnan(n_resp_end) == 0

        subplot(3,3,1)
        plot(t_resp_end,teta(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_resp_end,stim_angle_vel(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_resp_end,stim_angle_accel(n_resp_end),'ok','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_resp_end,roll(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_resp_end,slip(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_resp_end,pitch(n_resp_end),'ok','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_resp_end,V(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_resp_end,(An_hor(n_resp_end)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_resp_end,At_hor(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
    end


    % steady end
    if isnan(n_steady_end) == 0

        subplot(3,3,1)
        plot(t_steady_end,teta(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,4)
        plot(t_steady_end,stim_angle_vel(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,7)
        plot(t_steady_end,stim_angle_accel(n_steady_end),'ok','linewidth',2,'MarkerSize',5)

        subplot(3,3,2)
        plot(t_steady_end,roll(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,5)
        plot(t_steady_end,slip(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,8)
        plot(t_steady_end,pitch(n_steady_end),'ok','linewidth',2,'MarkerSize',5)

        subplot(3,3,3)
        plot(t_steady_end,V(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,6)
        plot(t_steady_end,(An_hor(n_steady_end)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,3,9)
        plot(t_steady_end,At_hor(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
    end
