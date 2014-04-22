
clc
clear
close all

% addpath('/home/florian/Dropbox/WORK/toolbox')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_9clusters.mat')

t = pathDB.t;
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t = t(trigger_frame:end,:);

for i = 1:size(pathDB.pos,2)
% for i = 4
    
    x = pathDB.pos(trigger_frame:end,i,1);
    y = pathDB.pos(trigger_frame:end,i,2);
    z = pathDB.pos(trigger_frame:end,i,3);

    u = pathDB.vel(trigger_frame:end,i,1);
    v = pathDB.vel(trigger_frame:end,i,2);
    w = pathDB.vel(trigger_frame:end,i,3);

    ax = pathDB.accel(trigger_frame:end,i,1);
    ay = pathDB.accel(trigger_frame:end,i,2);
    az = pathDB.accel(trigger_frame:end,i,3);


    stim_angle = pathDB.stim_angle(trigger_frame:end,i,1);
    IDX = pathDB.IDX(trigger_frame:end,i,1);
    
    cmap_k = settings.cmap_k;
    cmap_360 = settings.cmap_360;

    % vel & accel
    V = sqrt(u.^2 + v.^2 + w.^2);
    A = sqrt(ax.^2 + ay.^2 + az.^2);
    A_hor = sqrt(ax.^2 + ay.^2);

    % A
    An = pathDB.An(trigger_frame:end,i);
    At = pathDB.At(trigger_frame:end,i);
    An_hor = pathDB.An_hor(trigger_frame:end,i);
    At_hor = pathDB.At_hor(trigger_frame:end,i);
    alpha_dot_hor = pathDB.alpha_dot_hor(trigger_frame:end,i);
    
    % plot An-At
    
    
    % plot timelines
    figure(2)
    subplot(3,1,1)
    hold off
    subplot(3,1,2)
    hold off
    subplot(3,1,3)
    hold off
    for j = 1:10:length(IDX)
        j;
        if isnan(IDX(j))==0
            subplot(3,1,1)
            plot(t(j),stim_angle(j),'.','color',cmap_k(IDX(j),:))
            ylabel('heading')
            hold on

            subplot(3,1,2)
            plot(t(j),An_hor(j),'.','color',cmap_k(IDX(j),:))
            ylabel('An')
            hold on

            subplot(3,1,3)
            plot(t(j),At_hor(j),'.','color',cmap_k(IDX(j),:))
            ylabel('At')
            hold on
        end
    end
    
    
    % cluster shifts
    c=1;
    t_shift(c,1) = t(1);
    IDX_shift(c,1) = IDX(1);
    for j = 2:length(IDX)
        if isnan(IDX(j)) == 0 && IDX(j) ~= IDX(j-1)
            c=c+1;
            t_shift(c,1) = t(j);
            IDX_shift(c,1) = IDX(j);
        end
    end
    
    
    if IDX_shift(1) == 5
        
        % response time
        t_resp = t_shift(2);
        n_resp = find(t==t_resp);
        IDX_resp = IDX_shift(2);
        stim_angle_pre_resp = stim_angle(n_resp);
        V_pre_resp = V(n_resp);
        
        subplot(3,1,1)
        plot(t_resp,stim_angle(n_resp),'ok')
        subplot(3,1,2)
        plot(t_resp,An_hor(n_resp),'ok')
        subplot(3,1,3)
        plot(t_resp,At_hor(n_resp),'ok')

        
        % start of turn
        if min(IDX_shift)<4 || max(IDX_shift)>6
            
            turn_start_nr = min([find(IDX_shift>6);find(IDX_shift<4)]);
            IDX_turn_start = IDX_shift(turn_start_nr);
            t_turn_start = t_shift(turn_start_nr);
            n_turn_start = find(t==t_turn_start);
            stim_angle_pre_turn = stim_angle(n_turn_start);
            
            subplot(3,1,1)
            plot(t_turn_start,stim_angle(n_turn_start),'xk')
            subplot(3,1,2)
            plot(t_turn_start,An_hor(n_turn_start),'xk')
            
            % end of turn
            % left turn
            if IDX_turn_start < 4
                turn_stop = 0;
                n = find(t==t_turn_start);
                while turn_stop == 0
                    n = n+1;
                    if n > length(t)
                        n_turn_stop = nan;
                        t_turn_stop = nan;
                        stim_angle_post_turn = nan;
                        turn_stop = 1;
                    elseif IDX(n) > 3
                        n_turn_stop = n;
                        t_turn_stop = t(n);
                        stim_angle_post_turn = stim_angle(n);
                        turn_stop = 1;
                        
                        subplot(3,1,1)
                        plot(t_turn_stop,stim_angle(n_turn_stop),'xk')
                        subplot(3,1,2)
                        plot(t_turn_stop,An_hor(n_turn_stop),'xk')
                    end
                end
                
                dt_turn = t_turn_stop - t_turn_start;
                
                % max neg turn accel
                if isnan(dt_turn)==0
                    An_hor_min = min(An_hor(n_turn_start:n_turn_stop));
                else
                    An_hor_min = min(An_hor(n_turn_start:end));
                end
                n_An_hor_min = find(An_hor==An_hor_min);
                t_An_hor_min = t(n_An_hor_min);

                subplot(3,1,1)
                plot(t_An_hor_min,stim_angle(n_An_hor_min),'xk')
                subplot(3,1,2)
                plot(t_An_hor_min,An_hor_min,'xk')
                
            % right turn
            elseif IDX_turn_start > 6
                turn_stop = 0;
                n = find(t==t_turn_start);
                while turn_stop == 0
                    n = n+1;
                    if n > length(t)
                        n_turn_stop = nan;
                        t_turn_stop = nan;
                        stim_angle_post_turn = nan;
                        turn_stop = 1;
                    elseif IDX(n) < 7
                        n_turn_stop = n;
                        t_turn_stop = t(n);
                        stim_angle_post_turn = stim_angle(n);
                        turn_stop = 1;
                        
                        subplot(3,1,1)
                        plot(t_turn_stop,stim_angle(n_turn_stop),'xk')
                        subplot(3,1,2)
                        plot(t_turn_stop,An_hor(n_turn_stop),'xk')
                    end
                end
                
                dt_turn = t_turn_stop - t_turn_start;
                
                % max pos turn accel
                if isnan(dt_turn)==0
                    An_hor_max = max(An_hor(n_turn_start:n_turn_stop));
                else
                    An_hor_max = max(An_hor(n_turn_start:end));
                end
                n_An_hor_max = find(An_hor==An_hor_max);
                t_An_hor_max = t(n_An_hor_max);

                subplot(3,1,1)
                plot(t_An_hor_max,stim_angle(n_An_hor_max),'xk')
                subplot(3,1,2)
                plot(t_An_hor_max,An_hor(n_An_hor_max),'xk')
            end
        end

         
        % start of accelleration
        if isempty(find(IDX_shift==3, 1))==0 || isempty(find(IDX_shift==6, 1))==0 || isempty(find(IDX_shift==9, 1))==0 
            
            accel_start_nr = min([find(IDX_shift==3); find(IDX_shift==6); find(IDX_shift==9)]);
            IDX_accel_start = IDX_shift(accel_start_nr);
            t_accel_start = t_shift(accel_start_nr);
            n_accel_start = find(t==t_accel_start);
            V_pre_accel = V(n_accel_start);
            
            subplot(3,1,1)
            plot(t_accel_start,stim_angle(n_accel_start),'+k')
            subplot(3,1,3)
            plot(t_accel_start,At_hor(n_accel_start),'+k')

            % end of acceleration
            accel_stop = 0;
            n = find(t==t_accel_start);
            while accel_stop == 0
                n = n+1;
                if n > length(t)
                    t_accel_stop = nan;
                    accel_stop = 1;
                elseif IDX(n) == 3 || IDX(n) == 6 || IDX(n) == 9
                else
                    t_accel_stop = t(n);
                    n_accel_stop = n;
                    V_post_accel = V(n);
                    accel_stop = 1;
                    
                    subplot(3,1,1)
                    plot(t_accel_stop,stim_angle(n_accel_stop),'+k')
                    subplot(3,1,3)
                    plot(t_accel_stop,At_hor(n_accel_stop),'+k')
                end
            end
            
            dt_accel = t_accel_stop - t_accel_start;
            
            if isnan(dt_accel)==0
                At_hor_max = max(At_hor(n_accel_start:n_accel_stop));
                n_At_hor_max = find(At_hor==At_hor_max);
                t_At_hor_max = t(n_At_hor_max);
                
                subplot(3,1,1)
                plot(t_At_hor_max,stim_angle(n_At_hor_max),'+k')
                subplot(3,1,3)
                plot(t_At_hor_max,At_hor(n_At_hor_max),'+k')
            else
                At_hor_max = max(At_hor(n_accel_start:end));
                n_At_hor_max = find(At_hor==At_hor_max);
                t_At_hor_max = t(n_At_hor_max);
                
                subplot(3,1,1)
                plot(t_At_hor_max,stim_angle(n_At_hor_max),'+k')
                subplot(3,1,3)
                plot(t_At_hor_max,At_hor(n_At_hor_max),'+k')
            end
        end                
         
        % start of deceleration
        if isempty(find(IDX_shift==1, 1))==0 || isempty(find(IDX_shift==4, 1))==0 || isempty(find(IDX_shift==7, 1))==0 
            
            decel_start_nr = min([find(IDX_shift==1); find(IDX_shift==4); find(IDX_shift==7)]);
            IDX_decel_start = IDX_shift(decel_start_nr);
            t_decel_start = t_shift(decel_start_nr);
            n_decel_start = find(t==t_decel_start);
            V_pre_decel = V(n_decel_start);

            subplot(3,1,1)
            plot(t_decel_start,stim_angle(n_decel_start),'+k')
            subplot(3,1,3)
            plot(t_decel_start,At_hor(n_decel_start),'+k')
            
            % end of deceleration
            decel_stop = 0;
            n = find(t==t_decel_start);
            while decel_stop == 0
                n = n+1;
                if n > length(t)
                    t_decel_stop = nan;
                    n_decel_stop = nan;
                    V_post_decel = nan;
                    decel_stop = 1;
                elseif IDX(n) == 1 || IDX(n) == 4 || IDX(n) == 7
                else
                    t_decel_stop = t(n);
                    n_decel_stop = n;
                    V_post_decel = V(n);
                    decel_stop = 1;
                    
                    subplot(3,1,1)
                    plot(t_decel_stop,stim_angle(n_decel_stop),'+k')
                    subplot(3,1,3)
                    plot(t_decel_stop,At_hor(n_decel_stop),'+k')
                end
            end
            
            dt_decel = t_decel_stop - t_decel_start;
            
            if isnan(dt_decel)==0
                At_hor_min = min(At_hor(n_decel_start:n_decel_stop));
                n_At_hor_min = find(At_hor==At_hor_min);
                t_At_hor_min = t(n_At_hor_min);
                
                subplot(3,1,1)
                plot(t_At_hor_min,stim_angle(n_At_hor_min),'+k')
                subplot(3,1,3)
                plot(t_At_hor_min,At_hor(n_At_hor_min),'+k')
            else
                At_hor_min = min(At_hor(n_decel_start:end));
                n_At_hor_min = find(At_hor==At_hor_min);
                t_At_hor_min = t(n_At_hor_min);
                
                subplot(3,1,1)
                plot(t_At_hor_min,stim_angle(n_At_hor_min),'+k')
                subplot(3,1,3)
                plot(t_At_hor_min,At_hor(n_At_hor_min),'+k')
            end
        end
    end
    pause
end


    