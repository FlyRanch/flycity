function [IDX_shift,t_shift,t_resp,...
    t_turn_start,t_turn_stop,t_turn_max,...
    t_accel_start,t_accel_stop,t_accel_max,...
    t_decel_start,t_decel_stop,t_decel_max,...
    dt_turn,dt_accel,dt_decel,...
    An_hor_max,At_hor_max,At_hor_min,...
    stim_angle_trigger,stim_angle_pre_resp,stim_angle_trig2resp,stim_angle_pre_turn,...
    stim_angle_post_turn,stim_angle_post_resp,stim_angle_postresp2end,...
    V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
    V_post_resp,V_postresp2end,V_post_accel] = response_times_v2_hVA(pathDB,settings,i,toplot)

IDX_shift = nan;
t_shift = nan;

t_resp = nan;

t_turn_start = nan;
t_turn_stop = nan;
t_turn_max = nan;

t_accel_start = nan;
t_accel_stop = nan;
t_accel_max = nan; 

t_decel_start = nan;
t_decel_stop = nan;
t_decel_max = nan;

dt_turn = nan;
dt_accel = nan;
dt_decel = nan;

An_hor_max = nan;
At_hor_max = nan;
At_hor_min = nan;

stim_angle_trigger = nan;
stim_angle_pre_resp = nan;
stim_angle_trig2resp = nan;
stim_angle_pre_turn = nan;

stim_angle_post_turn = nan;
stim_angle_post_resp = nan;   
stim_angle_postresp2end = nan;   

V_trigger = nan;   
V_pre_resp = nan;   
V_trig2resp = nan;   
V_pre_accel = nan;   

V_post_resp = nan;   
V_postresp2end = nan;   
V_post_accel = nan;
    
t = pathDB.t;
trigger_frame = find(t == min(abs(t)));
dt = t(2) - t(1);
t = t(trigger_frame:end,:);

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

% plot timelines
if toplot == 1
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
            plot(t(j),V(j),'.','color',cmap_k(IDX(j),:))
            ylabel('V')
            hold on

            subplot(3,1,3)
            plot(t(j),An_hor(j),'.','color',cmap_k(IDX(j),:))
            hold on
            plot(t(j),At_hor(j),'.','color',cmap_k(IDX(j),:))
            ylabel('A')
            hold on
        end
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


if IDX_shift(1) == 5 && length(IDX_shift) > 1
    
    % @ trigger
    stim_angle_trigger = stim_angle(1);
    V_trigger = V(1);

    % response time
    t_resp = t_shift(2);
    n_resp = find(t==t_resp);
    IDX_resp = IDX_shift(2);
    stim_angle_pre_resp = mean(stim_angle(n_resp));
    V_pre_resp = V(n_resp);
    
    if toplot == 1
        subplot(3,1,1)
        plot(t_resp,stim_angle(n_resp),'ok','linewidth',2)
        subplot(3,1,2)
        plot(t_resp,V(n_resp),'ok','linewidth',2)
        subplot(3,1,3)
        plot(t_resp,An_hor(n_resp),'ok','linewidth',2)
        plot(t_resp,At_hor(n_resp),'ok','linewidth',2)
    end
    
    
    %% turn
    if min(IDX_shift)<4 || max(IDX_shift)>6

        % start of turn
        turn_start_nr = min([find(IDX_shift>6);find(IDX_shift<4)]);
        IDX_turn_start = IDX_shift(turn_start_nr);
        t_turn_start = t_shift(turn_start_nr);
        n_turn_start = find(t==t_turn_start);
        stim_angle_pre_turn = stim_angle(n_turn_start);

        if toplot == 1
            subplot(3,1,1)
            plot(t_turn_start,stim_angle(n_turn_start),'xk','linewidth',2)
            subplot(3,1,2)
            plot(t_turn_start,V(n_turn_start),'xk','linewidth',2)
            subplot(3,1,3)
            plot(t_turn_start,An_hor(n_turn_start),'xk','linewidth',2)
        end

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
                    
                    if toplot == 1
                        subplot(3,1,1)
                        plot(t_turn_stop,stim_angle(n_turn_stop),'xk','linewidth',2)
                        subplot(3,1,2)
                        plot(t_turn_stop,V(n_turn_stop),'xk','linewidth',2)
                        subplot(3,1,3)
                        plot(t_turn_stop,An_hor(n_turn_stop),'xk','linewidth',2)
                    end
                end
            end
            
            % turn time
            dt_turn = t_turn_stop - t_turn_start;

            % max NEGATIVE turn accel
            if isnan(dt_turn)==0
                An_hor_max = min(An_hor(n_turn_start:n_turn_stop));
            else
                An_hor_max = min(An_hor(n_turn_start:end));
            end
            n_turn_max = find(An_hor==An_hor_max);
            t_turn_max = t(n_turn_max);

            if toplot == 1
                subplot(3,1,1)
                plot(t_turn_max,stim_angle(n_turn_max),'xk','linewidth',2)
                subplot(3,1,2)
                plot(t_turn_max,V,'xk','linewidth',2)
                subplot(3,1,3)
                plot(t_turn_max,An_hor_max,'xk','linewidth',2)
            end

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


                    if toplot == 1
                        subplot(3,1,1)
                        plot(t_turn_stop,stim_angle(n_turn_stop),'xk','linewidth',2)
                        subplot(3,1,2)
                        plot(t_turn_stop,V(n_turn_stop),'xk','linewidth',2)
                        subplot(3,1,3)
                        plot(t_turn_stop,An_hor(n_turn_stop),'xk','linewidth',2)
                    end
                end
            end
            
            % turn time
            dt_turn = t_turn_stop - t_turn_start;

            % max POSITIVE turn accel
            if isnan(dt_turn)==0
                An_hor_max = max(An_hor(n_turn_start:n_turn_stop));
            else
                An_hor_max = max(An_hor(n_turn_start:end));
            end
            n_turn_max = find(An_hor==An_hor_max);
            t_turn_max = t(n_turn_max);

            if toplot == 1
                subplot(3,1,1)
                plot(t_turn_max,stim_angle(n_turn_max),'xk','linewidth',2)
                subplot(3,1,2)
                plot(t_turn_max,V(n_turn_max),'xk','linewidth',2)
                subplot(3,1,3)
                plot(t_turn_max,An_hor_max,'xk','linewidth',2)
            end
        end
    end
    
    %% Accelleration
    if isempty(find(IDX_shift==3, 1))==0 || isempty(find(IDX_shift==6, 1))==0 || isempty(find(IDX_shift==9, 1))==0 

        % start of accelleration
        accel_start_nr = min([find(IDX_shift==3); find(IDX_shift==6); find(IDX_shift==9)]);
        IDX_accel_start = IDX_shift(accel_start_nr);
        t_accel_start = t_shift(accel_start_nr);
        n_accel_start = find(t==t_accel_start);
        V_pre_accel = V(n_accel_start);

        if toplot == 1
            subplot(3,1,1)
            plot(t_accel_start,stim_angle(n_accel_start),'+k','linewidth',2)
            subplot(3,1,2)
            plot(t_accel_start,V(n_accel_start),'+k','linewidth',2)
            subplot(3,1,3)
            plot(t_accel_start,At_hor(n_accel_start),'+k','linewidth',2)
        end

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

                if toplot == 1
                    subplot(3,1,1)
                    plot(t_accel_stop,stim_angle(n_accel_stop),'+k','linewidth',2)
                    subplot(3,1,2)
                    plot(t_accel_stop,V(n_accel_stop),'+k','linewidth',2)
                    subplot(3,1,3)
                    plot(t_accel_stop,At_hor(n_accel_stop),'+k','linewidth',2)
                end
            end
        end
        
        % accel time
        dt_accel = t_accel_stop - t_accel_start;
        
        %% max accel
        if isnan(dt_accel)==0
            At_hor_max = max(At_hor(n_accel_start:n_accel_stop));
            n_accel_max = find(At_hor==At_hor_max);
            t_accel_max = t(n_accel_max);

            if toplot == 1
                subplot(3,1,1)
                plot(t_accel_max,stim_angle(n_accel_max),'+k','linewidth',2)
                subplot(3,1,2)
                plot(t_accel_max,V(n_accel_max),'+k','linewidth',2)
                subplot(3,1,3)
                plot(t_accel_max,At_hor_max,'+k','linewidth',2)
            end
        else
            At_hor_max = max(At_hor(n_accel_start:end));
            n_accel_max = find(At_hor==At_hor_max);
            t_accel_max = t(n_accel_max);

            if toplot == 1
                subplot(3,1,1)
                plot(t_accel_max,stim_angle(n_accel_max),'+k','linewidth',2)
                subplot(3,1,2)
                plot(t_accel_max,V(n_accel_max),'+k','linewidth',2)
                subplot(3,1,3)
                plot(t_accel_max,At_hor_max,'+k','linewidth',2)
            end
        end
    end                
    
    %% Deceleration
    if isempty(find(IDX_shift==1, 1))==0 || isempty(find(IDX_shift==4, 1))==0 || isempty(find(IDX_shift==7, 1))==0 

        % start of deceleration
        decel_start_nr = min([find(IDX_shift==1); find(IDX_shift==4); find(IDX_shift==7)]);
        
        if isnan(t_accel_start)==1 || decel_start_nr < accel_start_nr % decel before accel
            IDX_decel_start = IDX_shift(decel_start_nr);
            t_decel_start = t_shift(decel_start_nr);
            n_decel_start = find(t==t_decel_start);
            V_pre_decel = V(n_decel_start);

            if toplot == 1
                subplot(3,1,1)
                plot(t_decel_start,stim_angle(n_decel_start),'+k','linewidth',2)
                subplot(3,1,2)
                plot(t_decel_start,V(n_decel_start),'+k','linewidth',2)
                subplot(3,1,3)
                plot(t_decel_start,At_hor(n_decel_start),'+k','linewidth',2)
            end
            
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

                    if toplot == 1
                        subplot(3,1,1)
                        plot(t_decel_stop,stim_angle(n_decel_stop),'+k','linewidth',2)
                        subplot(3,1,2)
                        plot(t_decel_stop,V(n_decel_stop),'+k','linewidth',2)
                        subplot(3,1,3)
                        plot(t_decel_stop,At_hor(n_decel_stop),'+k','linewidth',2)
                    end
                end
            end

            dt_decel = t_decel_stop - t_decel_start;

            if isnan(dt_decel)==0
                At_hor_min = min(At_hor(n_decel_start:n_decel_stop));
                n_decel_min = find(At_hor==At_hor_min);
                t_decel_min = t(n_decel_min);

                if toplot == 1
                    subplot(3,1,1)
                    plot(t_decel_min,stim_angle(n_decel_min),'+k','linewidth',2)
                    subplot(3,1,2)
                    plot(t_decel_min,V(n_decel_min),'+k','linewidth',2)
                    subplot(3,1,3)
                    plot(t_decel_min,At_hor_min,'+k','linewidth',2)
                end
            else
                At_hor_min = min(At_hor(n_decel_start:end));
                n_decel_min = find(At_hor==At_hor_min);
                t_decel_min = t(n_decel_min);

                if toplot == 1
                    subplot(3,1,1)
                    plot(t_decel_min,stim_angle(n_decel_min),'+k','linewidth',2)
                    subplot(3,1,2)
                    plot(t_decel_min,V(n_decel_min),'+k','linewidth',2)
                    subplot(3,1,3)
                    plot(t_decel_min,At_hor_min,'+k','linewidth',2)
                end
            end
        end
    end
    
    
    %% mean values before and after response
    % before response
    if max(stim_angle(1:n_resp)) - min(stim_angle(1:n_resp)) < 90
        stim_angle_trig2resp = mean(stim_angle(1:n_resp));
    end
    V_trig2resp = mean(V(1:n_resp));
    
    % end of response
    nr_steady = find(IDX_shift==5);
    if length(nr_steady) > 1
        nr_resp_end = nr_steady(2);
        
        t_resp_end = t_shift(nr_resp_end);
        n_resp_end = find(t==t_resp_end);
        stim_angle_post_resp = stim_angle(n_resp_end);
        V_post_resp = V(n_resp_end);

        if toplot == 1
            subplot(3,1,1)
            plot(t_resp_end,stim_angle(n_resp_end),'ok','linewidth',2)
            subplot(3,1,2)
            plot(t_resp_end,V(n_resp_end),'ok','linewidth',2)
            subplot(3,1,3)
            plot(t_resp_end,An_hor(n_resp_end),'ok','linewidth',2)
            plot(t_resp_end,At_hor(n_resp_end),'ok','linewidth',2)
        end
        
        % mean values after response
        if length(IDX_shift) > nr_resp_end
            t_steady_end = t_shift(nr_resp_end+1);
            n_steady_end = find(t==t_steady_end);

            stim_angle_postresp2end = mean(stim_angle(n_resp_end:n_steady_end));
            V_postresp2end = mean(V(n_resp_end:n_steady_end));
            
            if toplot == 1
                subplot(3,1,1)
                plot(t_steady_end,stim_angle(n_steady_end),'ok','linewidth',2)
                subplot(3,1,2)
                plot(t_steady_end,V(n_steady_end),'ok','linewidth',2)
                subplot(3,1,3)
                plot(t_steady_end,An_hor(n_steady_end),'ok','linewidth',2)
                plot(t_steady_end,At_hor(n_steady_end),'ok','linewidth',2)
            end
        else
            stim_angle_postresp2end = nanmean(stim_angle(n_resp_end:end));
            V_postresp2end = nanmean(V(n_resp_end:end));
        end
    end
end


    