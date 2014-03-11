function     [IDX_shift,t_shift,t_resp,...
    t_turn_start,t_turn_stop,t_turn_max,...
    t_accel_start,t_accel_stop,t_accel_max,...
    t_decel_start,t_decel_stop,t_decel_max,...
    dt_turn,dt_accel,dt_decel,...
    An_hor_max,At_hor_max,At_hor_min,...
    stim_angle_vel_trigger,stim_angle_vel_pre_resp,stim_angle_vel_trig2resp,stim_angle_vel_pre_turn,...
    stim_angle_vel_post_turn,stim_angle_vel_post_resp,stim_angle_vel_postresp2end,...
    stim_angle_yaw_trigger,stim_angle_yaw_pre_resp,stim_angle_yaw_trig2resp,stim_angle_yaw_pre_turn,...
    stim_angle_yaw_post_turn,stim_angle_yaw_post_resp,stim_angle_yaw_postresp2end,...
    slip_trigger,slip_pre_resp,slip_trig2resp,slip_pre_turn,...
    slip_post_turn,slip_post_resp,slip_postresp2end,...
    V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
    V_post_resp,V_postresp2end,V_post_accel,...
    teta_pre_resp,teta_post_resp,teta_pre_turn,teta_post_turn,teta_pre_accel,teta_post_accel] =...
    calc_response_data_hVA_2plots(pathDB,patternDB,settings,i,toplot,skip,nstart);

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

stim_angle_vel_trigger = nan;
stim_angle_vel_pre_resp = nan;
stim_angle_vel_trig2resp = nan;
stim_angle_vel_pre_turn = nan;

stim_angle_vel_post_turn = nan;
stim_angle_vel_post_resp = nan;   
stim_angle_vel_postresp2end = nan;   

stim_angle_yaw_trigger = nan;
stim_angle_yaw_pre_resp = nan;
stim_angle_yaw_trig2resp = nan;
stim_angle_yaw_pre_turn = nan;

stim_angle_yaw_post_turn = nan;
stim_angle_yaw_post_resp = nan;   
stim_angle_yaw_postresp2end = nan;   

slip_trigger = nan;
slip_pre_resp = nan;
slip_trig2resp = nan;
slip_pre_turn = nan;

slip_post_turn = nan;
slip_post_resp = nan;   
slip_postresp2end = nan;   

V_trigger = nan;   
V_pre_resp = nan;   
V_trig2resp = nan;   
V_pre_accel = nan;   

V_post_resp = nan;   
V_postresp2end = nan;   
V_post_accel = nan;
    
teta_pre_resp = nan;   
teta_pre_turn = nan;
teta_pre_accel = nan;   

teta_post_resp = nan;   
teta_post_turn = nan;
teta_post_accel = nan;


t = pathDB.t;
dt = t(2) - t(1);

trigger_frame = find(t == min(abs(t)));
start_frame = trigger_frame + nstart;

t = t(start_frame:end,:);

x = pathDB.pos(start_frame:end,i,1);
y = pathDB.pos(start_frame:end,i,2);
z = pathDB.pos(start_frame:end,i,3);

u = pathDB.vel(start_frame:end,i,1);
v = pathDB.vel(start_frame:end,i,2);
w = pathDB.vel(start_frame:end,i,3);

ax = pathDB.accel(start_frame:end,i,1);
ay = pathDB.accel(start_frame:end,i,2);
az = pathDB.accel(start_frame:end,i,3);

IDX = pathDB.IDX(start_frame:end,i,1);
stim_angle_vel = pathDB.stim_angle_vel(start_frame:end,i,1);
stim_angle_yaw = pathDB.stim_angle_yaw(start_frame:end,i,1);
heading = pathDB.heading(start_frame:end,i,1);
slip = pathDB.slip(start_frame:end,i,1);

V = pathDB.V(start_frame:end,i,1);
A = pathDB.A(start_frame:end,i,1);
A_hor = pathDB.A_hor(start_frame:end,i,1);

An = pathDB.An(start_frame:end,i);
At = pathDB.At(start_frame:end,i);
An_hor = pathDB.An_hor(start_frame:end,i);
At_hor = pathDB.At_hor(start_frame:end,i);
alpha_dot_hor = pathDB.alpha_dot_hor(start_frame:end,i);

teta = patternDB.teta(start_frame:end,i);

cmap_k = settings.cmap_k;
% cmap_k_abs = settings.cmap_k_abs;
% cmap_360 = settings.cmap_360;

% plot timelines
if toplot == 1
    figure(1)
    subplot(3,1,1)
    hold off
    plot(t,teta,'-k')
    hold on
    axis([0 .13 0 180])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',0:90:180)
    xlabel('time')
    ylabel('optical angle')
    grid on
    
    subplot(3,1,2)
    hold off
    plot(t,stim_angle_vel,'-k')
    hold on
    axis([0 .13 -180 180])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',-180:180:180)
    xlabel('time')
    ylabel('heading angle')
    grid on

    subplot(3,1,3)
    hold off
    plot(t,stim_angle_yaw,'-k')
    hold on
    axis([0 .13 -180 180])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',-180:180:180)
    xlabel('time')
    ylabel('yaw angle')
    grid on

    figure(2)    
    subplot(3,1,1)
    hold off
    plot(t,V,'-k')
    hold on
    axis([0 .13 0 1])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',0:.5:1)
    xlabel('time')
    ylabel('V')
    grid on

    subplot(3,1,2)
    hold off
    plot(t,abs(An_hor),'-k')
    hold on
    axis([0 .13 0 15])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',0:5:15)
    xlabel('time')
    ylabel('An')
    grid on

    subplot(3,1,3)
    hold off
    plot(t,At_hor,'-k')
    hold on
    axis([0 .13 -5 15])
    set(gca,'XTick',0:.025:.125)
    set(gca,'YTick',[-5;0;15])
    xlabel('time')
    ylabel('At')
    grid on
    
    for j = 1:skip:length(IDX)
        length(IDX)-j;
        if isnan(IDX(j))==0
            figure(1)
            subplot(3,1,1)
            plot(t(j),teta(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            
            subplot(3,1,2)
            plot(t(j),stim_angle_vel(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)
            
            subplot(3,1,3)
            plot(t(j),stim_angle_yaw(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',5)

            figure(2)
            subplot(3,1,1)
            plot(t(j),V(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)

            subplot(3,1,2)
            plot(t(j),abs(An_hor(j)),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)
            
            subplot(3,1,3)
            plot(t(j),At_hor(j),'.','color',cmap_k(IDX(j),:),'MarkerSize',10)
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

% remove 1st nan
if isnan(IDX_shift(1)) == 1
    IDX_shift(1) = [];
    t_shift(1) = [];
end

if IDX_shift(1) == 5 && length(IDX_shift) > 1
    
    % @ trigger OR 1st frame
    t_first = t_shift(1);
    n_first = find(t==t_first);
    stim_angle_vel_trigger = stim_angle_vel(n_first);
    stim_angle_yaw_trigger = stim_angle_yaw(n_first);
    slip_trigger = slip(n_first);
    V_trigger = V(n_first);

    if toplot == 1
        figure(1)
        subplot(3,1,1)
        plot(t_first,teta(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,2)
        plot(t_first,stim_angle_vel(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,3)
        plot(t_first,stim_angle_yaw(n_first),'ok','linewidth',2,'MarkerSize',5)
        
        figure(2)
        subplot(3,1,1)
        plot(t_first,V(n_first),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,2)
        plot(t_first,abs(An_hor(n_first)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,3)
        plot(t_first,At_hor(n_first),'ok','linewidth',2,'MarkerSize',5)
    end
    
    % response time
    t_resp = t_shift(2);
    n_resp = find(t==t_resp);
    IDX_resp = IDX_shift(2);
    stim_angle_vel_pre_resp = mean(stim_angle_vel(n_resp));
    stim_angle_yaw_pre_resp = mean(stim_angle_yaw(n_resp));
    slip_pre_resp = mean(slip(n_resp));
    V_pre_resp = V(n_resp);
    teta_pre_resp = teta(n_resp);
    
    if toplot == 1
        figure(1)
        subplot(3,1,1)
        plot(t_resp,teta(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,2)
        plot(t_resp,stim_angle_vel(n_resp),'ok','linewidth',2,'MarkerSize',5)
         subplot(3,1,3)
       plot(t_resp,stim_angle_yaw(n_resp),'ok','linewidth',2,'MarkerSize',5)
        figure(2)
        subplot(3,1,1)
        plot(t_resp,V(n_resp),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,2)
        plot(t_resp,abs(An_hor(n_resp)),'ok','linewidth',2,'MarkerSize',5)
        subplot(3,1,3)
        plot(t_resp,At_hor(n_resp),'ok','linewidth',2,'MarkerSize',5)
    end
    
    
    %% turn
    if min(IDX_shift)<4 || max(IDX_shift)>6

        % start of turn
        turn_start_nr = min([find(IDX_shift>6);find(IDX_shift<4)]);
        IDX_turn_start = IDX_shift(turn_start_nr);
        t_turn_start = t_shift(turn_start_nr);
        n_turn_start = find(t==t_turn_start);
        stim_angle_vel_pre_turn = stim_angle_vel(n_turn_start);
        stim_angle_yaw_pre_turn = stim_angle_yaw(n_turn_start);
        slip_pre_turn = slip(n_turn_start);
        teta_pre_turn = teta(n_turn_start);

        if toplot == 1
            figure(1)
            subplot(3,1,1)
            plot(t_turn_start,teta(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
            subplot(3,1,2)
            plot(t_turn_start,stim_angle_vel(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
            subplot(3,1,3)
            plot(t_turn_start,stim_angle_yaw(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
            figure(2)
            subplot(3,1,1)
            plot(t_turn_start,V(n_turn_start),'xk','linewidth',2,'MarkerSize',5)
            subplot(3,1,2)
            plot(t_turn_start,abs(An_hor(n_turn_start)),'xk','linewidth',2,'MarkerSize',5)
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
                    stim_angle_vel_post_turn = nan;
                    stim_angle_yaw_post_turn = nan;
                    slip_post_turn = nan;
                    teta_post_turn = nan;
                    turn_stop = 1;
                elseif IDX(n) > 3
                    n_turn_stop = n;
                    t_turn_stop = t(n);
                    stim_angle_vel_post_turn = stim_angle_vel(n);
                    stim_angle_yaw_post_turn = stim_angle_yaw(n);
                    slip_post_turn = slip(n);
                    teta_post_turn = teta(n);
                    turn_stop = 1;
                    
                    if toplot == 1
                        figure(1)
                        subplot(3,1,1)
                        plot(t_turn_stop,teta(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,2)
                        plot(t_turn_stop,stim_angle_vel(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,3)
                        plot(t_turn_stop,stim_angle_yaw(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        figure(2)
                        subplot(3,1,1)
                        plot(t_turn_stop,V(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,2)
                        plot(t_turn_stop,abs(An_hor(n_turn_stop)),'xk','linewidth',2,'MarkerSize',5)
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
                figure(1)
                subplot(3,1,1)
                plot(t_turn_max,teta(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_turn_max,stim_angle_vel(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                 subplot(3,1,3)
               plot(t_turn_max,stim_angle_yaw(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_turn_max,V(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_turn_max,abs(An_hor_max),'xk','linewidth',2,'MarkerSize',5)
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
                    stim_angle_vel_post_turn = nan;
                    stim_angle_yaw_post_turn = nan;
                    slip_post_turn = nan;
                    teta_post_turn = nan;
                    turn_stop = 1;
                elseif IDX(n) < 7
                    n_turn_stop = n;
                    t_turn_stop = t(n);
                    stim_angle_vel_post_turn = stim_angle_vel(n);
                    stim_angle_yaw_post_turn = stim_angle_yaw(n);
                    slip_post_turn = slip(n);
                    teta_post_turn = teta(n);
                    turn_stop = 1;


                    if toplot == 1
                        figure(1)
                        subplot(3,1,1)
                        plot(t_turn_stop,teta(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,2)
                        plot(t_turn_stop,stim_angle_vel(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,3)
                        plot(t_turn_stop,stim_angle_yaw(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        figure(2)
                        subplot(3,1,1)
                        plot(t_turn_stop,V(n_turn_stop),'xk','linewidth',2,'MarkerSize',5)
                        subplot(3,1,2)
                        plot(t_turn_stop,abs(An_hor(n_turn_stop)),'xk','linewidth',2,'MarkerSize',5)
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
                figure(1)
                subplot(3,1,1)
                plot(t_turn_max,teta(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                 subplot(3,1,2)
               plot(t_turn_max,stim_angle_vel(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_turn_max,stim_angle_yaw(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_turn_max,V(n_turn_max),'xk','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_turn_max,abs(An_hor_max),'xk','linewidth',2,'MarkerSize',5)
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
        teta_pre_accel = teta(n_accel_start);

        if toplot == 1
            figure(1)
            subplot(3,1,1)
            plot(t_accel_start,teta(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
            subplot(3,1,2)
            plot(t_accel_start,stim_angle_vel(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
            subplot(3,1,3)
            plot(t_accel_start,stim_angle_yaw(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
            figure(2)
            subplot(3,1,1)
            plot(t_accel_start,V(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
            subplot(3,1,3)
            plot(t_accel_start,At_hor(n_accel_start),'+k','linewidth',2,'MarkerSize',5)
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
                teta_post_accel = teta(n);
                accel_stop = 1;

                if toplot == 1
                    figure(1)
                    subplot(3,1,1)
                    plot(t_accel_stop,teta(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
                     subplot(3,1,2)
                   plot(t_accel_stop,stim_angle_vel(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_accel_stop,stim_angle_yaw(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
                    figure(2)
                    subplot(3,1,1)
                    plot(t_accel_stop,V(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_accel_stop,At_hor(n_accel_stop),'+k','linewidth',2,'MarkerSize',5)
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
                figure(1)
                subplot(3,1,1)
                plot(t_accel_max,teta(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                 subplot(3,1,2)
               plot(t_accel_max,stim_angle_vel(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_accel_max,stim_angle_yaw(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_accel_max,V(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_accel_max,At_hor_max,'+k','linewidth',2,'MarkerSize',5)
            end
        else
            At_hor_max = max(At_hor(n_accel_start:end));
            n_accel_max = find(At_hor==At_hor_max);
            t_accel_max = t(n_accel_max);

            if toplot == 1
                figure(1)
                subplot(3,1,1)
                plot(t_accel_max,teta(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_accel_max,stim_angle_vel(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_accel_max,stim_angle_yaw(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_accel_max,V(n_accel_max),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_accel_max,At_hor_max,'+k','linewidth',2,'MarkerSize',5)
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
            teta_pre_decel = teta(n_decel_start);

            if toplot == 1
                figure(1)
                subplot(3,1,1)
                plot(t_decel_start,teta(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_decel_start,stim_angle_vel(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_decel_start,stim_angle_yaw(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_decel_start,V(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_decel_start,At_hor(n_decel_start),'+k','linewidth',2,'MarkerSize',5)
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
                    teta_post_decel = nan;
                    decel_stop = 1;
                elseif IDX(n) == 1 || IDX(n) == 4 || IDX(n) == 7
                else
                    t_decel_stop = t(n);
                    n_decel_stop = n;
                    V_post_decel = V(n);
                    teta_post_decel = teta(n);
                    decel_stop = 1;

                    if toplot == 1
                        figure(1)
                        subplot(3,1,1)
                        plot(t_decel_stop,teta(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
                        subplot(3,1,2)
                        plot(t_decel_stop,stim_angle_vel(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
                        subplot(3,1,3)
                        plot(t_decel_stop,stim_angle_yaw(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
                        figure(2)
                        subplot(3,1,1)
                        plot(t_decel_stop,V(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
                        subplot(3,1,3)
                        plot(t_decel_stop,At_hor(n_decel_stop),'+k','linewidth',2,'MarkerSize',5)
                    end
                end
            end

            dt_decel = t_decel_stop - t_decel_start;

            if isnan(dt_decel)==0
                At_hor_min = min(At_hor(n_decel_start:n_decel_stop));
                n_decel_min = find(At_hor==At_hor_min);
                t_decel_min = t(n_decel_min);

                if toplot == 1
                    figure(1)
                    subplot(3,1,1)
                    plot(t_decel_min,teta(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,2)
                    plot(t_decel_min,stim_angle_vel(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_decel_min,stim_angle_yaw(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    figure(2)
                    subplot(3,1,1)
                    plot(t_decel_min,V(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_decel_min,At_hor_min,'+k','linewidth',2,'MarkerSize',5)
                end
            else
                At_hor_min = min(At_hor(n_decel_start:end));
                n_decel_min = find(At_hor==At_hor_min);
                t_decel_min = t(n_decel_min);

                if toplot == 1
                    figure(1)
                    subplot(3,1,1)
                    plot(t_decel_min,teta(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,2)
                    plot(t_decel_min,stim_angle_vel(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_decel_min,stim_angle_yaw(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    figure(2)
                    subplot(3,1,1)
                    plot(t_decel_min,V(n_decel_min),'+k','linewidth',2,'MarkerSize',5)
                    subplot(3,1,3)
                    plot(t_decel_min,At_hor_min,'+k','linewidth',2,'MarkerSize',5)
                end
            end
        end
    end
    
    
    %% mean values before and after response
    % before response
    if max(stim_angle_vel(n_first:n_resp)) - min(stim_angle_vel(n_first:n_resp)) < 90
        stim_angle_vel_trig2resp = mean(stim_angle_vel(n_first:n_resp));
    end
    if max(stim_angle_yaw(n_first:n_resp)) - min(stim_angle_yaw(n_first:n_resp)) < 90
        stim_angle_yaw_trig2resp = mean(stim_angle_yaw(n_first:n_resp));
    end
    if max(slip(n_first:n_resp)) - min(slip(n_first:n_resp)) < 90
        slip_trig2resp = mean(slip(n_first:n_resp));
    end
    V_trig2resp = mean(V(n_first:n_resp));
    
    % end of response
    nr_steady = find(IDX_shift==5);
    if length(nr_steady) > 1
        nr_resp_end = nr_steady(2);
        
        t_resp_end = t_shift(nr_resp_end);
        n_resp_end = find(t==t_resp_end);
        stim_angle_vel_post_resp = stim_angle_vel(n_resp_end);
        stim_angle_yaw_post_resp = stim_angle_yaw(n_resp_end);
        slip_post_resp = slip(n_resp_end);
        V_post_resp = V(n_resp_end);
        teta_post_resp = teta(n_resp_end);

        if toplot == 1
            figure(1)
            subplot(3,1,1)
            plot(t_resp_end,teta(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
            subplot(3,1,2)
            plot(t_resp_end,stim_angle_vel(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
            subplot(3,1,3)
            plot(t_resp_end,stim_angle_yaw(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
            figure(2)
            subplot(3,1,1)
            plot(t_resp_end,V(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
            subplot(3,1,2)
            plot(t_resp_end,abs(An_hor(n_resp_end)),'ok','linewidth',2,'MarkerSize',5)
            subplot(3,1,3)
            plot(t_resp_end,At_hor(n_resp_end),'ok','linewidth',2,'MarkerSize',5)
        end
        
        % mean values after response
        if length(IDX_shift) > nr_resp_end
            t_steady_end = t_shift(nr_resp_end+1);
            n_steady_end = find(t==t_steady_end);

            stim_angle_vel_postresp2end = mean(stim_angle_vel(n_resp_end:n_steady_end));
            stim_angle_yaw_postresp2end = mean(stim_angle_yaw(n_resp_end:n_steady_end));
            slip_postresp2end = mean(slip(n_resp_end:n_steady_end));
            V_postresp2end = mean(V(n_resp_end:n_steady_end));
            
            if toplot == 1
                figure(1)
                subplot(3,1,1)
                plot(t_steady_end,teta(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_steady_end,stim_angle_vel(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_steady_end,stim_angle_yaw(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
                figure(2)
                subplot(3,1,1)
                plot(t_steady_end,V(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
                subplot(3,1,2)
                plot(t_steady_end,abs(An_hor(n_steady_end)),'ok','linewidth',2,'MarkerSize',5)
                subplot(3,1,3)
                plot(t_steady_end,At_hor(n_steady_end),'ok','linewidth',2,'MarkerSize',5)
            end
        else
            stim_angle_vel_postresp2end = nanmean(stim_angle_vel(n_resp_end:end));
            stim_angle_yaw_postresp2end = nanmean(stim_angle_yaw(n_resp_end:end));
            slip_postresp2end = nanmean(slip(n_resp_end:end));
            V_postresp2end = nanmean(V(n_resp_end:end));
        end
    end
end


    
