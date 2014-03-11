function         [IDX_shift,...
    t_shift,t_resp,t_resp_end,...
    t_turn_start,t_turn_stop,t_turn_max,...
    t_accel_start,t_accel_stop,t_accel_max,...
    t_decel_start,t_decel_stop,t_decel_max,...
    dt_turn,dt_accel,dt_decel,...
    A_hor_max,An_hor_max,At_hor_max,At_hor_min,...
    stim_angle_vel_trigger,stim_angle_vel_pre_resp,stim_angle_vel_trig2resp,stim_angle_vel_pre_turn,...
    stim_angle_vel_post_turn,stim_angle_vel_post_resp,stim_angle_vel_postresp2end,...
    stim_angle_accel_trigger,stim_angle_accel_pre_resp,stim_angle_accel_trig2resp,stim_angle_accel_pre_turn,...
    stim_angle_accel_resp_min,stim_angle_accel_resp_max,stim_angle_accel_resp_mean,stim_angle_accel_resp_std,...
    stim_angle_accel_post_turn,stim_angle_accel_post_resp,stim_angle_accel_postresp2end,...
    stim_angle_yaw_trigger,stim_angle_yaw_pre_resp,stim_angle_yaw_trig2resp,stim_angle_yaw_pre_turn,...
    stim_angle_yaw_post_turn,stim_angle_yaw_post_resp,stim_angle_yaw_postresp2end,...
    slip_trigger,slip_pre_resp,slip_trig2resp,slip_pre_turn,...
    slip_post_turn,slip_post_resp,slip_postresp2end,...
    pitch_trigger,pitch_pre_resp,pitch_trig2resp,pitch_pre_turn,...
    pitch_post_turn,pitch_post_resp,pitch_postresp2end,...
    roll_trigger,roll_pre_resp,roll_trig2resp,roll_pre_turn,...
    roll_post_turn,roll_post_resp,roll_postresp2end,...
    V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
    V_post_resp,V_postresp2end,V_post_accel,...
    teta_pre_resp,teta_post_resp,teta_pre_turn,teta_post_turn,teta_pre_accel,teta_post_accel] =...
    calc_response_data_pitch(pathDB,patternDB,settings,i,toplot,skip,nstart)


IDX_shift = nan;
t_shift = nan;
t_first = nan;

t_resp = nan;
t_resp_end = nan;

t_turn_start = nan;
t_turn_stop = nan;
t_turn_max = nan;

t_accel_start = nan;
t_accel_stop = nan;
t_accel_max = nan; 

t_decel_start = nan;
t_decel_stop = nan;
t_decel_max = nan;

t_steady_end = nan;

n_first = nan;

n_resp = nan;
n_resp_end = nan;

n_turn_start = nan;
n_turn_stop = nan;
n_turn_max = nan;

n_accel_start = nan;
n_accel_stop = nan;
n_accel_max = nan; 

n_decel_start = nan;
n_decel_stop = nan;
n_decel_min = nan;

n_steady_end = nan;

dt_turn = nan;
dt_accel = nan;
dt_decel = nan;

A_hor_max = nan;
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

stim_angle_accel_trigger = nan;
stim_angle_accel_pre_resp = nan;
stim_angle_accel_trig2resp = nan;
stim_angle_accel_pre_turn = nan;

stim_angle_accel_resp_min = nan;
stim_angle_accel_resp_max = nan;
stim_angle_accel_resp_mean = nan;
stim_angle_accel_resp_std = nan;

stim_angle_accel_post_turn = nan;
stim_angle_accel_post_resp = nan;   
stim_angle_accel_postresp2end = nan;   

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

pitch_trigger = nan;
pitch_pre_resp = nan;
pitch_trig2resp = nan;
pitch_pre_turn = nan;

pitch_post_turn = nan;
pitch_post_resp = nan;   
pitch_postresp2end = nan;   

roll_trigger = nan;
roll_pre_resp = nan;
roll_trig2resp = nan;
roll_pre_turn = nan;

roll_post_turn = nan;
roll_post_resp = nan;   
roll_postresp2end = nan;   

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
stim_angle_accel = pathDB.stim_angle_accel(start_frame:end,i,1);
stim_angle_yaw = pathDB.stim_angle_yaw(start_frame:end,i,1);

slip = pathDB.slip(start_frame:end,i,1);
pitch = pathDB.pitch(start_frame:end,i,1);
roll = pathDB.roll(start_frame:end,i,1);

V = pathDB.V(start_frame:end,i,1);
A = pathDB.A(start_frame:end,i,1);
A_hor = pathDB.A_hor(start_frame:end,i,1);

An = pathDB.An(start_frame:end,i);
At = pathDB.At(start_frame:end,i);
An_hor = pathDB.An_hor(start_frame:end,i);
At_hor = pathDB.At_hor(start_frame:end,i);
% alpha_dot_hor = pathDB.alpha_dot_hor(start_frame:end,i);

teta = patternDB.teta(start_frame:end,i);

cmap_k = settings.cmap_k;
% cmap_k_abs = settings.cmap_k_abs;
% cmap_360 = settings.cmap_360;


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
    stim_angle_accel_trigger = stim_angle_accel(n_first);
    stim_angle_yaw_trigger = stim_angle_yaw(n_first);
    
    slip_trigger = slip(n_first);
    pitch_trigger = pitch(n_first);
    roll_trigger = roll(n_first);
    
    V_trigger = V(n_first);
    
    % response time
    t_resp = t_shift(2);
    n_resp = find(t==t_resp);
    IDX_resp = IDX_shift(2);
    
    stim_angle_vel_pre_resp = mean(stim_angle_vel(n_resp));
    stim_angle_accel_pre_resp = mean(stim_angle_accel(n_resp));
    stim_angle_yaw_pre_resp = mean(stim_angle_yaw(n_resp));
    
    slip_pre_resp = mean(slip(n_resp));
    pitch_pre_resp = mean(pitch(n_resp));
    roll_pre_resp = mean(roll(n_resp));
    
    V_pre_resp = V(n_resp);
    teta_pre_resp = teta(n_resp);
    
    %% turn
    if min(IDX_shift)<4 || max(IDX_shift)>6

        % start of turn
        turn_start_nr = min([find(IDX_shift>6);find(IDX_shift<4)]);
        IDX_turn_start = IDX_shift(turn_start_nr);
        t_turn_start = t_shift(turn_start_nr);
        n_turn_start = find(t==t_turn_start);
        
        stim_angle_vel_pre_turn = stim_angle_vel(n_turn_start);
        stim_angle_accel_pre_turn = stim_angle_accel(n_turn_start);
        stim_angle_yaw_pre_turn = stim_angle_yaw(n_turn_start);
        
        slip_pre_turn = slip(n_turn_start);
        pitch_pre_turn = pitch(n_turn_start);
        roll_pre_turn = roll(n_turn_start);
        
        teta_pre_turn = teta(n_turn_start);


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
                    stim_angle_accel_post_turn = nan;
                    stim_angle_yaw_post_turn = nan;
                    
                    slip_post_turn = nan;
                    pitch_post_turn = nan;
                    roll_post_turn = nan;
                    
                    teta_post_turn = nan;
                    turn_stop = 1;
                elseif IDX(n) > 3
                    n_turn_stop = n;
                    t_turn_stop = t(n);
                    
                    stim_angle_vel_post_turn = stim_angle_vel(n);
                    stim_angle_accel_post_turn = stim_angle_accel(n);
                    stim_angle_yaw_post_turn = stim_angle_yaw(n);
                    
                    slip_post_turn = slip(n);
                    pitch_post_turn = pitch(n);
                    roll_post_turn = roll(n);
                    
                    teta_post_turn = teta(n);
                    turn_stop = 1;
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
                    stim_angle_accel_post_turn = nan;
                    stim_angle_yaw_post_turn = nan;
                    
                    slip_post_turn = nan;
                    pitch_post_turn = nan;
                    roll_post_turn = nan;
                    
                    teta_post_turn = nan;
                    turn_stop = 1;
                elseif IDX(n) < 7
                    n_turn_stop = n;
                    t_turn_stop = t(n);
                    
                    stim_angle_vel_post_turn = stim_angle_vel(n);
                    stim_angle_accel_post_turn = stim_angle_accel(n);
                    stim_angle_yaw_post_turn = stim_angle_yaw(n);
                    
                    slip_post_turn = slip(n);
                    pitch_post_turn = pitch(n);
                    roll_post_turn = roll(n);
                    
                    teta_post_turn = teta(n);
                    turn_stop = 1;
                end
            end
            
            % turn time
            dt_turn = t_turn_stop - t_turn_start;

            % max POSITIVE turn accel (ALSO NEG TURN WITH ABS(An_hor))
            if isnan(dt_turn)==0
                An_hor_max = An_hor (abs(An_hor) ==  max(abs(An_hor(n_turn_start:n_turn_stop)))  );
            else
                An_hor_max = An_hor (abs(An_hor) ==  max(abs(An_hor(n_turn_start:end)))  );
            end
            n_turn_max = find(An_hor==An_hor_max);
            t_turn_max = t(n_turn_max);
        end
    end
    
    %% Accelleration
    if isempty(find(IDX_shift==3, 1))==0 || isempty(find(IDX_shift==6, 1))==0 || isempty(find(IDX_shift==9, 1))==0 

        % start of accelleration
        accel_start_nr = min([find(IDX_shift==3); find(IDX_shift==6); find(IDX_shift==9)]);
        t_accel_start = t_shift(accel_start_nr);
        
        % only accel if within same maneuver as turn
        if isnan(t_turn_stop)==1 || t_accel_start < t_turn_stop
            
            n_accel_start = find(t==t_accel_start);
            IDX_accel_start = IDX_shift(accel_start_nr);
            V_pre_accel = V(n_accel_start);
            teta_pre_accel = teta(n_accel_start);

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
                end
            end

            % accel time
            dt_accel = t_accel_stop - t_accel_start;

            %% max accel
            if isnan(dt_accel)==0
                At_hor_max = max(At_hor(n_accel_start:n_accel_stop));
                n_accel_max = find(At_hor==At_hor_max);
                t_accel_max = t(n_accel_max);
            else
                At_hor_max = max(At_hor(n_accel_start:end));
                n_accel_max = find(At_hor==At_hor_max);
                t_accel_max = t(n_accel_max);
            end
        else
            t_accel_start = nan;
        end
        
    end                
    
    %% Deceleration
    if isempty(find(IDX_shift==1, 1))==0 || isempty(find(IDX_shift==4, 1))==0 || isempty(find(IDX_shift==7, 1))==0 

        % start of deceleration
        decel_start_nr = min([find(IDX_shift==1); find(IDX_shift==4); find(IDX_shift==7)]);
        
        if isnan(t_accel_start)==1 || decel_start_nr < accel_start_nr % decel before accel
            t_decel_start = t_shift(decel_start_nr);
            if isnan(t_turn_stop)==1 || t_decel_start < t_turn_stop % decel before end of turn

                IDX_decel_start = IDX_shift(decel_start_nr);
                n_decel_start = find(t==t_decel_start);
                V_pre_decel = V(n_decel_start);
                teta_pre_decel = teta(n_decel_start);

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
                    end
                end

                dt_decel = t_decel_stop - t_decel_start;

                if isnan(dt_decel)==0
                    At_hor_min = min(At_hor(n_decel_start:n_decel_stop));
                    n_decel_min = find(At_hor==At_hor_min);
                    t_decel_min = t(n_decel_min);
                else
                    At_hor_min = min(At_hor(n_decel_start:end));
                    n_decel_min = find(At_hor==At_hor_min);
                    t_decel_min = t(n_decel_min);
                end
            else
                t_decel_start = nan;
            end
        end
    end
    
    
    %% mean values before, after and during response
    % before response
    if max(stim_angle_vel(n_first:n_resp)) - min(stim_angle_vel(n_first:n_resp)) < 90
        stim_angle_vel_trig2resp = mean(stim_angle_vel(n_first:n_resp));
    end
    if max(stim_angle_accel(n_first:n_resp)) - min(stim_angle_accel(n_first:n_resp)) < 90
        stim_angle_accel_trig2resp = mean(stim_angle_accel(n_first:n_resp));
    end
    if max(stim_angle_yaw(n_first:n_resp)) - min(stim_angle_yaw(n_first:n_resp)) < 90
        stim_angle_yaw_trig2resp = mean(stim_angle_yaw(n_first:n_resp));
    end
    if max(slip(n_first:n_resp)) - min(slip(n_first:n_resp)) < 90
        slip_trig2resp = mean(slip(n_first:n_resp));
    end
    if max(pitch(n_first:n_resp)) - min(pitch(n_first:n_resp)) < 90
        pitch_trig2resp = mean(pitch(n_first:n_resp));
    end
    if max(roll(n_first:n_resp)) - min(roll(n_first:n_resp)) < 90
        roll_trig2resp = mean(roll(n_first:n_resp));
    end
    V_trig2resp = mean(V(n_first:n_resp));
    
    % end of response
    nr_steady = find(IDX_shift==5);
    if length(nr_steady) > 1
        nr_resp_end = nr_steady(2);
        
        t_resp_end = t_shift(nr_resp_end);
        n_resp_end = find(t==t_resp_end);
        
        stim_angle_vel_post_resp = stim_angle_vel(n_resp_end);
        stim_angle_accel_post_resp = stim_angle_accel(n_resp_end);
        stim_angle_yaw_post_resp = stim_angle_yaw(n_resp_end);
        
        slip_post_resp = slip(n_resp_end);
        pitch_post_resp = pitch(n_resp_end);
        roll_post_resp = roll(n_resp_end);
        
        V_post_resp = V(n_resp_end);
        teta_post_resp = teta(n_resp_end);
        
        % mean values after response
        if length(IDX_shift) > nr_resp_end
            t_steady_end = t_shift(nr_resp_end+1);
            n_steady_end = find(t==t_steady_end);

            stim_angle_vel_postresp2end = mean(stim_angle_vel(n_resp_end:n_steady_end));
            stim_angle_accel_postresp2end = mean(stim_angle_accel(n_resp_end:n_steady_end));
            stim_angle_yaw_postresp2end = mean(stim_angle_yaw(n_resp_end:n_steady_end));
            
            slip_postresp2end = mean(slip(n_resp_end:n_steady_end));
            pitch_postresp2end = mean(pitch(n_resp_end:n_steady_end));
            roll_postresp2end = mean(roll(n_resp_end:n_steady_end));
            
            V_postresp2end = mean(V(n_resp_end:n_steady_end));
        else
            stim_angle_vel_postresp2end = nanmean(stim_angle_vel(n_resp_end:end));
            stim_angle_accel_postresp2end = nanmean(stim_angle_accel(n_resp_end:end));
            stim_angle_yaw_postresp2end = nanmean(stim_angle_yaw(n_resp_end:end));
            
            slip_postresp2end = nanmean(slip(n_resp_end:end));
            pitch_postresp2end = nanmean(pitch(n_resp_end:end));
            roll_postresp2end = nanmean(roll(n_resp_end:end));
            
            V_postresp2end = nanmean(V(n_resp_end:end));
        end
        
        % mean values during response
        A_hor_max = max(A_hor(n_resp:n_resp_end));
        
        stim_angle_accel_resp_min = min(stim_angle_accel(n_resp:n_resp_end));
        stim_angle_accel_resp_max = max(stim_angle_accel(n_resp:n_resp_end));
        
        stim_angle_accel_resp_mean = circ_mean_deg_nonan(stim_angle_accel(n_resp:n_resp_end));
%         stim_angle_accel_resp_mean = nanmean(stim_angle_accel(n_resp:n_resp_end));
        
        [~, stim_angle_accel_resp_std] = circ_std_deg_nonan(stim_angle_accel(n_resp:n_resp_end));
%         stim_angle_accel_resp_std = nanstd(stim_angle_accel(n_resp:n_resp_end));
    else
        % mean values from start of response until end of seq
        A_hor_max = max(A_hor(n_resp:end));
        
        stim_angle_accel_resp_min = min(stim_angle_accel(n_resp:end));
        stim_angle_accel_resp_max = max(stim_angle_accel(n_resp:end));
        
        stim_angle_accel_resp_mean = circ_mean_deg_nonan(stim_angle_accel(n_resp:end));
%         stim_angle_accel_resp_mean = nanmean(stim_angle_accel(n_resp:end));
        
        [~, stim_angle_accel_resp_std] = circ_std_deg_nonan(stim_angle_accel(n_resp:end));
%         stim_angle_accel_resp_std = nanstd(stim_angle_accel(n_resp:end));
    end
end

if toplot == 1
    plot_flighttracks_clusters_separate
end

    
