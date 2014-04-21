IDX_sub = pathDB.IDX(start_frame:end,i,1);

IDX_shift = nan;
n_shift = nan;

n_first = nan;
n_resp = nan;
n_resp_end = nan;
n_steady_end = nan;

n_turn_start = nan;
n_turn_stop = nan;
n_turn_max = nan;

n_accel_start = nan;
n_accel_stop = nan;
n_accel_max = nan; 

n_decel_start = nan;
n_decel_stop = nan;
n_decel_min = nan;

turn_start_nr = nan;
turn_stop_nr = nan;

accel_start_nr = nan;
accel_stop_nr = nan;

decel_start_nr = nan;
decel_stop_nr = nan;

% cluster shifts
c=1;
n_shift(c,1) = start_frame;
IDX_shift(c,1) = IDX_sub(1);
for j = 2:length(IDX_sub)
    if isnan(IDX_sub(j)) == 0 && IDX_sub(j) ~= IDX_sub(j-1)
        c=c+1;
        n_shift(c,1) = start_frame + j-1;
        IDX_shift(c,1) = IDX_sub(j);
    end
end

% remove 1st nan
if isnan(IDX_shift(1)) == 1
    IDX_shift(1) = [];
    n_shift(1) = [];
end

if IDX_shift(1) == 5 && length(IDX_shift) > 1
    
    % @ trigger OR 1st frame
    n_first = n_shift(1);
    n_resp = n_shift(2);

    %% turn
    if min(IDX_shift)<4 || max(IDX_shift)>6

        % start of turn
        turn_start_nr = min([find(IDX_shift>6);find(IDX_shift<4)]);
        n_turn_start = n_shift(turn_start_nr);
        IDX_turn_start = IDX_shift(turn_start_nr);
        
        % end of turn
        % left turn
        if IDX_turn_start < 4
            turn_stop = 0;
            j = turn_start_nr;
            while turn_stop == 0
                j = j+1;
                if j > length(IDX_shift)
                    n_turn_stop = nan;
                    turn_stop = 1;
                elseif IDX_shift(j) > 3
                    turn_stop_nr = j;
                    n_turn_stop = n_shift(j);
                    turn_stop = 1;
                end
            end

        % right turn
        elseif IDX_turn_start > 6
            turn_stop = 0;
            j = turn_start_nr;
            while turn_stop == 0
                j = j+1;
                if j > length(IDX_shift)
                    n_turn_stop = nan;
                    turn_stop = 1;
                elseif IDX_shift(j) < 7
                    turn_stop_nr = j;
                    n_turn_stop = n_shift(j);
                    turn_stop = 1;
                end
            end
        end

        % max turn accel
        if isnan(n_turn_stop)==0
            n_turn_max = find(abs(An_hor) == max(abs(An_hor(n_turn_start:n_turn_stop))));
        else
            n_turn_max = find(abs(An_hor) == max(abs(An_hor(n_turn_start:end))));
        end
    end
    
    %% Accelleration
    if isempty(find(IDX_shift==3, 1))==0 || isempty(find(IDX_shift==6, 1))==0 || isempty(find(IDX_shift==9, 1))==0 

        % start of accelleration
        accel_start_nr = min([find(IDX_shift==3); find(IDX_shift==6); find(IDX_shift==9)]);
        
        % only accel if within same maneuver as turn
        if isnan(n_turn_stop)==1 || accel_start_nr < turn_stop_nr
            n_accel_start = n_shift(accel_start_nr);
            IDX_accel_start = IDX_shift(accel_start_nr);

            % end of acceleration
            accel_stop = 0;
            j = accel_start_nr;
            while accel_stop == 0
                j = j+1;
                if j > length(IDX_shift)
                    IDX_accel_stop = nan;
                    n_accel_stop = nan;
                    accel_stop = 1;
                elseif IDX_shift(j) == 3 || IDX_shift(j) == 6 || IDX_shift(j) == 9
                else
                    IDX_accel_stop = IDX_shift(j);
                    n_accel_stop = n_shift(j);
                    accel_stop_nr = j;
                    accel_stop = 1;
                end
            end

            %% max accel
            if isnan(n_accel_stop)==0
                At_hor_max = max(At_hor(n_accel_start:n_accel_stop));
                n_accel_max = find(At_hor==At_hor_max);
            else
                At_hor_max = max(At_hor(n_accel_start:end));
                n_accel_max = find(At_hor==At_hor_max);
            end
        end
        
    end                
    
    %% Deceleration
    if isempty(find(IDX_shift==1, 1))==0 || isempty(find(IDX_shift==4, 1))==0 || isempty(find(IDX_shift==7, 1))==0 

        % start of deceleration
        decel_start_nr = min([find(IDX_shift==1); find(IDX_shift==4); find(IDX_shift==7)]);
        
        if isnan(n_accel_start)==1 || decel_start_nr < accel_start_nr % decel before accel
            if isnan(n_turn_stop)==1 || decel_start_nr < turn_stop_nr % decel before end of turn
                
                n_decel_start = n_shift(decel_start_nr);
                IDX_decel_start = IDX_shift(decel_start_nr);

                % end of deceleration
                decel_stop = 0;
                j = decel_start_nr;
                while decel_stop == 0
                    j = j+1;
                    if j > length(IDX_shift)
                        IDX_decel_stop = nan;
                        n_decel_stop = nan;
                        decel_stop = 1;
                    elseif IDX_shift(j) == 1 || IDX_shift(j) == 4 || IDX_shift(j) == 7
                    else
                        IDX_decel_stop = IDX_shift(j);
                        n_decel_stop = n_shift(j);
                        decel_stop_nr = j;
                        decel_stop = 1;
                    end
                end

                if isnan(n_decel_stop)==0
                    At_hor_min = min(At_hor(n_decel_start:n_decel_stop));
                    n_decel_min = find(At_hor==At_hor_min);
                else
                    At_hor_min = min(At_hor(n_decel_start:end));
                    n_decel_min = find(At_hor==At_hor_min);
                end
            end
        end
    end
    
    % end of response
    nr_steady = find(IDX_shift==5);
    if length(nr_steady) > 1
        nr_resp_end = nr_steady(2);
        n_resp_end = n_shift(nr_resp_end);
        
        % end of steady seq after resp
        if length(IDX_shift) > nr_resp_end
            n_steady_end = n_shift(nr_resp_end+1);
        else
            n_steady_end = find(isnan(An_hor)==0, 1, 'last' );
        end
    end
    
end
