% calc variables at start & stop of response

%% nstart & nstop
n_turn_start = nan(size(responseDB.t_turn_start));
n_turn_stop = nan(size(responseDB.t_turn_start));
n_accel_start = nan(size(responseDB.t_turn_start));
n_accel_stop = nan(size(responseDB.t_turn_start));
n_decel_start = nan(size(responseDB.t_turn_start));
n_decel_stop = nan(size(responseDB.t_turn_start));

for i = 1:length(responseDB.t_turn_start)
    
    if isnan(responseDB.t_turn_start(i))==0
        n_turn_start(i,1) = find(t==responseDB.t_turn_start(i));
    end
    
    if isnan(responseDB.t_turn_stop(i))==0
        n_turn_stop(i,1) = find(t==responseDB.t_turn_stop(i));
    end
    
    if isnan(responseDB.t_accel_start(i))==0
        n_accel_start(i,1) = find(t==responseDB.t_accel_start(i));
    end
    
    if isnan(responseDB.t_accel_stop(i))==0
        n_accel_stop(i,1) = find(t==responseDB.t_accel_stop(i));
    end
    
    if isnan(responseDB.t_decel_start(i))==0
        n_decel_start(i,1) = find(t==responseDB.t_decel_start(i));
    end
    
    if isnan(responseDB.t_decel_stop(i))==0
        n_decel_stop(i,1) = find(t==responseDB.t_decel_stop(i));
    end
end
n_pre = min([n_turn_start n_accel_start n_decel_start]')';
n_post = max([n_turn_stop n_accel_stop n_decel_stop]')'-1;

%% tstart & tstop
t_pre = nan(size(n_pre));
t_post = nan(size(n_pre));
for i = 1:length(n_pre)
    if isnan(n_pre(i))==0
        t_pre(i,1) = t(n_pre(i));
    end
        
    if isnan(n_post(i))==0
        t_post(i,1) = t(n_post(i));
    end
end


