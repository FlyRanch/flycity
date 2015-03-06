% make pitch MOD DB

file_name = [save_name,'_rot_'];
list = dir([file_name,'*mat']);

% init vars with nans
t_pitch = nan(seq_length_max,length(list));
Fx_pitch = nan(seq_length_max,length(list));
Fy_pitch = nan(seq_length_max,length(list));
Fz_pitch = nan(seq_length_max,length(list));
Mx_pitch = nan(seq_length_max,length(list));
My_pitch = nan(seq_length_max,length(list));
Mz_pitch = nan(seq_length_max,length(list));
stroke_L_pitch = nan(seq_length_max,length(list));
pitch_L_pitch = nan(seq_length_max,length(list));
dev_L_pitch = nan(seq_length_max,length(list));
stroke_R_pitch = nan(seq_length_max,length(list));
pitch_R_pitch = nan(seq_length_max,length(list));
dev_R_pitch = nan(seq_length_max,length(list));

% load data
for i = 1:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    load(file_now)
%     load(list(i).name)
    
    run_type_pitch{i} = exp_type;
    mod_value_pitch(i,1) = mod_val;
    
    %% data NOW
    t_now = t;
    
    % forces in fly ref frame (different from borf)
    Fx_now =  F_robo2fly * ft(:,1);  % thrust, fwd pos
    Fy_now =  F_robo2fly * ft(:,3);  % sideways, right pos
    Fz_now = -F_robo2fly * ft(:,2); % vertical, down pos
    
    Mx_now =  M_robo2fly * ft(:,4);  % Roll, right pos
    My_now =  M_robo2fly * ft(:,6);  % pitch, nose up pos
    Mz_now = -M_robo2fly * ft(:,5); % yaw, right pos
    
    stroke_L_now = rad2deg(kine(:,5));
    pitch_L_now = rad2deg(kine(:,1));
    dev_L_now = rad2deg(kine(:,2));
    
    stroke_R_now = rad2deg(kine(:,6));
    pitch_R_now = rad2deg(kine(:,4));
    dev_R_now = rad2deg(kine(:,3));
    
    %% calc wb mean data
    N_wb_now = round(length(t_now)/Nwb);
    n_start = wb_start * N_wb_now;
    n_stop = wb_stop * N_wb_now;
    
    frobo_now = (1/max(t_now))*Nwb;
    f_now = f_robo2fly * frobo_now;
    freq_pitch(i,1) = f_now;
    
    Fx_mean_pitch(i,1) = nanmean(Fx_now(n_start:n_stop));
    Fy_mean_pitch(i,1) = nanmean(Fy_now(n_start:n_stop));
    Fz_mean_pitch(i,1) = nanmean(Fz_now(n_start:n_stop));
    
    Mx_mean_pitch(i,1) = nanmean(Mx_now(n_start:n_stop));
    My_mean_pitch(i,1) = nanmean(My_now(n_start:n_stop));
    Mz_mean_pitch(i,1) = nanmean(Mz_now(n_start:n_stop));
    
    Astroke_pitch(i,1) = max([stroke_R_now;stroke_L_now]) - min([stroke_R_now;stroke_L_now]);
    Apitch_pitch(i,1) = max([pitch_R_now;pitch_L_now]) - min([pitch_R_now;pitch_L_now]);
    Adev_pitch(i,1) = max([dev_R_now;dev_L_now]) - min([dev_R_now;dev_L_now]);
    
    %% store pitch data
    t_pitch(1:length(t),i) = t_now;
    
    Fx_pitch(1:length(t),i) = Fx_now;
    Fy_pitch(1:length(t),i) = Fy_now;
    Fz_pitch(1:length(t),i) = Fz_now;
    
    Mx_pitch(1:length(t),i) = Mx_now;
    My_pitch(1:length(t),i) = My_now;
    Mz_pitch(1:length(t),i) = Mz_now;
    
    stroke_L_pitch(1:length(t),i) = stroke_L_now;
    pitch_L_pitch(1:length(t),i) = pitch_L_now;
    dev_L_pitch(1:length(t),i) = dev_L_now;
    
    stroke_R_pitch(1:length(t),i) = stroke_R_now;
    pitch_R_pitch(1:length(t),i) = pitch_R_now;
    dev_R_pitch(1:length(t),i) = dev_R_now;
    
    
end
F_mean_pitch = sqrt(Fx_mean_pitch.^2+Fy_mean_pitch.^2+Fz_mean_pitch.^2);

if length(t_pitch) > seq_length_max
    seq_length_max
    'TOO SHORT'
end



















    