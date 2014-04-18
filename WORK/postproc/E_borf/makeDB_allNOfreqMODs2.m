% allNOfreq MODs data

file_name = [save_name,'_all_const_freq_'];
list = dir([file_name,'*mat']);

% % init vars with nans
% t_allNOfreq = nan(seq_length_max,length(list));
% Fx_allNOfreq = nan(seq_length_max,length(list));
% Fy_allNOfreq = nan(seq_length_max,length(list));
% Fz_allNOfreq = nan(seq_length_max,length(list));
% Mx_allNOfreq = nan(seq_length_max,length(list));
% My_allNOfreq = nan(seq_length_max,length(list));
% Mz_allNOfreq = nan(seq_length_max,length(list));
% stroke_L_allNOfreq = nan(seq_length_max,length(list));
% pitch_L_allNOfreq = nan(seq_length_max,length(list));
% dev_L_allNOfreq = nan(seq_length_max,length(list));
% stroke_R_allNOfreq = nan(seq_length_max,length(list));
% pitch_R_allNOfreq = nan(seq_length_max,length(list));
% dev_R_allNOfreq = nan(seq_length_max,length(list));

t_allNOfreq = [];
Fx_allNOfreq = [];
Fy_allNOfreq = [];
Fz_allNOfreq = [];
Mx_allNOfreq = [];
My_allNOfreq = [];
Mz_allNOfreq = [];
stroke_L_allNOfreq = [];
pitch_L_allNOfreq = [];
dev_L_allNOfreq = [];
stroke_R_allNOfreq = [];
pitch_R_allNOfreq = [];
dev_R_allNOfreq = [];

n=0;

% load data
for i = 1:2:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
    %     load(list(i).name)
    
        n=n+1;
        run_type_allNOfreq{n} = exp_type;
        mod_value_allNOfreq(n,1) = mod_val;

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
        freq_allNOfreq(n,1) = f_now;

        Fx_mean_allNOfreq(n,1) = nanmean(Fx_now(n_start:n_stop));
        Fy_mean_allNOfreq(n,1) = nanmean(Fy_now(n_start:n_stop));
        Fz_mean_allNOfreq(n,1) = nanmean(Fz_now(n_start:n_stop));

        Mx_mean_allNOfreq(n,1) = nanmean(Mx_now(n_start:n_stop));
        My_mean_allNOfreq(n,1) = nanmean(My_now(n_start:n_stop));
        Mz_mean_allNOfreq(n,1) = nanmean(Mz_now(n_start:n_stop));

        Astroke_allNOfreq(n,1) = max([stroke_R_now;stroke_L_now]) - min([stroke_R_now;stroke_L_now]);
        Apitch_allNOfreq(n,1) = max([pitch_R_now;pitch_L_now]) - min([pitch_R_now;pitch_L_now]);
        Adev_allNOfreq(n,1) = max([dev_R_now;dev_L_now]) - min([dev_R_now;dev_L_now]);

        %% store allNOfreqNOfreq data
        t_allNOfreq(1:length(t),n) = t_now;

        Fx_allNOfreq(1:length(t),n) = Fx_now;
        Fy_allNOfreq(1:length(t),n) = Fy_now;
        Fz_allNOfreq(1:length(t),n) = Fz_now;

        Mx_allNOfreq(1:length(t),n) = Mx_now;
        My_allNOfreq(1:length(t),n) = My_now;
        Mz_allNOfreq(1:length(t),n) = Mz_now;

        stroke_L_allNOfreq(1:length(t),n) = stroke_L_now;
        pitch_L_allNOfreq(1:length(t),n) = pitch_L_now;
        dev_L_allNOfreq(1:length(t),n) = dev_L_now;

        stroke_R_allNOfreq(1:length(t),n) = stroke_R_now;
        pitch_R_allNOfreq(1:length(t),n) = pitch_R_now;
        dev_R_allNOfreq(1:length(t),n) = dev_R_now;

    end
end
F_mean_allNOfreq = sqrt(Fx_mean_allNOfreq.^2+Fy_mean_allNOfreq.^2+Fz_mean_allNOfreq.^2);

if length(t_allNOfreq) > seq_length_max
    seq_length_max
    'TOO SHORT'
end




















