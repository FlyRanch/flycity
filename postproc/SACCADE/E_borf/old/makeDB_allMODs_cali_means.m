% all MODs data

file_name = [save_name,'_all_'];
list = dir([file_name,'*mat']);

% % init vars with nans
% t_all = nan(seq_length_max,length(list));
% Fx_all = nan(seq_length_max,length(list));
% Fy_all = nan(seq_length_max,length(list));
% Fz_all = nan(seq_length_max,length(list));
% Mx_all = nan(seq_length_max,length(list));
% My_all = nan(seq_length_max,length(list));
% Mz_all = nan(seq_length_max,length(list));
% stroke_L_all = nan(seq_length_max,length(list));
% pitch_L_all = nan(seq_length_max,length(list));
% dev_L_all = nan(seq_length_max,length(list));
% stroke_R_all = nan(seq_length_max,length(list));
% pitch_R_all = nan(seq_length_max,length(list));
% dev_R_all = nan(seq_length_max,length(list));

t_all = [];
Fx_all = [];
Fy_all = [];
Fz_all = [];
Mx_all = [];
My_all = [];
Mz_all = [];
stroke_L_all = [];
pitch_L_all = [];
dev_L_all = [];
stroke_R_all = [];
pitch_R_all = [];
dev_R_all = [];

n=0;

% load data
for i = 1:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
        %     load(list(i).name)
        
        n=n+1;
        run_type_all{n} = exp_type;
        mod_value_all(n,1) = mod_val;
    
        t_all(1:seq_length_max,end+1) = nan;
        Fx_all(1:seq_length_max,end+1) = nan;
        Fy_all(1:seq_length_max,end+1) = nan;
        Fz_all(1:seq_length_max,end+1) = nan;
        Mx_all(1:seq_length_max,end+1) = nan;
        My_all(1:seq_length_max,end+1) = nan;
        Mz_all(1:seq_length_max,end+1) = nan;
        stroke_L_all(1:seq_length_max,end+1) = nan;
        pitch_L_all(1:seq_length_max,end+1) = nan;
        dev_L_all(1:seq_length_max,end+1) = nan;
        stroke_R_all(1:seq_length_max,end+1) = nan;
        pitch_R_all(1:seq_length_max,end+1) = nan;
        dev_R_all(1:seq_length_max,end+1) = nan;


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
        freq_all(n,1) = f_now;

        Fx_mean_all(n,1) = nanmean(Fx_now(n_start:n_stop));
        Fy_mean_all(n,1) = nanmean(Fy_now(n_start:n_stop));
        Fz_mean_all(n,1) = nanmean(Fz_now(n_start:n_stop));

        Mx_mean_all(n,1) = nanmean(Mx_now(n_start:n_stop));
        My_mean_all(n,1) = nanmean(My_now(n_start:n_stop));
        Mz_mean_all(n,1) = nanmean(Mz_now(n_start:n_stop));

        Astroke_all(n,1) = max([stroke_R_now;stroke_L_now]) - min([stroke_R_now;stroke_L_now]);
        Apitch_all(n,1) = max([pitch_R_now;pitch_L_now]) - min([pitch_R_now;pitch_L_now]);
        Adev_all(n,1) = max([dev_R_now;dev_L_now]) - min([dev_R_now;dev_L_now]);

        %% store all data
        t_all(1:length(t),n) = t_now;

        Fx_all(1:length(t),n) = Fx_now;
        Fy_all(1:length(t),n) = Fy_now;
        Fz_all(1:length(t),n) = Fz_now;

        Mx_all(1:length(t),n) = Mx_now;
        My_all(1:length(t),n) = My_now;
        Mz_all(1:length(t),n) = Mz_now;

        stroke_L_all(1:length(t),n) = stroke_L_now;
        pitch_L_all(1:length(t),n) = pitch_L_now;
        dev_L_all(1:length(t),n) = dev_L_now;

        stroke_R_all(1:length(t),n) = stroke_R_now;
        pitch_R_all(1:length(t),n) = pitch_R_now;
        dev_R_all(1:length(t),n) = dev_R_now;
    end
end

F_mean_all = sqrt(Fx_mean_all.^2+Fy_mean_all.^2+Fz_mean_all.^2);

if length(t_all) > seq_length_max
    seq_length_max
    'TOO SHORT'
end




















