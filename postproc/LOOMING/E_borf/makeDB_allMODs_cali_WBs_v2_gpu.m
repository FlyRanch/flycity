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
        
        n_wb = length(kine)/Nwb_total;
        wb_kin_start = (wb_start-1)*n_wb+1;
        wb_kin_stop = wb_kin_start+n_wb-1;
        
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
        t_now = t(wb_kin_start:wb_kin_stop);

        % kinematics
        stroke_L_now = rad2deg(kine(wb_kin_start:wb_kin_stop,5));
        pitch_L_now = rad2deg(kine(wb_kin_start:wb_kin_stop,1));
        dev_L_now = rad2deg(kine(wb_kin_start:wb_kin_stop,2));

        stroke_R_now = rad2deg(kine(wb_kin_start:wb_kin_stop,6));
        pitch_R_now = rad2deg(kine(wb_kin_start:wb_kin_stop,4));
        dev_R_now = rad2deg(kine(wb_kin_start:wb_kin_stop,3));
        
        % forces in fly ref frame (different from borf)
        clear Fx_robo Fy_robo Fz_robo Mx_robo My_robo Mz_robo 
        for i=1:Nwb
            Fx_robo(:,i) =  ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,1);  % thrust, fwd pos
            Fy_robo(:,i) =  ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,3);  % sideways, right pos
            Fz_robo(:,i) = -ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,2); % vertical, down pos

            Mx_robo(:,i) =  ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,4);  % Roll, right pos
            My_robo(:,i) =  ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,6);  % pitch, nose up pos
            Mz_robo(:,i) = -ft(wb_kin_start+(i-1)*n_wb:wb_kin_start+i*n_wb-1,5); % yaw, right pos
        end
        
        %% calc cali data
        Fx_cali_L = nan(size(stroke_L_now));
        Fy_cali_L = nan(size(stroke_L_now));
        Fz_cali_L = nan(size(stroke_L_now));
        Fx_cali_R = nan(size(stroke_L_now));
        Fy_cali_R = nan(size(stroke_L_now));
        Fz_cali_R = nan(size(stroke_L_now));
        Mx_cali_L = nan(size(stroke_L_now));
        My_cali_L = nan(size(stroke_L_now));
        Mz_cali_L = nan(size(stroke_L_now));
        Mx_cali_R = nan(size(stroke_L_now));
        My_cali_R = nan(size(stroke_L_now));
        Mz_cali_R = nan(size(stroke_L_now));
        
        Fx_now = nan(size(Fx_robo));
        Fy_now = nan(size(Fx_robo));
        Fz_now = nan(size(Fx_robo));
        Mx_now = nan(size(Fx_robo));
        My_now = nan(size(Fx_robo));
        Mz_now = nan(size(Fx_robo));
        
        for j=1:length(stroke_L_now)
            length(stroke_L_now)-j
            Fx_cali_L(j,1) = interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fx_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Fy_cali_L(j,1) = interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fy_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Fz_cali_L(j,1) = -interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fz_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            
            Fx_cali_R(j,1) = interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fx_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Fy_cali_R(j,1) = interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fy_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Fz_cali_R(j,1) = -interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fz_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            
            Mx_cali_L(j,1) = interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mx_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            My_cali_L(j,1) = interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.My_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Mz_cali_L(j,1) = -interp3_gpu(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mz_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            
            Mx_cali_R(j,1) = interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mx_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            My_cali_R(j,1) = interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.My_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));
            Mz_cali_R(j,1) = -interp3_gpu(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mz_array,...
                stroke_L_now(j),pitch_L_now(j),dev_L_now(j));

            Fx_now(j,:) = Fx_robo(j,:) - Fx_cali_L(j) - Fx_cali_R(j);
            Fy_now(j,:) = Fy_robo(j,:) - Fy_cali_L(j) - Fy_cali_R(j);
            Fz_now(j,:) = Fz_robo(j,:) - Fz_cali_L(j) - Fz_cali_R(j);
            
            Mx_now(j,:) = Mx_robo(j,:) - Mx_cali_L(j) - Mx_cali_R(j);
            My_now(j,:) = My_robo(j,:) - My_cali_L(j) - My_cali_R(j);
            Mz_now(j,:) = Mz_robo(j,:) - Mz_cali_L(j) - Mz_cali_R(j);
        end

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




















