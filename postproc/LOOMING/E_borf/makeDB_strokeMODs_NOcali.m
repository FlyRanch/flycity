% stroke MODs data

file_name = [save_name,'_stroke_'];
list = dir([file_name,'*mat']);

% % init vars with nans
% t_stroke = nan(seq_length_max,length(list));
% Fx_stroke = nan(seq_length_max,length(list));
% Fy_stroke = nan(seq_length_max,length(list));
% Fz_stroke = nan(seq_length_max,length(list));
% Mx_stroke = nan(seq_length_max,length(list));
% My_stroke = nan(seq_length_max,length(list));
% Mz_stroke = nan(seq_length_max,length(list));
% stroke_L_stroke = nan(seq_length_max,length(list));
% pitch_L_stroke = nan(seq_length_max,length(list));
% dev_L_stroke = nan(seq_length_max,length(list));
% stroke_R_stroke = nan(seq_length_max,length(list));
% pitch_R_stroke = nan(seq_length_max,length(list));
% dev_R_stroke = nan(seq_length_max,length(list));

t_stroke = [];
Fx_stroke = [];
Fy_stroke = [];
Fz_stroke = [];
Mx_stroke = [];
My_stroke = [];
Mz_stroke = [];
stroke_L_stroke = [];
pitch_L_stroke = [];
dev_L_stroke = [];
stroke_R_stroke = [];
pitch_R_stroke = [];
dev_R_stroke = [];

n=0;

% load data
for i = 1:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
        %     load(list(i).name)
        
        n=n+1;
        run_type_stroke{n} = exp_type;
        mod_value_stroke(n,1) = mod_val;
    
        t_stroke(1:seq_length_max,end+1) = nan;
        Fx_stroke(1:seq_length_max,end+1) = nan;
        Fy_stroke(1:seq_length_max,end+1) = nan;
        Fz_stroke(1:seq_length_max,end+1) = nan;
        Mx_stroke(1:seq_length_max,end+1) = nan;
        My_stroke(1:seq_length_max,end+1) = nan;
        Mz_stroke(1:seq_length_max,end+1) = nan;
        stroke_L_stroke(1:seq_length_max,end+1) = nan;
        pitch_L_stroke(1:seq_length_max,end+1) = nan;
        dev_L_stroke(1:seq_length_max,end+1) = nan;
        stroke_R_stroke(1:seq_length_max,end+1) = nan;
        pitch_R_stroke(1:seq_length_max,end+1) = nan;
        dev_R_stroke(1:seq_length_max,end+1) = nan;


        %% data NOW
        t_now = t;

        % forces in fly ref frame (different from borf)
        Fx_robo =  ft(:,1);  % thrust, fwd pos
        Fy_robo =  ft(:,3);  % sideways, right pos
        Fz_robo = -ft(:,2); % vertical, down pos

        Mx_robo =  ft(:,4);  % Roll, right pos
        My_robo =  ft(:,6);  % pitch, nose up pos
        Mz_robo = -ft(:,5); % yaw, right pos

        stroke_L_now = rad2deg(kine(:,5));
        pitch_L_now = rad2deg(kine(:,1));
        dev_L_now = rad2deg(kine(:,2));

        stroke_R_now = rad2deg(kine(:,6));
        pitch_R_now = rad2deg(kine(:,4));
        dev_R_now = rad2deg(kine(:,3));
        
%         % calc cali data
%         Fx_cali_L = nan(size(Fx_robo));
%         Fy_cali_L = nan(size(Fx_robo));
%         Fz_cali_L = nan(size(Fx_robo));
%         Fx_cali_R = nan(size(Fx_robo));
%         Fy_cali_R = nan(size(Fx_robo));
%         Fz_cali_R = nan(size(Fx_robo));
%         Mx_cali_L = nan(size(Mx_robo));
%         My_cali_L = nan(size(Mx_robo));
%         Mz_cali_L = nan(size(Mx_robo));
%         Mx_cali_R = nan(size(Mx_robo));
%         My_cali_R = nan(size(Mx_robo));
%         Mz_cali_R = nan(size(Mx_robo));
%         
%         Fx_now = nan(size(Mx_robo));
%         Fy_now = nan(size(Mx_robo));
%         Fz_now = nan(size(Mx_robo));
%         Mx_now = nan(size(Mx_robo));
%         My_now = nan(size(Mx_robo));
%         Mz_now = nan(size(Mx_robo));
%         
%         for j=1:length(stroke_L_now)
%             length(stroke_L_now)-j
%             Fx_cali_L(j,1) = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fx_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Fy_cali_L(j,1) = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fy_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Fz_cali_L(j,1) = -interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Fz_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             
%             Fx_cali_R(j,1) = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fx_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Fy_cali_R(j,1) = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fy_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Fz_cali_R(j,1) = -interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Fz_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             
%             Mx_cali_L(j,1) = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mx_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             My_cali_L(j,1) = interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.My_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Mz_cali_L(j,1) = -interp3(cali_L.str_mesh,cali_L.rot_mesh,cali_L.dev_mesh,cali_L.Mz_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             
%             Mx_cali_R(j,1) = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mx_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             My_cali_R(j,1) = interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.My_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%             Mz_cali_R(j,1) = -interp3(cali_R.str_mesh,cali_R.rot_mesh,cali_R.dev_mesh,cali_R.Mz_array,...
%                 stroke_L_now(j),pitch_L_now(j),dev_L_now(j),'spline');
%         end
%         
%         Fx_now = (Fx_robo - Fx_cali_L - Fx_cali_R) * F_robo2fly;
%         Fy_now = (Fy_robo - Fy_cali_L - Fy_cali_R) * F_robo2fly;
%         Fz_now = (Fz_robo - Fz_cali_L - Fz_cali_R) * F_robo2fly;
% 
%         Mx_now = (Mx_robo - Mx_cali_L - Mx_cali_R) * M_robo2fly;
%         My_now = (My_robo - My_cali_L - My_cali_R) * M_robo2fly;
%         Mz_now = (Mz_robo - Mz_cali_L - Mz_cali_R) * M_robo2fly;

        %% NO cali
        Fx_now = Fx_robo * F_robo2fly;
        Fy_now = Fy_robo * F_robo2fly;
        Fz_now = Fz_robo * F_robo2fly;

        Mx_now = Mx_robo * M_robo2fly;
        My_now = My_robo * M_robo2fly;
        Mz_now = Mz_robo * M_robo2fly;

        %% calc wb mean data
        N_wb_now = round(length(t_now)/Nwb);
        n_start = wb_start * N_wb_now;
        n_stop = wb_stop * N_wb_now;

        frobo_now = (1/max(t_now))*Nwb;
        f_now = f_robo2fly * frobo_now;
        freq_stroke(n,1) = f_now;

        Fx_mean_stroke(n,1) = nanmean(Fx_now(n_start:n_stop));
        Fy_mean_stroke(n,1) = nanmean(Fy_now(n_start:n_stop));
        Fz_mean_stroke(n,1) = nanmean(Fz_now(n_start:n_stop));

        Mx_mean_stroke(n,1) = nanmean(Mx_now(n_start:n_stop));
        My_mean_stroke(n,1) = nanmean(My_now(n_start:n_stop));
        Mz_mean_stroke(n,1) = nanmean(Mz_now(n_start:n_stop));

        Astroke_stroke(n,1) = max([stroke_R_now;stroke_L_now]) - min([stroke_R_now;stroke_L_now]);
        Apitch_stroke(n,1) = max([pitch_R_now;pitch_L_now]) - min([pitch_R_now;pitch_L_now]);
        Adev_stroke(n,1) = max([dev_R_now;dev_L_now]) - min([dev_R_now;dev_L_now]);

        %% store stroke data
        t_stroke(1:length(t),n) = t_now;

        Fx_stroke(1:length(t),n) = Fx_now;
        Fy_stroke(1:length(t),n) = Fy_now;
        Fz_stroke(1:length(t),n) = Fz_now;

        Mx_stroke(1:length(t),n) = Mx_now;
        My_stroke(1:length(t),n) = My_now;
        Mz_stroke(1:length(t),n) = Mz_now;

        stroke_L_stroke(1:length(t),n) = stroke_L_now;
        pitch_L_stroke(1:length(t),n) = pitch_L_now;
        dev_L_stroke(1:length(t),n) = dev_L_now;

        stroke_R_stroke(1:length(t),n) = stroke_R_now;
        pitch_R_stroke(1:length(t),n) = pitch_R_now;
        dev_R_stroke(1:length(t),n) = dev_R_now;
    end
end

F_mean_stroke = sqrt(Fx_mean_stroke.^2+Fy_mean_stroke.^2+Fz_mean_stroke.^2);

if length(t_stroke) > seq_length_max
    seq_length_max
    'TOO SHORT'
end




















