% allNOfreq MODs data

% file_name = [save_name,'_all_const_freq_'];
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

stroke_L_allNOfreq = [];
pitch_L_allNOfreq = [];
dev_L_allNOfreq = [];

stroke_R_allNOfreq = [];
pitch_R_allNOfreq = [];
dev_R_allNOfreq = [];

Fx_allNOfreq = [];
Fy_allNOfreq = [];
Fz_allNOfreq = [];
Mx_allNOfreq = [];
My_allNOfreq = [];
Mz_allNOfreq = [];
n=0;

% load data
for i = 1:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
        
        %% FILTER RAW DATA
        ft(:,1) = filter_borf_raw(ft(:,1),butter_cut,butter_n,Nwb_total);
        ft(:,2) = filter_borf_raw(ft(:,2),butter_cut,butter_n,Nwb_total);
        ft(:,3) = filter_borf_raw(ft(:,3),butter_cut,butter_n,Nwb_total);
        ft(:,4) = filter_borf_raw(ft(:,4),butter_cut,butter_n,Nwb_total);
        ft(:,5) = filter_borf_raw(ft(:,5),butter_cut,butter_n,Nwb_total);
        ft(:,6) = filter_borf_raw(ft(:,6),butter_cut,butter_n,Nwb_total);

        %     load(list(i).name)
        
        n_wb = length(kine)/Nwb_total;
        wb_kin_start = (wb_start-1)*n_wb+1;
        wb_kin_stop = wb_kin_start+n_wb-1;
        
        n=n+1;
        run_type_allNOfreq{n} = exp_type;
        mod_value_allNOfreq(n,1) = mod_val;
    
        t_allNOfreq(1:seq_length_max,end+1) = nan;
        
        stroke_L_allNOfreq(1:seq_length_max,end+1) = nan;
        pitch_L_allNOfreq(1:seq_length_max,end+1) = nan;
        dev_L_allNOfreq(1:seq_length_max,end+1) = nan;
        
        stroke_R_allNOfreq(1:seq_length_max,end+1) = nan;
        pitch_R_allNOfreq(1:seq_length_max,end+1) = nan;
        dev_R_allNOfreq(1:seq_length_max,end+1) = nan;
        
        Fx_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;
        Fy_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;
        Fz_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;
        Mx_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;
        My_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;
        Mz_allNOfreq(1:seq_length_max,1:Nwb,end+1) = nan;


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
%         Fx_cali_L = nan(size(stroke_L_now));
%         Fy_cali_L = nan(size(stroke_L_now));
%         Fz_cali_L = nan(size(stroke_L_now));
%         Fx_cali_R = nan(size(stroke_L_now));
%         Fy_cali_R = nan(size(stroke_L_now));
%         Fz_cali_R = nan(size(stroke_L_now));
%         Mx_cali_L = nan(size(stroke_L_now));
%         My_cali_L = nan(size(stroke_L_now));
%         Mz_cali_L = nan(size(stroke_L_now));
%         Mx_cali_R = nan(size(stroke_L_now));
%         My_cali_R = nan(size(stroke_L_now));
%         Mz_cali_R = nan(size(stroke_L_now));
%         
% %         Fx_now = nan(size(Fx_robo));
% %         Fy_now = nan(size(Fx_robo));
% %         Fz_now = nan(size(Fx_robo));
% %         Mx_now = nan(size(Fx_robo));
% %         My_now = nan(size(Fx_robo));
% %         Mz_now = nan(size(Fx_robo));
%         
%         for j=1:length(stroke_L_now)
% %             length(stroke_L_now)-j
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
% 
% %             Fx_now(j,:) = Fx_robo(j,:) - Fx_cali_L(j) - Fx_cali_R(j);
% %             Fy_now(j,:) = Fy_robo(j,:) - Fy_cali_L(j) - Fy_cali_R(j);
% %             Fz_now(j,:) = Fz_robo(j,:) - Fz_cali_L(j) - Fz_cali_R(j);
% %             Mx_now(j,:) = Mx_robo(j,:) - Mx_cali_L(j) - Mx_cali_R(j);
% %             My_now(j,:) = My_robo(j,:) - My_cali_L(j) - My_cali_R(j);
% %             Mz_now(j,:) = Mz_robo(j,:) - Mz_cali_L(j) - Mz_cali_R(j);
%         end
%         
%         
%         %% remove weight effect
%         Fx_cali_L_array = [];
%         Fy_cali_L_array = [];
%         Fz_cali_L_array = [];
%         Mx_cali_L_array = [];
%         My_cali_L_array = [];
%         Mz_cali_L_array = [];
%         Fx_cali_R_array = [];
%         Fy_cali_R_array = [];
%         Fz_cali_R_array = [];
%         Mx_cali_R_array = [];
%         My_cali_R_array = [];
%         Mz_cali_R_array = [];
%         
%         for i=1:Nwb
%             Fx_cali_L_array = [Fx_cali_L_array Fx_cali_L];
%             Fy_cali_L_array = [Fy_cali_L_array Fy_cali_L];
%             Fz_cali_L_array = [Fz_cali_L_array Fz_cali_L];
%             Mx_cali_L_array = [Mx_cali_L_array Mx_cali_L];
%             My_cali_L_array = [My_cali_L_array My_cali_L];
%             Mz_cali_L_array = [Mz_cali_L_array Mz_cali_L];
%             Fx_cali_R_array = [Fx_cali_R_array Fx_cali_R];
%             Fy_cali_R_array = [Fy_cali_R_array Fy_cali_R];
%             Fz_cali_R_array = [Fz_cali_R_array Fz_cali_R];
%             Mx_cali_R_array = [Mx_cali_R_array Mx_cali_R];
%             My_cali_R_array = [My_cali_R_array My_cali_R];
%             Mz_cali_R_array = [Mz_cali_R_array Mz_cali_R];
%         end
% 
%         Fx_now = Fx_robo - Fx_cali_L_array - Fx_cali_R_array;
%         Fy_now = Fy_robo - Fy_cali_L_array - Fy_cali_R_array;
%         Fz_now = Fz_robo - Fz_cali_L_array - Fz_cali_R_array;
%         Mx_now = Mx_robo - Mx_cali_L_array - Mx_cali_R_array;
%         My_now = My_robo - My_cali_L_array - My_cali_R_array;
%         Mz_now = Mz_robo - Mz_cali_L_array - Mz_cali_R_array;

        %% NOcali
        Fx_now = Fx_robo;
        Fy_now = Fy_robo;
        Fz_now = Fz_robo;
        Mx_now = Mx_robo;
        My_now = My_robo;
        Mz_now = Mz_robo;

        % robo2fly
        Fx_fly = Fx_now * F_robo2fly;
        Fy_fly = Fy_now * F_robo2fly;
        Fz_fly = Fz_now * F_robo2fly;

        Mx_fly = Mx_now * M_robo2fly;
        My_fly = My_now * M_robo2fly;
        Mz_fly = Mz_now * M_robo2fly;

        %% calc wb mean data
        Fx_mean_allNOfreq(n,1) = nanmean(Fx_fly(:));
        Fy_mean_allNOfreq(n,1) = nanmean(Fy_fly(:));
        Fz_mean_allNOfreq(n,1) = nanmean(Fz_fly(:));

        Mx_mean_allNOfreq(n,1) = nanmean(Mx_fly(:));
        My_mean_allNOfreq(n,1) = nanmean(My_fly(:));
        Mz_mean_allNOfreq(n,1) = nanmean(Mz_fly(:));

        %% store allNOfreq data
        t_allNOfreq(1:length(t_now),n) = t_now;

        stroke_L_allNOfreq(1:length(t_now),n) = stroke_L_now;
        pitch_L_allNOfreq(1:length(t_now),n) = pitch_L_now;
        dev_L_allNOfreq(1:length(t_now),n) = dev_L_now;

        stroke_R_allNOfreq(1:length(t_now),n) = stroke_R_now;
        pitch_R_allNOfreq(1:length(t_now),n) = pitch_R_now;
        dev_R_allNOfreq(1:length(t_now),n) = dev_R_now;

        Fx_allNOfreq(1:length(t_now),:,n) = Fx_fly;
        Fy_allNOfreq(1:length(t_now),:,n) = Fy_fly;
        Fz_allNOfreq(1:length(t_now),:,n) = Fz_fly;

        Mx_allNOfreq(1:length(t_now),:,n) = Mx_fly;
        My_allNOfreq(1:length(t_now),:,n) = My_fly;
        Mz_allNOfreq(1:length(t_now),:,n) = Mz_fly;
    end
end

if length(t_allNOfreq) > seq_length_max
    seq_length_max
    'TOO SHORT'
end




















