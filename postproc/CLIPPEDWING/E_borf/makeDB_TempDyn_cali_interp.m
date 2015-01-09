% all MODs data

% file_name = [save_name,'_all_'];
list = dir([file_name,'*mat']);

% load wb data & make cali db & determine wingbeat start
    file_now = [file_name,num2str(0),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
        
        n_wb = length(kine)/Nwb;
        
        stroke_L_now = rad2deg(kine(:,5));
        pitch_L_now = rad2deg(kine(:,1));
        dev_L_now = rad2deg(kine(:,2));

        stroke_R_now = rad2deg(kine(:,6));
        pitch_R_now = rad2deg(kine(:,4));
        dev_R_now = rad2deg(kine(:,3));
        
        %% wingbeat start
%         wb_start(1) = 1;
%         for j=2:Nwb
%             
%             stroke_now = stroke_L_now(n_wb*(j-1.5):n_wb*(j-.5));
%             wb_start(j,1) = find(stroke_now==min(stroke_now)) + n_wb*(j-1.5)-1;
%             wb_stop(j-1,1) = wb_start(j);
%         end
%         wb_stop(Nwb) = length(kine);
% 
%         t_start = t(round(wb_start));
%         t_stop = t(round(wb_stop));
%         t_mean = mean([t_start t_stop]')';

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
        
        for j=1:length(stroke_L_now)
            
            n_str_L = find(str_steps_interp == round(stroke_L_now(j)));
            n_rot_L = find(rot_steps_interp == round(pitch_L_now(j)));
            n_dev_L = find(dev_steps_interp == round(dev_L_now(j)));
            
            n_str_R = find(str_steps_interp == round(stroke_R_now(j)));
            n_rot_R = find(rot_steps_interp == round(pitch_R_now(j)));
            n_dev_R = find(dev_steps_interp == round(dev_R_now(j)));
            
            Fx_cali_L(j,1) = cali_L.Fx_interp(n_str_L,n_rot_L,n_dev_L);
            Fy_cali_L(j,1) = cali_L.Fy_interp(n_str_L,n_rot_L,n_dev_L);
            Fz_cali_L(j,1) = cali_L.Fz_interp(n_str_L,n_rot_L,n_dev_L);
            Mx_cali_L(j,1) = cali_L.Mx_interp(n_str_L,n_rot_L,n_dev_L);
            My_cali_L(j,1) = cali_L.My_interp(n_str_L,n_rot_L,n_dev_L);
            Mz_cali_L(j,1) = cali_L.Mz_interp(n_str_L,n_rot_L,n_dev_L);
            
            Fx_cali_R(j,1) = cali_R.Fx_interp(n_str_R,n_rot_R,n_dev_R);
            Fy_cali_R(j,1) = cali_R.Fy_interp(n_str_R,n_rot_R,n_dev_R);
            Fz_cali_R(j,1) = cali_R.Fz_interp(n_str_R,n_rot_R,n_dev_R);
            Mx_cali_R(j,1) = cali_R.Mx_interp(n_str_R,n_rot_R,n_dev_R);
            My_cali_R(j,1) = cali_R.My_interp(n_str_R,n_rot_R,n_dev_R);
            Mz_cali_R(j,1) = cali_R.Mz_interp(n_str_R,n_rot_R,n_dev_R);
            
        end
    end
        

% load force data
for i = 1:length(list)
    file_now = [file_name,num2str(i-1),'.mat'];
    
    if exist(file_now) == 2
        load(file_now)
        
        %% data NOW
        Fx_robo_now =  ft(:,1);  % thrust, fwd pos
        Fy_robo_now =  ft(:,3);  % sideways, right pos
        Fz_robo_now = -ft(:,2); % vertical, down pos

        Mx_robo_now =  ft(:,4);  % Roll, right pos
        My_robo_now =  ft(:,6);  % pitch, nose up pos
        Mz_robo_now = -ft(:,5); % yaw, right pos
        
        %% remove weight effect
        Fx_now = Fx_robo_now - Fx_cali_L - Fx_cali_R;
        Fy_now = Fy_robo_now - Fy_cali_L - Fy_cali_R;
        Fz_now = Fz_robo_now - Fz_cali_L - Fz_cali_R;
        Mx_now = Mx_robo_now - Mx_cali_L - Mx_cali_R;
        My_now = My_robo_now - My_cali_L - My_cali_R;
        Mz_now = Mz_robo_now - Mz_cali_L - Mz_cali_R;
        
        % robo2fly
        Fx_fly(:,i) = Fx_now * F_robo2fly;
        Fy_fly(:,i) = Fy_now * F_robo2fly;
        Fz_fly(:,i) = Fz_now * F_robo2fly;

        Mx_fly(:,i) = Mx_now * M_robo2fly;
        My_fly(:,i) = My_now * M_robo2fly;
        Mz_fly(:,i) = Mz_now * M_robo2fly;
        
    end
end

        %% calc wb mean data
        Fx_mean_all = nanmean(Fx_fly')';
        Fy_mean_all = nanmean(Fy_fly')';
        Fz_mean_all = nanmean(Fz_fly')';

        Mx_mean_all = nanmean(Mx_fly')';
        My_mean_all = nanmean(My_fly')';
        Mz_mean_all = nanmean(Mz_fly')';
        
        %% WB mean data
        for i=1:Nwb
            Fx_mean_mean(i,1) = nanmean(Fx_mean_all(wb_start(i):wb_stop(i)));
            Fy_mean_mean(i,1) = nanmean(Fy_mean_all(wb_start(i):wb_stop(i)));
            Fz_mean_mean(i,1) = nanmean(Fz_mean_all(wb_start(i):wb_stop(i)));
            
            Mx_mean_mean(i,1) = nanmean(Mx_mean_all(wb_start(i):wb_stop(i)));
            My_mean_mean(i,1) = nanmean(My_mean_all(wb_start(i):wb_stop(i)));
            Mz_mean_mean(i,1) = nanmean(Mz_mean_all(wb_start(i):wb_stop(i)));
        end

        %% store all data
        t_all = t;

        stroke_L_all = stroke_L_now;
        pitch_L_all = pitch_L_now;
        dev_L_all = dev_L_now;

        stroke_R_all = stroke_R_now;
        pitch_R_all = pitch_R_now;
        dev_R_all = dev_R_now;

        Fx_all = Fx_fly;
        Fy_all = Fy_fly;
        Fz_all = Fz_fly;

        Mx_all = Mx_fly;
        My_all = My_fly;
        Mz_all = Mz_fly;


