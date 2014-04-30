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
        wb_start(1) = 1;
        for j=2:Nwb
            
            stroke_now = stroke_L_now(n_wb*(j-1.5):n_wb*(j-.5));
            wb_start(j,1) = find(stroke_now==min(stroke_now)) + n_wb*(j-1.5)-1;
            wb_stop(j-1,1) = wb_start(j);
        end
        wb_stop(Nwb) = length(kine);

        t_start = t(round(wb_start));
        t_stop = t(round(wb_stop));
        t_mean = mean([t_start t_stop]')';
    end