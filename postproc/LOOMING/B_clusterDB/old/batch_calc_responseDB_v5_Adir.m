
clc
clear
close all
% addpath('/home/florian/Dropbox/WORK/toolbox')

% filenames
% inputfile = 'flightpathDB_2flies_9clusters_Atmin.mat';
% skip = 1
inputfile = 'flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3.mat';
% inputfile = 'flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3.mat';
load(inputfile)

skip = 10
nstart = 175;

outputfile = [inputfile(1:end-4),'_response.mat'];
figdir = 'accel_plots_2.5n-3n2.5';
mkdir(figdir)

toplot=1;
toplot=0;


responseDB.IDX_shift = nan(10,size(pathDB.pos,2));
responseDB.t_shift = nan(10,size(pathDB.pos,2));

for i = 1:size(pathDB.pos,2)
    size(pathDB.pos,2)-i

%     [IDX_shift,...
%     t_shift,t_resp,t_resp_end,...
%     t_turn_start,t_turn_stop,t_turn_max,...
%     t_accel_start,t_accel_stop,t_accel_max,...
%     t_decel_start,t_decel_stop,t_decel_max,...
%     dt_turn,dt_accel,dt_decel,...
%     A_hor_pre_resp,A_hor_pre_turn,A_hor_pre_accel,...
%     A_hor_post_resp,A_hor_post_turn,A_hor_post_accel,...
%     A_hor_max,A_hor_mean,A_max,A_mean,...
%     An_hor_max,At_hor_max,At_hor_min,...
%     t_Ahormax,V_Ahormax,An_hor_Ahormax,At_hor_Ahormax,...
%     stim_angle_vel_Ahormax,stim_angle_accel_Ahormax,stim_angle_yaw_Ahormax,...
%     accel_angle_hor_vel_Ahormax,accel_angle_hor_body_Ahormax,...
%     slip_Ahormax,pitch_Ahormax,roll_Ahormax,...
%     t_Amax,V_Amax,An_hor_Amax,At_hor_Amax,...
%     stim_angle_vel_Amax,stim_angle_accel_Amax,stim_angle_yaw_Amax,...
%     accel_angle_hor_vel_Amax,accel_angle_hor_body_Amax,...
%     slip_Amax,pitch_Amax,roll_Amax,...
%     stim_angle_vel_trigger,stim_angle_vel_pre_resp,stim_angle_vel_trig2resp,stim_angle_vel_pre_turn,...
%     stim_angle_vel_post_turn,stim_angle_vel_post_resp,stim_angle_vel_postresp2end,...
%     stim_angle_accel_trigger,stim_angle_accel_pre_resp,stim_angle_accel_trig2resp,stim_angle_accel_pre_turn,...
%     stim_angle_accel_post_turn,stim_angle_accel_post_resp,stim_angle_accel_postresp2end,...
%     stim_angle_accel_resp_min,stim_angle_accel_resp_max,stim_angle_accel_resp_mean,stim_angle_accel_resp_std,...
%     stim_angle_yaw_trigger,stim_angle_yaw_pre_resp,stim_angle_yaw_trig2resp,stim_angle_yaw_pre_turn,...
%     stim_angle_yaw_post_turn,stim_angle_yaw_post_resp,stim_angle_yaw_postresp2end,...
%     accel_angle_hor_vel_pre_resp,accel_angle_hor_vel_pre_turn,accel_angle_hor_vel_pre_accel,...
%     accel_angle_hor_vel_post_turn,accel_angle_hor_vel_post_resp,accel_angle_hor_vel_post_accel,...
%     accel_angle_hor_body_pre_resp,accel_angle_hor_body_pre_turn,accel_angle_hor_body_pre_accel,...
%     accel_angle_hor_body_post_turn,accel_angle_hor_body_post_resp,accel_angle_hor_body_post_accel,...
%     accel_angle_hor_vel_resp_mean,accel_angle_hor_body_resp_mean,...
%     slip_trigger,slip_pre_resp,slip_trig2resp,slip_pre_turn,...
%     slip_post_turn,slip_post_resp,slip_postresp2end,...
%     pitch_trigger,pitch_pre_resp,pitch_trig2resp,pitch_pre_turn,...
%     pitch_post_turn,pitch_post_resp,pitch_postresp2end,...
%     roll_trigger,roll_pre_resp,roll_trig2resp,roll_pre_turn,...
%     roll_post_turn,roll_post_resp,roll_postresp2end,...
%     V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
%     V_post_resp,V_postresp2end,V_post_accel,...
%     teta_pre_resp,teta_post_resp,teta_pre_turn,teta_post_turn,teta_pre_accel,teta_post_accel] =...
%     calc_response_data_Adir(pathDB,patternDB,settings,i,toplot,skip,nstart);

    calc_response_data_F(pathDB,patternDB,settings,i,toplot,skip,nstart);


    responseDB.IDX_shift(1:length(t_shift),i) = IDX_shift;
    responseDB.t_shift(1:length(t_shift),i) = t_shift;

    responseDB.t_resp(i,1) = t_resp;
    responseDB.t_resp_end(i,1) = t_resp_end;

    responseDB.t_turn_start(i,1) = t_turn_start;
    responseDB.t_turn_stop(i,1) = t_turn_stop;
    responseDB.t_turn_max(i,1) = t_turn_max;

    responseDB.t_accel_start(i,1) = t_accel_start;
    responseDB.t_accel_stop(i,1) = t_accel_stop;
    responseDB.t_accel_max(i,1) = t_accel_max; 

    responseDB.t_decel_start(i,1) = t_decel_start;
    responseDB.t_decel_stop(i,1) = t_decel_stop;
    responseDB.t_decel_max(i,1) = t_decel_max;

    responseDB.dt_turn(i,1) = dt_turn;
    responseDB.dt_accel(i,1) = dt_accel;
    responseDB.dt_decel(i,1) = dt_decel;

    responseDB.A_hor_pre_resp(i,1) = A_hor_pre_resp;   
    responseDB.A_hor_pre_turn(i,1) = A_hor_pre_turn;   
    responseDB.A_hor_pre_accel(i,1) = A_hor_pre_accel;   

    responseDB.A_hor_post_resp(i,1) = A_hor_post_resp;   
    responseDB.A_hor_post_turn(i,1) = A_hor_post_turn;   
    responseDB.A_hor_post_accel(i,1) = A_hor_post_accel;   

    responseDB.A_hor_max(i,1) = A_hor_max;
    responseDB.A_hor_mean(i,1) = A_hor_mean;

    responseDB.A_max(i,1) = A_max;
    responseDB.A_mean(i,1) = A_mean;

    responseDB.An_hor_max(i,1) = An_hor_max;
    responseDB.At_hor_max(i,1) = At_hor_max;
    responseDB.At_hor_min(i,1) = At_hor_min;

    responseDB.t_Ahormax(i,1) = nanmean(t_Ahormax);
    responseDB.V_Ahormax(i,1) = nanmean(V_Ahormax);
    responseDB.An_hor_Ahormax(i,1) = nanmean(An_hor_Ahormax);
    responseDB.At_hor_Ahormax(i,1) = nanmean(At_hor_Ahormax);

    responseDB.stim_angle_vel_Ahormax(i,1) = nanmean(stim_angle_vel_Ahormax);
    responseDB.stim_angle_accel_Ahormax(i,1) = nanmean(stim_angle_accel_Ahormax);
    responseDB.stim_angle_yaw_Ahormax(i,1) = nanmean(stim_angle_yaw_Ahormax);

    responseDB.accel_angle_hor_vel_Ahormax(i,1) = nanmean(accel_angle_hor_vel_Ahormax);
    responseDB.accel_angle_hor_body_Ahormax(i,1) = nanmean(accel_angle_hor_body_Ahormax);

    responseDB.slip_Ahormax(i,1) = nanmean(slip_Ahormax);
    responseDB.pitch_Ahormax(i,1) = nanmean(pitch_Ahormax);
    responseDB.roll_Ahormax(i,1) = nanmean(roll_Ahormax);

    responseDB.t_Amax(i,1) = nanmean(t_Amax);
    responseDB.V_Amax(i,1) = nanmean(V_Amax);
    responseDB.An_hor_Amax(i,1) = nanmean(An_hor_Amax);
    responseDB.At_hor_Amax(i,1) = nanmean(At_hor_Amax);

    responseDB.stim_angle_vel_Amax(i,1) = nanmean(stim_angle_vel_Amax);
    responseDB.stim_angle_accel_Amax(i,1) = nanmean(stim_angle_accel_Amax);
    responseDB.stim_angle_yaw_Amax(i,1) = nanmean(stim_angle_yaw_Amax);

    responseDB.accel_angle_hor_vel_Amax(i,1) = nanmean(accel_angle_hor_vel_Amax);
    responseDB.accel_angle_hor_body_Amax(i,1) = nanmean(accel_angle_hor_body_Amax);

    responseDB.slip_Amax(i,1) = nanmean(slip_Amax);
    responseDB.pitch_Amax(i,1) = nanmean(pitch_Amax);
    responseDB.roll_Amax(i,1) = nanmean(roll_Amax);

    responseDB.F_Amax(i,1) = nanmean(slip_Amax);
    responseDB.Fsp_pitch_Amax(i,1) = nanmean(Fsp_pitch_Amax);
    responseDB.Fsp_roll_Amax(i,1) = nanmean(Fsp_roll_Amax);

    responseDB.stim_angle_vel_trigger(i,1) = stim_angle_vel_trigger;
    responseDB.stim_angle_vel_pre_resp(i,1) = stim_angle_vel_pre_resp;
    responseDB.stim_angle_vel_trig2resp(i,1) = stim_angle_vel_trig2resp;
    responseDB.stim_angle_vel_pre_turn(i,1) = stim_angle_vel_pre_turn;

    responseDB.stim_angle_vel_post_turn(i,1) = stim_angle_vel_post_turn;
    responseDB.stim_angle_vel_post_resp(i,1) = stim_angle_vel_post_resp;   
    responseDB.stim_angle_vel_postresp2end(i,1) = stim_angle_vel_postresp2end;   

    responseDB.stim_angle_accel_trigger(i,1) = stim_angle_accel_trigger;
    responseDB.stim_angle_accel_pre_resp(i,1) = stim_angle_accel_pre_resp;
    responseDB.stim_angle_accel_trig2resp(i,1) = stim_angle_accel_trig2resp;
    responseDB.stim_angle_accel_pre_turn(i,1) = stim_angle_accel_pre_turn;

    responseDB.stim_angle_accel_post_turn(i,1) = stim_angle_accel_post_turn;
    responseDB.stim_angle_accel_post_resp(i,1) = stim_angle_accel_post_resp;   
    responseDB.stim_angle_accel_postresp2end(i,1) = stim_angle_accel_postresp2end;   
    
    responseDB.stim_angle_accel_resp_min(i,1) = stim_angle_accel_resp_min;
    responseDB.stim_angle_accel_resp_max(i,1) = stim_angle_accel_resp_max;
    responseDB.stim_angle_accel_resp_mean(i,1) = stim_angle_accel_resp_mean;
    responseDB.stim_angle_accel_resp_std(i,1) = stim_angle_accel_resp_std;

    responseDB.stim_angle_yaw_trigger(i,1) = stim_angle_yaw_trigger;
    responseDB.stim_angle_yaw_pre_resp(i,1) = stim_angle_yaw_pre_resp;
    responseDB.stim_angle_yaw_trig2resp(i,1) = stim_angle_yaw_trig2resp;
    responseDB.stim_angle_yaw_pre_turn(i,1) = stim_angle_yaw_pre_turn;

    responseDB.stim_angle_yaw_post_turn(i,1) = stim_angle_yaw_post_turn;
    responseDB.stim_angle_yaw_post_resp(i,1) = stim_angle_yaw_post_resp;   
    responseDB.stim_angle_yaw_postresp2end(i,1) = stim_angle_yaw_postresp2end;   
    
    responseDB.accel_angle_hor_vel_pre_resp(i,1) = accel_angle_hor_vel_pre_resp;
    responseDB.accel_angle_hor_vel_pre_turn(i,1) = accel_angle_hor_vel_pre_turn;
    responseDB.accel_angle_hor_vel_pre_accel(i,1) = accel_angle_hor_vel_pre_accel;
        
    responseDB.accel_angle_hor_vel_post_turn(i,1) = accel_angle_hor_vel_post_turn;
    responseDB.accel_angle_hor_vel_post_resp(i,1) = accel_angle_hor_vel_post_resp;
    responseDB.accel_angle_hor_vel_post_accel(i,1) = accel_angle_hor_vel_post_accel;
        
    responseDB.accel_angle_hor_body_pre_resp(i,1) = accel_angle_hor_body_pre_resp;
    responseDB.accel_angle_hor_body_pre_turn(i,1) = accel_angle_hor_body_pre_turn;
    responseDB.accel_angle_hor_body_pre_accel(i,1) = accel_angle_hor_body_pre_accel;
        
    responseDB.accel_angle_hor_body_post_turn(i,1) = accel_angle_hor_body_post_turn;
    responseDB.accel_angle_hor_body_post_resp(i,1) = accel_angle_hor_body_post_resp;
    responseDB.accel_angle_hor_body_post_accel(i,1) = accel_angle_hor_body_post_accel;


    responseDB.accel_angle_hor_vel_resp_mean(i,1) = (accel_angle_hor_vel_resp_mean);
    responseDB.accel_angle_hor_body_resp_mean(i,1) = (accel_angle_hor_body_resp_mean);

    responseDB.slip_trigger(i,1) = slip_trigger;
    responseDB.slip_pre_resp(i,1) = slip_pre_resp;
    responseDB.slip_trig2resp(i,1) = slip_trig2resp;
    responseDB.slip_pre_turn(i,1) = slip_pre_turn;

    responseDB.slip_post_turn(i,1) = slip_post_turn;
    responseDB.slip_post_resp(i,1) = slip_post_resp;   
    responseDB.slip_postresp2end(i,1) = slip_postresp2end;   

    responseDB.pitch_trigger(i,1) = pitch_trigger;
    responseDB.pitch_pre_resp(i,1) = pitch_pre_resp;
    responseDB.pitch_trig2resp(i,1) = pitch_trig2resp;
    responseDB.pitch_pre_turn(i,1) = pitch_pre_turn;

    responseDB.pitch_post_turn(i,1) = pitch_post_turn;
    responseDB.pitch_post_resp(i,1) = pitch_post_resp;   
    responseDB.pitch_postresp2end(i,1) = pitch_postresp2end;   

    responseDB.roll_trigger(i,1) = roll_trigger;
    responseDB.roll_pre_resp(i,1) = roll_pre_resp;
    responseDB.roll_trig2resp(i,1) = roll_trig2resp;
    responseDB.roll_pre_turn(i,1) = roll_pre_turn;

    responseDB.roll_post_turn(i,1) = roll_post_turn;
    responseDB.roll_post_resp(i,1) = roll_post_resp;   
    responseDB.roll_postresp2end(i,1) = roll_postresp2end;   

    responseDB.V_trigger(i,1) = V_trigger;   
    responseDB.V_pre_resp(i,1) = V_pre_resp;   
    responseDB.V_trig2resp(i,1) = V_trig2resp;   
    responseDB.V_pre_accel(i,1) = V_pre_accel;   

    responseDB.V_post_resp(i,1) = V_post_resp;   
    responseDB.V_postresp2end(i,1) = V_postresp2end;   
    responseDB.V_post_accel(i,1) = V_post_accel;

    if toplot == 1
        cd(figdir)
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.fig'])
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.png'])
        cd ..
        % pause
    end
end

save(outputfile,'pathDB','patternDB','settings','responseDB');

    
