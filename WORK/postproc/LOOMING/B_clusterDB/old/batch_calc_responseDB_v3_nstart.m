
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

    [IDX_shift,t_shift,t_resp,...
    t_turn_start,t_turn_stop,t_turn_max,...
    t_accel_start,t_accel_stop,t_accel_max,...
    t_decel_start,t_decel_stop,t_decel_max,...
    dt_turn,dt_accel,dt_decel,...
    An_hor_max,At_hor_max,At_hor_min,...
    stim_angle_vel_trigger,stim_angle_vel_pre_resp,stim_angle_vel_trig2resp,stim_angle_vel_pre_turn,...
    stim_angle_vel_post_turn,stim_angle_vel_post_resp,stim_angle_vel_postresp2end,...
    stim_angle_yaw_trigger,stim_angle_yaw_pre_resp,stim_angle_yaw_trig2resp,stim_angle_yaw_pre_turn,...
    stim_angle_yaw_post_turn,stim_angle_yaw_post_resp,stim_angle_yaw_postresp2end,...
    slip_trigger,slip_pre_resp,slip_trig2resp,slip_pre_turn,...
    slip_post_turn,slip_post_resp,slip_postresp2end,...
    V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
    V_post_resp,V_postresp2end,V_post_accel,...
    teta_pre_resp,teta_post_resp,teta_pre_turn,teta_post_turn,teta_pre_accel,teta_post_accel] =...
    calc_response_data_hVA_absAn(pathDB,patternDB,settings,i,toplot,skip,nstart);

    responseDB.IDX_shift(1:length(t_shift),i) = IDX_shift;
    responseDB.t_shift(1:length(t_shift),i) = t_shift;

    responseDB.t_resp(i,1) = t_resp;

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

    responseDB.An_hor_max(i,1) = An_hor_max;
    responseDB.At_hor_max(i,1) = At_hor_max;
    responseDB.At_hor_min(i,1) = At_hor_min;

    responseDB.stim_angle_vel_trigger(i,1) = stim_angle_vel_trigger;
    responseDB.stim_angle_vel_pre_resp(i,1) = stim_angle_vel_pre_resp;
    responseDB.stim_angle_vel_trig2resp(i,1) = stim_angle_vel_trig2resp;
    responseDB.stim_angle_vel_pre_turn(i,1) = stim_angle_vel_pre_turn;

    responseDB.stim_angle_vel_post_turn(i,1) = stim_angle_vel_post_turn;
    responseDB.stim_angle_vel_post_resp(i,1) = stim_angle_vel_post_resp;   
    responseDB.stim_angle_vel_postresp2end(i,1) = stim_angle_vel_postresp2end;   

    responseDB.stim_angle_yaw_trigger(i,1) = stim_angle_yaw_trigger;
    responseDB.stim_angle_yaw_pre_resp(i,1) = stim_angle_yaw_pre_resp;
    responseDB.stim_angle_yaw_trig2resp(i,1) = stim_angle_yaw_trig2resp;
    responseDB.stim_angle_yaw_pre_turn(i,1) = stim_angle_yaw_pre_turn;

    responseDB.stim_angle_yaw_post_turn(i,1) = stim_angle_yaw_post_turn;
    responseDB.stim_angle_yaw_post_resp(i,1) = stim_angle_yaw_post_resp;   
    responseDB.stim_angle_yaw_postresp2end(i,1) = stim_angle_yaw_postresp2end;   

    responseDB.slip_trigger(i,1) = slip_trigger;
    responseDB.slip_pre_resp(i,1) = slip_pre_resp;
    responseDB.slip_trig2resp(i,1) = slip_trig2resp;
    responseDB.slip_pre_turn(i,1) = slip_pre_turn;

    responseDB.slip_post_turn(i,1) = slip_post_turn;
    responseDB.slip_post_resp(i,1) = slip_post_resp;   
    responseDB.slip_postresp2end(i,1) = slip_postresp2end;   

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

    
