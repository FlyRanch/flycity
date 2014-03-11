
clc
clear
close all
% addpath('/home/florian/Dropbox/WORK/toolbox')

% filenames
inputfile = 'flightpathDB_pos_qbodyEKF_9clusters.mat';
outputfile = 'flightpathDB_pos_qbodyEKF_9clusters_responseDB.mat';
figdir = 'accel_plots';
mkdir(figdir)

toplot=1;
% toplot=0;

load(inputfile)

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
    stim_angle_trigger,stim_angle_pre_resp,stim_angle_trig2resp,stim_angle_pre_turn,...
    stim_angle_post_turn,stim_angle_post_resp,stim_angle_postresp2end,...
    V_trigger,V_pre_resp,V_trig2resp,V_pre_accel,...
    V_post_resp,V_postresp2end,V_post_accel] = response_times_v3_hVA(pathDB,settings,i,toplot);    


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

responseDB.stim_angle_trigger(i,1) = stim_angle_trigger;
responseDB.stim_angle_pre_resp(i,1) = stim_angle_pre_resp;
responseDB.stim_angle_trig2resp(i,1) = stim_angle_trig2resp;
responseDB.stim_angle_pre_turn(i,1) = stim_angle_pre_turn;

responseDB.stim_angle_post_turn(i,1) = stim_angle_post_turn;
responseDB.stim_angle_post_resp(i,1) = stim_angle_post_resp;   
responseDB.stim_angle_postresp2end(i,1) = stim_angle_postresp2end;   

responseDB.V_trigger(i,1) = V_trigger;   
responseDB.V_pre_resp(i,1) = V_pre_resp;   
responseDB.V_trig2resp(i,1) = V_trig2resp;   
responseDB.V_pre_accel(i,1) = V_pre_accel;   

responseDB.V_post_resp(i,1) = V_post_resp;   
responseDB.V_postresp2end(i,1) = V_postresp2end;   
responseDB.V_post_accel(i,1) = V_post_accel;

cd(figdir)
saveas(gcf,['flightpath_',num2str(i),'.png'])
cd ..

% pause
end

save(outputfile,'pathDB','settings','responseDB');

    