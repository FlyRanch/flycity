
% clc
% clear
% close all
% % addpath('/home/florian/Dropbox/WORK/toolbox')
% 
% % filenames
% % inputfile = 'flightpathDB_2flies_9clusters_Atmin.mat';
% % skip = 1
% inputfile = 'flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3.mat';
% % inputfile = 'flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3.mat';
% outputfile = [inputfile(1:end-4),'_response.mat'];
% load(inputfile)
% 
% nstart = 175;
% trigger_frame = find(pathDB.t == min(abs(pathDB.t)));
% start_frame = trigger_frame + nstart;
% 
% toplot=1;
% toplot=0;
% if toplot == 1
%     skip = 10
%     cmap_k = settings.cmap_k;
%     figdir = 'accel_plots_2.5n-3n2.5';
%     mkdir(figdir)
% end

responseDB.IDX_shift = nan(10,size(pathDB.pos,2));
responseDB.t_shift = nan(10,size(pathDB.pos,2));

for i = 1:size(pathDB.pos,2)
    size(pathDB.pos,2)-i
    
t = pathDB.t;

% teta = patternDB.teta(:,i);
stim_angle_vel = pathDB.stim_angle_vel(:,i);
stim_angle_accel = pathDB.stim_angle_accel(:,i);

slip = pathDB.slip(:,i);
pitch = pathDB.pitch(:,i);
roll = pathDB.roll(:,i);

V = pathDB.V(:,i);
An_hor = pathDB.An_hor(:,i);
At_hor = pathDB.At_hor(:,i);

IDX = pathDB.IDX(:,i);

%     calc_response_data_F;
    calc_response_data_n_CLIP;
    
    responseDB.IDX_shift(1:length(n_shift),i) = IDX_shift;
    responseDB.n_shift(1:length(n_shift),i) = n_shift;

    responseDB.n_first(i,1) = n_first;
    responseDB.n_resp(i,1) = n_resp;
    responseDB.n_resp_end(i,1) = n_resp_end;
    responseDB.n_steady_end(i,1) = n_steady_end;

    responseDB.n_turn_start(i,1) = n_turn_start;
    responseDB.n_turn_stop(i,1) = n_turn_stop;
    responseDB.n_turn_max(i,1) = n_turn_max;

    responseDB.n_accel_start(i,1) = n_accel_start;
    responseDB.n_accel_stop(i,1) = n_accel_stop;
    responseDB.n_accel_max(i,1) = n_accel_max; 

    responseDB.n_decel_start(i,1) = n_decel_start;
    responseDB.n_decel_stop(i,1) = n_decel_stop;
    responseDB.n_decel_min(i,1) = n_decel_min;

    if toplot == 1
        plot_flighttracks_clusters_separate_n_CLIP
        cd(figdir)
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.fig'])
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.png'])
        cd ..
%         pause
    end
end
% 
% save(outputfile,'pathDB','patternDB','settings','responseDB');

    
