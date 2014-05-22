clc
clear
close all

%% const
% % fly data
% Lwing = .00299; %m
% ARwing_fly  = 3.1579; % wing aspect ratio of fly (based on Elzinga's R/c)
% 
% % 8 fem & 7 male, starved for 24h, by FTM. 50%male/female
% Mfly = 1.8e-6;
% Mg_fly = Mfly*9.81;
% 
% % Yaw data from Dickinson et al 2010 & converted to Dh scale
% % roll and pitch damping scaled from yaw acordoing to Cheng et al 2009 (table 2)
% 
% % inertial moment coeff [Nm/s^2] !!! RAD based !!!
% Iyaw   = 1.87e-12;
% Iroll  = Iyaw;
% Ipitch = 2*Iyaw;
% 
% % damping coeff [Nm/s] !!! RAD based !!!
% Cyaw   = -1.14e-10;
% Croll  = 0.55 * Cyaw;
% Cpitch = 0.19 * Cyaw;
% 
% save('flyVars.mat','Lwing','ARwing_fly','Mfly','Mg_fly','Iyaw','Iroll','Ipitch','Cyaw','Croll','Cpitch')
% 

%% load data

% load('flyVars_addedmassSphere.mat')
% load('flyVars_addedmassDisc.mat')
% load('flyVars_addedmassDiscSphere.mat')
var_name = dir('flyVar*');
var_name = var_name.name;
load(var_name);

% loadname = 'flightpathDB_pos_qbodyEKF_INCroll_9clusters_1.4n-1.7n1.9_strokeplane47.5deg_startframe100.mat';
loadname = 'flightpathDB_pos_qbodyEKF_INCroll_strokeplane47.5deg_rotaxis8deg_startframe100.mat';
load(loadname)

savename = ['kin',loadname];

% load  kin
load('Wing_kinematics_data_strokeplane47.5deg.mat')
% load('Wing_kinematics_data.mat')

% load seq info
load('seq_infoDB.mat')

% seq lists
seqs_path = settings.seq;
seqs_kin = nan([length(seq_names),2]);
for i = 1:length(seqs_kin)
    seq_now = seq_names{i};
    seqs_kin(i,1) = str2num(seq_now(end-13:end-6));
    seqs_kin(i,2) = str2num(seq_now(end-3:end));
end

% U
Uwing_L = sqrt(u_wing_L.^2+v_wing_L.^2+w_wing_L.^2);
Uwing_R = sqrt(u_wing_R.^2+v_wing_R.^2+w_wing_R.^2);

% build DBs
settings.seq_kin = nan(size(settings.seq));

kinDB.n_down_start_L = nan(size(down_time_L,1),length(seqs_path));
kinDB.n_up_start_L = nan(size(down_time_L,1),length(seqs_path));

kinDB.wingtip_path_L = nan(size(wingtip_path_L,1),size(wingtip_path_L,2),length(seqs_path));
kinDB.joint_pos_L = nan(size(joint_pos_L,1),length(seqs_path));

kinDB.pitch_L = nan(size(eta_L,1),length(seqs_path));
kinDB.dev_L = nan(size(theta_L,1),length(seqs_path));
kinDB.stroke_L = nan(size(phi_L,1),length(seqs_path));

kinDB.Uwing_L = nan(size(Uwing_L,1),size(Uwing_L,2),length(seqs_path));
kinDB.aoa_L = nan(size(Uwing_L,1),size(Uwing_L,2),length(seqs_path));

kinDB.n_down_start_R = nan(size(down_time_R,1),length(seqs_path));
kinDB.n_up_start_R = nan(size(down_time_R,1),length(seqs_path));

kinDB.wingtip_path_R = nan(size(wingtip_path_R,1),size(wingtip_path_R,2),length(seqs_path));
kinDB.joint_pos_R = nan(size(joint_pos_R,1),length(seqs_path));

kinDB.pitch_R = nan(size(eta_R,1),length(seqs_path));
kinDB.dev_R = nan(size(theta_R,1),length(seqs_path));
kinDB.stroke_R = nan(size(phi_R,1),length(seqs_path));

kinDB.Uwing_R = nan(size(Uwing_L,1),size(Uwing_L,2),length(seqs_path));
kinDB.aoa_R = nan(size(Uwing_L,1),size(Uwing_L,2),length(seqs_path));

kinDB.wing_length = nan(size(wing_length,1),length(seqs_path))';

trigger_frame_global = find(pathDB.t==0);
for seqp = 1:length(kinDB.wing_length)
    length(kinDB.wing_length) - seqp
    if (length(kinDB.wing_length) - seqp) == 61
        er=1
    end
    
    date_now = seqs_path(seqp,1);
    seq_now = seqs_path(seqp,2);
    
    % update trigger frame
    seq_nr = find(seq_info.date == date_now & seq_info.seq == seq_now);
    if isempty(seq_nr) == 0
        trigger_frame_now = seq_info.trigger_frame(seq_nr(1));
    else
        trigger_frame_now = trigger_frame_global;
        'NO TRIGGERFRAME MATCH'
        settings.seq(seqp,:)
%                 pause
    end
    
    trigger_frame_now = trigger_frame_global;
    
    
    % kin seq nr
    seqk = find(seqs_kin(:,1) == date_now & seqs_kin(:,2) == seq_now);
    if isempty(seqk)==0
        settings.seq_kin(seqp,:) = seqs_kin(seqk,:);

        % start of seq
        n_start = find(isnan(eta_L(:,seqk))==0, 1 );
        dn = trigger_frame_global - trigger_frame_now;

        kinDB.n_down_start_L(:,seqp) = down_time_L(:,1,seqk) + n_start-1 + dn;
        kinDB.n_up_start_L(:,seqp) = up_time_L(:,1,seqk) + n_start-1 + dn;

        kinDB.n_down_start_R(:,seqp) = down_time_R(:,1,seqk) + n_start-1 + dn;
        kinDB.n_up_start_R(:,seqp) = up_time_R(:,1,seqk) + n_start-1 + dn;

        kinDB.joint_pos_L(:,seqp) = joint_pos_L(:,seqk);
        kinDB.joint_pos_R(:,seqp) = joint_pos_R(:,seqk);
        kinDB.wing_length(seqp,1) = wing_length(:,seqk);

        for fr = 1:size(eta_L,1)
            if isnan(wingtip_path_L(fr,1,seqk))==0

                kinDB.wingtip_path_L(fr+dn,:,seqp) = wingtip_path_L(fr,:,seqk);
                kinDB.wingtip_path_R(fr+dn,:,seqp) = wingtip_path_R(fr,:,seqk);

                kinDB.pitch_L(fr+dn,seqp) = eta_L(fr,seqk);
                kinDB.dev_L(fr+dn,seqp) = theta_L(fr,seqk);
                kinDB.stroke_L(fr+dn,seqp) = phi_L(fr,seqk);

                kinDB.pitch_R(fr+dn,seqp) = eta_R(fr,seqk);
                kinDB.dev_R(fr+dn,seqp) = theta_R(fr,seqk);
                kinDB.stroke_R(fr+dn,seqp) = phi_R(fr,seqk);
                
                kinDB.Uwing_L(fr+dn,:,seqp) = Uwing_L(fr,:,seqk);
                kinDB.Uwing_R(fr+dn,:,seqp) = Uwing_R(fr,:,seqk);
                
                kinDB.aoa_L(fr+dn,:,seqp) = alfa_L(fr,:,seqk);
                kinDB.aoa_R(fr+dn,:,seqp) = alfa_R(fr,:,seqk);
            end
        end
    end
end


save(savename,'kinDB','pathDB','patternDB','responseDB','settings');
