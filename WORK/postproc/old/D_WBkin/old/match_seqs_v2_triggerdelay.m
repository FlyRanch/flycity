clc
clear
close all

% load body kin
loadname = 'flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3.mat';
if exist(loadname) ~= 2
    loadname = 'flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3.mat';
end
load(loadname)

savename = ['kin',loadname];

% load  kin
load('Wing_kinematics_data_Johan.mat')

% load seq info
load('seq_infoDB.mat')

% seq lists
seqs_path = settings.seq;
seqs_kin = nan([length(seq_names),2]);
for i = 1:length(seqs_kin)
    seq_now = seq_names{i};
    seqs_kin(i,1) = str2num(seq_now(1:8));
    seqs_kin(i,2) = str2num(seq_now(11:14));
end

% build DBs
settings.seq_kin = nan(size(settings.seq));

kinDB.n_down_start_L = nan(size(down_time_L,1),length(seqs_path));
kinDB.n_up_start_L = nan(size(down_time_L,1),length(seqs_path));

kinDB.wingtip_path_L = nan(size(wingtip_path_L,1),size(wingtip_path_L,2),length(seqs_path));
kinDB.joint_pos_L = nan(size(joint_pos_L,1),length(seqs_path));

kinDB.pitch_L = nan(size(eta_L,1),length(seqs_path));
kinDB.dev_L = nan(size(theta_L,1),length(seqs_path));
kinDB.stroke_L = nan(size(phi_L,1),length(seqs_path));

kinDB.n_down_start_R = nan(size(down_time_R,1),length(seqs_path));
kinDB.n_up_start_R = nan(size(down_time_R,1),length(seqs_path));

kinDB.wingtip_path_R = nan(size(wingtip_path_R,1),size(wingtip_path_R,2),length(seqs_path));
kinDB.joint_pos_R = nan(size(joint_pos_R,1),length(seqs_path));

kinDB.pitch_R = nan(size(eta_R,1),length(seqs_path));
kinDB.dev_R = nan(size(theta_R,1),length(seqs_path));
kinDB.stroke_R = nan(size(phi_R,1),length(seqs_path));

kinDB.wing_length = nan(size(wing_length,1),length(seqs_path))';

trigger_frame_global = find(pathDB.t==0);
for seqp = 1:length(kinDB.wing_length)
    length(kinDB.wing_length) - seqp
    
    date_now = seqs_path(seqp,1);
    seq_now = seqs_path(seqp,2);
    
    % update trigger frame
    seq_nr = find(seq_info.date == date_now & seq_info.seq == seq_now);
    if isempty(seq_nr) == 0
        trigger_frame_now = seq_info.trigger_frame(seq_nr(1));
    else
        trigger_frame_now = trigger_frame_global;
        'NO MATCH'
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
            end
        end
    end
end


save(savename,'kinDB','pathDB','patternDB','responseDB','settings');