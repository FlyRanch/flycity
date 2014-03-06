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

seqs_path = settings.seq;

seqs_kin = nan([length(seq_names),2]);
for i = 1:length(seqs_kin)
    seq_now = seq_names{i};
    seqs_kin(i,1) = str2num(seq_now(1:8));
    seqs_kin(i,2) = str2num(seq_now(11:14));
end

settings.seq_kin = nan(size(settings.seq));

kinDB.down_time_L = nan(size(down_time_L,1),size(down_time_L,2),length(seqs_path));
kinDB.up_time_L = nan(size(down_time_L,1),size(down_time_L,2),length(seqs_path));
kinDB.wingtip_path_L = nan(size(wingtip_path_L,1),size(wingtip_path_L,2),length(seqs_path));

kinDB.joint_pos_L = nan(size(joint_pos_L,1),length(seqs_path));
kinDB.pitch_L = nan(size(eta_L,1),length(seqs_path));
kinDB.dev_L = nan(size(theta_L,1),length(seqs_path));
kinDB.stroke_L = nan(size(phi_L,1),length(seqs_path));

kinDB.down_time_R = nan(size(down_time_R,1),size(down_time_R,2),length(seqs_path));
kinDB.up_time_R = nan(size(down_time_R,1),size(down_time_R,2),length(seqs_path));
kinDB.wingtip_path_R = nan(size(wingtip_path_R,1),size(wingtip_path_R,2),length(seqs_path));

kinDB.joint_pos_R = nan(size(joint_pos_R,1),length(seqs_path));
kinDB.pitch_R = nan(size(eta_R,1),length(seqs_path));
kinDB.dev_R = nan(size(theta_R,1),length(seqs_path));
kinDB.stroke_R = nan(size(phi_R,1),length(seqs_path));

kinDB.wing_length = nan(size(wing_length,1),length(seqs_path));
kinDB.wingbeat_time = nan(size(wingbeat_time,1),size(wingbeat_time,2),length(seqs_path));


for i = 1:length(seqs_path)
    
    date_now = seqs_path(i,1);
    seq_now = seqs_path(i,2);
    
    for j = 1:length(seqs_kin)
        if seqs_kin(j,1) == date_now && seqs_kin(j,2) == seq_now
            
            settings.seq_kin(i,:) = seqs_kin(j,:);
            
            kinDB.down_time_L(:,:,i) = down_time_L(:,:,j);
            kinDB.up_time_L(:,:,i) = up_time_L(:,:,j);
            kinDB.wingtip_path_L(:,:,i) = wingtip_path_L(:,:,j);
            
            kinDB.joint_pos_L(:,i) = joint_pos_L(:,j);
            kinDB.pitch_L(:,i) = eta_L(:,j);
            kinDB.dev_L(:,i) = theta_L(:,j);
            kinDB.stroke_L(:,i) = phi_L(:,j);

            kinDB.down_time_R(:,:,i) = down_time_R(:,:,j);
            kinDB.up_time_R(:,:,i) = up_time_R(:,:,j);
            kinDB.wingtip_path_R(:,:,i) = wingtip_path_R(:,:,j);
            
            kinDB.joint_pos_R(:,i) = joint_pos_R(:,j);
            kinDB.pitch_R(:,i) = eta_R(:,j);
            kinDB.dev_R(:,i) = theta_R(:,j);
            kinDB.stroke_R(:,i) = phi_R(:,j);

            kinDB.wing_length(:,i) = wing_length(:,j);
            kinDB.wingbeat_time(:,:,i) = wingbeat_time(:,:,j);
        end
    end
end

% N down&upstroke
n_start = nan(size(n_pre));

n_down_start_L = nan(size(down_time_L,3),size(down_time_L,1));
n_up_start_L = n_down_start_L;
n_down_start_R = n_down_start_L;
n_up_start_R = n_down_start_L;

for i = 1:length(n_pre)
    if isnan(min(dev_L(:,i)))==0
        n_start(i,1) = find(isnan(dev_L(:,i))==0, 1 );
        n_start2(i,1) = find(isnan(roll(:,i))==0, 1 );

        n_stop(i,1) = find(isnan(dev_L(:,i))==0, 1, 'last' );
        n_stop2(i,1) = find(isnan(roll(:,i))==0, 1, 'last' );

        n_down_start_L(i,:) = down_time_L(:,1,i) + n_start(i,1);
        n_up_start_L(i,:) = up_time_L(:,1,i) + n_start(i,1);

        n_down_start_R(i,:) = down_time_R(:,1,i) + n_start(i,1);
        n_up_start_R(i,:) = up_time_R(:,1,i) + n_start(i,1);
    end
        
end




save(savename,'kinDB','pathDB','patternDB','responseDB','settings');
