%% Determine parameters from wingbeat data

% Grab parameters of interst, *note: phi_R is already defined in wingbeats.m

% phi_L = pathDB_grabber(pathDB.phi_L); % comment out if smooth peaks is
% enabled
theta_R = pathDB_grabber(pathDB.theta_R);
theta_L = pathDB_grabber(pathDB.theta_L);
eta_R = pathDB_grabber(pathDB.eta_R);
eta_L = pathDB_grabber(pathDB.eta_L);

% Grab parameters of interest and normalize to stimulus start

% phi_R_stim = nan(length(frames),seqs);
% phi_L_stim = nan(length(frames),seqs);
% theta_R_stim = nan(length(frames),seqs);
% theta_L_stim = nan(length(frames),seqs);
% eta_R_stim= nan(length(frames),seqs);
% eta_L_stim = nan(length(frames),seqs);

% for i = 1:seqs
% %     phi_R_stim(1:length(pathDB.phi_R(trims(i,6):end,i)),i) = pathDB.phi_R(trims(i,6):end,i);
% %     phi_L_stim(1:length(pathDB.phi_L(trims(i,6):end,i)),i) = pathDB.phi_L(trims(i,6):end,i);
%     theta_R_stim(1:length(pathDB.theta_R(trims(i,6):end,i)),i) = pathDB.theta_R(trims(i,6):end,i);
%     theta_L_stim(1:length(pathDB.theta_L(trims(i,6):end,i)),i) = pathDB.theta_L(trims(i,6):end,i);
%     eta_R_stim(1:length(pathDB.eta_R(trims(i,6):end,i)),i) = pathDB.eta_R(trims(i,6):end,i);
%     eta_L_stim(1:length(pathDB.eta_L(trims(i,6):end,i)),i) = pathDB.eta_L(trims(i,6):end,i);
% end

% compare
% figure
% hold on
% plot(phi_R(:,1))
% plot(linspace(200-offset(1),200-offset(1)+length(phi_R_stim(:,1)),length(phi_R_stim(:,1))),phi_R_stim(:,1),'-r')


% trim variable to wingbeat lengths

stroke_wb_R_MATT = rad2deg(trim2wingbeat_wb(phi_R));
stroke_wb_L_MATT = rad2deg(trim2wingbeat_wb(phi_L));
dev_wb_R_MATT = rad2deg(trim2wingbeat_wb(theta_R));
dev_wb_L_MATT = rad2deg(trim2wingbeat_wb(theta_L));
pitch_wb_R_MATT = rad2deg(trim2wingbeat_wb(eta_R));
pitch_wb_L_MATT = rad2deg(trim2wingbeat_wb(eta_L));

stroke_wb_R_MATT(51:end,:,:) = [];
stroke_wb_L_MATT(51:end,:,:) = [];
dev_wb_R_MATT(51:end,:,:) = [];
dev_wb_L_MATT(51:end,:,:) = [];
pitch_wb_R_MATT(51:end,:,:) = [];
pitch_wb_L_MATT(51:end,:,:) = [];

% trim variable to wingbeat lengths, with offset for stimulus
% stroke_wb_R_MATT_stim = rad2deg(trim2wingbeat_stim(phi_R_stim));
% stroke_wb_L_MATT_stim = rad2deg(trim2wingbeat_stim(phi_L_stim));
% dev_wb_R_MATT_stim = rad2deg(trim2wingbeat_stim(theta_R_stim));
% dev_wb_L_MATT_stim = rad2deg(trim2wingbeat_stim(theta_L_stim));
% pitch_wb_R_MATT_stim = rad2deg(trim2wingbeat_stim(eta_R_stim));
% pitch_wb_L_MATT_stim = rad2deg(trim2wingbeat_stim(eta_L_stim));

% % for k = 1:length(dev_wb_L_MATT_stim_bins(1,1,:))
% for i = 1:12 % length(dev_wb_L_MATT_stim_bins(1,:,1))
%     
%     B = rand(1);
%     C = rand(1);
%     D = rand(1);
%     
%     
%     
%     for j = 1:81 % length(dev_wb_L_MATT_stim_bins(:,i,1))
%         
%         A = dev_wb_L_MATT_stim_bins(:,i,j);
%         
%     figure(1)
%     hold on
%     plot(A,'-','color',[B,C,D])
%     
%     end
%     
%     pause
% 
% end