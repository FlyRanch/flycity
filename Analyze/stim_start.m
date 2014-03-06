clear all
load('loadvariables_6_3_2013.mat');
load('trims.mat');

seqs = length(pathDB.x(1,:));
seqs_plot = fliplr(1:81);

% data normalized to stimulus start
phi_R_stim = nan(length(frames),seqs);
phi_L_stim = nan(length(frames),seqs);
theta_R_stim = nan(length(frames),seqs);
theta_L_stim = nan(length(frames),seqs);
eta_R_stim= nan(length(frames),seqs);
eta_L_stim = nan(length(frames),seqs);

for i = 1:seqs
    phi_R_stim(1:length(pathDB.phi_R(trims(i,6):end,i)),i) = pathDB.phi_R(trims(i,6):end,i);
    phi_L_stim(1:length(pathDB.phi_L(trims(i,6):end,i)),i) = pathDB.phi_L(trims(i,6):end,i);
    theta_R_stim(1:length(pathDB.theta_R(trims(i,6):end,i)),i) = pathDB.theta_R(trims(i,6):end,i);
    theta_L_stim(1:length(pathDB.theta_L(trims(i,6):end,i)),i) = pathDB.theta_L(trims(i,6):end,i);
    eta_R_stim(1:length(pathDB.eta_R(trims(i,6):end,i)),i) = pathDB.eta_R(trims(i,6):end,i);
    eta_L_stim(1:length(pathDB.eta_L(trims(i,6):end,i)),i) = pathDB.eta_L(trims(i,6):end,i);
end


wingbeats_stim_start

global i

phi_R_stim_trim = nan(frames,20,seqs);


for i = 1
    A = phi_R_stim_start(:,i);
%     A = A(isnan(A)==0);
    
    B = trim2wingbeat_stim(A);
%     phi_R_stim_trim(1:length(B(:,1)),1:length(B(1,:)),i) = B; 
    
end

h = 1;
g = 41;
for i = 1:10
    figure(1)
    hold on
    x = [h:g];
    plot(x,B(:,i))
    h = h+length(B(:,i));
    g = g+length(B(:,i));
end
    

phi_R_stim_trim(:,:,1);





% bins = 200;
% 
% for i = 1:seqs
%     
% phi_R_stim(1:length(,j,i) = bin(phi_R_stim_start(:,i),bins);
% 
% end




phi_L_stim_start = nan(length(frames),seqs);
theta_R_stim_start = nan(length(frames),seqs);
theta_L_stim_start = nan(length(frames),seqs);
eta_R_stim_start = nan(length(frames),seqs);
eta_L_stim_start = nan(length(frames),seqs);
