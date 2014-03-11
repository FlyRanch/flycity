%% Bin Generator
warning('off')

varlist = {
'stroke_wb_MATT_bins_up_turn'
'stroke_wb_MATT_bins_dn_turn'
'dev_wb_MATT_bins_up_turn'
'dev_wb_MATT_bins_dn_turn'
'pitch_wb_MATT_bins_up_turn'
'pitch_wb_MATT_bins_dn_turn'
'stroke_wb_R_MATT_bins_steady'
'stroke_wb_L_MATT_bins_steady'
'dev_wb_R_MATT_bins_steady'
'dev_wb_L_MATT_bins_steady'
'pitch_wb_R_MATT_bins_steady'
'pitch_wb_L_MATT_bins_steady'
'stroke_wb_MATT_bins_up_turn'
'stroke_wb_MATT_bins_dn_turn'
'dev_wb_MATT_bins_up_turn'
'dev_wb_MATT_bins_dn_turn'
'pitch_wb_MATT_bins_up_turn'
'pitch_wb_MATT_bins_dn_turn'
'stroke_wb_R_MATT_bins_steady'
'stroke_wb_L_MATT_bins_steady'
'dev_wb_R_MATT_bins_steady'
'dev_wb_L_MATT_bins_steady'
'pitch_wb_R_MATT_bins_steady'
'pitch_wb_L_MATT_bins_steady'
'stroke_wb_R_MATT_bins_mean_steady'
'stroke_wb_L_MATT_bins_mean_steady'
'dev_wb_MATT_bins_mean_up_turn'
'dev_wb_MATT_bins_mean_dn_turn'
'dev_wb_R_MATT_bins_mean_steady'
'dev_wb_L_MATT_bins_mean_steady'
'pitch_wb_MATT_bins_mean_up_turn'
'pitch_wb_MATT_bins_mean_dn_turn'
'pitch_wb_R_MATT_bins_mean_steady'
'pitch_wb_L_MATT_bins_mean_steady'}';

clear(varlist{:})

stim_exp_90 = [30:32,48:50,51:53,64:66,74:76];
stim_exp_270 = [36:41,58:60,70:73];
stim_exp_0 = [1:29,42:44,54:56,61:63,77:79];


%% Bin data and interpolate

bins = 200;
% seqs_max = 83;

% % Allocate space
% stroke_wb_L_MATT_bins = [];
% stroke_wb_R_MATT_bins = [];
% dev_wb_R_MATT_bins = [];
% dev_wb_L_MATT_bins = [];
% pitch_wb_R_MATT_bins = [];
% pitch_wb_L_MATT_bins = [];

% Hack into bins of length 'bins' and interpolate
% for i = 1:seqs
% stroke_wb_L_MATT_bins(:,:,i) = bin(stroke_wb_L_MATT(:,:,i),bins); % input, number of bins
% stroke_wb_R_MATT_bins(:,:,i) = bin(stroke_wb_R_MATT(:,:,i),bins);
% dev_wb_R_MATT_bins(:,:,i) = bin(dev_wb_R_MATT(:,:,i),bins);
% dev_wb_L_MATT_bins(:,:,i) = bin(dev_wb_L_MATT(:,:,i),bins);
% pitch_wb_R_MATT_bins(:,:,i) = bin(pitch_wb_R_MATT(:,:,i),bins);
% pitch_wb_L_MATT_bins(:,:,i) = bin(pitch_wb_L_MATT(:,:,i),bins);
% end

% Determine if turning or not, then which wing it turning, then seperate
for i = [stim_exp_90, stim_exp_270]
    for j = 1:length(stroke_wb_L_MATT(1,:,i))
        if turning_def(j,7,i) > 0.15
           
            if find(stim_exp_90 == i) > 0
            
                if strcmp(wing{j,1,i},'right') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = bin(stroke_wb_R_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_up_turn(:,j,i) = bin(dev_wb_R_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = bin(pitch_wb_R_MATT(:,j,i),bins);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = bin(stroke_wb_L_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = bin(dev_wb_L_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = bin(pitch_wb_L_MATT(:,j,i),bins);
                
                elseif strcmp(wing{j,1,i},'left') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = bin(stroke_wb_L_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_up_turn(:,j,i) = bin(dev_wb_L_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = bin(pitch_wb_L_MATT(:,j,i),bins);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = bin(stroke_wb_R_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = bin(dev_wb_R_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = bin(pitch_wb_R_MATT(:,j,i),bins);
                end
                
            elseif find(stim_exp_270 == i) > 0    
                    
                if strcmp(wing{j,1,i},'left') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = bin(stroke_wb_R_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_up_turn(:,j,i) = bin(dev_wb_R_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = bin(pitch_wb_R_MATT(:,j,i),bins);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = bin(stroke_wb_L_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = bin(dev_wb_L_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = bin(pitch_wb_L_MATT(:,j,i),bins);
                
                elseif strcmp(wing{j,1,i},'right') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = bin(stroke_wb_L_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_up_turn(:,j,i) = bin(dev_wb_L_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = bin(pitch_wb_L_MATT(:,j,i),bins);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = bin(stroke_wb_R_MATT(:,j,i),bins);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = bin(dev_wb_R_MATT(:,j,i),bins);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = bin(pitch_wb_R_MATT(:,j,i),bins);
                end
            end
            
        elseif turning_def(j,7,i) < 0.1
            
            stroke_wb_R_MATT_bins_steady(:,j,i) = bin(stroke_wb_R_MATT(:,j,i),bins);
            stroke_wb_L_MATT_bins_steady(:,j,i) = bin(stroke_wb_L_MATT(:,j,i),bins);
            
            dev_wb_R_MATT_bins_steady(:,j,i) = bin(dev_wb_R_MATT(:,j,i),bins);
            dev_wb_L_MATT_bins_steady(:,j,i) = bin(dev_wb_L_MATT(:,j,i),bins);
            
            pitch_wb_R_MATT_bins_steady(:,j,i) = bin(pitch_wb_R_MATT(:,j,i),bins);
            pitch_wb_L_MATT_bins_steady(:,j,i) = bin(pitch_wb_L_MATT(:,j,i),bins);
            
        end
    end
end

%% Zero to NaN
%%%%%%%%%%%%%%%%%%%%
stroke_wb_MATT_bins_up_turn(stroke_wb_MATT_bins_up_turn == 0) = nan;
stroke_wb_MATT_bins_dn_turn(stroke_wb_MATT_bins_dn_turn == 0) = nan;

dev_wb_MATT_bins_up_turn(dev_wb_MATT_bins_up_turn == 0) = nan;
dev_wb_MATT_bins_dn_turn(dev_wb_MATT_bins_dn_turn == 0) = nan;

pitch_wb_MATT_bins_up_turn(pitch_wb_MATT_bins_up_turn == 0) = nan;
pitch_wb_MATT_bins_dn_turn(pitch_wb_MATT_bins_dn_turn == 0) = nan;
%%%%%%%%%%%%%%%%%%%%%
stroke_wb_R_MATT_bins_steady(stroke_wb_R_MATT_bins_steady == 0) = nan;
stroke_wb_L_MATT_bins_steady(stroke_wb_L_MATT_bins_steady == 0) = nan;

dev_wb_R_MATT_bins_steady(dev_wb_R_MATT_bins_steady == 0) = nan;
dev_wb_L_MATT_bins_steady(dev_wb_L_MATT_bins_steady == 0) = nan;

pitch_wb_R_MATT_bins_steady(pitch_wb_R_MATT_bins_steady == 0) = nan;
pitch_wb_L_MATT_bins_steady(pitch_wb_L_MATT_bins_steady == 0) = nan;
%%%%%%%%%%%%%%%%%%%%%

%% Remove outliers
tol = 2;
stroke_wb_MATT_bins_up_turn = outlier(stroke_wb_MATT_bins_up_turn,tol);
stroke_wb_MATT_bins_dn_turn = outlier(stroke_wb_MATT_bins_dn_turn,tol);

dev_wb_MATT_bins_up_turn = outlier(dev_wb_MATT_bins_up_turn,tol);
dev_wb_MATT_bins_dn_turn = outlier(dev_wb_MATT_bins_dn_turn,tol);

pitch_wb_MATT_bins_up_turn = outlier(pitch_wb_MATT_bins_up_turn,tol);
pitch_wb_MATT_bins_dn_turn = outlier(pitch_wb_MATT_bins_dn_turn,tol);

stroke_wb_R_MATT_bins_steady = outlier(stroke_wb_R_MATT_bins_steady,tol);
stroke_wb_L_MATT_bins_steady = outlier(stroke_wb_L_MATT_bins_steady,tol);

dev_wb_R_MATT_bins_steady = outlier(dev_wb_R_MATT_bins_steady,tol);
dev_wb_L_MATT_bins_steady = outlier(dev_wb_L_MATT_bins_steady,tol);

pitch_wb_R_MATT_bins_steady = outlier(pitch_wb_R_MATT_bins_steady,tol);
pitch_wb_L_MATT_bins_steady = outlier(pitch_wb_L_MATT_bins_steady,tol);

%% RESHAPE
% stroke_wb_MATT_bins_up_turn = reshape(stroke_wb_MATT_bins_up_turn,length(stroke_wb_MATT_bins_up_turn(:,1,1)),(length(stroke_wb_MATT_bins_up_turn(1,:,1))*length(stroke_wb_MATT_bins_up_turn(1,1,:))));
% stroke_wb_MATT_bins_dn_turn = reshape(stroke_wb_MATT_bins_dn_turn,length(stroke_wb_MATT_bins_dn_turn(:,1,1)),(length(stroke_wb_MATT_bins_dn_turn(1,:,1))*length(stroke_wb_MATT_bins_dn_turn(1,1,:))));
stroke_wb_MATT_bins_up_turn_num = num_winbeats(stroke_wb_MATT_bins_up_turn);
stroke_wb_MATT_bins_dn_turn_num = num_winbeats(stroke_wb_MATT_bins_dn_turn);

% stroke_wb_R_MATT_bins_steady = reshape(stroke_wb_R_MATT_bins_steady,length(stroke_wb_R_MATT_bins_steady(:,1,1)),(length(stroke_wb_R_MATT_bins_steady(1,:,1))*length(stroke_wb_R_MATT_bins_steady(1,1,:))));
% stroke_wb_L_MATT_bins_steady = reshape(stroke_wb_L_MATT_bins_steady,length(stroke_wb_L_MATT_bins_steady(:,1,1)),(length(stroke_wb_L_MATT_bins_steady(1,:,1))*length(stroke_wb_L_MATT_bins_steady(1,1,:))));
stroke_wb_MATT_bins_steady = [stroke_wb_R_MATT_bins_steady,stroke_wb_L_MATT_bins_steady];
stroke_wb_MATT_bins_steady_num = num_winbeats(stroke_wb_MATT_bins_steady);

% dev_wb_MATT_bins_up_turn = reshape(dev_wb_MATT_bins_up_turn,length(dev_wb_MATT_bins_up_turn(:,1,1)),(length(dev_wb_MATT_bins_up_turn(1,:,1))*length(dev_wb_MATT_bins_up_turn(1,1,:))));
% dev_wb_MATT_bins_dn_turn = reshape(dev_wb_MATT_bins_dn_turn,length(dev_wb_MATT_bins_dn_turn(:,1,1)),(length(dev_wb_MATT_bins_dn_turn(1,:,1))*length(dev_wb_MATT_bins_dn_turn(1,1,:))));
dev_wb_MATT_bins_up_turn_num = num_winbeats(dev_wb_MATT_bins_up_turn);
dev_wb_MATT_bins_dn_turn_num = num_winbeats(dev_wb_MATT_bins_dn_turn);

% dev_wb_R_MATT_bins_steady = reshape(dev_wb_R_MATT_bins_steady,length(dev_wb_R_MATT_bins_steady(:,1,1)),(length(dev_wb_R_MATT_bins_steady(1,:,1))*length(dev_wb_R_MATT_bins_steady(1,1,:))));
% dev_wb_L_MATT_bins_steady = reshape(dev_wb_L_MATT_bins_steady,length(dev_wb_L_MATT_bins_steady(:,1,1)),(length(dev_wb_L_MATT_bins_steady(1,:,1))*length(dev_wb_L_MATT_bins_steady(1,1,:))));
dev_wb_MATT_bins_steady = [dev_wb_R_MATT_bins_steady,dev_wb_L_MATT_bins_steady];
dev_wb_MATT_bins_steady_num = num_winbeats(dev_wb_MATT_bins_steady);

% pitch_wb_MATT_bins_up_turn = reshape(pitch_wb_MATT_bins_up_turn,length(pitch_wb_MATT_bins_up_turn(:,1,1)),(length(pitch_wb_MATT_bins_up_turn(1,:,1))*length(pitch_wb_MATT_bins_up_turn(1,1,:))));
% pitch_wb_MATT_bins_dn_turn = reshape(pitch_wb_MATT_bins_dn_turn,length(pitch_wb_MATT_bins_dn_turn(:,1,1)),(length(pitch_wb_MATT_bins_dn_turn(1,:,1))*length(pitch_wb_MATT_bins_dn_turn(1,1,:))));
pitch_wb_MATT_bins_up_turn_num = num_winbeats(pitch_wb_MATT_bins_up_turn);
pitch_wb_MATT_bins_dn_turn_num = num_winbeats(pitch_wb_MATT_bins_dn_turn);

% pitch_wb_R_MATT_bins_steady = reshape(pitch_wb_R_MATT_bins_steady,length(pitch_wb_R_MATT_bins_steady(:,1,1)),(length(pitch_wb_R_MATT_bins_steady(1,:,1))*length(pitch_wb_R_MATT_bins_steady(1,1,:))));
% pitch_wb_L_MATT_bins_steady = reshape(pitch_wb_L_MATT_bins_steady,length(pitch_wb_L_MATT_bins_steady(:,1,1)),(length(pitch_wb_L_MATT_bins_steady(1,:,1))*length(pitch_wb_L_MATT_bins_steady(1,1,:))));
pitch_wb_MATT_bins_steady = [pitch_wb_R_MATT_bins_steady,pitch_wb_L_MATT_bins_steady];
pitch_wb_MATT_bins_steady_num = num_winbeats(pitch_wb_MATT_bins_steady);

%% MEAN
for i = 1:bins
% turning
stroke_wb_MATT_bins_mean_up_turn(i,1:3) = circ_mean_deg_nonan(stroke_wb_MATT_bins_up_turn(i,:)');
stroke_wb_MATT_bins_mean_dn_turn(i,1:3) = circ_mean_deg_nonan(stroke_wb_MATT_bins_dn_turn(i,:)');

% steady
stroke_wb_MATT_bins_mean_steady(i,1:3) = circ_mean_deg_nonan(stroke_wb_MATT_bins_steady(i,:)');

% turning
dev_wb_MATT_bins_mean_up_turn(i,1:3) = circ_mean_deg_nonan(dev_wb_MATT_bins_up_turn(i,:)');
dev_wb_MATT_bins_mean_dn_turn(i,1:3) = circ_mean_deg_nonan(dev_wb_MATT_bins_dn_turn(i,:)');

% steady
dev_wb_MATT_bins_mean_steady(i,1:3) = circ_mean_deg_nonan(dev_wb_R_MATT_bins_steady(i,:)');

% turning
pitch_wb_MATT_bins_mean_up_turn(i,1:3) = circ_mean_deg_nonan(pitch_wb_MATT_bins_up_turn(i,:)');
pitch_wb_MATT_bins_mean_dn_turn(i,1:3) = circ_mean_deg_nonan(pitch_wb_MATT_bins_dn_turn(i,:)');

% steady
pitch_wb_MATT_bins_mean_steady(i,1:3) = circ_mean_deg_nonan(pitch_wb_R_MATT_bins_steady(i,:)');
end


% figure(1)
% hold on
% plot(stroke_wb_MATT_bins_dn_turn,'-r')
% plot(stroke_wb_MATT_bins_up_turn,'-b')

figure(1)
subplot(1,3,1)
hold on
% plot(stroke_wb_L_MATT_bins_mean_steady(:,1:3),'-','color','blue')
% plot(stroke_wb_R_MATT_bins_mean_steady(:,1:3),'--','color','blue')
plot(stroke_wb_MATT_bins_mean_up_turn(:,1),'-','color','red','linewidth',2)
plot(stroke_wb_MATT_bins_mean_dn_turn(:,1),'-','color','blue','linewidth',2)
plot(stroke_wb_MATT_bins_mean_steady(:,1:3),'-','color','black')
% plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','black')

plot(stroke_wb_MATT_bins_mean_up_turn(:,2:3),'--','color','red')

plot(stroke_wb_MATT_bins_mean_dn_turn(:,1),'-','color','blue','linewidth',2)
plot(stroke_wb_MATT_bins_mean_dn_turn(:,2:3),'--','color','blue')
% plot(stroke_wb_R_MATT_bins_mean_steady(:,1:3),'-','color','blue')
% plot(stroke_wb_L_MATT_bins_mean_steady(:,1:3),'-','color','blue')
% plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','black')

legend('up wing','down wing', 'steady - free');

hold off

subplot(1,3,2)
hold on
% plot(dev_wb_L_MATT_bins_mean_steady(:,1:3),'-','color','blue')
% plot(dev_wb_R_MATT_bins_mean_steady(:,1:3),'--','color','blue')
plot(dev_wb_MATT_bins_mean_up_turn(:,1),'-','color','red','linewidth',2)
plot(dev_wb_MATT_bins_mean_up_turn(:,2:3),'--','color','red')

plot(dev_wb_MATT_bins_mean_dn_turn(:,1),'-','color','blue','linewidth',2)
plot(dev_wb_MATT_bins_mean_dn_turn(:,2:3),'--','color','blue')

plot(dev_wb_MATT_bins_mean_steady(:,1:3),'-','color','black')
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','black')
hold off

subplot(1,3,3)
hold on
% plot(pitch_wb_L_MATT_bins_mean_steady(:,1:3),'-','color','blue')
% plot(pitch_wb_R_MATT_bins_mean_steady(:,1:3),'--','color','blue')
plot(pitch_wb_MATT_bins_mean_up_turn(:,1),'-','color','red','linewidth',2)
plot(pitch_wb_MATT_bins_mean_up_turn(:,2:3),'--','color','red')

plot(pitch_wb_MATT_bins_mean_dn_turn(:,1),'-','color','blue','linewidth',2)
plot(pitch_wb_MATT_bins_mean_dn_turn(:,2:3),'--','color','blue')

plot(pitch_wb_MATT_bins_mean_steady(:,1:3),'-','color','black')
% plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-','color','black')
hold off


