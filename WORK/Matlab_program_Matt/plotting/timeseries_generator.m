%% Construct Bins by Individual

stim_exp_90 = [30:32,48:50,51:53,64:66,74:76];
stim_exp_270 = [36:41,58:60,70:73];
stim_exp_0 = [11:29,42:44,54:56,61:63,77:79];
stim_exp_180 = [33:35,45:47,57,67:70];

stim_drum_1hz = [85,93];
stim_drum_5hz = [87,95,97];
stim_drum_10hz = [89,91,99,102,103,104,111,112];


%% Bin data and interpolate

bins = 200;
% seqs_max = 83;

% Hack into bins of length 'bins' and interpolate
for i = 1:seqs;
    stroke_wb_L_MATT_bins(:,:,i) = bin(stroke_wb_L_MATT(:,:,i),bins); % input, number of bins
    stroke_wb_R_MATT_bins(:,:,i) = bin(stroke_wb_R_MATT(:,:,i),bins);
    dev_wb_R_MATT_bins(:,:,i) = bin(dev_wb_R_MATT(:,:,i),bins);
    dev_wb_L_MATT_bins(:,:,i) = bin(dev_wb_L_MATT(:,:,i),bins);
    pitch_wb_R_MATT_bins(:,:,i) = bin(pitch_wb_R_MATT(:,:,i),bins);
    pitch_wb_L_MATT_bins(:,:,i) = bin(pitch_wb_L_MATT(:,:,i),bins);
end

stroke_wb_L_MATT_bins(stroke_wb_L_MATT_bins == 0) = nan;
stroke_wb_R_MATT_bins(stroke_wb_R_MATT_bins == 0) = nan;
dev_wb_R_MATT_bins(dev_wb_R_MATT_bins == 0) = nan;
dev_wb_L_MATT_bins(dev_wb_L_MATT_bins == 0) = nan;
pitch_wb_R_MATT_bins(pitch_wb_R_MATT_bins == 0) = nan;
pitch_wb_L_MATT_bins(pitch_wb_L_MATT_bins == 0) = nan;


% for i = 1:10
%     figure(i)
%     hold on
% for j = 1:10
% plot(dev_wb_R_MATT_bins(:,j,i),'-r')
% plot(dev_wb_L_MATT_bins(:,j,i),'-b')
% end
% end

for i = stim_drum_10hz;
    for j = 1:length(stroke_wb_L_MATT_bins(1,:,i))
        if turning_def(j,7,i) > 0.15
            
                if strcmp(wing{j,1,i},'right') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_L_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_L_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_L_MATT_bins(:,j,i);
                
                elseif strcmp(wing{j,1,i},'left') == 1
                    stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_L_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_L_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_L_MATT_bins(:,j,i);
                
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
                end
                
        elseif turning_def(j,7,i) < 0.15
            
                    stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
                    
                    stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
                    dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
                    pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
                    
        end
    end
end

stroke_wb_MATT_bins_up_turn(stroke_wb_MATT_bins_up_turn == 0) = nan;
dev_wb_MATT_bins_up_turn(dev_wb_MATT_bins_up_turn == 0) = nan;
pitch_wb_MATT_bins_up_turn(pitch_wb_MATT_bins_up_turn == 0) = nan;
stroke_wb_MATT_bins_dn_turn(stroke_wb_MATT_bins_dn_turn == 0) = nan;
dev_wb_MATT_bins_dn_turn(dev_wb_MATT_bins_dn_turn == 0) = nan;
pitch_wb_MATT_bins_dn_turn(pitch_wb_MATT_bins_dn_turn == 0) = nan;



%% plot all parameters in seqs
clear A B C D E F



for i = 1:length(stroke_wb_MATT_bins_up_turn(1,:,1)) % wingbeats

for j = 1:length(stroke_wb_MATT_bins_up_turn(1,1,:)) % individual
                    
    A(:,j) = stroke_wb_MATT_bins_up_turn(:,i,j);
    B(:,j) = dev_wb_MATT_bins_up_turn(:,i,j);
    C(:,j) = pitch_wb_MATT_bins_up_turn(:,i,j);
    
    D(:,j) = stroke_wb_MATT_bins_dn_turn(:,i,j);
    E(:,j) = dev_wb_MATT_bins_dn_turn(:,i,j);
    F(:,j) = pitch_wb_MATT_bins_dn_turn(:,i,j);
    
end

count = 0;
for k = 1:length(A(:,1));
    count = count+1;

    stroke_up(count,1) = nanmean(A(k,:));
    stroke_up(count,2) = (stroke_up(count,1)+((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
    stroke_up(count,3) = (stroke_up(count,1)-((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
    
    dev_up(count,1) = nanmean(B(k,:));
    dev_up(count,2) = (dev_up(count,1)+((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
    dev_up(count,3) = (dev_up(count,1)-((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
    
    pitch_up(count,1) = nanmean(C(k,:));
    pitch_up(count,2) = (pitch_up(count,1)+((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
    pitch_up(count,3) = (pitch_up(count,1)-((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
    
    stroke_dn(count,1) = nanmean(D(k,:));
    stroke_dn(count,2) = (stroke_dn(count,1)+((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
    stroke_dn(count,3) = (stroke_dn(count,1)-((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
    
    dev_dn(count,1) = nanmean(E(k,:));
    dev_dn(count,2) = (dev_dn(count,1)+((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
    dev_dn(count,3) = (dev_dn(count,1)-((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
    
    pitch_dn(count,1) = nanmean(F(k,:));
    pitch_dn(count,2) = (pitch_dn(count,1)+((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
    pitch_dn(count,3) = (pitch_dn(count,1)-((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
    
end

stroke_up_wingbeats(:,1:3,i) = stroke_up;
stroke_dn_wingbeats(:,1:3,i) = stroke_dn;
dev_up_wingbeats(:,1:3,i) = dev_up;
dev_dn_wingbeats(:,1:3,i) = dev_dn;
pitch_up_wingbeats(:,1:3,i) = pitch_up;
pitch_dn_wingbeats(:,1:3,i) = pitch_dn;

end

