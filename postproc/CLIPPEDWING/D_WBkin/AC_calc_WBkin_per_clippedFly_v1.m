clc
clear
close all

loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

%% steady wb
freq_steady = f_wb_steady_meanCIstd;
stroke_steady = stroke_wb_steady_bins_meanCIstd;
Astroke_steady = max(stroke_steady(:,1)) - min(stroke_steady(:,1));

            
%% loop through SecondMomentRatio's
SecondMomentRatio_list = unique(SecondMomentRatio);


%% clipped fly WBs
% make nan DB arrays
S2_ratio_all = nan(150,length(SecondMomentRatio_list));
S3_ratio_all = nan(150,length(SecondMomentRatio_list));

freq_all = nan(150,length(SecondMomentRatio_list));

Astroke_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_intact_all = nan(150,length(SecondMomentRatio_list));

Astroke_ratio_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_intact_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_intact_clip_all = nan(150,length(SecondMomentRatio_list));


for seq_now = 1:length(SecondMomentRatio_list)
    
    S2_now = SecondMomentRatio_list(seq_now)
    wb = find(SecondMomentRatio==S2_now & steady_nr_mean_wb==1);
    
                % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);
    
%             if mean(seq_nr_now) ~= min(seq_nr_now);
%                 ERROR = 1
%                 break
%             end
            
            %% clip data
            clip_side_now = clip_side(wb);
            clip_type_now = clip_type(wb);

            FirstMomentNormClipped_now = FirstMomentNormClipped(wb);
            FirstMomentNormIntact_now = FirstMomentNormIntact(wb);
            FirstMomentRatio_now = FirstMomentRatio(wb);

            SecondMomentNormClipped_now = SecondMomentNormClipped(wb);
            SecondMomentNormIntact_now = SecondMomentNormIntact(wb);
            SecondMomentRatio_now = SecondMomentRatio(wb);

            ThirdMomentNormClipped_now = ThirdMomentNormClipped(wb);
            ThirdMomentNormIntact_now = ThirdMomentNormIntact(wb);
            ThirdMomentRatio_now = ThirdMomentRatio(wb);

            AreaNormClipped_now = AreaNormClipped(wb);
            AreaNormIntact_now = AreaNormIntact(wb);
            AreaRatio_now = AreaRatio(wb);

            LengthClipped_pixels_now = LengthClipped_pixels(wb);
            LengthIntact_pixels_now = LengthIntact_pixels(wb);
            LengthRatio_now = LengthRatio(wb);

            % body kin data
            vel_now = V_mean_wb(wb);
            
            slip_now = slip_mean_wb(wb);
            roll_now = roll_mean_wb(wb);
            pitch_now = pitch_mean_wb(wb);
            
            Fsp_pitch_now = Fsp_pitch_mean_wb(wb);
            Fsp_roll_now = Fsp_roll_mean_wb(wb);
            
            % wingbeat kin data
            f_wb_now = mean([f_wb_L(wb),f_wb_R(wb)],2);
            
            stroke_clip_now = stroke_wb_R_bins(:,wb);
            stroke_intact_now = stroke_wb_L_bins(:,wb);

            pitch_wb_clip_now = pitch_wb_R_bins(:,wb);
            pitch_wb_intact_now = pitch_wb_L_bins(:,wb);

            dev_wb_clip_now = dev_wb_R_bins(:,wb);
            dev_wb_intact_now = dev_wb_L_bins(:,wb);

            U_wb_clip_now = U_wb_R_bins(:,wb);
            U_wb_intact_now = U_wb_L_bins(:,wb);

            aoa_wb_clip_now = aoa_wb_R_bins(:,wb);
            aoa_wb_intact_now = aoa_wb_L_bins(:,wb);

            Dstroke = Dstroke_wb_bins(:,wb);
            Dpitch = Dpitch_wb_bins(:,wb);
            Ddev = Ddev_wb_bins(:,wb);
            DU = DU_wb_bins(:,wb);
            Daoa = Daoa_wb_bins(:,wb);
            
            %% calc amplitudes
            Astroke_clip_now = (max(stroke_clip_now)-min(stroke_clip_now))';
            Astroke_intact_now = (max(stroke_intact_now)-min(stroke_intact_now))';
            
            Astroke_ratio_clip_now = Astroke_clip_now ./ Astroke_steady;
            Astroke_ratio_intact_now = Astroke_intact_now ./ Astroke_steady;
            Astroke_ratio_intact_clip_now = Astroke_intact_now ./ Astroke_clip_now;

            %% store data in DB
            S2_ratio_all(1:length(SecondMomentRatio_now),seq_now) = SecondMomentRatio_now;
            S3_ratio_all(1:length(ThirdMomentRatio_now),seq_now) = ThirdMomentRatio_now;
            
            freq_all(1:length(f_wb_now),seq_now) = f_wb_now;
   
            Astroke_clip_all(1:length(Astroke_clip_now),seq_now) = Astroke_clip_now;
            Astroke_intact_all(1:length(Astroke_intact_now),seq_now) = Astroke_intact_now;
   
            Astroke_ratio_clip_all(1:length(Astroke_ratio_clip_now),seq_now) = Astroke_ratio_clip_now;
            Astroke_ratio_intact_all(1:length(Astroke_ratio_intact_now),seq_now) = Astroke_ratio_intact_now;
            Astroke_ratio_intact_clip_all(1:length(Astroke_ratio_intact_clip_now),seq_now) = Astroke_ratio_intact_clip_now;
            
            plot(stroke_clip_now,'r')
            hold on
            plot(stroke_intact_now,'b')
            hold off
            seq_nr_now
            seq_now
            mean(Astroke_ratio_intact_clip_now)
            S2_now
%             pause
end

S2_ratio_mean = nanmean(S2_ratio_all)';
S3_ratio_mean = nanmean(S3_ratio_all)';

freq_mean = nanmean(freq_all)';

Astroke_clip_mean = nanmean(Astroke_clip_all)';
Astroke_intact_mean = nanmean(Astroke_intact_all)';

Astroke_ratio_clip_mean = nanmean(Astroke_ratio_clip_all)';
Astroke_ratio_intact_mean = nanmean(Astroke_ratio_intact_all)';
Astroke_ratio_intact_clip_mean = nanmean(Astroke_ratio_intact_clip_all)';


%% plot
% % datapoints with color
% for i = 1:length(Astroke_ratio_intact_clip_mean)
%     color_nr = round(98*(Astroke_ratio_intact_clip_mean(i)-.5)+1);
%     plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor'cmap(color_nr,:),'markersize',5)
% end


subplot(2,1,1)
plot(S2_ratio_mean,Astroke_ratio_intact_clip_mean,'ob')

subplot(2,1,2)
plot(S3_ratio_mean,Astroke_ratio_intact_clip_mean,'or')

save('WBdataset_ClipNintact_wingbeat_kin.mat','S2_ratio_mean','S3_ratio_mean','freq_mean',...
    'Astroke_clip_mean','Astroke_intact_mean',...
    'Astroke_ratio_clip_mean','Astroke_ratio_intact_mean','Astroke_ratio_intact_clip_mean');





