clc
clear
close all
warning off

plot_on = 0;
seq_nr = 0;
cd('FlatWingFrames')

list = dir;
for ls = 3:length(list)
    
    
    if list(ls).isdir == 1
        
        length(list) - ls
        seq_nr = seq_nr+1;
    
        dir_now = list(ls).name
        cd(dir_now)
        
        [data_now] = MomentofAreafromFlatWingFrames_4batch_txtNxlsNimages_S2S3(plot_on);
        
        %% save data
        WingClipData.date(seq_nr,1) = str2num(dir_now(1:8));
        WingClipData.seq(seq_nr,1) = str2num(dir_now(11:end));
        
        % wing area
        WingClipData.WingAreaRatio(seq_nr,1) = data_now.WingAreaRatio;
        WingClipData.WingAreaC_pixels(seq_nr,1) = data_now.WingAreaC;
        WingClipData.WingAreaI_pixels(seq_nr,1) = data_now.WingAreaI;
        
        % wing length
        WingClipData.WingLengthRatio(seq_nr,1) = data_now.WingLengthRatio;
        WingClipData.WingLengthC_pixels(seq_nr,1) = data_now.WingLengthC;
        WingClipData.WingLengthI_pixels(seq_nr,1) = data_now.WingLengthI;

        % S2&S3: Hinge->Centroid method
        WingClipData.SecondMomentRatioCentroid(seq_nr,1) = data_now.SecondMomentRatioCentroid;
        WingClipData.SecondMomentCC_pixels(seq_nr,1) = data_now.SecondMomentCC;
        WingClipData.SecondMomentCI_pixels(seq_nr,1) = data_now.SecondMomentCI;
        WingClipData.SecondMomentCC_norm(seq_nr,1) = data_now.SecondMomentCC_norm;
        WingClipData.SecondMomentCI_norm(seq_nr,1) = data_now.SecondMomentCI_norm;

        WingClipData.ThirdMomentRatioCentroid(seq_nr,1) = data_now.ThirdMomentRatioCentroid;
        WingClipData.ThirdMomentCC_pixels(seq_nr,1) = data_now.ThirdMomentCC;
        WingClipData.ThirdMomentCI_pixels(seq_nr,1) = data_now.ThirdMomentCI;
        WingClipData.ThirdMomentCC_norm(seq_nr,1) = data_now.ThirdMomentCC_norm;
        WingClipData.ThirdMomentCI_norm(seq_nr,1) = data_now.ThirdMomentCI_norm;

        % S2&S3: Hinge->Tip method
        WingClipData.SecondMomentRatioTip(seq_nr,1) = data_now.SecondMomentRatioTip;
        WingClipData.SecondMomentTC_pixels(seq_nr,1) = data_now.SecondMomentTC;
        WingClipData.SecondMomentTI_pixels(seq_nr,1) = data_now.SecondMomentTI;
        WingClipData.SecondMomentTC_norm(seq_nr,1) = data_now.SecondMomentTC_norm;
        WingClipData.SecondMomentTI_norm(seq_nr,1) = data_now.SecondMomentTI_norm;

        WingClipData.ThirdMomentRatioTip(seq_nr,1) = data_now.ThirdMomentRatioTip;
        WingClipData.ThirdMomentTC_pixels(seq_nr,1) = data_now.ThirdMomentTC;
        WingClipData.ThirdMomentTI_pixels(seq_nr,1) = data_now.ThirdMomentTI;
        WingClipData.ThirdMomentTC_norm(seq_nr,1) = data_now.ThirdMomentTC_norm;
        WingClipData.ThirdMomentTI_norm(seq_nr,1) = data_now.ThirdMomentTI_norm;

        cd ..
    end
end

cd ..
save('WingClipDatabase.mat','WingClipData')