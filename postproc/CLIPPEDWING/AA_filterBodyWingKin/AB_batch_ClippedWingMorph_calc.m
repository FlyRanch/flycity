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
        
        [data_now] = MomentofAreafromFlatWingFrames_4batch_txtNxls(plot_on);
        
        %% save data
        WingClipData.date(seq_nr,1) = str2num(dir_now(1:8));
        WingClipData.seq(seq_nr,1) = str2num(dir_now(11:end));
        
        WingClipData.WingAreaRatio(seq_nr,1) = data_now.WingAreaRatio;
        WingClipData.SecondMomentRatioCentroid(seq_nr,1) = data_now.SecondMomentRatioCentroid;

        WingClipData.WingAreaC(seq_nr,1) = data_now.WingAreaC;
        WingClipData.WingAreaI(seq_nr,1) = data_now.WingAreaI;


        WingClipData.WingLengthC(seq_nr,1) = data_now.WingLengthC;
        WingClipData.WingLengthI(seq_nr,1) = data_now.WingLengthI;
        WingClipData.WingLengthRatio(seq_nr,1) = data_now.WingLengthRatio;

        %Hinge->Centroid method
        WingClipData.SecondMomentCC(seq_nr,1) = data_now.SecondMomentCC;
        WingClipData.SecondMomentCI(seq_nr,1) = data_now.SecondMomentCI;


        %Hinge->Tip method
        WingClipData.SecondMomentTC(seq_nr,1) = data_now.SecondMomentTC;
        WingClipData.SecondMomentTI(seq_nr,1) = data_now.SecondMomentTI;
        WingClipData.SecondMomentRatioTip(seq_nr,1) = data_now.SecondMomentRatioTip;

        cd ..
    end
end

cd ..
save('WingClipDatabase.mat','WingClipData')