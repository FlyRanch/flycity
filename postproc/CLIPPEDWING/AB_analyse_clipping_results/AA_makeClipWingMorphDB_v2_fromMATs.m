clear all; close all; clc;

list = dir('*.mat')

for i=1:length(list)
    
    file_now = list(i).name;
    load(file_now);
    
    WingClipData.date(i,1) = str2num(file_now(1:8));
    WingClipData.seq(i,1) = str2num(file_now(11:14));
    
    WingClipData.SecondMomentCC(i,1) = SecondMomentCC;
    WingClipData.SecondMomentCI(i,1) = SecondMomentCI;
    WingClipData.SecondMomentTC(i,1) = SecondMomentTC;
    WingClipData.SecondMomentTI(i,1) = SecondMomentTI;
    WingClipData.SecondMomentRatioCentroid(i,1) = SecondMomentRatioCentroid;
    WingClipData.SecondMomentRatioTip(i,1) = SecondMomentRatioTip;
    
    WingClipData.WingAreaC(i,1) = WingAreaC;
    WingClipData.WingAreaI(i,1) = WingAreaI;
    WingClipData.WingAreaRatio(i,1) = WingAreaRatio;
    
    WingClipData.WingLengthC(i,1) = WingLengthC;
    WingClipData.WingLengthI(i,1) = WingLengthI;
    WingClipData.WingLengthRatio(i,1) = WingLengthRatio;
    
end

% save data
cd ..
save('WingClipData.mat','WingClipData')