function [WingClipData] = MomentofAreafromFlatWingFrames_4batch_txtNxlsNimages_CORRECTED(plot_on)

% clear all; close all; clc;

%% hinge location (x,y) FLIPPED

%% load clipped data
new_name = dir('Results_Clipped*txt');

% rename xls to txt
if isempty(new_name)
    filename = dir('Results_Clipped*xls');
    filename = filename.name;
    [current_dir, current_name, current_ext] = fileparts(filename);
    old_name = strcat(current_name, current_ext);
    new_ext ='.txt';
    new_name=strcat(current_name,new_ext);
    copyfile(old_name, new_name);
else
    new_name = new_name.name;
end

% open txt file
[ndata] = importdata(new_name);
HingePt=ndata.data;
HingeC=[HingePt(2,7), HingePt(2,6)];



%% load intact data
new_name = dir('Results_Intact*txt');

% rename xls to txt
if isempty(new_name)
    filename = dir('Results_Intact*xls');
    filename = filename.name;
    [current_dir, current_name, current_ext] = fileparts(filename);
    old_name = strcat(current_name, current_ext);
    new_name=strcat(current_name,new_ext);
    copyfile(old_name, new_name);
else
    new_name = new_name.name;
end

% open txt file
[ndata] = importdata(new_name);
HingePt=ndata.data;
HingeI=[HingePt(2,7), HingePt(2,6)];

%% load wing images
filename = dir('Mask_Clipped*tif');
filename = filename.name;
im_clipped=imread(filename);

filename = dir('Mask_Int*tif');
filename = filename.name;
im_intact=imread(filename);

%% Converts im_clipped and im_intact into matrix, 0= empty, 1= wing area (pixel)
im_clipped=double(im_clipped);
im_clipped=im_clipped/255;
im_clipped=im_clipped-1;
im_clipped=-im_clipped;
WingAreaC=sum(im_clipped(:));

im_intact=double(im_intact);
im_intact=im_intact/255;
im_intact=im_intact-1;
im_intact=-im_intact;
WingAreaI=sum(im_intact(:));

WingAreaRatio= WingAreaC/WingAreaI; %Temporarily suppressed output

%% Centroid Computation
[AX, AY] = size(im_clipped);
IntAC=[0, 0];
IntAI=[0, 0];
DistC=zeros(AX,AY);
DistI=zeros(AX,AY);
for i=1:AX;
    for j=1:AY;
        if im_clipped(i,j)==1;
            IntAC=IntAC+[i, j];
            DistC(i,j)= hypot(j-HingeC(2),i-HingeC(1));
        end
        if im_intact(i,j)==1;
            IntAI=IntAI+[i, j];
            DistI(i,j)= hypot(j-HingeI(2),i-HingeI(1));
        end
        
    end
   
end
WingLengthC=max(DistC(:));
WingTipC=[0,0];
[WingTipC(1),WingTipC(2)] = find(DistC==max(DistC(:)));

WingLengthI=max(DistI(:));
WingTipI=[0,0];
[WingTipI(1),WingTipI(2)] = find(DistI==max(DistI(:)));

WingLengthRatio= WingLengthC/WingLengthI; %Temporarily suppressed output

CentroidC= (1/WingAreaC)*IntAC;
CentroidI= (1/WingAreaI)*IntAI;

%% primary vectors
HCC=CentroidC-HingeC; %Hinge-Centroid Vector (Clipped)
HCI=CentroidI-HingeI; %Hinge-Centroid Vector (Intact)

HTC=WingTipC-HingeC; %Hinge-Tip Vector (Clipped)
HTI=WingTipI-HingeI; %Hinge-Tip Vector (Intact)

% HCC = HCC/norm(HCC);
% HCI = HCI/norm(HCI);
% 
% HTC = HTC/norm(HTC);
% HTI = HTI/norm(HTI);

%% orthogonal vectors to primary vectors
AXISCC = null(HCC)';
AXISCI = null(HCI)';

AXISTC = null(HTC)';
AXISTI = null(HTI)';

%% Plot Centroids onto Wing Area, compare axis

if plot_on == 1
    figure(1)
    imagesc(DistC);figure(gcf);
    hold on
    axis equal
    axis tight
    axis([HingeC(2)-WingLengthC HingeC(2)+WingLengthC HingeC(1)-WingLengthC HingeC(1)+WingLengthC])
    % colormap gray

    plot(HingeC(2),HingeC(1), '.r')
    plot(CentroidC(2),CentroidC(1),'.r')
    quiver(HingeC(2),HingeC(1),1*HCC(2),1*HCC(1),0,'r') %Plot Hinge point + AxisC
    quiver(HingeC(2),HingeC(1),norm(HCC)*AXISCC(2),norm(HCC)*AXISCC(1),0,'r') %Plot Hinge point + AxisC

    plot(HingeC(2),HingeC(1), 'og')
    plot(WingTipC(2),WingTipC(1),'og')
    quiver(HingeC(2),HingeC(1),1*HTC(2),1*HTC(1),0,'g') %Plot Hinge point + AxisC
    quiver(HingeC(2),HingeC(1),norm(HTC)*AXISTC(2),norm(HTC)*AXISTC(1),0,'g') %Plot Hinge point + AxisC

    figure(2)
    imagesc(DistI);figure(gcf);
    hold on
    axis equal
    axis tight
    axis([HingeI(2)-WingLengthI HingeI(2)+WingLengthI HingeI(1)-WingLengthI HingeI(1)+WingLengthI])
    % colormap gray

    plot(HingeI(2),HingeI(1), '.r')
    plot(CentroidI(2),CentroidI(1),'.r')
    quiver(HingeI(2),HingeI(1),1*HCI(2),1*HCI(1),0,'r') %Plot Hinge point + AxisI
    quiver(HingeI(2),HingeI(1),norm(HCI)*AXISCI(2),norm(HCI)*AXISCI(1),0,'r') %Plot Hinge point + AxisI

    plot(HingeI(2),HingeI(1), 'og')
    plot(WingTipI(2),WingTipI(1),'og')
    quiver(HingeI(2),HingeI(1),1*HTI(2),1*HTI(1),0,'g') %Plot Hinge point + AxisI
    quiver(HingeI(2),HingeI(1),norm(HTI)*AXISTI(2),norm(HTI)*AXISTI(1),0,'g') %Plot Hinge point + AxisI
end

%% Integral - Point-axis distance
rhoCC=zeros(AX,AY);
rhoCI=zeros(AX,AY);

rhoTC=zeros(AX,AY);
rhoTI=zeros(AX,AY);

for i=1:AX;
    for j=1:AY;
        if im_clipped(i,j)==1;
            
            DC=[(HingeC(1)-i),(HingeC(2)-j)];
            
            rhoCC(i,j)= norm(cross([AXISCC, 0],[DC, 0])) / norm(AXISCC);
            rhoTC(i,j)= norm(cross([AXISTC, 0],[DC, 0])) / norm(AXISTC);

        end
        if im_intact(i,j)==1;
            
            DI= [(HingeI(1)-i),(HingeI(2)-j)];
            
            rhoCI(i,j)= norm(cross([AXISCI, 0],[DI, 0])) / norm(AXISCI);
            rhoTI(i,j)= norm(cross([AXISTI, 0],[DI, 0])) / norm(AXISTI);
        end
    end
end

%% First Moment of Area Calculation
rho1CC= rhoCC;
rho1TC= rhoTC;
rho1CI= rhoCI;
rho1TI= rhoTI;

FirstMomentCC= sum(rho1CC(:));
FirstMomentTC= sum(rho1TC(:));

FirstMomentCI= sum(rho1CI(:));
FirstMomentTI= sum(rho1TI(:));

% FirstMomentCC_norm = FirstMomentCC/WingAreaC/WingLengthC^1;
% FirstMomentTC_norm = FirstMomentTC/WingAreaC/WingLengthC^1;
% 
% FirstMomentCI_norm = FirstMomentCI/WingAreaI/WingLengthI^1;
% FirstMomentTI_norm = FirstMomentTI/WingAreaI/WingLengthI^1;

FirstMomentCC_norm = FirstMomentCC/WingLengthC^2/WingLengthC^1;
FirstMomentTC_norm = FirstMomentTC/WingLengthC^2/WingLengthC^1;

FirstMomentCI_norm = FirstMomentCI/WingLengthI^2/WingLengthI^1;
FirstMomentTI_norm = FirstMomentTI/WingLengthI^2/WingLengthI^1;

FirstMomentRatioCentroid = FirstMomentCC/FirstMomentCI;
FirstMomentRatioTip = FirstMomentTC/FirstMomentTI;

%% Second Moment of Area Calculation
rho2CC= rhoCC.^2;
rho2TC= rhoTC.^2;
rho2CI= rhoCI.^2;
rho2TI= rhoTI.^2;

SecondMomentCC= sum(rho2CC(:));
SecondMomentTC= sum(rho2TC(:));

SecondMomentCI= sum(rho2CI(:));
SecondMomentTI= sum(rho2TI(:));

% SecondMomentCC_norm = SecondMomentCC/WingAreaC/WingLengthC^2;
% SecondMomentTC_norm = SecondMomentTC/WingAreaC/WingLengthC^2;
% 
% SecondMomentCI_norm = SecondMomentCI/WingAreaI/WingLengthI^2;
% SecondMomentTI_norm = SecondMomentTI/WingAreaI/WingLengthI^2;

SecondMomentCC_norm = SecondMomentCC/WingLengthC^2/WingLengthC^2;
SecondMomentTC_norm = SecondMomentTC/WingLengthC^2/WingLengthC^2;

SecondMomentCI_norm = SecondMomentCI/WingLengthI^2/WingLengthI^2;
SecondMomentTI_norm = SecondMomentTI/WingLengthI^2/WingLengthI^2;

SecondMomentRatioCentroid = SecondMomentCC/SecondMomentCI;
SecondMomentRatioTip = SecondMomentTC/SecondMomentTI;


%% Third Moment of Area Calculation
rho3CC= rhoCC.^3;
rho3TC= rhoTC.^3;
rho3CI= rhoCI.^3;
rho3TI= rhoTI.^3;

ThirdMomentCC= sum(rho3CC(:));
ThirdMomentTC= sum(rho3TC(:));

ThirdMomentCI= sum(rho3CI(:));
ThirdMomentTI= sum(rho3TI(:));

% ThirdMomentCC_norm = ThirdMomentCC/WingAreaC/WingLengthC^3;
% ThirdMomentTC_norm = ThirdMomentTC/WingAreaC/WingLengthC^3;
% 
% ThirdMomentCI_norm = ThirdMomentCI/WingAreaI/WingLengthI^3;
% ThirdMomentTI_norm = ThirdMomentTI/WingAreaI/WingLengthI^3;

ThirdMomentCC_norm = ThirdMomentCC/WingLengthC^2/WingLengthC^3;
ThirdMomentTC_norm = ThirdMomentTC/WingLengthC^2/WingLengthC^3;

ThirdMomentCI_norm = ThirdMomentCI/WingLengthI^2/WingLengthI^3;
ThirdMomentTI_norm = ThirdMomentTI/WingLengthI^2/WingLengthI^3;

ThirdMomentRatioCentroid = ThirdMomentCC/ThirdMomentCI;
ThirdMomentRatioTip = ThirdMomentTC/ThirdMomentTI;

%% Centroid second moment

if plot_on == 1
figure(3)
imagesc(rhoCC);figure(gcf);
hold on
axis equal
axis tight
axis([HingeC(2)-WingLengthC HingeC(2)+WingLengthC HingeC(1)-WingLengthC HingeC(1)+WingLengthC])
% colormap gray

plot(HingeC(2),HingeC(1), '.r')
plot(CentroidC(2),CentroidC(1),'.r')
quiver(HingeC(2),HingeC(1),1*HCC(2),1*HCC(1),0,'r') %Plot Hinge point + AxisC
quiver(HingeC(2),HingeC(1),norm(HCC)*AXISCC(2),norm(HCC)*AXISCC(1),0,'r') %Plot Hinge point + AxisC

plot(HingeC(2),HingeC(1), 'og')
plot(WingTipC(2),WingTipC(1),'og')
quiver(HingeC(2),HingeC(1),1*HTC(2),1*HTC(1),0,'g') %Plot Hinge point + AxisC
quiver(HingeC(2),HingeC(1),norm(HTC)*AXISTC(2),norm(HTC)*AXISTC(1),0,'g') %Plot Hinge point + AxisC

figure(4)
imagesc(rhoCI);figure(gcf);
hold on
axis equal
axis tight
axis([HingeI(2)-WingLengthI HingeI(2)+WingLengthI HingeI(1)-WingLengthI HingeI(1)+WingLengthI])
% colormap gray

plot(HingeI(2),HingeI(1), '.r')
plot(CentroidI(2),CentroidI(1),'.r')
quiver(HingeI(2),HingeI(1),1*HCI(2),1*HCI(1),0,'r') %Plot Hinge point + AxisI
quiver(HingeI(2),HingeI(1),norm(HCI)*AXISCI(2),norm(HCI)*AXISCI(1),0,'r') %Plot Hinge point + AxisI

plot(HingeI(2),HingeI(1), 'og')
plot(WingTipI(2),WingTipI(1),'og')
quiver(HingeI(2),HingeI(1),1*HTI(2),1*HTI(1),0,'g') %Plot Hinge point + AxisI
quiver(HingeI(2),HingeI(1),norm(HTI)*AXISTI(2),norm(HTI)*AXISTI(1),0,'g') %Plot Hinge point + AxisI

%% Tip second moment
figure(5)
imagesc(rhoTC);figure(gcf);
hold on
axis equal
axis tight
axis([HingeC(2)-WingLengthC HingeC(2)+WingLengthC HingeC(1)-WingLengthC HingeC(1)+WingLengthC])
% colormap gray

plot(HingeC(2),HingeC(1), '.r')
plot(CentroidC(2),CentroidC(1),'.r')
quiver(HingeC(2),HingeC(1),1*HCC(2),1*HCC(1),0,'r') %Plot Hinge point + AxisC
quiver(HingeC(2),HingeC(1),norm(HCC)*AXISCC(2),norm(HCC)*AXISCC(1),0,'r') %Plot Hinge point + AxisC

plot(HingeC(2),HingeC(1), 'og')
plot(WingTipC(2),WingTipC(1),'og')
quiver(HingeC(2),HingeC(1),1*HTC(2),1*HTC(1),0,'g') %Plot Hinge point + AxisC
quiver(HingeC(2),HingeC(1),norm(HTC)*AXISTC(2),norm(HTC)*AXISTC(1),0,'g') %Plot Hinge point + AxisC

figure(6)
imagesc(rhoTI);figure(gcf);
hold on
axis equal
axis tight
axis([HingeI(2)-WingLengthI HingeI(2)+WingLengthI HingeI(1)-WingLengthI HingeI(1)+WingLengthI])
% colormap gray

plot(HingeI(2),HingeI(1), '.r')
plot(CentroidI(2),CentroidI(1),'.r')
quiver(HingeI(2),HingeI(1),1*HCI(2),1*HCI(1),0,'r') %Plot Hinge point + AxisI
quiver(HingeI(2),HingeI(1),norm(HCI)*AXISCI(2),norm(HCI)*AXISCI(1),0,'r') %Plot Hinge point + AxisI

plot(HingeI(2),HingeI(1), 'og')
plot(WingTipI(2),WingTipI(1),'og')
quiver(HingeI(2),HingeI(1),1*HTI(2),1*HTI(1),0,'g') %Plot Hinge point + AxisI
quiver(HingeI(2),HingeI(1),norm(HTI)*AXISTI(2),norm(HTI)*AXISTI(1),0,'g') %Plot Hinge point + AxisI

end

%% correct S2 & S3 for erronious hinge loc (!!!!!!!!!!!!ONLY FOR TRAILING EDGE CUT WINGS!!!!!!!!!!!)

% DO NOT CORRECT WING AREA (STALK NEAR HINGE IS SMALL)
WingAreaC_corr = WingAreaC;
WingAreaI_corr = WingAreaI;
WingAreaRatio_corr = WingAreaRatio;

% WING LENGTH IS MEAN LENGTH
WingLength_mean = mean([WingLengthC WingLengthI]);

WingLengthC_corr = WingLength_mean;
WingLengthI_corr = WingLength_mean;
WingLengthRatio_corr = WingLengthC_corr/WingLengthI_corr;

% CG POSITION
dCC = FirstMomentCC / WingAreaC;
dCI = FirstMomentCI / WingAreaI;

% TRANSLATE MEAN HINGE WRT CG
dCC_corr = dCC + (WingLength_mean - WingLengthC);
dCI_corr = dCI + (WingLength_mean - WingLengthI);

% calc corrected values
FirstMomentCC_corr = FirstMomentCC + WingAreaC* (dCC_corr^1-dCC^1);
FirstMomentCI_corr = FirstMomentCI + WingAreaC* (dCI_corr^1-dCI^1);
FirstMomentRatioCentroid_corr = FirstMomentCC_corr/FirstMomentCI_corr;

SecondMomentCC_corr = SecondMomentCC + WingAreaC* (dCC_corr^2-dCC^2);
SecondMomentCI_corr = SecondMomentCI + WingAreaC* (dCI_corr^2-dCI^2);
SecondMomentRatioCentroid_corr = SecondMomentCC_corr/SecondMomentCI_corr;

ThirdMomentCC_corr = ThirdMomentCC + WingAreaC* (dCC_corr^3-dCC^3);
ThirdMomentCI_corr = ThirdMomentCI + WingAreaC* (dCI_corr^3-dCI^3);
ThirdMomentRatioCentroid_corr = ThirdMomentCC_corr/ThirdMomentCI_corr;

% norms
% FirstMomentCC_norm_corr = FirstMomentCC_corr/WingAreaC_corr/WingLengthC_corr^1;
% FirstMomentCI_norm_corr = FirstMomentCI_corr/WingAreaI_corr/WingLengthI_corr^1;
% 
% SecondMomentCC_norm_corr = SecondMomentCC_corr/WingAreaC_corr/WingLengthC_corr^2;
% SecondMomentCI_norm_corr = SecondMomentCI_corr/WingAreaI_corr/WingLengthI_corr^2;
% 
% ThirdMomentCC_norm_corr = ThirdMomentCC_corr/WingAreaC_corr/WingLengthC_corr^3;
% ThirdMomentCI_norm_corr = ThirdMomentCI_corr/WingAreaI_corr/WingLengthI_corr^3;

FirstMomentCC_norm_corr = FirstMomentCC_corr/WingLengthC_corr^2/WingLengthC_corr^1;
FirstMomentCI_norm_corr = FirstMomentCI_corr/WingLengthI_corr^2/WingLengthI_corr^1;

SecondMomentCC_norm_corr = SecondMomentCC_corr/WingLengthC_corr^2/WingLengthC_corr^2;
SecondMomentCI_norm_corr = SecondMomentCI_corr/WingLengthI_corr^2/WingLengthI_corr^2;

ThirdMomentCC_norm_corr = ThirdMomentCC_corr/WingLengthC_corr^2/WingLengthC_corr^3;
ThirdMomentCI_norm_corr = ThirdMomentCI_corr/WingLengthI_corr^2/WingLengthI_corr^3;

%% store CORRECTED data
% wing length
WingClipData.WingLengthRatio = WingLengthRatio_corr;
WingClipData.WingLengthC = WingLengthC_corr;
WingClipData.WingLengthI = WingLengthI_corr;

% wing area
WingClipData.WingAreaRatio = WingAreaRatio_corr;
WingClipData.WingAreaC = WingAreaC_corr;
WingClipData.WingAreaI = WingAreaI_corr;
WingClipData.WingAreaC_norm = WingAreaC_corr/WingLengthC_corr^2;
WingClipData.WingAreaI_norm = WingAreaI_corr/WingLengthI_corr^2;

%Hinge->Centroid method
WingClipData.FirstMomentRatioCentroid = FirstMomentRatioCentroid_corr;
WingClipData.FirstMomentCC = FirstMomentCC_corr;
WingClipData.FirstMomentCI = FirstMomentCI_corr;
WingClipData.FirstMomentCC_norm = FirstMomentCC_norm_corr;
WingClipData.FirstMomentCI_norm = FirstMomentCI_norm_corr;

WingClipData.SecondMomentRatioCentroid = SecondMomentRatioCentroid_corr;
WingClipData.SecondMomentCC = SecondMomentCC_corr;
WingClipData.SecondMomentCI = SecondMomentCI_corr;
WingClipData.SecondMomentCC_norm = SecondMomentCC_norm_corr;
WingClipData.SecondMomentCI_norm = SecondMomentCI_norm_corr;

WingClipData.ThirdMomentRatioCentroid = ThirdMomentRatioCentroid_corr;
WingClipData.ThirdMomentCC = ThirdMomentCC_corr;
WingClipData.ThirdMomentCI = ThirdMomentCI_corr;
WingClipData.ThirdMomentCC_norm = ThirdMomentCC_norm_corr;
WingClipData.ThirdMomentCI_norm = ThirdMomentCI_norm_corr;
% 
% %Hinge->Tip method
% WingClipData.FirstMomentRatioTip = FirstMomentRatioTip_corr;
% WingClipData.FirstMomentTC = FirstMomentTC_corr;
% WingClipData.FirstMomentTI = FirstMomentTI_corr;
% WingClipData.FirstMomentTC_norm = FirstMomentTC_norm_corr;
% WingClipData.FirstMomentTI_norm = FirstMomentTI_norm_corr;
% 
% WingClipData.SecondMomentRatioTip = SecondMomentRatioTip_corr;
% WingClipData.SecondMomentTC = SecondMomentTC_corr;
% WingClipData.SecondMomentTI = SecondMomentTI_corr;
% WingClipData.SecondMomentTC_norm = SecondMomentTC_norm_corr;
% WingClipData.SecondMomentTI_norm = SecondMomentTI_norm_corr;
% 
% WingClipData.ThirdMomentRatioTip = ThirdMomentRatioTip_corr;
% WingClipData.ThirdMomentTC = ThirdMomentTC_corr;
% WingClipData.ThirdMomentTI = ThirdMomentTI_corr;
% WingClipData.ThirdMomentTC_norm = ThirdMomentTC_norm_corr;
% WingClipData.ThirdMomentTI_norm = ThirdMomentTI_norm_corr;

%% store data
% % wing length
% WingClipData.WingLengthRatio = WingLengthRatio;
% WingClipData.WingLengthC = WingLengthC;
% WingClipData.WingLengthI = WingLengthI;
% 
% % wing area
% WingClipData.WingAreaRatio = WingAreaRatio;
% WingClipData.WingAreaC = WingAreaC;
% WingClipData.WingAreaI = WingAreaI;
% WingClipData.WingAreaC_norm = WingAreaC/WingLengthC^2;
% WingClipData.WingAreaI_norm = WingAreaI/WingLengthI^2;
% 
% %Hinge->Centroid method
% WingClipData.FirstMomentRatioCentroid = FirstMomentRatioCentroid;
% WingClipData.FirstMomentCC = FirstMomentCC;
% WingClipData.FirstMomentCI = FirstMomentCI;
% WingClipData.FirstMomentCC_norm = FirstMomentCC_norm;
% WingClipData.FirstMomentCI_norm = FirstMomentCI_norm;
% 
% WingClipData.SecondMomentRatioCentroid = SecondMomentRatioCentroid;
% WingClipData.SecondMomentCC = SecondMomentCC;
% WingClipData.SecondMomentCI = SecondMomentCI;
% WingClipData.SecondMomentCC_norm = SecondMomentCC_norm;
% WingClipData.SecondMomentCI_norm = SecondMomentCI_norm;
% 
% WingClipData.ThirdMomentRatioCentroid = ThirdMomentRatioCentroid;
% WingClipData.ThirdMomentCC = ThirdMomentCC;
% WingClipData.ThirdMomentCI = ThirdMomentCI;
% WingClipData.ThirdMomentCC_norm = ThirdMomentCC_norm;
% WingClipData.ThirdMomentCI_norm = ThirdMomentCI_norm;
% 
% %Hinge->Tip method
% WingClipData.FirstMomentRatioTip = FirstMomentRatioTip;
% WingClipData.FirstMomentTC = FirstMomentTC;
% WingClipData.FirstMomentTI = FirstMomentTI;
% WingClipData.FirstMomentTC_norm = FirstMomentTC_norm;
% WingClipData.FirstMomentTI_norm = FirstMomentTI_norm;
% 
% WingClipData.SecondMomentRatioTip = SecondMomentRatioTip;
% WingClipData.SecondMomentTC = SecondMomentTC;
% WingClipData.SecondMomentTI = SecondMomentTI;
% WingClipData.SecondMomentTC_norm = SecondMomentTC_norm;
% WingClipData.SecondMomentTI_norm = SecondMomentTI_norm;
% 
% WingClipData.ThirdMomentRatioTip = ThirdMomentRatioTip;
% WingClipData.ThirdMomentTC = ThirdMomentTC;
% WingClipData.ThirdMomentTI = ThirdMomentTI;
% WingClipData.ThirdMomentTC_norm = ThirdMomentTC_norm;
% WingClipData.ThirdMomentTI_norm = ThirdMomentTI_norm;