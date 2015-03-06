function [WingClipData] = MomentofAreafromFlatWingFrames_4batch_txtNxlsNimages_S2S3(plot_on)

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

FirstMomentCC_norm = FirstMomentCC/WingAreaC/WingLengthC^1;
FirstMomentTC_norm = FirstMomentTC/WingAreaC/WingLengthC^1;

FirstMomentCI_norm = FirstMomentCI/WingAreaI/WingLengthI^1;
FirstMomentTI_norm = FirstMomentTI/WingAreaI/WingLengthI^1;

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

SecondMomentCC_norm = SecondMomentCC/WingAreaC/WingLengthC^2;
SecondMomentTC_norm = SecondMomentTC/WingAreaC/WingLengthC^2;

SecondMomentCI_norm = SecondMomentCI/WingAreaI/WingLengthI^2;
SecondMomentTI_norm = SecondMomentTI/WingAreaI/WingLengthI^2;

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

ThirdMomentCC_norm = ThirdMomentCC/WingAreaC/WingLengthC^3;
ThirdMomentTC_norm = ThirdMomentTC/WingAreaC/WingLengthC^3;

ThirdMomentCI_norm = ThirdMomentCI/WingAreaI/WingLengthI^3;
ThirdMomentTI_norm = ThirdMomentTI/WingAreaI/WingLengthI^3;

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

%% store data
% wing length
WingClipData.WingLengthRatio = WingLengthRatio;
WingClipData.WingLengthC = WingLengthC;
WingClipData.WingLengthI = WingLengthI;

% wing area
WingClipData.WingAreaRatio = WingAreaRatio;
WingClipData.WingAreaC = WingAreaC;
WingClipData.WingAreaI = WingAreaI;
WingClipData.WingAreaC_norm = WingAreaC/WingLengthC^2;
WingClipData.WingAreaI_norm = WingAreaI/WingLengthI^2;

%Hinge->Centroid method
WingClipData.FirstMomentRatioCentroid = FirstMomentRatioCentroid;
WingClipData.FirstMomentCC = FirstMomentCC;
WingClipData.FirstMomentCI = FirstMomentCI;
WingClipData.FirstMomentCC_norm = FirstMomentCC_norm;
WingClipData.FirstMomentCI_norm = FirstMomentCI_norm;

WingClipData.SecondMomentRatioCentroid = SecondMomentRatioCentroid;
WingClipData.SecondMomentCC = SecondMomentCC;
WingClipData.SecondMomentCI = SecondMomentCI;
WingClipData.SecondMomentCC_norm = SecondMomentCC_norm;
WingClipData.SecondMomentCI_norm = SecondMomentCI_norm;

WingClipData.ThirdMomentRatioCentroid = ThirdMomentRatioCentroid;
WingClipData.ThirdMomentCC = ThirdMomentCC;
WingClipData.ThirdMomentCI = ThirdMomentCI;
WingClipData.ThirdMomentCC_norm = ThirdMomentCC_norm;
WingClipData.ThirdMomentCI_norm = ThirdMomentCI_norm;

%Hinge->Tip method
WingClipData.FirstMomentRatioTip = FirstMomentRatioTip;
WingClipData.FirstMomentTC = FirstMomentTC;
WingClipData.FirstMomentTI = FirstMomentTI;
WingClipData.FirstMomentTC_norm = FirstMomentTC_norm;
WingClipData.FirstMomentTI_norm = FirstMomentTI_norm;

WingClipData.SecondMomentRatioTip = SecondMomentRatioTip;
WingClipData.SecondMomentTC = SecondMomentTC;
WingClipData.SecondMomentTI = SecondMomentTI;
WingClipData.SecondMomentTC_norm = SecondMomentTC_norm;
WingClipData.SecondMomentTI_norm = SecondMomentTI_norm;

WingClipData.ThirdMomentRatioTip = ThirdMomentRatioTip;
WingClipData.ThirdMomentTC = ThirdMomentTC;
WingClipData.ThirdMomentTI = ThirdMomentTI;
WingClipData.ThirdMomentTC_norm = ThirdMomentTC_norm;
WingClipData.ThirdMomentTI_norm = ThirdMomentTI_norm;
