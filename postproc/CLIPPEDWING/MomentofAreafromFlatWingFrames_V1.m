clear all; close all; clc;

%ARBITRARILY ASSIGNING HINGE POINTS UNTIL LOAD DATA IS FIXED:
%  filename = dir('Results_Clipped*txt');
%  filename = filename.name;
% [ndata, headertext] = importdata(filename)
% xlsfinfo(filename)
% hc= load(filename);

HingeI=[557,508];
HingeC=[589,509];

filename = dir('Mask_Clipped*tif');
filename = filename.name;
a=imread(filename);

filename = dir('Mask_Int*tif');
filename = filename.name;
b=imread(filename);

%% Converts a and b into matrix, 0= empty, 1= wing area (pixel)
a=double(a);
a=a/255;
a=a-1;
a=-a;
A_Clipped=sum(a(:));

b=double(b);
b=b/255;
b=b-1;
b=-b;
A_Intact=sum(b(:));
WingAreaRatio= A_Clipped/A_Intact

%% Centroid Computation
[AX, AY] = size(a);
IntAC=[0, 0];
IntAI=[0, 0];
DistC=zeros(AX,AY);
DistI=zeros(AX,AY);
for i=1:AX;
    for j=1:AY;
        if a(i,j)==1;
            IntAC=IntAC+[i, j];
            DistC(i,j)= hypot(j-HingeC(1),i-HingeC(2));
        end
        if b(i,j)==1;
            IntAI=IntAI+[i, j];
            DistI(i,j)= hypot(j-HingeI(1),i-HingeI(2));
        end
        
    end
   
end
WingLengthC=max(DistC(:));
WingTipC=[0,0];
[WingTipC(1),WingTipC(2)] = find(DistC==max(DistC(:)));

WingLengthI=max(DistI(:));
WingTipI=[0,0];
[WingTipI(1),WingTipI(2)] = find(DistI==max(DistI(:)));

CentroidC= (1/A_Clipped)*IntAC;
HCC=CentroidC-HingeC; %Hinge-Centroid Vector (Clipped)
CentroidI= (1/A_Intact)*IntAI;
HCI=CentroidI-HingeI; %Hinge-Centroid Vector (Intact)
AXISC=null(CentroidC-HingeC)
AXISI=null(CentroidI-HingeI)


%% Axis perpendicular to Hinge->Centroid vector
HC=zeros(1,3);
HI=zeros(1,3);
HC=[HCC(1); HCC(2); 0];
HI=[HCI(1); HCI(2); 0];
Z=[0 0 1];
AxisC1= cross(HC,Z);
AxisC= AxisC1./norm(AxisC1); %unit vector
AxisC1=[AxisC1(1), AxisC1(2)]
AxisI1= cross(HI,Z);
AxisI= AxisI1./norm(AxisI1); %unit vector

%Plot Centroids onto Wing Area, compare axis
figure(1)
imagesc(a);figure(gcf);
hold on
axis equal
plot(CentroidC(2),CentroidC(1),'*w')
plot(WingTipC(2), WingTipC(1),'.w')
plot(HingeC(1),HingeC(2), '.w')
quiver(HingeC(1),HingeC(2),AxisC1(1),AxisC1(2),'-w') %Plot Hinge point + AxisC


figure(2)
imagesc(b);figure(gcf);
hold on
axis equal
plot(CentroidI(2),CentroidI(1),'*w')
plot(WingTipI(2), WingTipI(1),'.w')
plot(HingeI(1),HingeI(2), '.w')
quiver(HingeI(1),HingeI(2),AxisI1(1),AxisI1(2),'-w') %Plot Hinge point + AxisI


%% Integral - Point-axis distance
rhoC=zeros(AX,AY);
rhoI=zeros(AX,AY);
distance=zeros(AX,AY);
for i=1:AX;
    for j=1:AY;
        if a(i,j)==1;
%             pt = [i, j,0];
%             HingeC=[HingeC(1), HingeC(2), 0] 
%             HAC= HingeC+ AxisC;
%             distanceC(i,j)=point_to_line(pt,HingeC,HAC);
            DC=[abs(j-HingeC(1)); abs(i-HingeC(2)); 0];
            rhoC(i,j)= norm(cross(AxisC',DC)) / norm(AxisC');
        end
        if b(i,j)==1;
%             pt = [i, j,0];
%             HingeI=[HingeI(1), HingeI(2), 0] 
%             HAI= HingeI+ AxisI;
%             distanceI(i,j)=point_to_line(pt,HingeI,HAI);
            DI= [abs(j-HingeI(1)); abs(i-HingeI(2)); 0];
            rhoI(i,j)= norm(cross(AxisI',DI)) / norm(AxisI');
        end
    end
end

%% Second Moment of Area Calculation
MAC= rhoC.^2;
MAI= rhoI.^2;
Second_MomentC= sum(MAC(:));
Second_MomentI= sum(MAI(:));
Second_MomentC_norm = Second_MomentC/A_Clipped/WingLengthC
Second_MomentI_norm = Second_MomentI/A_Intact/WingLengthI
SMA_ratio= Second_MomentC/Second_MomentI
