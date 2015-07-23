clc
clear
close all
warning off

plot_on = 0;
% plot_on = 1;

x_start = 381;
name_start = 6;

%% calc correction factor
% robo data
Lrobo = .23;
Lmodel = .1945;
Lratio_robomodel = Lrobo / Lmodel;

% load intact wing
name_intact = dir('*00.jpg');
name_intact = name_intact.name;

% load & process image
im = imread(name_intact);
im = double(im(:,:,1));
im = im/max(im(:));
im = abs(im-1);

if plot_on == 1
    % test cut
    x_start = 381;
    im(:,1:x_start-1)=0.5;
    imshow(im)
    colormap gray
    pause
end

im(:,1:x_start-1)=0;
im(im==0)=nan;

% calc intact wing winglength in pixels
[y_wing,x_wing] = find(im==1);
l_model_pixels = max(x_wing)-min(x_wing);
l_wing_pixels = Lratio_robomodel * l_model_pixels;

Dl_pixels = l_wing_pixels - l_model_pixels;

%% run calc routine
list = dir('*.jpg');
for ls = 1:length(list)
    
    im_name = list(ls).name;
    
    % cut amount
    cut_perc = str2num(im_name(name_start:name_start+1));
    if cut_perc == 0
        cut_perc = 100;
    end
    
    % cut type
    if im_name(name_start+2)=='.'       % NO cut
        cut_type = 0;
    elseif im_name(name_start+2)=='d'   % distal cut (tip)
        cut_type = 1;
    elseif im_name(name_start+2)=='t'   % trailing edge cut
        cut_type = 2;
    end
    
    % load & process image
    im = imread(im_name);
    im = double(im(:,:,1));
    im = im/max(im(:));
    im = abs(im-1);
    
    if plot_on == 1
        % test cut
        im(:,1:x_start-1)=0.5;
        imshow(im)
        colormap gray
        pause
    end
    
    im(:,1:x_start-1)=0;
    im(im==0)=nan;
    
    %% !!!!!!!!!! weighted image: start at Dl_pixels !!!!!!!!!!!!!!!!
%     n=0;
    n = Dl_pixels;
    im_weight = im;
    for i = x_start:size(im,2)
        n=n+1;
        im_weight(:,i) = im_weight(:,i)*n;
    end
    
    if plot_on == 1
        % test im_weight
        imagesc(im_weight)
        colormap jet
        pause
    end

    
    % calc morph vars
    [y_wing,x_wing] = find(im==1);
    
%     l_pixels = max(x_wing)-min(x_wing);
    l_pixels = max(x_wing)-min(x_wing) + Dl_pixels; % !!! axis correction !!!

    A_pixels = nansum(im(:));
    
    S1_pixels = nansum(im_weight(:));
    S2_pixels = nansum(im_weight(:).^2);
    S3_pixels = nansum(im_weight(:).^3);

    AR = l_pixels^2/A_pixels;
    c_pixels = A_pixels/l_pixels;
    CoA_pixels = S1_pixels/A_pixels;
    
    % correct for difference in rotation point
    
    
    % norms based on wing length
    A_norm = A_pixels/l_pixels^2;
    S1_norm = S1_pixels/l_pixels^3;
    S2_norm = S2_pixels/l_pixels^4;
    S3_norm = S3_pixels/l_pixels^5;
    c_norm = c_pixels/l_pixels;
    CoA_norm = CoA_pixels/l_pixels;
    
    % norms based on wing area
    l_normA = l_pixels/A_pixels^(1/2);
    S1_normA = S1_pixels/A_pixels^(3/2);
    S2_normA = S2_pixels/A_pixels^(4/2);
    S3_normA = S3_pixels/A_pixels^(5/2);
    c_normA = c_pixels/A_pixels^(1/2);
    CoA_normA = CoA_pixels/A_pixels^(1/2);
    
    %% morph at no cut
    if cut_perc==100
        l_pixels_nocut = l_pixels;
        A_pixels_nocut = A_pixels;

        S1_pixels_nocut = S1_pixels;
        S2_pixels_nocut = S2_pixels;
        S3_pixels_nocut = S3_pixels;
        
        AR_nocut = AR;
        c_pixels_nocut = c_pixels;
        CoA_pixels_nocut = CoA_pixels;
    end

    %% save data
    BorfMorphCutData.cut_perc(ls,1) = cut_perc;
    BorfMorphCutData.cut_type(ls,1) = cut_type;

    % wing length
    BorfMorphCutData.WingLength_pixels(ls,1) = l_pixels;
    BorfMorphCutData.WingLength_normA(ls,1) = l_normA;

    % wing area
    BorfMorphCutData.WingArea_pixels(ls,1) = A_pixels;
    BorfMorphCutData.WingArea_norm(ls,1) = A_norm;

    % S1&S2&S3
    BorfMorphCutData.FirstMoment_pixels(ls,1) = S1_pixels;
    BorfMorphCutData.FirstMoment_norm(ls,1) = S1_norm;
    BorfMorphCutData.FirstMoment_normA(ls,1) = S1_normA;

    BorfMorphCutData.SecondMoment_pixels(ls,1) = S2_pixels;
    BorfMorphCutData.SecondMoment_norm(ls,1) = S2_norm;
    BorfMorphCutData.SecondMoment_normA(ls,1) = S2_normA;

    BorfMorphCutData.ThirdMoment_pixels(ls,1) = S3_pixels;
    BorfMorphCutData.ThirdMoment_norm(ls,1) = S3_norm;
    BorfMorphCutData.ThirdMoment_normA(ls,1) = S3_normA;
    
    % aditional
    BorfMorphCutData.AR(ls,1) = AR;
    
    BorfMorphCutData.c_pixels(ls,1) = c_pixels;
    BorfMorphCutData.c_norm(ls,1) = c_norm;
    BorfMorphCutData.c_normA(ls,1) = c_normA;
    
    BorfMorphCutData.CoA_pixels(ls,1) = CoA_pixels;
    BorfMorphCutData.CoA_norm(ls,1) = CoA_norm;
    BorfMorphCutData.CoA_normA(ls,1) = CoA_normA;
end

%% calc morph ratio (cut/NOcut)
BorfMorphCutData.WingLength_ratio   = BorfMorphCutData.WingLength_pixels    /l_pixels_nocut;
BorfMorphCutData.WingArea_ratio     = BorfMorphCutData.WingArea_pixels      /A_pixels_nocut;

BorfMorphCutData.FirstMoment_ratio  = BorfMorphCutData.FirstMoment_pixels   /S1_pixels_nocut;
BorfMorphCutData.SecondMoment_ratio = BorfMorphCutData.SecondMoment_pixels  /S2_pixels_nocut;
BorfMorphCutData.ThirdMoment_ratio  = BorfMorphCutData.ThirdMoment_pixels   /S3_pixels_nocut;

BorfMorphCutData.AR_ratio  = BorfMorphCutData.AR   /AR_nocut;
BorfMorphCutData.c_ratio  = BorfMorphCutData.c_pixels   /c_pixels_nocut;
BorfMorphCutData.CoA_ratio  = BorfMorphCutData.CoA_pixels   /CoA_pixels_nocut;

BorfMorphCutData.Lrobo = Lrobo;
BorfMorphCutData.Lmodel = Lmodel;

%% save
% cd ..
save('BorfMorphCutDatabase.mat','BorfMorphCutData')

%% test
if plot_on == 1
    load('borf_clipped_wing_geometry.mat')
    figure
    axis equal
    hold on
    axis equal

    plot(A_ratio,BorfMorphCutData.WingArea_ratio,'*r')
    plot(CoA_ratio,BorfMorphCutData.CoA_ratio,'*g')
    plot(S1_ratio,BorfMorphCutData.FirstMoment_ratio,'*')
    plot(S2_ratio,BorfMorphCutData.SecondMoment_ratio,'*c')

    legend('Area','CoA','S1','S2')
    plot([0 1],[0 1],'k')
end