function [Features,IMBW] = feat_detect(n,Features,PAR)
% n is the frame number
%  Calculate segmented images
IMBW = FlySegment(0,n,n,PAR.BG,PAR);

%% FTMmod 20120623 only cuurent frame
% %% IMBW changed from 'frame' to 'frame-PAR.startframe+1' FTMmod 20120603
% %Only take the images from the current frame
% IMBW = squeeze(IMBW(:,:,n-PAR.startframe+1,:));
% 

%% Iterate over each camera
for mm = 1:PAR.numcam
    data = IMBW(:,:,mm);    
    
    % Background image generation
    Background = PAR.BG(:,:,mm); 

    %load the raw grayscale image
%     %Number of digit places for the total number of frames
%     digits = length(num2str(PAR.numframes));
%         
%     input_filename = sprintf(['%s/cam%03d/%s%0' num2str(digits) 'd%s'], ...
%         PAR.imagepath, mm, PAR.stub, n,PAR.image_filter(2:end));
    
    %% photron file struc FTMmod
    input_filename = ([PAR.imagepath, PAR.stub(1:3), num2str(mm), PAR.stub(5:end), '/', PAR.stub(1:3), num2str(mm), PAR.stub(5:end), sprintf(['%0',num2str(PAR.digits),'d'], n),PAR.image_filter(2:end)]);
    IMgray = imread(input_filename);
    %%
    
%     %% remove bg FTMmod
%         IMgray = medfilt2(IMgray, [3 3]);
%         Background = medfilt2(Background, [3 3]);
% 
%         IMgray = imadjust(IMgray);
%         Background = imadjust(Background);   
% 
%         IMgray = imsubtract(Background, IMgray);
%     %%

    flypix = IMgray(data == 1);
    %THLmod hack to fix problem with dead pixels - maybe come up with
    %something better than this.
    
    flypix = flypix(flypix < 150);
%     close all
%     figure, imagesc(IMgray), colormap gray
%     saveas(gcf,['flyim_raw',num2str(mm),'.jpg'])
%     figure, imagesc(IMgray.*data), colormap gray
%     saveas(gcf,['flyim_fly',num2str(mm),'.jpg'])
    
    %define initial conditions for EM algorithm
    init.W = [.5 .5];
        %% FTMmod image intensity dependent initial condition
    init.M = [mode(double(flypix)) mean(double(flypix))];
%     init.M = [mode(double(flypix)) 55];
    %%
    init.V(:,:,1) = 3;
    init.V(:,:,2) = 100;
    [W,M,V,L,Q] = EM_GM_fast(double(flypix),2,[],[],1,init);
    %Calculate the local minimum between the two gaussian peaks.
    xidx = find(min(M) < Q(:,1) & Q(:,1) < max(M));
    [vel,pidx] = min(Q(xidx,2));
    tmpx = Q(xidx,1);
    thresh = tmpx(pidx);
    
    %% FTMmod
%     thresh = thresh * PAR.WingBodyRatio;
    
    %% added for ill defined thresholds FTMmod 20120608
    if isempty(thresh)==1
        thresh = nanmean(M);
        disp('here')
    else
%     figure, plot(Q(:,1),Q(:,2))
%     saveas(gcf,['flyim_gausian',num2str(mm),'.jpg'])
    end
    
%     %% changed thresh to not include wings in IMbody FTMmod 20120608
%         thresh = mean(M);

    IMbody = zeros(size(data));
    IMbody(min(flypix) <= IMgray & IMgray <= thresh) = 1;

    %% FTMmod remove body stuff outside data
    IMbody=IMbody.*double(data);
    %%
%     close all
%      figure,imshow(IMbody)
%      pause
% %% BGsub body is light instead of dark FTMmod 20120618
% %     IMbody = zeros(size(data));    
% %     IMbody(IMgray >= thresh) = 1;
% %%
    %Convert to 256 grayscale image
    %data(logical(IMbody)) = 0;
    
    %subplot(2,2,1), imshow(255*data)
    
    % remove small areas (smaller than 50 pixels)
    data = bwareaopen(data, 50);
    
    %Fill in stray pixels
    data = bwmorph(data,'bridge',inf);
    data = bwmorph(data,'fill',inf);
    data = bwmorph(data,'clean',inf);
    data = imfill(data,'holes');
    
    IMbody = bwmorph(IMbody,'bridge',inf);
    IMbody = bwmorph(IMbody,'fill',inf);
    IMbody = bwmorph(IMbody,'clean',inf);
    IMbody = imfill(IMbody,'holes');
    
    % smooth edges
    SE = strel('disk',2); 
%     SE = strel('square',2)
    IMbody = imerode(IMbody,SE);
    IMbody = imdilate(IMbody,SE);
    data = imerode(data,SE);
    data = imdilate(data,SE);
    
%    figure(1),hold off
%    imshow(IMbody)
%    figure(2),hold off
%    imshow(data)
%    pause
    
%    % remove areas larger than IMfull (data)
%    CC = bwconncomp(data);
%    numPixeldata = cellfun(@numel,CC.PixelIdxList);
%    numPixelMax = max(numPixeldata);
%    CC = bwconncomp(IMbody);
%    numPixels = cellfun(@numel,CC.PixelIdxList);
%    [toobig,idx] = find(numPixels>numPixelMax);
%    for i=1:size(idx,2)
%        IMbody(CC.PixelIdxList{idx(i)}) = 0;
%    end
%    
%     remove small areas (smaller than 1/10 of fly or 200pixels)
%    if size(numPixeldata,2)==1
%     %   IMbody = bwareaopen(IMbody, round(numPixelMax/10));
%     else
%        IMbody = bwareaopen(IMbody, 200);
%     end
%     
% figure(1), imshow((IMbody + data)/2), colormap gray
% pause(.1)
% saveas(gcf,['flyim_bodywing',num2str(mm),'.jpg'])
% close all
    %% save
    Features(n).IMbodyfull{mm} = IMbody.*255;
    
    
    % Add artificial occlusion
    if ~isempty(PAR.OccludeShape{mm})
        %keyboard
        fisp('feat detect')
        tmp = data;
        BW = roipoly(tmp,PAR.OccludeShape{mm}(:,1),PAR.OccludeShape{mm}(:,2));
        tmp(BW) = 0;
        data = tmp;
    end

    Features(n).IMfull{mm} = data.*255;
    
%     %% save images FTMmod 20120607
%     imwrite(Features(n).IMfull{mm},    [PAR.solutionpath 'Images_' PAR.solutiondirname '/flyimage_' sprintf(['%0',num2str(PAR.digits),'d'], n) '_cam' num2str(mm) '_full.bmp']);
%     imwrite(Features(n).IMbodyfull{mm},[PAR.solutionpath 'Images_' PAR.solutiondirname '/flyimage_' sprintf(['%0',num2str(PAR.digits),'d'], n) '_cam' num2str(mm) '_body.bmp']);
    

end % iterate over 'm', repeat for each camera