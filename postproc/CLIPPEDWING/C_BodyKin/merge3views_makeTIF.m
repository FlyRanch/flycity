%close all;
clear
clc;
close all
warning off

% skip = 3
skip = 1
fps = 7500/skip;

% start = 1
% stop = 5588
start = 1957
stop = 3100

% % crop settings: pooing fruitfly (640x360)
% crop = true
% x_min = 250
% x_max = x_min + 640 -1
% y_min = 500
% y_max = y_min + 360 -1
crop = false

PathName = [cd,'/'];
cams = dir('*006');

for i=1:length(cams)
    cd(cams(i).name);
    ims = dir('*.tif');
    ims_all(i,1:length(ims)) = ims;
    cd(PathName)
end
    
cd(PathName)
% cd ..
% 
% movfilename = ['3projections_',num2str(fps),'fps.avi'];
% mov = avifile(movfilename, 'compression', 'none');

solname = [cd,'/merged'];
mkdir(solname);

% im_nr = 0;
for i=start:skip:stop
    counter = stop-i
    for cam=1:length(cams)
        cd(PathName)
        cd(cams(cam).name);
        
        im = imread(ims_all(cam,i).name);
        
        if crop == true
            im = im(y_min:y_max,x_min:x_max);
        end

        im = double(im);
        im = im/max(im(:));
        im = uint8(255*im);
%         im = imadjust(im);
%         imshow(im);
%         colormap gray
        
        im_merge(:,(cam-1)*size(im,2)+1:cam*size(im,2),:) = im;
    end
%         imshow(im_merge);
%         colormap greay
    
    
%         cd(solname)
%         tic
%         mov = addframe(mov,gcf);
%         toc
% 
%         im_nr = im_nr +1;
%         imwrite(im_merge,[solname 'merged/merged' num2str(im_nr) '.tif']);
%         imwrite(im_merge8bit,[solname 'merged/merged' num2str(im_nr) '.tif']);
        imwrite(im_merge,[solname '/' ims(i).name]);
%         cd(PathName)
end

% cd(PathName)
% cd ..
% tic
% mov = close(mov);
% toc
% cd(PathName)

% end