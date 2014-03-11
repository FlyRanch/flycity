%close all;
clear
clc;
close all
warning off

skip = 1

% skip = 3
% start = 2794
% stop = 4063
% 
% % scale = 0.5 
% x_min = 69
% y_max = 852

% scale = 0.55
x_min = 1
y_max = 937

PathName = [cd,'/'];
cams = dir('cam*');

for i=1:length(cams)
    cd(cams(i).name);
    ims = dir('*.tif');
    ims_all(i,1:length(ims)) = ims;
    cd(PathName)
end
    
% cd(cams(1).name);
% ims = dir('*.tif');
% im_templ = imread(ims(1).name);
% im_x = size(im_templ,1);
% im_y = size(im_templ,2);
% im_merge = [im_templ im_templ im_templ];
% im_merge = uint8(im_merge);
% cd ..

fps = 7500/skip;

cd(PathName)
% cd ..
% 
% movfilename = ['3projections_',num2str(fps),'fps.avi'];
% mov = avifile(movfilename, 'compression', 'none');

solname = [cd,'/'];
% create projection folders
mkdir([solname,'/merged']);

im_nr = 0;
for i=1:skip:length(ims)
    counter = length(ims)-i
    
    for cam=1:length(cams)
        cd(PathName)
        cd(cams(cam).name);
        
        im = imread(ims_all(cam,i).name);
        im = im(x_min:end,1:y_max,:);
%         im = double(im);
%         im = im/max(im(:));
        im_merge(:,(cam-1)*size(im,2)+1:cam*size(im,2),:,:) = im;
    end
    
%     im_merge8bit = uint8(255*im_merge);
%         imshow(im_merge);
% 
%         cd(solname)
%         tic
%         mov = addframe(mov,gcf);
%         toc
% 
%         save(figure(1),[solname 'merged/' ims(i).name])
% im_nr = im_nr +1;
%         imwrite(im_merge,[solname 'merged/merged' num2str(im_nr) '.tif']);
        imwrite(im_merge,[solname 'merged/' ims_all(cam,i).name]);
%         imwrite(im_merge,[solname 'merged/' ims(i).name]);
        cd(PathName)
end

% cd(PathName)
% cd ..
% tic
% mov = close(mov);
% toc
% cd(PathName)

% end