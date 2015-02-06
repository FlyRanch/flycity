%close all;
clear
clc;
close all
warning off

skip = 3
start = 2794
stop = 4063
stop = 2810

PathName = [cd,'/'];
cams = dir('*002');

for i=1:length(cams)
    cd(cams(i).name);
    ims = dir('*.tif');
    ims_all(i,1:length(ims)) = ims;
    cd(PathName)
end
    
cd(cams(1).name);
ims = dir('*.tif');
im_templ = imread(ims(1).name);
im_x = size(im_templ,1);
im_y = size(im_templ,2);
im_merge = [im_templ im_templ im_templ];
cd ..

fps = 7500/skip;

cd(PathName)
% cd ..
% 
movfilename = ['3views_',num2str(fps),'fps.avi'];
mov = avifile(movfilename, 'compression', 'none');
solname = [cd,'/'];

% create projection folders
mkdir([solname,'/merged']);

im_nr = 0;
for i=start:skip:stop
    for cam=1:length(cams)
        cd(PathName)
        cd(cams(cam).name);
        
        im = imread(ims_all(cam,i).name);
        im_merge(:,(cam-1)*im_y+1:cam*im_y,:) = im;
    end
    
    im_merge8bit = uint8(255*double(im_merge)/65535);
    for j=1:3
    im_movie(:,:,j) = im_merge8bit;
    end

%     
%         imshow(im_merge, [min(im_merge(:)) max(im_merge(:))]);
% 
%         cd(solname)
%         tic
%         mov = addframe(mov,gcf);
%         toc

%         save(figure(1),[solname 'merged/' ims(i).name])
im_nr = im_nr +1;
M(im_nr) = im2frame(uint8(255*double(im_merge)/65535));
%         imwrite(im_merge,[solname 'merged/merged' num2str(im_nr) '.tif']);
%         imwrite(im_merge,[solname 'merged/' ims(i).name]);
        cd(PathName)
end

cd(PathName)
tic
mov = close(mov);
toc
