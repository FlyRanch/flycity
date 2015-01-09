%close all;
clear
clc;
close all
warning off

cutoff_x = 69
cutoff_y = 852

PathName = [cd,'/'];
cams = dir('cam*');

cd(cams(1).name);
ims = dir('fly*');
im_templ = imread(ims(1).name);
im_templ = im_templ(cutoff_x:end,1:cutoff_y,:);
im_x = size(im_templ,1);
im_y = size(im_templ,2);
im_merge = [im_templ im_templ im_templ];
cd ..

if exist('skip')==0
    skip = 3
end
fps = 7500/skip;

cd(PathName)
cd ..

movfilename = ['3projections_',num2str(fps),'fps.avi'];
mov = avifile(movfilename, 'compression', 'none');
solname = [cd,'/'];


for i=1:length(ims)
    for cam=1:length(cams)

        cd(PathName)
        cd(cams(cam).name);
        im = imread(ims(i).name);
        im = im(cutoff_x:end,1:cutoff_y,:);
        
        im_merge(:,(cam-1)*im_y+1:cam*im_y,:) = im;
    end
        
        figure(1)
        hold off
        imshow(im_merge);
        pause(0.001)
    
        cd(solname)
        tic
        mov = addframe(mov,gcf);
        toc
        cd(PathName)
end

cd(PathName)
cd ..
tic
mov = close(mov);
toc
cd(PathName)

% end