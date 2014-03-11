%close all;
clear
clc;
close all
warning off

skip = 1
fps = 7500/skip;

start = 500
stop = 1600

% crop settings: pooing fruitfly (640x360)
x_min = 250
x_max = x_min + 640 -1
y_min = 500
y_max = y_min + 360 -1

PathName = [cd,'/'];
ims = dir('*.tif');


cd(PathName)
solname = [cd,'/subset'];
mkdir(solname);

% im_nr = 0;
for i=start:skip:stop
    COUNTER = stop -i
    
    im = imread(ims(i).name);
    im = im(y_min:y_max,x_min:x_max,1);                            % grayscale movie: remove matrix 2&3

%     im = double(im);
%     im = im/max(im(:));
%     im = uint8(255*im);
    
%     im = imadjust(im);
%     imshow(im);
% 
%     cd(solname)
%     tic
%     mov = addframe(mov,gcf);
%     toc
% 
    imwrite(im,[solname '/' ims(i).name]);
%         cd(PathName)
end

% cd(PathName)
% cd ..
% tic
% mov = close(mov);
% toc
% cd(PathName)
