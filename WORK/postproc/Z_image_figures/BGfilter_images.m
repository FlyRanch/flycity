close all
clc
clear


% settings
save_ims = 0
save_ims = 1

fig_on = 0
% fig_on = 1

crop_on = 0
% crop_on = 1

% name
% fig_name = 'top_filteredREV_60fps'
fig_name = 'top_filteredREV'

dir_now = cd;
dir_name = [fig_name,'_',dir_now(end-21:end-14),dir_now(end-4:end)];
cd ..
mkdir(dir_name);
cd(dir_name)
dir_save = cd;
cd(dir_now)

% im #s
bg_nr = 5080;

% start = 1;
start = 1817;
% start = 2795;

% stop = 5588;
stop = 2795;

skip = 1;
skip = 3;
% skip = 125;
% skip = 10;

x_min = 400;
x_max = 900;
y_min = 400;
y_max = 850;


ls=dir('*.tif')
im_bg = double(imread(ls(bg_nr).name));

% if fig_on == 1
    im = im_bg;
    im = im-min(im(:));
    im = im./max(im(:));
    figure(1),imshow(im)
    pause
% end

im_sum = zeros(size(im_bg));
if crop_on == 1
    im_sum = im_sum(x_min:x_max,y_min:y_max);
end

% for i=1:length(ls)
for i=start:skip:stop
    im = double(imread(ls(i).name));
    
    if crop_on == 1
        im = im(x_min:x_max,y_min:y_max);
    end

    im = im_bg - im;
    im = im-min(im(:));
    im = im./max(im(:));
%     im = imadjust(im,[0.3 1],[]);
    im = imadjust(im,[0.3 .9],[]);
    im = 1-im;
    
    im_sum = im_sum + im;
    
    if fig_on == 1
        figure(1),imshow(im)
        pause
    end
    
    if save_ims == 1
        cd(dir_save)
        nr_now = num2str(i+1000000);
        nr_now = nr_now(2:end);
        imwrite(im,[fig_name,nr_now,'.tif'])
        cd(dir_now)
    end
    
end

% % im_sum = im_sum-min(im_sum(:));
% % im_sum = im_sum./max(im_sum(:));
% % im_sum = imadjust(im_sum,[0.3 1],[]);
% im_sum = imadjust(im_sum,[0 .9],[]);
% figure,imshow(im_sum)
% 
%         cd(dir_save)
% imwrite(im_sum,[fig_name,'_compiled.tif'])
%         cd(dir_now)

