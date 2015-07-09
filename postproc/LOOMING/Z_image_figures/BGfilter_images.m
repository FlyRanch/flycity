close all
clc
clear


% settings
save_ims = 0
save_ims = 1

fig_on = 0
% fig_on = 1

crop_on = 0
crop_on = 1

% name
dir_now = cd;

camnr = dir_now(end-9);
fig_name = ['cam' camnr '_filteredREV']
dir_name = [fig_name,'_',dir_now(end-21:end-14),dir_now(end-4:end)];

cd ..
mkdir(dir_name);
cd(dir_name)
dir_save = cd;
cd(dir_now)

% im #s
bg_nr = 1;
% bg_nr = 5588;
% bg_nr = 5589;

start = 1;
% start = 1000;
start = 795;

stop = 5588;
stop = 3000;

skip = 1;
% skip = 3;
% skip = 125;
% skip = 100;

%% clipped
if str2num(camnr) == 1
    x_min = 101;
    x_max = 900;
    y_min = 201;
    y_max = 800;

elseif str2num(camnr) == 2
    x_min = 201;
    x_max = 600;
    y_min = 601;
    y_max = 900;

elseif str2num(camnr) == 3
    x_min = 1;
    x_max = 1000;
    y_min = 1;
    y_max = 900;
end


ls=dir('*.tif')
im_bg = double(imread(ls(bg_nr).name));
if crop_on == 1
    im_bg = im_bg(y_min:y_max,x_min:x_max);
end

% if fig_on == 1
    im = im_bg;
    im = im-min(im(:));
    im = im./max(im(:));
    figure(1),imshow(im)
    pause
% end

im_sum = zeros(size(im_bg));
% if crop_on == 1
%     im_sum = im_sum(y_min:y_max,x_min:x_max);
% end

% for i=1:length(ls)
for i=start:skip:stop
    
    stop - i
    
    im = double(imread(ls(i).name));
    
    if crop_on == 1
        im = im(y_min:y_max,x_min:x_max);
    end

    im = im_bg - im;
    im = im-min(im(:));
    im = im./max(im(:));
%     im = imadjust(im,[.3 1],[]);
%     im = imadjust(im,[.3 .9],[]);
    im = 1-im;
    
    im = imadjust(im,[.2 .7],[]);
    
%     im = imadjust(im,[.2 .8],[]);
%     im = medfilt2(im,[2 2]);
    
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

