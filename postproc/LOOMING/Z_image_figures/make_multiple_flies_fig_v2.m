% close all
clc
clear


crop_on = 0
% crop_on = 1

fig_on = 0

% name
fig_name = 'top_compiled'
dir_now = cd;
fig_name = [fig_name,'_',dir_now(end-11:end)];

% % im #s
% trig = 2795
% 
% % % 20121120S01 #4
% % frames = [trig-400 trig-220 trig+5 trig+315 trig+505 trig+600 trig+678 trig+740];
% 
% % % 20121120S04 #6
% % frames = [2192 2411 2621 2795 2996 3203 3401 3518 3614 3707];
% 
% % % 20121120S05 #XXX
% % frames = [trig-115 trig-3 trig+115 trig+235 trig+400 trig+525 trig+600 trig+675 trig+750];
% 
% % 20121129S03 #8
% frames = [trig-1770 trig-1590 trig-1270 trig-693 trig+350 trig+500 trig+630 trig+774];
% 
% % % 20121205S01 #12
% % frames = [1895 2175 2485 2797 3085 3275 3365 3455];
% 
% % % 20121205S03 #14
% % frames = [trig-300 trig-150 trig trig+155 trig+290 trig+465 trig+600 trig+720 trig+844];
% 
% % % 20121205S06 #37
% % frames = [trig-800 trig-400 trig trig+410 trig+660 trig+820 trig+950];
% 
% crop pixels
x_min = 400;
x_max = 900;
y_min = 400;
y_max = 850;

ls=dir('*.tif');
im = double(imread(ls(1).name));
im_sum = zeros(size(im));

for i=1:length(ls)
    im = double(imread(ls(i).name));
% 
% for i=1:length(frames)
%     im = double(imread(ls(frames(i)).name));
    
    im_sum = im_sum + im;
    
    if fig_on == 1
        figure(1),imagesc(im),colormap(gray)
        pause
    end
end

if crop_on == 1
    im_sum = im_sum(x_min:x_max,y_min:y_max);
end

% figure
hold off
imagesc(im_sum)
colormap(gray)
axis equal
pause

im_sum = im_sum-min(im_sum(:));
im_sum = im_sum./max(im_sum(:));
im_sum2 = imadjust(im_sum,[0.4 1],[]);
% im_sum = imadjust(im_sum,[0 .9],[]);
im_sum2= medfilt2(im_sum2,[2 2]);
figure,imshow(im_sum2)

cd ..
imwrite(im_sum,[fig_name,'.tif'])
save([fig_name,'_frame_numbers'])
cd(dir_now)

