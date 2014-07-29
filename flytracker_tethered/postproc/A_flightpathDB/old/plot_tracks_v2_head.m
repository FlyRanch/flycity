clear
clc
% close all

trigger_frame_end_delay = 1850; % !!!NOT CONFIRMED!!!
trigger_frame_center = 2795;

trigger_frame = trigger_frame_end_delay;
trigger_frame = trigger_frame_center;
start_dir = 3;
start_frame = 2;
skip_frame=10;

alldirs=dir;
% for i=3:size(alldirs)-2
for i=start_dir:size(alldirs)
    if alldirs(i).isdir==1
        cd(alldirs(i).name);
        
        if exist('flytracks')==7
            cd('flytracks')
            allfiles=dir('*.mat');

            for n=start_frame:skip_frame:length(allfiles)
                load(allfiles(n).name)
                
                path(n-1,:) = xh(1:3)';
                nr(n-1,:) = str2num(allfiles(n).name(4:9));
                
                hold on
                
                if nr(n-1,:) < trigger_frame
                    plot3(xh(1),xh(2),xh(3),'.b')
                else
                    plot3(xh(1),xh(2),xh(3),'.r')
                end
            end
            
            clear path
            cd ..
        end
        cd ..
    end
end
axis equal
grid on