% clear
clc
close all

savefile = 'flightpathDB_noskip.mat';
file_name = 'fly';
dirname = 'flytracks';

digits = 6;
% skip_frame=10;
skip_frame=1;
fps = 7500;
frame_end = 5588;
frames = [1:skip_frame:frame_end]';

start_dir = 3;
start_file = 4; % exclude ManualFit
% skip_file=10;
skip_file=1;

trigger_frame_center = 2795;
frame_center = round(trigger_frame_center/skip_frame);
dN = 0;
trigger_frame = trigger_frame_center;
DBt = (frames-trigger_frame)/fps;

DBx=nan(size(frames));
DBy=nan(size(frames));
DBz=nan(size(frames));

seq = 0;

%% trigger ~ end - 0.5sec

% trigger_frame_end_delay = 1850; % !!!NOT CONFIRMED!!!
% dN = trigger_frame_center - trigger_frame_end_delay;
% dN_closest = (round(dN/skip_frame))*skip_frame;
% frames = frames - dN_closest;

%% construct DB

alldirs=dir;

for i=start_dir:length(alldirs)
    length(alldirs)-i
    if alldirs(i).isdir==1
        cd(alldirs(i).name);
        
        if exist(dirname)==7
            cd(dirname)
            
            seq=seq+1;
            DBseq(seq,1) = str2num(alldirs(i).name(1:8));
            DBseq(seq,2) = str2num(alldirs(i).name(11:14));
            
            DBx(:,seq)=nan;
            DBy(:,seq)=nan;
            DBz(:,seq)=nan;
            
            allfiles=dir;
            for n=start_file:skip_file:length(allfiles)
                name_now = allfiles(n).name;
                if exist(name_now)==2
                    load(name_now)

                    frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                    closest_frame = (round(frame_now/skip_frame))*skip_frame;
                    loc = find(frames==closest_frame);
                    
                    DBx(loc,seq) = xh(1);
                    DBy(loc,seq) = xh(2);
                    DBz(loc,seq) = xh(3);
                    
%                     DBx(loc,seq) = InternalVariablesDS.xh_(1);
%                     DBy(loc,seq) = InternalVariablesDS.xh_(2);
%                     DBz(loc,seq) = InternalVariablesDS.xh_(3);
                end
            end
            cd ..
        end
        cd ..
    end
end

save(savefile,'DBseq','DBt','DBx','DBy','DBz')