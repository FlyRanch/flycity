function [pathDB,settings] = make_flightpathDB_v8(pathDB,settings,expansion,incq)

file_name = 'fly';
dirname = 'flytracks';

digits = 6;
start_dir = 3;
start_file = 2; % exclude ManualFit

if exist('pathDB') == 1
    t = pathDB.t;
    pos = pathDB.pos_obs;
    
    if incq == 1
        qbody = pathDB.qbody_obs;
        qwingL = pathDB.qwingL_obs;
        qwingR = pathDB.qwingR_obs;
    end
end

%% construct pathDB
seq = size(pos,2);
alldirs=dir;
for i=start_dir:length(alldirs)
    length(alldirs)-i
    if alldirs(i).isdir==1
        cd(alldirs(i).name);
        
        if exist(dirname)==7
            cd(dirname)
            
            seq=seq+1;
            settings.seq(seq,1) = str2num(alldirs(i).name(1:8));
            settings.seq(seq,2) = str2num(alldirs(i).name(11:14));
            settings.expansion.speed(seq,1) = expansion.speed;
            settings.expansion.VerPos(seq,1) = expansion.VerPos;
            settings.expansion.HorPos(seq,1) = expansion.HorPos;
            settings.expansion.stepwise(seq,1) = expansion.stepwise;
            settings.expansion.maxangle(seq,1) = expansion.maxangle;
            settings.expansion.arista(seq,1) = expansion.arista;
            settings.expansion.OFF(seq,1) = expansion.OFF;
            
            pos(1:length(t),seq,1:3)=nan;
            
            if incq == 1
                qbody(1:length(t),seq,1:4)=nan;
                qwingL(1:length(t),seq,1:4)=nan;
                qwingR(1:length(t),seq,1:4)=nan;
            end

            allfiles=dir('*.mat');
            for n=start_file:length(allfiles)
                name_now = allfiles(n).name;
                if exist(name_now)==2
                    load(name_now)

                    frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                    
                    pos(frame_now,seq,1:3) = xh(1:3);
                    
                    if incq == 1
                        qbody(frame_now,seq,1:4) = xh(4:7);
                        qwingL(frame_now,seq,1:4) = xh(8:11);
                        qwingR(frame_now,seq,1:4) = xh(12:15);
                    end
                end
            end
            cd ..
        end
        cd ..
    end
end

pathDB.pos_obs = pos;
if incq == 1
    pathDB.qbody_obs = qbody;
    pathDB.qwingL_obs = qwingL;
    pathDB.qwingR_obs = qwingR;
end
