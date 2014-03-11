function [settings,pathDB] = make_flightpathDB_v3_subsets(settings,pathDB,expansion)

file_name = 'fly';
dirname = 'flytracks';

digits = 6;
start_dir = 3;
start_file = 4; % exclude ManualFit

t=pathDB.t;
x=pathDB.x_obs;
y=pathDB.y_obs;
z=pathDB.z_obs;


%% construct pathDB
seq = size(x,2);
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
            settings.expansion.pos(seq,1) = expansion.pos;
            
            x(1:length(t),seq)=nan;
            y(1:length(t),seq)=nan;
            z(1:length(t),seq)=nan;
            
            allfiles=dir;
            for n=start_file:length(allfiles)
                name_now = allfiles(n).name;
                if exist(name_now)==2
                    load(name_now)

                    frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                    
                    x(frame_now,seq) = xh(1);
                    y(frame_now,seq) = xh(2);
                    z(frame_now,seq) = xh(3);
%                     
%                     x(frame_now,seq) = InternalVariablesDS.xh_(1);
%                     y(frame_now,seq) = InternalVariablesDS.xh_(2);
%                     z(frame_now,seq) = InternalVariablesDS.xh_(3);
                end
            end
            cd ..
        end
        cd ..
    end
end

pathDB.x_obs = x;
pathDB.y_obs = y;
pathDB.z_obs = z;
