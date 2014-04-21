function make_pathDB1(settings,pathDB)

savefile = 'pathDB1.mat';
file_name = 'fly';
dirname = 'flytracks';

digits = 6;
start_dir = 1;
start_file = 4; % exclude ManualFit


if exist('pathDB') == 1
    t=pathDB.t;
    x=pathDB.x;
    y=pathDB.y;
    z=pathDB.z;
    bodyQ1=pathDB.qb1;
    bodyQ2=pathDB.qb2;
    bodyQ3=pathDB.qb3;
    bodyQ4=pathDB.qb4;
    wingQL1=pathDB.qL1;
    wingQL2=pathDB.qL2;
    wingQL3=pathDB.qL3;
    wingQL4=pathDB.qL4;
    wingQR1=pathDB.qR1;
    wingQR2=pathDB.qR2;
    wingQR3=pathDB.qR3;
    wingQR4=pathDB.qR4;
end

%% construct pathDB
seq = size(x,2);
% alldirs=dir;
alldirs=settings.sequence_names;
for i=start_dir:size(alldirs,1)
    size(alldirs,1)-i
    
%     if alldirs(i).isdir==1
%         cd(alldirs(i).name);
    dir_now = alldirs{i};
    cd(dir_now);
        
        if exist(dirname)==7
            cd(dirname)
            
            seq=seq+1;
%             settings.seq(seq,1) = str2num(alldirs(i).name(1:8));
%             settings.seq(seq,2) = str2num(alldirs(i).name(11:14));
            settings.seq(seq,1) = str2num(dir_now(end-13:end-6));
            settings.seq(seq,2) = str2num(dir_now(end-3:end));
            
            x(1:length(t),seq)=nan;
            y(1:length(t),seq)=nan;
            z(1:length(t),seq)=nan;

            bodyQ1(1:length(t),seq) = nan;
            bodyQ2(1:length(t),seq) = nan;
            bodyQ3(1:length(t),seq) = nan;
            bodyQ4(1:length(t),seq) = nan;
            
            wingQL1(1:length(t),seq)=nan;
            wingQL2(1:length(t),seq)=nan;
            wingQL3(1:length(t),seq)=nan;
            wingQL4(1:length(t),seq)=nan;
            
            wingQR1(1:length(t),seq)=nan;
            wingQR2(1:length(t),seq)=nan;
            wingQR3(1:length(t),seq)=nan;
            wingQR4(1:length(t),seq)=nan;

%             allfiles=dir;
%             for n=start_file:length(allfiles)
            allfiles=dir('fly*');
            for n=1:length(allfiles)
                name_now = allfiles(n).name;
                if exist(name_now)==2
                    load(name_now)

                    frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                    
                    x(frame_now,seq) = xh(1);
                    y(frame_now,seq) = xh(2);
                    z(frame_now,seq) = xh(3);
                    
                    bodyQ1(frame_now,seq) = xh(4);
                    bodyQ2(frame_now,seq) = xh(5);
                    bodyQ3(frame_now,seq) = xh(6);
                    bodyQ4(frame_now,seq) = xh(7);
                    
                    wingQL1(frame_now,seq) = xh(8);
                    wingQL2(frame_now,seq) = xh(9);
                    wingQL3(frame_now,seq) = xh(10);
                    wingQL4(frame_now,seq) = xh(11);
                    
                    wingQR1(frame_now,seq) = xh(12);
                    wingQR2(frame_now,seq) = xh(13);
                    wingQR3(frame_now,seq) = xh(14);
                    wingQR4(frame_now,seq) = xh(15);
                   
                end
            end
            cd ..
        end
%         cd ..
%         cd ..
%     end
end


save(savefile,'t','x','y','z','bodyQ1','bodyQ2','bodyQ3','bodyQ4','wingQL1','wingQL2','wingQL3','wingQL4','wingQR1','wingQR2','wingQR3','wingQR4')



