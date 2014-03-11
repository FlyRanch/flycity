function [settings,pathDB] = make_flightpathDB_v4_orientations(settings,pathDB,expansion)

file_name = 'fly';
dirname = 'flytracks';

digits = 6;
start_dir = 3;
start_file = 4; % exclude ManualFit

if exist('pathDB') == 1
    t=pathDB.t;
    x=pathDB.x_obs;
    y=pathDB.y_obs;
    z=pathDB.z_obs;

    roll=pathDB.roll_obs;
    pitch=pathDB.pitch_obs;
    yaw=pathDB.yaw_obs;

    wingL1=pathDB.wingL1_obs;
    wingL2=pathDB.wingL2_obs;
    wingL3=pathDB.wingL3_obs;

    wingR1=pathDB.wingR1_obs;
    wingR2=pathDB.wingR2_obs;
    wingR3=pathDB.wingR3_obs;
end

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

            roll(1:length(t),seq)=nan;
            pitch(1:length(t),seq)=nan;
            yaw(1:length(t),seq)=nan;

            wingL1(1:length(t),seq)=nan;
            wingL2(1:length(t),seq)=nan;
            wingL3(1:length(t),seq)=nan;

            wingR1(1:length(t),seq)=nan;
            wingR2(1:length(t),seq)=nan;
            wingR3(1:length(t),seq)=nan;

            allfiles=dir;
            for n=start_file:length(allfiles)
                name_now = allfiles(n).name;
                if exist(name_now)==2
                    load(name_now)

                    frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                    
                    x(frame_now,seq) = xh(1);
                    y(frame_now,seq) = xh(2);
                    z(frame_now,seq) = xh(3);
                    
                    [roll_now,pitch_now,yaw_now] = quat2angle(xh(4:7)');
                    roll(frame_now,seq) = rad2deg(roll_now);
                    pitch(frame_now,seq) = rad2deg(pitch_now);
                    yaw(frame_now,seq) = rad2deg(yaw_now);
                    
                    [wingL1_now,wingL2_now,wingL3_now] = quat2angle(xh(8:11)');
                    wingL1(frame_now,seq) = rad2deg(wingL1_now);
                    wingL2(frame_now,seq) = rad2deg(wingL2_now);
                    wingL3(frame_now,seq) = rad2deg(wingL3_now);
                    
                    [wingR1_now,wingR2_now,wingR3_now] = quat2angle(xh(12:15)');
                    wingR1(frame_now,seq) = rad2deg(wingR1_now);
                    wingR2(frame_now,seq) = rad2deg(wingR2_now);
                    wingR3(frame_now,seq) = rad2deg(wingR3_now);
                    
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

pathDB.roll_obs = roll;
pathDB.pitch_obs = pitch;
pathDB.yaw_obs = yaw;

pathDB.wingL1_obs = wingL1;
pathDB.wingL2_obs = wingL2;
pathDB.wingL3_obs = wingL3;

pathDB.wingR1_obs = wingR1;
pathDB.wingR2_obs = wingR2;
pathDB.wingR3_obs = wingR3;

