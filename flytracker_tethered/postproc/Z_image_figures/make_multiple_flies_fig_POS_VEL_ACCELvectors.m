% make vel & accel vectors for top projection images

close all
clear
clc

% settings
scale_vel = 1/200;
scale_accel = 1/5000;

color_vel = [0,1,1];
color_accel = [1,0,0];

vec_width = 1;


% with or without roll
roll_on = 1
% roll_on = 0

if roll_on == 1
    dbname = dir('*INCroll*');
else
    dbname = dir('*NOroll*');
end
dbname = dbname.name;
load(dbname)

list = dir('top*mat')
for i = 1:length(list)
    seq_name = list(i).name;
    load(seq_name)
    date = str2num(dir_now(end-12:end-5));
    seq = str2num(dir_now(end-3:end));
    
    seq_nr = find(settings.seq(:,1)==date & settings.seq(:,2)==seq);
    
    hold off
    for j = 1:length(frames)
        frame_nr = frames(j);
        
%         quiver(pathDB.pos(frame_nr,seq_nr,1),pathDB.pos(frame_nr,seq_nr,2),...
%             pathDB.vel(frame_nr,seq_nr,1),pathDB.vel(frame_nr,seq_nr,2),scale_vel,'color',color_vel,'linewidth',vec_width)
%         hold on
%         quiver(pathDB.pos(frame_nr,seq_nr,1),pathDB.pos(frame_nr,seq_nr,2),...
%             pathDB.accel(frame_nr,seq_nr,1),pathDB.accel(frame_nr,seq_nr,2),scale_accel,'color',color_accel,'linewidth',vec_width)
        
        plot([pathDB.pos(frame_nr,seq_nr,1),pathDB.pos(frame_nr,seq_nr,1)+pathDB.vel(frame_nr,seq_nr,1)*scale_vel],...
            [pathDB.pos(frame_nr,seq_nr,2),pathDB.pos(frame_nr,seq_nr,2)+pathDB.vel(frame_nr,seq_nr,2)*scale_vel],...
            'color',color_vel,'linewidth',vec_width)
        hold on
        plot([pathDB.pos(frame_nr,seq_nr,1),pathDB.pos(frame_nr,seq_nr,1)+pathDB.accel(frame_nr,seq_nr,1)*scale_accel],...
            [pathDB.pos(frame_nr,seq_nr,2),pathDB.pos(frame_nr,seq_nr,2)+pathDB.accel(frame_nr,seq_nr,2)*scale_accel],...
            'color',color_accel,'linewidth',vec_width)
    end
%     quiver(0,0,1,0,scale_vel,'color',color_vel,'linewidth',vec_width)
%     quiver(0,0,0,1,scale_accel,'color',color_accel,'linewidth',vec_width)
    plot([0,.005],[0,0],'k','linewidth',vec_width)
    plot([0,1*scale_vel],[0,0],'color',color_vel,'linewidth',vec_width)
    plot([0,25*scale_accel],[0,0],'color',color_accel,'linewidth',vec_width)
    axis equal
    axis([-.01 .04 -.01 .04])
    axis off
    
    saveas(gca, [seq_name(1:end-4),'_vectors.fig'])
    saveas(gca, [seq_name(1:end-4),'_vectors.tif'])
    plot2svg([seq_name(1:end-4),'_vectors.svg'])
end
