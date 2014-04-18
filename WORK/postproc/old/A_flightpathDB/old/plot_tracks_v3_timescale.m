% clear
% clc
% % close all

file_name = 'fly';
digits = 6;
skip_frame=10;
fps = 7500;
frame_end = 5588;
frames = [1:skip_frame:frame_end]';

trigger_frame_center = 2795;
dN = 0;
trigger_frame = trigger_frame_center;
DBt = (frames-trigger_frame)/fps;

%% trigger ~ end - 0.5sec
trigger_frame_end_delay = 1850; % !!!NOT CONFIRMED!!!
dN = trigger_frame_center - trigger_frame_end_delay;
frames = frames - dN
% trigger_frame = trigger_frame_end_delay;
seq = 6
%% construct DB

start_dir = 3;
alldirs=dir;

% DBx=nan(size(frames));
% DBy=nan(size(frames));
% DBz=nan(size(frames));
% 
% seq = 0;
for i=start_dir:size(alldirs)
    if alldirs(i).isdir==1
        cd(alldirs(i).name);
        
        if exist('flytracks')==7
            cd('flytracks')
            
            seq=seq+1;
            
            DBx(:,seq)=nan;
            DBy(:,seq)=nan;
            DBz(:,seq)=nan;
            
            for n=1:length(frames)
                
                frame_now = frames(n);
                name_now = [file_name,sprintf(['%0',num2str(digits),'d'], frame_now),'.mat'];
                
                if exist(name_now) == 2
                    
                    load(name_now)
                    
                    DBx(n,seq) = xh(1);
                    DBy(n,seq) = xh(2);
                    DBz(n,seq) = xh(3);
                end
            end
            cd ..
        end
        cd ..
    end
end

plot3(DBx,DBy,DBz)
plot3(DBx(280,:),DBy(280,:),DBz(280,:),'r*')
hold on
plot3(DBx,DBy,DBz)
axis equal
grid on

%% make movie

movfilename = '3Dflightpaths.avi';
mov = avifile(movfilename, 'compression', 'none');


start = min(find(nansum(DBx,2)~=0));
stop =  max(find(nansum(DBx,2)~=0));

figure
hold on
axis equal
grid on
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
xlabel('x')
ylabel('y')
zlabel('z')

for k = start:stop
    
    if DBt(k)<0
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
    else
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
    end
    
    mov = addframe(mov,gcf);
end

mov = close(mov);


%% make 4 views movie

movfilename = '3Dflightpaths_4views.avi';
mov = avifile(movfilename, 'compression', 'none');


start = min(find(nansum(DBx,2)~=0));
stop =  max(find(nansum(DBx,2)~=0));

figure

subplot(2,2,1)
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
view(1);
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,2,2)
view(2);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);

subplot(2,2,3)
view(2);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);

subplot(2,2,4)
view(2);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);

hold on
axis equal
grid on

for k = start:stop
    
    if DBt(k)<0
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
    else
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
    end
    
    mov = addframe(mov,gcf);
end

mov = close(mov);


