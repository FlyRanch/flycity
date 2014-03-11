clear
clc
close all

file_name = 'fly';
digits = 6;
skip_frame=1;
movie_step=10;
fps = 7500;
frame_end = 5588;
frames = [1:skip_frame:frame_end]';

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
% frames = frames - dN;
% % trigger_frame = trigger_frame_end_delay;
% % seq = 6

%% construct DB

% start_dir = 3;
% alldirs=dir;
% 
% for i=start_dir:size(alldirs)
%     if alldirs(i).isdir==1
%         cd(alldirs(i).name);
        
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
%         cd ..
%     end
% end
% 
% plot3(DBx,DBy,DBz)
% plot3(DBx(frame_center,:),DBy(frame_center,:),DBz(frame_center,:),'r*')
% hold on
% plot3(DBx,DBy,DBz,'.')
% axis equal
% grid on

%% make movie

movfilename = '3Dflightpaths_internal.avi';
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

for k = start:movie_step:stop
    
    if DBt(k)<0
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
    else
        plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
    end
    
    mov = addframe(mov,gcf);
end

mov = close(mov);


% %% make 4 views movie
% 
% movfilename = '3Dflightpaths_4views.avi';
% mov = avifile(movfilename, 'compression', 'none');
% 
% 
% start = min(find(nansum(DBx,2)~=0));
% stop =  max(find(nansum(DBx,2)~=0));
% 
% figure
% 
% subplot(2,2,1)
% view(2);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% axis equal
% grid on
% set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
% 
% subplot(2,2,2)
% view(90,0);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% axis equal
% grid on
% set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
% 
% subplot(2,2,3)
% view(0,0);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% axis equal
% grid on
% set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
% 
% subplot(2,2,4)
% view(3);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% axis equal
% grid on
% set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
% 
% 
% for k = start:stop
%     
%     if DBt(k)<0
%         for s = 1:4
%             subplot(2,2,s)
%             plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
%         end
%     else
%         for s = 1:4
%             subplot(2,2,s)
%             plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
%         end
%     end
%     
%     mov = addframe(mov,gcf);
% end
% 
% mov = close(mov);
% 
% 
