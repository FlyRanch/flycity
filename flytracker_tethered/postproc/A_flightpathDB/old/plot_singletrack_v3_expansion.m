% clear
% clc
close all

% green ring
center = [x(find(t==min(abs(t)))) y(find(t==min(abs(t)))) z(find(t==min(abs(t))))];
teta = 0:1:360;
r = 0.120; % [m] panel-fly distance

ring_x = r * cosd(teta);
ring_y = r * sind(teta);

ring_x = center(1)+ring_x;
ring_y = center(2)+ring_y;

DN_trig = 400
DN_resp = 200

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
% trigger_frame = find(DBt == min(abs(DBt)));

t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp * fps)
DN_escape = round(t_escape * fps)

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



% make movie

movfilename = '3Dflightpaths_expansion.avi';
mov = avifile(movfilename, 'compression', 'none');

start = min(find(nansum(DBx,2)~=0));
stop =  max(find(nansum(DBx,2)~=0));
seqs = size(DBx,2);

figure
hold on
axis equal
grid on
set(gca,'xlim',[min(DBx(:)) max(DBx(:))],'ylim',[min(DBy(:)) max(DBy(:))],'zlim',[min(DBz(:)) max(DBz(:))]);
xlabel('x')
ylabel('y')
zlabel('z')

for k = start:movie_step:stop
    
%     if DBt(k)<0
%         plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
%     else
%         plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
%     end
%     
    for l = 1:seqs
        
        last_frame = max(find(isnan(DBy(:,l))==0));
        if k<trigger_frame
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'b.')
        elseif k<trigger_frame + DN_trig
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'r.')
        elseif k<last_frame - DN_resp
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'y.')
        else
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'g.')
        end
    end

    tic
    mov = addframe(mov,gcf);
    toc
end

mov = close(mov);


%% make 4 views movie

% movfilename = '3Dflightpaths_4views.avi';
% mov = avifile(movfilename, 'compression', 'none');
% 
% 
% start = min(find(nansum(DBx,2)~=0));
% stop =  max(find(nansum(DBx,2)~=0));
% seqs = size(DBx,2);
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
% for k = start:movie_step:stop
% 
%     for l = 1:seqs
%         last_frame = max(find(isnan(DBy(:,l))==0));
% 
%         if k<trigger_frame
%             for s = 1:4
%                 subplot(2,2,s)
%                 plot3(DBx(k,l),DBy(k,l),DBz(k,l),'b.')
%             end
%         elseif k<trigger_frame + DN_trig
%             for s = 1:4
%                 subplot(2,2,s)
%                 plot3(DBx(k,l),DBy(k,l),DBz(k,l),'r.')
%             end
%         elseif k<last_frame - DN_resp
%             for s = 1:4
%                 subplot(2,2,s)
%                 plot3(DBx(k,l),DBy(k,l),DBz(k,l),'y.')
%             end
%         else
%             for s = 1:4
%                 subplot(2,2,s)
%                 plot3(DBx(k,l),DBy(k,l),DBz(k,l),'g.')
%             end
%         end
%     end
% %         for s = 1:4
% %             subplot(2,2,s)
% %             plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
% %         end
% %     else
% %         for s = 1:4
% %             subplot(2,2,s)
% %             plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
% %         end
% %     end
%     
%     tic
%     mov = addframe(mov,gcf);
%     toc
% end
% 
% 
% mov = close(mov);
% 
% 
