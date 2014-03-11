% clc
% clear
% close all
% load('flightpathDB_10steps_all.mat')
% 
trigger_frame = find(DBt == min(abs(DBt)));

dt = DBt(2) - DBt(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)


%% make top view movie

movie_step = 1;
movfilename = '3Dflightpaths_topview.avi';
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

    for l = 1:seqs
        last_frame = max(find(isnan(DBy(:,l))==0));
        if k<trigger_frame
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'b.')
        elseif k<trigger_frame + DN_noresp
            plot3(DBx(k,l),DBy(k,l),DBz(k,l),'r.')
        elseif k<last_frame - DN_escape
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

% movie_step = 1;
% movfilename = '3Dflightpaths_4views.avi';
% mov = avifile(movfilename, 'compression', 'none');
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
%         elseif k<trigger_frame + DN_noresp
%             for s = 1:4
%                 subplot(2,2,s)
%                 plot3(DBx(k,l),DBy(k,l),DBz(k,l),'r.')
%             end
%         elseif k<last_frame - DN_escape
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
