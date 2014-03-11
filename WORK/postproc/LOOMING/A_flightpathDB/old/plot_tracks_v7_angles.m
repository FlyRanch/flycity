% make movie
DN_trig = 40
DN_resp = 20
trigger_frame = 280

movfilename = '3Dflightpaths_angles.avi';
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

for k = start:stop
    
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

%     if DBt(k)<0
%         plot3(DBx(k,:),DBy(k,:),DBz(k,:),'b.')
%     else
%         plot3(DBx(k,:),DBy(k,:),DBz(k,:),'r.')
%     end
    
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
