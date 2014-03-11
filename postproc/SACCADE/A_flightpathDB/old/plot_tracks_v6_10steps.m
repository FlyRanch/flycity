% make movie

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
