clc
clear
close all

load('flightpathDB.mat')
% load('expansion_steps_DOUBLESPEED.mat')
load('expansion_steps.mat')

t = pathDB.t;
x = nan(size(t));
y = nan(size(t));
z = nan(size(t));

%% reverse tracks for expansion direction = 180 deg
j=0;
seqs = length(settings.expansion.speed);
for i = 1:seqs
    if settings.expansion.speed(i) == 0 && settings.expansion.HorPos(i) == 0
        j = j +1;
        speed(j,1) = settings.expansion.speed(i);
        HorPos(j,1) = settings.expansion.HorPos(i);
        VerPos(j,1) = settings.expansion.VerPos(i);
        
%         x(:,j) = pathDB.x_obs(:,i)/1000;
%         y(:,j) = pathDB.y_obs(:,i)/1000;
%         z(:,j) = pathDB.z_obs(:,i)/1000;
        
        x(:,j) = pathDB.x(:,i);
        y(:,j) = pathDB.y(:,i);
        z(:,j) = pathDB.z(:,i);
    end
end

trigger_frame = find(t == min(abs(t)));

dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

%% panel ring
% center = [x(find(t==min(abs(t)))) y(find(t==min(abs(t)))) z(find(t==min(abs(t))))];
center = [nanmean(x(trigger_frame,:)) nanmean(y(trigger_frame,:)) nanmean(z(trigger_frame,:))];
teta_ring = 0:1:360;
% r = 0.120; % [m] panel-fly distance
r = 0.028; % closely around flies

x_ring = r * cosd(teta_ring);
y_ring = r * sind(teta_ring);

x_ring = center(1)+x_ring;
y_ring = center(2)+y_ring;

% expansion step frames
exp_steps(:,5)= round(exp_steps(:,4)*7500 + trigger_frame);

%% make movies
sideview = true
fourviews = false

%% make side view movie

if sideview == true
    movie_step = 15;
    fps = 7500/15;
    movfilename = ['3Dflightpaths_sideview_slow_vertOrigin_',num2str(fps),'fps.avi'];
%     movfilename = ['3Dflightpaths_sideview_slow_NOrevexp_',num2str(fps),'fps.avi'];
    mov = avifile(movfilename, 'compression', 'none');

    start = find(nansum(x,2)~=0, 1 );
    stop =  find(nansum(x,2)~=0, 1, 'last' );
    seqs = size(x,2);

    figure
    hold on
    axis equal
    view(90,0);
% %     grid on
%     plot(x_ring,y_ring,'-g','linewidth',5)
%     set(gca,'xlim',[min(x_ring(:)) max(x_ring(:))],'ylim',[min(y_ring(:)) max(y_ring(:))]);
% %     xlabel('x')
% %     ylabel('y')
% %     zlabel('z')
    hold on
    axis equal
    grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);
    
 set(gca, 'XTickLabel', []);
 set(gca, 'YTickLabel', []);
 set(gca, 'ZTickLabel', []);

 set(gca, 'XTick', []);
 set(gca, 'YTick', []);
 set(gca, 'ZTick', []);

    exp_frame = 1;
    for k = start:movie_step:stop
        counter = stop-k
% 
%         if k<trigger_frame
%             plot3(x(k,:),y(k,:),z(k,:),'bo','MarkerSize',3,'MarkerFaceColor','b')
%         else
%             plot3(x(k,:),y(k,:),z(k,:),'ro','MarkerSize',3,'MarkerFaceColor','r')
%         end
        
        for i = 1:size(x,2)
            if VerPos(i)==-1
                plot3(x(k,i),y(k,i),z(k,i),'bo','MarkerSize',3,'MarkerFaceColor','b')
            elseif VerPos(i)==0
                plot3(x(k,i),y(k,i),z(k,i),'ro','MarkerSize',3,'MarkerFaceColor','r')
            elseif VerPos(i)==1
                plot3(x(k,i),y(k,i),z(k,i),'go','MarkerSize',3,'MarkerFaceColor','g')
            end
        end
        
%         for l = 1:seqs
%             last_frame = find(isnan(y(:,l))==0, 1, 'last' );
%             if k<trigger_frame
%                 plot3(x(k,l),y(k,l),z(k,l),'b.')
%             elseif k<trigger_frame + DN_noresp
%                 plot3(x(k,l),y(k,l),z(k,l),'r.')
%             elseif k<last_frame - DN_escape
%                 plot3(x(k,l),y(k,l),z(k,l),'y.')
%             else
%                 plot3(x(k,l),y(k,l),z(k,l),'r.')
%             end
%         end
%             
%             %% exp pattern
%             if exp_frame <= length(exp_steps)
%                 if k >= exp_steps(exp_frame,5)
%                     teta_exp = -exp_steps(exp_frame,3)/2:exp_steps(exp_frame,3)/2;
%                     
%                     % exp origin at 0 deg
%                     teta_exp = teta_exp - 90;
% %                     % exp origin at 180 deg
% %                     teta_exp = teta_exp + 90;
% 
%                     x_exp = r * cosd(teta_exp);
%                     y_exp = r * sind(teta_exp);
% 
%                     x_exp = center(1)+x_exp;
%                     y_exp = center(2)+y_exp;
% 
%                     plot(x_exp,y_exp,'-k','linewidth',5)
% 
%                     exp_frame = exp_frame +1;
%                 end
%             end
% 
%         tic
%         mov = addframe(mov,gcf);
%         toc;
        mov = addframe(mov,gcf);
    end

    mov = close(mov);
end

%% make 4 views movie

if fourviews == true
    movie_step = 1;
    movfilename = '3Dflightpaths_4views_slow_NOrevexp.avi';
    mov = avifile(movfilename, 'compression', 'none');

    start = find(nansum(x,2)~=0, 1 );
    stop =  find(nansum(x,2)~=0, 1, 'last' );
    seqs = size(x,2);

    figure

    subplot(2,2,1)
    view(2);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    axis equal
    grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);

    subplot(2,2,2)
    view(90,0);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    axis equal
    grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);

    subplot(2,2,3)
    view(0,0);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    axis equal
    grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);

    subplot(2,2,4)
    view(3);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    axis equal
    grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);


    for k = start:movie_step:stop

        for l = 1:seqs
            last_frame = max(find(isnan(y(:,l))==0));

            if k<trigger_frame
                for s = 1:4
                    subplot(2,2,s)
                    plot3(x(k,l),y(k,l),z(k,l),'b.')
                end
            elseif k<trigger_frame + DN_noresp
                for s = 1:4
                    subplot(2,2,s)
                    plot3(x(k,l),y(k,l),z(k,l),'r.')
                end
            elseif k<last_frame - DN_escape
                for s = 1:4
                    subplot(2,2,s)
                    plot3(x(k,l),y(k,l),z(k,l),'y.')
                end
            else
                for s = 1:4
                    subplot(2,2,s)
                    plot3(x(k,l),y(k,l),z(k,l),'g.')
                end
            end
        end
    %         for s = 1:4
    %             subplot(2,2,s)
    %             plot3(x(k,:),y(k,:),z(k,:),'b.')
    %         end
    %     else
    %         for s = 1:4
    %             subplot(2,2,s)
    %             plot3(x(k,:),y(k,:),z(k,:),'r.')
    %         end
    %     end

        tic
        mov = addframe(mov,gcf);
        toc;
    end


    mov = close(mov);

end
