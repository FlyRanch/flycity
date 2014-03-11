clc
% clear
close all
% load('flightpath_10steps_all.mat')
% 
trigger_frame = find(t == min(abs(t)));

dt = t(2) - t(1);
t_noresp = 0.05;
t_escape = 0.025;

DN_noresp = round(t_noresp / dt)
DN_escape = round(t_escape / dt)

topview = true
fourviews = true
fourviews = false

%% make top view movie

if topview == true
    movie_step = 15;
    movfilename = '3Dflightpaths_topview_kalman_500fps.avi';
    mov = avifile(movfilename, 'compression', 'none');

    start = min(find(nansum(x,2)~=0));
    stop =  max(find(nansum(x,2)~=0));
    seqs = size(x,2);

    figure
    hold on
    axis equal
%     grid on
    set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'zlim',[min(z(:)) max(z(:))]);
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
 set(gca, 'XTick', []);
 set(gca, 'YTick', []);
 set(gca, 'ZTick', []);
 
 
    for k = start:movie_step:stop

        for l = 1:seqs
            last_frame = max(find(isnan(y(:,l))==0));
            if k<trigger_frame
                plot3(x(k,l),y(k,l),z(k,l),'b.')
            elseif k<trigger_frame + DN_noresp
                plot3(x(k,l),y(k,l),z(k,l),'r.')
            elseif k<last_frame - DN_escape
                plot3(x(k,l),y(k,l),z(k,l),'y.')
            else
                plot3(x(k,l),y(k,l),z(k,l),'g.')
            end
        end

        tic
        mov = addframe(mov,gcf);
        toc
    end

    mov = close(mov);
end

%% make 4 views movie

if fourviews == true
    movie_step = 1;
    movfilename = '3Dflightpaths_4views_kalman.avi';
    mov = avifile(movfilename, 'compression', 'none');

    start = min(find(nansum(x,2)~=0));
    stop =  max(find(nansum(x,2)~=0));
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
        toc
    end


    mov = close(mov);

end
