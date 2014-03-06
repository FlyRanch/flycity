%% stroke
figure(1)
title('Stroke Angle')
box on; grid off;
fr = 1/7500;
xstep = 400+[0 750 1500 2250 3000 3750];
xlab = fr.*(xstep-400);
xlabel('length of time after stimulus (s)');
ylabel('wing position (deg)');
set(gca,'XTick',xstep,'XTickLabel',{'pt.expand',xlab(2:end)},'TickLength',[0 0],'LineWidth', 1.2);
xlim([400 4800]); ylim([-60 95]);
x = [1:200]';
count = 0;
for i = 1:length(stroke_dn_wingbeats(1,1,:))
    hold on
    z = x';

    y1 = stroke_up_wingbeats(:,3,i)';
    y2 = stroke_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, stroke_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, stroke_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = stroke_dn_wingbeats(:,1)';
    y1 = stroke_dn_wingbeats(:,3,i)';
    y2 = stroke_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, stroke_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, stroke_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

%     title('stroke - up(red) and down(blue)')

%     legend([h1 h2],'up wing','down wing','box','off','EdgeColor','white','Position','best','linewidth',0)
    
    if max(x) > 400 && max(x) < 4800
        count = count+1;
%         if isinteger(count/2) == 1
        if count < 10 
            text((max(x)-120),90,num2str(count))
        else
            text((max(x)-150),90,num2str(count))
        end
%         end
        text((max(x)-200),90,'|')
        text((max(x)),90,'|')
    end
        hold off
    x = [1:200]' + max(x);
end

print -dpng -r600 rotatingdrum_stroke_timeseries.png
    
%% dev
figure(2)
title('Stroke Deviation')
box on; grid off;
fr = 1/7500;
xstep = 400+[0 750 1500 2250 3000 3750];
xlab = fr.*(xstep-400);
xlabel('length of time after stimulus (s)');
ylabel('wing position (deg)');
set(gca,'XTick',xstep,'XTickLabel',{'pt.expand',xlab(2:end)},'TickLength',[0 0],'LineWidth', 1.2);
xlim([400 4800]); ylim([-24 24]);
x = [1:200]';
count = 0;
for i = 1:length(dev_dn_wingbeats(1,1,:))
    hold on
    z = x';

    y1 = dev_up_wingbeats(:,3,i)';
    y2 = dev_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, dev_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, dev_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = dev_dn_wingbeats(:,1)';
    y1 = dev_dn_wingbeats(:,3,i)';
    y2 = dev_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, dev_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, dev_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

%     title('dev - up(red) and down(blue)')

%     legend([h1 h2],'up wing','down wing','box','off','EdgeColor','white','Position','best','linewidth',0)

        if max(x) > 400 && max(x) < 4800
        count = count+1;
%         if isinteger(count/2) == 1
        if count < 10
            text((max(x)-120),21,num2str(count))
        else
            text((max(x)-150),21,num2str(count))
        end
%         end
        text((max(x)-200),21,'|')
        text((max(x)),21,'|')
    end
    
        hold off
    x = [1:200]' + max(x);
end

print -dpng -r600 rotatingdrum_dev_timeseries.png
    
%% pitch
figure(3)
title('Wing Pitch')
box on; grid off;
fr = 1/7500;
xstep = 400+[0 750 1500 2250 3000 3750];
xlab = fr.*(xstep-400);
xlabel('length of time after stimulus (s)');
ylabel('wing position (deg)');
set(gca,'XTick',xstep,'XTickLabel',{'pt.expand',xlab(2:end)},'TickLength',[0 0],'LineWidth', 1.2);
xlim([400 4800]); ylim([15 157]);
x = [1:200]';
count = 0;
for i = 1:length(pitch_dn_wingbeats(1,1,:))
    hold on
    z = x';

    y1 = pitch_up_wingbeats(:,3,i)';
    y2 = pitch_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, pitch_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, pitch_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = pitch_dn_wingbeats(:,1)';
    y1 = pitch_dn_wingbeats(:,3,i)';
    y2 = pitch_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, pitch_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, pitch_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

%     title('pitch - up(red) and down(blue)')
%     legend([h1 h2],'up wing','down wing','box','off','EdgeColor','white','Position','best','linewidth',0)

        if max(x) > 400 && max(x) < 4800
        count = count+1;
%         if isinteger(count/2) == 1
        if count < 10
            text((max(x)-120),152,num2str(count))
        else
            text((max(x)-150),152,num2str(count))
        end
%         end
        text((max(x)-200),152,'|')
        text((max(x)),152,'|')
    end
    
    
    hold off
    
    x = [1:200]' + max(x);
end

print -dpng -r600 rotatingdrum_pitch_timeseries.png