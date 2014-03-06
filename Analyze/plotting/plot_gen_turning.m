figure(1)
%% stroke angle
subplot(1,3,1)
grid off; box on; hold on
title('Stroke Angle')
ylabel('wing position (deg)')
xlim([1 200]); ylim([-70 90]);
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(stroke_wb_MATT_bins_up_turn,'color','blue')
% plot(stroke_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

x = [1:200];
mid = stroke_wb_MATT_bins_mean_up_turn(:,1)';
y1 = stroke_wb_MATT_bins_mean_up_turn(:,3)';
y2 = stroke_wb_MATT_bins_mean_up_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% h2 = plot(y1,'color',1/255*[250,107,107]);
% h3 = plot(y2,'color',1/255*[250,107,107]);

x = [1:200];
mid = stroke_wb_MATT_bins_mean_dn_turn(:,1)';
y1 = stroke_wb_MATT_bins_mean_dn_turn(:,3)';
y2 = stroke_wb_MATT_bins_mean_dn_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h4 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color',1/255*[107,145,250]);
% h6 = plot(y2,'color',1/255*[107,145,250]);

% plot(stroke_wb_MATT_bins_mean_steady(:,1:3),'--','color','black','linewidth',1.5)
plot(stroke_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',1.5)
% plot(stroke_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5)
% plot(stroke_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5)
% plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')
hold off

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y

%% stroke deviation
subplot(1,3,2)
grid off; box on
hold on
title('Stroke Deviation')
xlabel('stroke cycle')
axis([0 200 -20 25])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(dev_wb_MATT_bins_up_turn,'color','blue')
% plot(dev_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

x = [1:200];
mid = dev_wb_MATT_bins_mean_up_turn(:,1)';
y1 = dev_wb_MATT_bins_mean_up_turn(:,3)';
y2 = dev_wb_MATT_bins_mean_up_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% h2 = plot(y1,'color',1/255*[250,107,107]);
% h3 = plot(y2,'color',1/255*[250,107,107]);

x = [1:200];
mid = dev_wb_MATT_bins_mean_dn_turn(:,1)';
y1 = dev_wb_MATT_bins_mean_dn_turn(:,3)';
y2 = dev_wb_MATT_bins_mean_dn_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h4 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color',1/255*[107,145,250]);
% h6 = plot(y2,'color',1/255*[107,145,250]);

plot(dev_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',1.5)
% plot(dev_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',1.5)
% plot(dev_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5)
% plot(dev_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5)
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')
hold off

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y

%% wing pitch
subplot(1,3,3)
grid off; box on
hold on
title('Wing Pitch')
axis([0 200 35 150])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(pitch_wb_MATT_bins_up_turn,'color','blue')
% plot(pitch_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

x = [1:200];
mid = pitch_wb_MATT_bins_mean_up_turn(:,1)';
y1 = pitch_wb_MATT_bins_mean_up_turn(:,3)';
y2 = pitch_wb_MATT_bins_mean_up_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% h2 = plot(y1,'color',1/255*[250,107,107]);
% h3 = plot(y2,'color',1/255*[250,107,107]);

x = [1:200];
mid = pitch_wb_MATT_bins_mean_dn_turn(:,1)';
y1 = pitch_wb_MATT_bins_mean_dn_turn(:,3)';
y2 = pitch_wb_MATT_bins_mean_dn_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h4 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color',1/255*[107,145,250]);
% h6 = plot(y2,'color',1/255*[107,145,250]);

% h7 = plot(pitch_wb_MATT_bins_mean_steady(:,1:3),'--','color','black','linewidth',1.5);
h7 = plot(pitch_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',1.5);
% plot(pitch_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5);
% plot(pitch_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5);
% plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-','color','red');
hold off

% legend([h1 h4 h7],'up wing','down wing','steady flight','box','off','EdgeColor','white','Position','best','linewidth',0)

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y h7

print -dpng -r600 rotating_drum_10hz_6.13.2013.png
% x = [1:200]';
% [ha hb hc] = shadeplot(x,pitch_wb_MATT_bins_mean_dn_turn(:,2),stroke_wb_MATT_bins_mean_dn_turn(:,1),'r',0.8,1)