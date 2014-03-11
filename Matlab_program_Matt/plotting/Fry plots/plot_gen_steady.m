%% Load Variables

% addpath('/home/matt/Dropbox/Analyze/variables/');
% 
% load('FRYcomp.mat');
% 
% stroke_fry_free_steady = FRYcomp.stroke_fry_free_steady;
% stroke_fry_tethered_steady = FRYcomp.stroke_fry_tethered_steady;
% dev_fry_free_steady = FRYcomp.dev_fry_free_steady;
% dev_fry_tethered_steady = FRYcomp.dev_fry_tethered_steady
% pitch_fry_free_steady = FRYcomp.pitch_fry_free_steady;
% pitch_fry_tethered_steady = FRYcomp.pitch_fry_tethered_steady
% 
% stroke_steady_tethered_1016wb = FRYcomp.stroke_steady_tethered_1016wb;
% dev_steady_tethered_1016wb = FRYcomp.dev_steady_tethered_1016wb;
% pitch_steady_tethered_1016wb = FRYcomp.pitch_steady_tethered_1016wb;
% 
% stroke_wb_steady_bins_meanCIstd = FRYcomp.stroke_wb_steady_bins_meanCIstd;
% dev_wb_steady_bins_meanCIstd = FRYcomp.dev_wb_steady_bins_meanCIstd;
% pitch_wb_steady_bins_meanCIstd = FRYcomp.pitch_wb_steady_bins_meanCIstd;


figure(1)
%% stroke angle
subplot(1,3,1)
grid off; box on; hold on
title('Stroke Angle')
ylabel('wing position (deg)')
xlim([1 200]); ylim([-70 90]);
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

x = [1:200];
mid = stroke_steady_tethered_1016wb(:,1)';
y1 = stroke_steady_tethered_1016wb(:,3)';
y2 = stroke_steady_tethered_1016wb(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

hold on
% plot(stroke_fry_tethered_steady(:,1),stroke_fry_tethered_steady(:,2),'--','color',1/255*[51,255,51],'LineWidth', 3);
% plot(stroke_fry_free_steady(:,1),stroke_fry_free_steady(:,2),'--','color',1/255*[51,51,255],'LineWidth', 3);
plot(stroke_fry_tethered_steady(:,1),stroke_fry_tethered_steady(:,2),'--','color',1/255*[51,255,51],'LineWidth', 3);
plot(stroke_fry_free_steady(:,1),stroke_fry_free_steady(:,2),'--','color',1/255*[51,51,255],'LineWidth', 3);



% plot(MOD.stroke_wb_down_90deg_rollmax)

% h3 = fill(X,Y,1/255*[255,0,0],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% h2 = plot(y1,'color',1/255*[255,0,0]);
% h3 = plot(y2,'color',1/255*[255,0,0]);

x = [1:200];
mid = stroke_wb_steady_bins_meanCIstd(:,1)';
y1 = stroke_wb_steady_bins_meanCIstd(:,3)';
y2 = stroke_wb_steady_bins_meanCIstd(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

% h4 = fill(X,Y,'black','EdgeColor','none');
% alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color','black','linewidth',2);
% h6 = plot(y2,'color','black','linewidth',2);

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y

%% stroke deviation
subplot(1,3,2)
grid off; box on
hold on
title('Stroke Deviation')
xlabel('stroke cycle')
% axis([0 200 -20 25])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

x = [1:200];
mid = dev_steady_tethered_1016wb(:,1)';
y1 = dev_steady_tethered_1016wb(:,3)';
y2 = dev_steady_tethered_1016wb(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

hold on
plot(dev_fry_tethered_steady(:,1),dev_fry_tethered_steady(:,2),'--','color',1/255*[51,255,51],'LineWidth', 3);
plot(dev_fry_free_steady(:,1),dev_fry_free_steady(:,2),'--','color',1/255*[51,51,255],'LineWidth', 3);

% figure(1)
% hold on
% for i = 1:length(MOD.MODval)
%     MOD.MODval(i)
%     display('fwd')
%     plot(MOD.dev_wb_fwd_90deg_yawmax(:,i),'-g')
%     pause
%     MOD.MODval(i)
%     display('rwd')
%     plot(MOD.dev_wb_rwd_90deg_yawmax(:,i),'-y')
%     pause
% end

% h3 = fill(X,Y,1/255*[255,0,0],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% h2 = plot(y1,'color',1/255*[255,0,0]);
% h3 = plot(y2,'color',1/255*[255,0,0]);

x = [1:200];
mid = dev_wb_steady_bins_meanCIstd(:,1)';
y1 = dev_wb_steady_bins_meanCIstd(:,3)';
y2 = dev_wb_steady_bins_meanCIstd(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

% h4 = fill(X,Y,'black','EdgeColor','none');
% alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color','black','linewidth',2);
% h6 = plot(y2,'color','black','linewidth',2);

% plot(dev_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',1.5)
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','black')
hold off

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y

%% wing pitch
subplot(1,3,3)
grid off; box on
hold on
title('Wing Pitch')
axis([0 200 -5 155])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

x = [1:200];
mid = pitch_steady_tethered_1016wb(:,1)';
y1 = pitch_steady_tethered_1016wb(:,3)';
y2 = pitch_steady_tethered_1016wb(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

hold on
h1 = plot(pitch_fry_tethered_steady(:,1),pitch_fry_tethered_steady(:,2),'--','color',1/255*[51,255,51],'LineWidth', 3);
h2 = plot(pitch_fry_free_steady(:,1),pitch_fry_free_steady(:,2),'--','color',1/255*[51,51,255],'LineWidth', 3);

% h3 = fill(X,Y,1/255*[255,0,0],'EdgeColor','none');
% alpha(0.7)
% plot(mid,'color',[0.5 0.5 0.5]);
% plot(y1,'color',1/255*[255,0,0]);
% plot(y2,'color',1/255*[255,0,0]);

x = [1:200];
mid = pitch_wb_steady_bins_meanCIstd(:,1)'; % florian
y1 = pitch_wb_steady_bins_meanCIstd(:,3)';
y2 = pitch_wb_steady_bins_meanCIstd(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

% h4 = fill(X,Y,'black','EdgeColor','none'); % florian
% alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% plot(y1,'color','black','linewidth',2);
% plot(y2,'color','black','linewidth',2);

hold off

legend([h2 h1],'Free flight: Fry','Tethered flight: Fry','box','off','EdgeColor','white','Position','best','linewidth',0)
% legend([h2 h4 h1 h3],'Free flight: Fry','Free flight: Florian','Tethered flight: Fry','Tethered flight: George','box','off','EdgeColor','white','Position','best','linewidth',0)

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y h7

%% Save Graph
% print -dpng -r600 steady_comparison.png
% print -dpng -r600 steady_comparison_legend.png
% print -dpng -r600 fry_steady.png
% print -dpng -r600 fry_steady_legend.png

%% Save variables

% FRYcomp = {};
% 
% FRYcomp.stroke_fry_free_steady = stroke_fry_free_steady;
% FRYcomp.stroke_fry_tethered_steady = stroke_fry_tethered_steady;
% FRYcomp.dev_fry_free_steady = dev_fry_free_steady;
% FRYcomp.dev_fry_tethered_steady = dev_fry_tethered_steady
% FRYcomp.pitch_fry_free_steady = pitch_fry_free_steady;
% FRYcomp.pitch_fry_tethered_steady = pitch_fry_tethered_steady
% 
% 
% FRYcomp.stroke_steady_tethered_1016wb = stroke_steady_tethered_1016wb;
% FRYcomp.dev_steady_tethered_1016wb = dev_steady_tethered_1016wb;
% FRYcomp.pitch_steady_tethered_1016wb = pitch_steady_tethered_1016wb;
% 
% FRYcomp.stroke_wb_steady_bins_meanCIstd = stroke_wb_steady_bins_meanCIstd;
% FRYcomp.dev_wb_steady_bins_meanCIstd = dev_wb_steady_bins_meanCIstd;
% FRYcomp.pitch_wb_steady_bins_meanCIstd = pitch_wb_steady_bins_meanCIstd;
% 
% save('FRYcomp.mat','FRYcomp');

