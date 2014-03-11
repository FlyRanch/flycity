%%%%%%%%%%%%%%%%%%%
rollMOD = find(MOD.MODval == 0.5);
yawMOD = find(MOD.MODval == -0.5);
pitchMOD = find(MOD.MODval == 1);
forceMOD = find(MOD.MODval == 1);

stroke_MOD_rollmax_up_florian = MOD.stroke_wb_up_90deg_rollmax(:,rollMOD);
stroke_MOD_rollmax_dn_florian = MOD.stroke_wb_down_90deg_rollmax(:,rollMOD);
stroke_MOD_yawmax_fwd_florian = MOD.stroke_wb_fwd_90deg_yawmax(:,yawMOD);
stroke_MOD_yawmax_rwd_florian = MOD.stroke_wb_rwd_90deg_yawmax(:,yawMOD);
stroke_MOD_pitchmax_florian = MOD.stroke_wb_90deg_pitchmax(:,pitchMOD);
stroke_MOD_forcemax_florian = MOD.stroke_wb_90deg_forcemax(:,forceMOD);

dev_MOD_rollmax_up_florian = MOD.dev_wb_up_90deg_rollmax(:,rollMOD);
dev_MOD_rollmax_dn_florian = MOD.dev_wb_down_90deg_rollmax(:,rollMOD);
dev_MOD_yawmax_fwd_florian = MOD.dev_wb_fwd_90deg_yawmax(:,yawMOD);
dev_MOD_yawmax_rwd_florian = MOD.dev_wb_rwd_90deg_yawmax(:,yawMOD);
dev_MOD_pitchmax_florian = MOD.dev_wb_90deg_pitchmax(:,pitchMOD);
dev_MOD_forcemax_florian = MOD.dev_wb_90deg_forcemax(:,forceMOD);

pitch_MOD_rollmax_up_florian = MOD.pitch_wb_up_90deg_rollmax(:,rollMOD);
pitch_MOD_rollmax_dn_florian = MOD.pitch_wb_down_90deg_rollmax(:,rollMOD);
pitch_MOD_yawmax_fwd_florian = MOD.pitch_wb_fwd_90deg_yawmax(:,yawMOD);
pitch_MOD_yawmax_rwd_florian = MOD.pitch_wb_rwd_90deg_yawmax(:,yawMOD);
pitch_MOD_pitchmax_florian = MOD.pitch_wb_90deg_pitchmax(:,pitchMOD);
pitch_MOD_forcemax_florian = MOD.pitch_wb_90deg_forcemax(:,forceMOD);

%%%%%%%%%%%%%
for i = 1:length(MOD.MODval)
stroke_wb_down_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
stroke_wb_up_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
stroke_wb_fwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
stroke_wb_rwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
stroke_wb_90deg_pitchmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1);
stroke_wb_90deg_forcemax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * strokeMOD_wb_Fenhance_bins_meanCIstd(:,1);

dev_wb_down_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
dev_wb_up_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
dev_wb_fwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
dev_wb_rwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
dev_wb_90deg_pitchmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_PitchAccel_bins_meanCIstd(:,1);
dev_wb_90deg_forcemax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * devMOD_wb_Fenhance_bins_meanCIstd(:,1);

pitch_wb_down_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
pitch_wb_up_90deg_rollmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
pitch_wb_fwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
pitch_wb_rwd_90deg_yawmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
pitch_wb_90deg_pitchmax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1);
pitch_wb_90deg_forcemax(:,i) = stroke_wb_MATT_bins_mean_steady(:,1) + MOD.MODval(i) * pitchMOD_wb_Fenhance_bins_meanCIstd(:,1);
end


% for i = 1:length(MOD.stroke_wb_up_90deg_rollmax(1,:))
%     
%     diff_roll_dn(i) = mean(abs(MOD.stroke_wb_down_90deg_rollmax(:,i)-stroke_wb_MATT_bins_mean_dn_turn(:,1)));
%     diff_roll_up(i) = mean(abs(MOD.stroke_wb_up_90deg_rollmax(:,i)-stroke_wb_MATT_bins_mean_up_turn(:,1)));
%     diff_yaw_fwd(i) = mean(abs(MOD.stroke_wb_fwd_90deg_yawmax(:,i)-stroke_wb_MATT_bins_mean_dn_turn(:,1)));
%     diff_yaw_rwd(i) = mean(abs(MOD.stroke_wb_rwd_90deg_yawmax(:,i)-stroke_wb_MATT_bins_mean_up_turn(:,1)));
%     diff_pitch(i) = mean(abs(MOD.stroke_wb_90deg_pitchmax(:,i)-stroke_wb_MATT_bins_mean_dn_turn(:,1)));
%     diff_force(i) = mean(abs(MOD.stroke_wb_90deg_forcemax(:,i)-stroke_wb_MATT_bins_mean_up_turn(:,1)));
% 
% end
% 
% roll_dn = (find(diff_roll_dn == min(diff_roll_dn)));
% roll_up = (find(diff_roll_up == min(diff_roll_up)));
% yaw_fwd = (find(diff_yaw_fwd == min(diff_yaw_fwd)));
% yaw_rwd = (find(diff_yaw_rwd == min(diff_yaw_rwd)));
% pitch = (find(diff_pitch == min(diff_pitch)));
% force =(find(diff_force == min(diff_force)));
% 
% MOD_roll_dn = MOD.MODval(roll_dn)
% MOD_roll_dn = MOD.MODval(roll_up)
% 
% MOD_yaw_fwd = MOD.MODval(yaw_fwd)
% MOD_yaw_rwd = MOD.MODval(yaw_rwd)
% 
% MOD_pitch = MOD.MODval(pitch)
% MOD_force = MOD.MODval(force)

% hold on
% plot(MOD.stroke_wb_up_90deg_rollmax(:,roll_dn),'--')
% plot(MOD.stroke_wb_down_90deg_rollmax(:,roll_up),'--')
% plot(MOD.stroke_wb_fwd_90deg_yawmax(:,yaw_fwd),'--')
% plot(MOD.stroke_wb_rwd_90deg_yawmax(:,yaw_rwd),'--')
% plot(MOD.stroke_wb_90deg_pitchmax(:,pitch),'--')
% plot(MOD.stroke_wb_90deg_forcemax(:,force),'--')

figure(1)
%% stroke angle
% subplot(1,3,1)
grid off; box on; hold on
title('Stroke Angle')
ylabel('wing position (deg)')
xlim([1 200]); ylim([-70 90]);
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(stroke_wb_MATT_bins_up_turn,'color','blue')
% plot(stroke_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','black','linewidth',2)
plot(stroke_MOD_rollmax_up_florian,'--','color','red','linewidth',2)
plot(stroke_MOD_rollmax_dn_florian,'-','color','red','linewidth',2)
% plot(stroke_MOD_yawmax_fwd_florian,'--','color','green','linewidth',2)
% plot(stroke_MOD_yawmax_rwd_florian,'-','color','green','linewidth',2)
% plot(stroke_MOD_forcemax_florian,'-','color','blue','linewidth',2)
print -dpng -r600 MOD_rollmax.png

text(20,min(stroke_wb_steady_bins_meanCIstd(:,1)), ...
['free steady: ',num2str(min(stroke_wb_steady_bins_meanCIstd(:,1)))])

text(20,min(stroke_MOD_rollmax_up_florian), ...
['free rollmax: ',num2str(min(stroke_MOD_rollmax_up_florian))])

text(20,min(stroke_MOD_forcemax_florian), ...
['free forcemax: ',num2str(min(stroke_MOD_forcemax_florian))])


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

% plot(stroke_steady(:,1),'--','color','black','linewidth',1.5)
plot(stroke_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',2)
% plot(stroke_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5)
% plot(stroke_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5)
% plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')

text(20,min(stroke_wb_MATT_bins_mean_steady(:,1))+1, ...
['tethered steady: ',num2str(min(stroke_wb_MATT_bins_mean_steady(:,1)))])

text(20,min(stroke_wb_MATT_bins_mean_dn_turn(:,1)), ...
['tethered turning: ',num2str(min(stroke_wb_MATT_bins_mean_dn_turn(:,1)))])

hold off

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y
print -dpng -r600 ptexpand_stroke_left&right_MODvalues.png

figure(2)
%% stroke deviation
% subplot(1,3,2)
grid off; box on
hold on
title('Stroke Deviation')
xlabel('stroke cycle')
axis([0 200 -20 20])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(dev_wb_MATT_bins_up_turn,'color','blue')
% plot(dev_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

plot(dev_MOD_rollmax_up_florian,'--','color','red','linewidth',2)
plot(dev_MOD_rollmax_dn_florian,'-','color','red','linewidth',2)
plot(dev_MOD_yawmax_fwd_florian,'-','color','green','linewidth',2)
plot(dev_MOD_yawmax_rwd_florian,'-','color','green','linewidth',2)
plot(dev_MOD_forcemax_florian,'-','color','blue','linewidth',2)
plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','black','linewidth',2)

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

% plot(dev_steady(:,1),'--','color','black','linewidth',1.5)
plot(dev_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',2)
% plot(dev_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5)
% plot(dev_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5)
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')
hold off

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y

% print -dpng -r600 ptexpand_dev_left&right_MODvalues.png

%% wing pitch
figure(3)
% subplot(1,3,3)
grid off; box on
hold on
title('Wing Pitch')
axis([0 200 5 155])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

% plot(pitch_wb_MATT_bins_up_turn,'color','blue')
% plot(pitch_wb_MATT_bins_dn_turn,'color',[0.5 0.5 0.5])

h1 = plot(pitch_MOD_rollmax_up_florian,'--','color','red','linewidth',2);
plot(pitch_MOD_rollmax_dn_florian,'-','color','red','linewidth',2)
h2 = plot(pitch_MOD_yawmax_fwd_florian,'-','color','green','linewidth',2);
plot(pitch_MOD_yawmax_rwd_florian,'-','color','green','linewidth',2)
h3 = plot(pitch_MOD_forcemax_florian,'-','color','blue','linewidth',2);
h4 = plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-','color','black','linewidth',2);

x = [1:200];
mid = pitch_wb_MATT_bins_mean_up_turn(:,1)';
y1 = pitch_wb_MATT_bins_mean_up_turn(:,3)';
y2 = pitch_wb_MATT_bins_mean_up_turn(:,2)';
X = [x,fliplr(x)];
Y = [y1,fliplr(y2)];

h5 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
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

h6 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
alpha(0.8)
% plot(mid,'color',[0.5 0.5 0.5]);
% h5 = plot(y1,'color',1/255*[107,145,250]);
% h6 = plot(y2,'color',1/255*[107,145,250]);

% h7 = plot(pitch_steady(:,1),'--','color','black','linewidth',1.5);
h7 = plot(pitch_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',2);
% plot(pitch_wb_MATT_bins_mean_steady(:,2),'--','color','black','linewidth',1.5);
% plot(pitch_wb_MATT_bins_mean_steady(:,3),'--','color','black','linewidth',1.5);
% plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-','color','red');
hold off

% legend([h1 h2 h3 h4(1) h5 h6 h7],'free: max roll','free: max yaw','free: max pitch','free: steady','tethered: up wing','tethered: down wing','tethered: steady','box','off','EdgeColor','white','Position','best','linewidth',0)
% legend([h1 h2 h3 h4(1) h5],'free: max roll','free: max yaw','free: max pitch','box','off','EdgeColor','white','Position','best','linewidth',0)

clear h1 h2 h3 h4 h5 h6 x y1 y2 X Y h7

% print -dpng -r600 ptexpand_pitch_left&right_MODvalues.png


hold on
plot(dev_wb_MATT_bins_mean_steady(:,1:3),'-','color','black')
plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')

plot(stroke_wb_MATT_bins_mean_steady(:,1:3),'-','color','black')
plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-','color','red')

for j = [1:40];
for i = 1:200
A(i,j) = dev_wb_MATT_bins_mean_steady(i,1)'-(j*cos(i/200*(2*pi)));
end
end

for i = 1:length(A(1,:))
plot(A(:,i))
end


% minimum = find(dev_wb_MATT_bins_mean_steady == min(dev_wb_MATT_bins_mean_steady(:,1)));
minimum = find(stroke_wb_MATT_bins_mean_steady == min(stroke_wb_MATT_bins_mean_steady(:,1)));

minimum2 = find(stroke_wb_steady_bins_meanCIstd(:,1) == min(stroke_wb_steady_bins_meanCIstd(:,1)));

A = dev_wb_MATT_bins_mean_steady(1:minimum,1)';
B = dev_wb_MATT_bins_mean_steady(minimum+1:end,1)';

C = linspace(1,100,length(A));
D = linspace(101,200,length(B));

E(2,:) = [A,B];
E(1,:) = [C,D];

J(:,1) = dev_wb_steady_bins_meanCIstd(1:108,1);
J(:,2) = linspace(1,100,108);
H(:,1) = dev_wb_steady_bins_meanCIstd(109:200,1);
H(:,2) = linspace(101,200,(200-108));

I = [J(:,1)',H(:,1)'];
K = [J(:,2)',H(:,2)'];


p2 = polyfit(K,I,25);
vals2 = polyval(p2, [1:200]);

figure(1)
box on
title('Stroke Deviation')
xlabel('stroke cycle')
axis([0 200 -20 20])
set(gca,'XTick',[0 50 100 150 200],'XTickLabel',{' ';'down';' ';'up';' '},'TickLength',[0 0],'LineWidth', 1.2);

hold on
% h1 = plot(dev_wb_steady_bins_meanCIstd(:,1),'--','color','black','linewidth',2);
% h2 = plot(E(1,:),E(2,:),'--','color','red','linewidth',2);
% h2 = plot(dev_wb_MATT_bins_mean_steady(:,1),'--','color','black','linewidth',2);
p = polyfit(E(1,:),E(2,:),25);
vals = polyval(p, [1:200]);
% plot([1:200],polyval(p, [1:200]))
h3 = plot(-(vals2'-vals'),'--','color','blue','linewidth',2);


x = [1:200];
h4 = plot(x,16*cos(x/200*2*pi),'--','color','green','linewidth',3);

% legend([h3 h4],'delta deviation', '16*cos(x)','box','off','EdgeColor','white','Position','best','linewidth',0)

print -dpng -r600 delta_corrected_dev.png





% MODplot = {};
% 
% MODplot.stroke_MOD_rollmax_up_florian = stroke_MOD_rollmax_up_florian;
% MODplot.stroke_MOD_rollmax_dn_florian = stroke_MOD_rollmax_dn_florian;
% MODplot.stroke_MOD_yawmax_fwd_florian = stroke_MOD_yawmax_fwd_florian;
% MODplot.stroke_MOD_yawmax_rwd_florian = stroke_MOD_yawmax_rwd_florian;
% MODplot.stroke_MOD_forcemax_florian = stroke_MOD_forcemax_florian;
% MODplot.stroke_wb_steady_bins_meanCIstd = stroke_wb_steady_bins_meanCIstd;
% 
% MODplot.dev_MOD_rollmax_up_florian = dev_MOD_rollmax_up_florian;
% MODplot.dev_MOD_rollmax_dn_florian = dev_MOD_rollmax_dn_florian;
% MODplot.dev_MOD_yawmax_fwd_florian = dev_MOD_yawmax_fwd_florian;
% MODplot.dev_MOD_yawmax_rwd_florian = dev_MOD_yawmax_rwd_florian;
% MODplot.dev_MOD_forcemax_florian = dev_MOD_forcemax_florian;
% MODplot.dev_wb_steady_bins_meanCIstd = dev_wb_steady_bins_meanCIstd;
% 
% MODplot.pitch_MOD_rollmax_up_florian = pitch_MOD_rollmax_up_florian;
% MODplot.pitch_MOD_rollmax_dn_florian = pitch_MOD_rollmax_dn_florian;
% MODplot.pitch_MOD_yawmax_fwd_florian = pitch_MOD_yawmax_fwd_florian;
% MODplot.pitch_MOD_yawmax_rwd_florian = pitch_MOD_yawmax_rwd_florian;
% MODplot.pitch_MOD_forcemax_florian = pitch_MOD_forcemax_florian;
% MODplot.pitch_wb_steady_bins_meanCIstd = pitch_wb_steady_bins_meanCIstd;
% 
% MODplot.stroke_wb_MATT_bins_mean_up_turn = stroke_wb_MATT_bins_mean_up_turn;
% MODplot.dev_wb_MATT_bins_mean_up_turn = dev_wb_MATT_bins_mean_up_turn;
% MODplot.pitch_wb_MATT_bins_mean_up_turn = pitch_wb_MATT_bins_mean_up_turn;
% 
% MODplot.stroke_wb_MATT_bins_mean_dn_turn = stroke_wb_MATT_bins_mean_dn_turn;
% MODplot.dev_wb_MATT_bins_mean_dn_turn = dev_wb_MATT_bins_mean_dn_turn;
% MODplot.pitch_wb_MATT_bins_mean_dn_turn = pitch_wb_MATT_bins_mean_dn_turn;
% 
% MODplot.wb_turn = 255;
% 
% save('MODplot.mat','MODplot');






