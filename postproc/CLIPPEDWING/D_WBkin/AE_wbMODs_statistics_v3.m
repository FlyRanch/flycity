clear
clc

load('WBdataset_all_1252WBs_S2nS3.mat')
load('WBdataset_steady_1603WBs.mat')
f_wb = mean([f_wb_L,f_wb_R]')';

RS2_WBsteady = SecondMomentRatio(steady_nr_mean_wb==1);
RS3_WBsteady = ThirdMomentRatio(steady_nr_mean_wb==1);
f_WBsteady = f_wb(steady_nr_mean_wb==1);
pitch_WBsteady = pitch_mean_wb(steady_nr_mean_wb==1);
roll_WBsteady = roll_mean_wb(steady_nr_mean_wb==1);
slip_WBsteady = slip_mean_wb(steady_nr_mean_wb==1);

% wb based
f_all = [f_wb_steady;f_WBsteady];
f_groups(1:length(f_wb_steady),1)=1;
f_groups(length(f_wb_steady)+1:length(f_all),1)=2;
[p,tbl,stats] = kruskalwallis(f_all,f_groups)

[pval med P]=circ_cmtest(deg2rad(pitch_global_steady),deg2rad(pitch_WBsteady))
figure
rose(deg2rad(pitch_global_steady),360)
hold on
rose(deg2rad(pitch_WBsteady),360)

% indiv based
seq_intact_unique = unique(seq_nr_steady);
for i = 1:length(seq_intact_unique)
    n_now = find(seq_nr_steady==seq_intact_unique(i));
    
    f_intact_mean(i,1) = mean(f_wb_steady(n_now));
    pitch_intact_mean(i,1) = mean(pitch_global_steady(n_now));
end

RS2_WBsteady_unique = unique(RS2_WBsteady);
for i = 1:length(RS2_WBsteady_unique)
    n_now = find(RS2_WBsteady==RS2_WBsteady_unique(i));
    
    RS2_mean(i,1) = mean(RS2_WBsteady(n_now));
    RS3_mean(i,1) = mean(RS3_WBsteady(n_now));
    
    f_mean(i,1) = mean(f_WBsteady(n_now));
    pitch_mean(i,1) = mean(pitch_WBsteady(n_now));
    roll_mean(i,1) = mean(roll_WBsteady(n_now));
    slip_mean(i,1) = mean(slip_WBsteady(n_now));
end

% freq
f_mean_all = [f_intact_mean;f_mean];
f_mean_groups(1:length(f_intact_mean),1)=1;
f_mean_groups(length(f_intact_mean)+1:length(f_mean_all),1)=2;
[p,tbl,stats] = kruskalwallis(f_mean_all,f_mean_groups)

% pitch
[pitch_intact_mean_mean pitch_intact_ul pitch_intact_ll]=circ_mean_deg_nonan(pitch_intact_mean)
[pitch_mean_mean pitch_mean_ul pitch_mean_ll]=circ_mean_deg_nonan(pitch_mean)

[pitch_intact_mean_med]=rad2deg(circ_median(deg2rad(pitch_intact_mean)))
[pitch_mean_med]=rad2deg(circ_median(deg2rad(pitch_mean)))

[pval med P]=circ_cmtest(deg2rad(pitch_intact_mean),deg2rad(pitch_mean))

% roll
[roll_mean_mean roll_mean_ul roll_mean_ll]=circ_mean_deg_nonan(roll_mean)
pval=circ_medtest(deg2rad(roll_mean),0)
[h mu ul ll] = circ_mtest(deg2rad(roll_mean),0)
figure
circ_plot(deg2rad(roll_mean))

% slip
[slip_mean_mean slip_mean_ul slip_mean_ll]=circ_mean_deg_nonan(slip_mean)
pval=circ_medtest(deg2rad(slip_mean),0)
[h mu ul ll] = circ_mtest(deg2rad(slip_mean),0)

%% plot intact pitch
figure

angles = deg2rad(pitch_intact_mean);
angle_mean = deg2rad(pitch_intact_mean_mean);
angle_ll = deg2rad(pitch_intact_ll);
angle_ul = deg2rad(pitch_intact_ul);

c_now = [0 0 1];

subplot(2,2,1)
h=rose(angles,180)
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'b') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot damaged pitch
angles = deg2rad(pitch_mean);
angle_mean = deg2rad(pitch_mean_mean);
angle_ll = deg2rad(pitch_mean_ll);
angle_ul = deg2rad(pitch_mean_ul);

c_now = [1 0 0];

subplot(2,2,1)
h=rose(angles,180)
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot roll
angles = deg2rad(roll_mean);
angle_mean = deg2rad(roll_mean_mean);
angle_ll = deg2rad(roll_mean_ll);
angle_ul = deg2rad(roll_mean_ul);

c_now = [1 0 0];

subplot(2,2,2)
h=rose(angles,90)
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     

%% plot slip
angles = deg2rad(slip_mean);
angle_mean = deg2rad(slip_mean_mean);
angle_ll = deg2rad(slip_mean_ll);
angle_ul = deg2rad(slip_mean_ul);

c_now = [1 0 0];

subplot(2,2,3)
h=rose(angles,90)
set(h,'color','k')
x = get(h, 'XData') ;
y = get(h, 'YData') ;
p = patch(x, y, 'r') ;
set(p,'Facealpha',.3)
hold on
axis image

%// Plot mean
a = axis;
a = a(2); %// size of axis
plot([0 cos(angle_mean)*a], [0 sin(angle_mean)*a],'color',c_now,'linewidth',2)
plot([0 cos(angle_ll)*a], [0 sin(angle_ll)*a],'color','k','linewidth',.5)
plot([0 cos(angle_ul)*a], [0 sin(angle_ul)*a],'color','k','linewidth',.5)

% %// Plot error as many shaded triangles that compose a circular wedge
% t = linspace(angle_ll,angle_ul,100); %// increase "100" if needed
% for k = 1:numel(t)-1
%     h = patch([0 cos(t(k))*a cos(t(k+1))*a 0], ...
%         [0 sin(t(k))*a sin(t(k+1))*a 0], .5*c_now, 'edgecolor', 'none');
%         %// change color [.5 0 0] to something else if desired. Note also alpha
%     set(h,'Facealpha',.3) %// make transparent
% end     










