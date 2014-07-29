
%% heatmap hist roll vs yaw
x_hist = yaw_hist;
y_hist = roll_hist;

subplot(2,2,1)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('yaw','fontsize',10) 
ylabel('roll','fontsize',10) 
% colorbar

%% heatmap hist roll vs pitch
x_hist = pitch_hist;
y_hist = roll_hist;

subplot(2,2,2)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('pitch','fontsize',10) 
ylabel('roll','fontsize',10) 
% colorbar

%% heatmap hist pitch vs yaw
x_hist = pitch_hist;
y_hist = yaw_hist;

subplot(2,2,4)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('pitch','fontsize',10) 
ylabel('yaw','fontsize',10) 
% colorbar
