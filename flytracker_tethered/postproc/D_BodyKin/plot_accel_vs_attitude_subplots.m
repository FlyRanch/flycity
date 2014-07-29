
%% heatmap hist A vs yaw
x_hist = yaw_hist;
y_hist = A_hist;

subplot(3,3,1)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('yaw','fontsize',10) 
ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs yaw
x_hist = yaw_hist;
y_hist = An_hist;

subplot(3,3,4)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('yaw','fontsize',10) 
ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs yaw
x_hist = yaw_hist;
y_hist = At_hist;

subplot(3,3,7)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('yaw','fontsize',10) 
ylabel('At','fontsize',10) 
% colorbar

%% heatmap hist A vs roll
x_hist = roll_hist;
y_hist = A_hist;

subplot(3,3,2)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('roll','fontsize',10) 
% ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs roll
x_hist = roll_hist;
y_hist = An_hist;

subplot(3,3,5)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('roll','fontsize',10) 
% ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs roll
x_hist = roll_hist;
y_hist = At_hist;

subplot(3,3,8)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('roll','fontsize',10) 
% ylabel('At','fontsize',10) 
% colorbar

%% heatmap hist A vs pitch
x_hist = pitch_hist;
y_hist = A_hist;

subplot(3,3,3)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('pitch','fontsize',10) 
% ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs pitch
x_hist = pitch_hist;
y_hist = An_hist;

subplot(3,3,6)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
% xlabel('pitch','fontsize',10) 
% ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs pitch
x_hist = pitch_hist;
y_hist = At_hist;

subplot(3,3,9)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

if normy_on==1
    for i = 1:length(biny)
        yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
    end
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([min(binx) max(binx) -max(biny) -min(biny)])
set(gca,'XTick',min(binx): (max(binx)-min(binx))/4 :max(binx)) 
set(gca,'YTick',[-max(biny): (max(biny)-min(biny))/4 :-min(biny)])
set(gca,'YTicklabel',-[-max(biny): (max(biny)-min(biny))/4 :-min(biny)],'fontsize',8)
xlabel('pitch','fontsize',10) 
% ylabel('At','fontsize',10) 
% colorbar