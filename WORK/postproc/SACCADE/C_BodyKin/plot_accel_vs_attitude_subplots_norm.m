
An_hist = An_hor(:);
At_hist = At_hor(:);
% An_hist = An_hor(IDX~=5);
% At_hist = At_hor(IDX~=5);
A_hist = sqrt(An_hist.^2 + At_hist.^2);


%% heatmap hist A vs yaw
x_hist = yaw(:);
y_hist = A_hist;

subplot(3,3,1)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 0])
axis([min(binx) max(binx) min(biny) 0])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
set(gca,'YTick',-20:10:0,'YTicklabel',20:-10:0,'fontsize',8)
% xlabel('yaw','fontsize',10) 
ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs yaw
x_hist = yaw(:);
y_hist = An_hist;

subplot(3,3,4)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 20])
axis([min(binx) max(binx) min(biny) max(biny)])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
% xlabel('yaw','fontsize',10) 
ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs yaw
x_hist = yaw(:);
y_hist = At_hist;

subplot(3,3,7)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 10])
axis([min(binx) max(binx) min(biny) max(biny)])
set(gca,'XTick',-180:90:180) 
set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
xlabel('yaw','fontsize',10) 
ylabel('At','fontsize',10) 
% colorbar

%% heatmap hist A vs roll
x_hist = roll(:);
y_hist = A_hist;

subplot(3,3,2)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 0])
axis([min(binx) max(binx) min(biny) 0])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
% set(gca,'YTick',-20:10:0,'YTicklabel',20:-10:0,'fontsize',8)
set(gca,'YTick',-20:10:0,'YTicklabel',[],'fontsize',8)
% xlabel('roll','fontsize',10) 
% ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs roll
x_hist = roll(:);
y_hist = An_hist;

subplot(3,3,5)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 20])
axis([min(binx) max(binx) min(biny) max(biny)])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
% set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
set(gca,'YTick',-20:10:20,'YTicklabel',[],'fontsize',8)
% xlabel('roll','fontsize',10) 
% ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs roll
x_hist = roll(:);
y_hist = At_hist;

subplot(3,3,8)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 10])
axis([min(binx) max(binx) min(biny) max(biny)])
set(gca,'XTick',-180:90:180) 
% set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
set(gca,'YTick',-20:10:20,'YTicklabel',[],'fontsize',8)
xlabel('roll','fontsize',10) 
% ylabel('At','fontsize',10) 
% colorbar

%% heatmap hist A vs pitch
x_hist = pitch(:);
y_hist = A_hist;

subplot(3,3,3)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-45 135 -20 0])
axis([-180 180 -20 0])
axis([min(binx) max(binx) min(biny) 0])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
% set(gca,'YTick',-20:10:0,'YTicklabel',20:-10:0,'fontsize',8)
set(gca,'YTick',-20:10:0,'YTicklabel',[],'fontsize',8)
% xlabel('pitch','fontsize',10) 
% ylabel('A','fontsize',10) 
% colorbar

%% heatmap hist An vs pitch
x_hist = pitch(:);
y_hist = An_hist;

subplot(3,3,6)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-45 135 -20 20])
axis([-180 180 -20 20])
axis([min(binx) max(binx) min(biny) max(biny)])
% set(gca,'XTick',-180:45:180) 
set(gca,'XTick',-180:45:180,'XTicklabel',[]) 
% set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
set(gca,'YTick',-20:10:20,'YTicklabel',[],'fontsize',8)
% xlabel('pitch','fontsize',10) 
% ylabel('An','fontsize',10) 
% colorbar

%% heatmap hist At vs pitch
x_hist = pitch(:);
y_hist = At_hist;

subplot(3,3,9)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);

% for i = 1:length(binx)
%     yx_hist_log(:,i) = yx_hist_log(:,i) / max(yx_hist_log(:,i)); % normalize per time bin
% end

for i = 1:length(biny)
    yx_hist_log(i,:) = yx_hist_log(i,:) / max(yx_hist_log(i,:)); % normalize per time bin
end

imagesc(binx,-biny,yx_hist_log)
% axis equal
axis([-180 180 -20 10])
axis([min(binx) max(binx) min(biny) max(biny)])
set(gca,'XTick',-180:90:180) 
% set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
set(gca,'YTick',-20:10:20,'YTicklabel',[],'fontsize',8)
xlabel('pitch','fontsize',10) 
% ylabel('At','fontsize',10) 
% colorbar