% kmean flighttracks

clc
clear
close all

% load('flightpathDB_pos_INCq.mat')
% load('flightpathDB_pos_qbodyEKF.mat')
% load('flightpathDB_pos_qbodyEKF_INCroll_2clusters_Ahor2.75mps2_strokeplane47.5deg_startframe1.mat')
load('flightpathDB_pos_qbodyEKF_INCroll_9clusters_1.4n-1.7n1.9_strokeplane47.5deg_startframe100.mat')

mkdir('EMcluster')
cd('EMcluster')

% kcluster db
tk=pathDB.t;
for i=1:size(pathDB.V,2)-1
    tk = [tk;pathDB.t];
end
t = tk;

stim_angle_velk = pathDB.stim_angle_vel(:);
% r_hork = pathDB.r_hor(:);
Vk = pathDB.V(:);

Ak=pathDB.A(:);
Ank=pathDB.An(:);
Atk=pathDB.At(:);
A_hork=pathDB.A_hor(:);
At_hork=pathDB.At_hor(:);
An_hork=pathDB.An_hor(:);

% alpha_dot_hork = pathDB.alpha_dot_hor(:);

% all data unclustered
figure
skip = 1

plot(An_hork(1:skip:end),At_hork(1:skip:end),'.k','markersize',5)

grid on
axis equal
xlabel('An','fontsize',18)
ylabel('At','fontsize',18)
set(gca,'xlim',[-20 20],'ylim',[-20 20])
set(gca,'XTick',[-20:10:20],'fontsize',12)
set(gca,'YTick',[-20:10:20],'fontsize',12)


saveas(gca,['AnAt_unclustered_skip',num2str(skip),'.fig'])
saveas(gca,['AnAt_unclustered_skip',num2str(skip),'.png'])
% plot2svg(['AnAt_unclustered_skip',num2str(skip),'.svg'])

%% cluster variables
X_accel = [tk An_hork At_hork];
X_accel = [An_hork At_hork];
% X_accel = [tk An_hork At_hork stim_angle_velk Vk];
% X_accel = [Ak];

%% EM cluster
skip = 1
% skip = 100

options = statset('Display','final');
gm = gmdistribution.fit(X_accel,2,'Options',options);

plot(X_accel(1:skip:end,1),X_accel(1:skip:end,2),'.k','markersize',2)
hold on
ezcontour(@(x,y)pdf(gm,[x y]));
hold off
axis equal
xlabel('An','fontsize',18)
ylabel('At','fontsize',18)
set(gca,'xlim',[-6 6],'ylim',[-6 6])
set(gca,'XTick',[-21:3:21],'fontsize',12)
set(gca,'YTick',[-21:3:21],'fontsize',12)

saveas(gca,['EMmean_clusters_contour_skip',num2str(skip),'.fig'])
saveas(gca,['EMmean_clusters_contour_skip',num2str(skip),'.png'])
% plot2svg(['EMmean_clusters_contour_skip',num2str(skip),'.svg'])

%% EMmean_2clusters
skip = 1
% skip = 100

idx = cluster(gm,X_accel);
cluster1 = (idx == 1);
cluster2 = (idx == 2);
X1 = X_accel(cluster1,:);
X2 = X_accel(cluster2,:);

% scatter(X_accel(cluster1,1),X_accel(cluster1,2),10,'r+');
% hold on
% scatter(X_accel(cluster2(1:skip:end),1),X_accel(cluster2(1:skip:end),2),10,'bo');
% hold off
% legend('Cluster 1','Cluster 2','Location','NW')
plot(X1(1:skip:end,1),X1(1:skip:end,2),'.r','markersize',5);
hold on
plot(X2(1:skip:end,1),X2(1:skip:end,2),'.','color',[.5 .5 .5],'markersize',5);
hold off
% legend('Cluster 1','Cluster 2','Location','NW')
axis equal
xlabel('An','fontsize',18)
ylabel('At','fontsize',18)
set(gca,'xlim',[-6 6],'ylim',[-6 6])
set(gca,'XTick',[-21:3:21],'fontsize',12)
set(gca,'YTick',[-21:3:21],'fontsize',12)

saveas(gca,['EMmean_2clusters_skip',num2str(skip),'.fig'])
saveas(gca,['EMmean_2clusters_skip',num2str(skip),'.png'])
% plot2svg(['EMmean_2clusters_skip',num2str(skip),'.svg'])

%% EMmean_2clusters_probdistr_skip
skip = 1
skip = 100

P = posterior(gm,X_accel);
P1 = P(cluster1,:);
P2 = P(cluster2,:);

scatter(X1(1:skip:end,1),X1(1:skip:end,2),1,P1(1:skip:end,1),'.')
hold on
scatter(X2(1:skip:end,1),X2(1:skip:end,2),1,P2(1:skip:end,1),'.')
hold off
legend('Cluster 1','Cluster 2','Location','NW')
clrmap = jet(80); colormap(clrmap(9:72,:))
ylabel(colorbar,'Component 1 Posterior Probability')
axis equal
xlabel('An','fontsize',18)
ylabel('At','fontsize',18)
set(gca,'xlim',[-6 6],'ylim',[-6 6])
set(gca,'XTick',[-21:3:21],'fontsize',12)
set(gca,'YTick',[-21:3:21],'fontsize',12)

% saveas(gca,['EMmean_2clusters_probdistr_skip',num2str(skip),'.fig'])
% saveas(gca,['EMmean_2clusters_probdistr_skip',num2str(skip),'.png'])
saveas(gca,['EMmean_2clusters_probdistr_skip',num2str(skip),'.jpg'])
% plot2svg(['EMmean_2clusters_probdistr_skip',num2str(skip),'.svg'])

%% EMmean_2clusters_probdistr_score
skip = 1
% skip = 100

[~,order] = sort(P(:,1));
plot(1:size(X_accel,1),P(order,1),'r-',1:size(X_accel,1),P(order,2),'b-','linewidth',2);
legend({'Cluster 1 Score' 'Cluster 2 Score'},'location','NW');
ylabel('Cluster Membership Score','fontsize',12)
xlabel('Point Ranking','fontsize',12)
% axis equal
set(gca,'xlim',[0 2e5],'ylim',[0 1])
set(gca,'XTick',[0:1e5:2e5],'fontsize',12)
set(gca,'YTick',[0:.5:1],'fontsize',12)

saveas(gca,['EMmean_2clusters_probdistr_score_skip',num2str(skip),'.fig'])
saveas(gca,['EMmean_2clusters_probdistr_score_skip',num2str(skip),'.png'])
plot2svg(['EMmean_2clusters_probdistr_score_skip',num2str(skip),'.svg'])


%% k-mean cluster
% 
% % % inc seeds
% % k=9;
% % seeds = [0.1 0 0; 0.1 10 0; 0.1 -5 0; 0.1 0 10; 0.1 0 -5; 0.1 10 10; 0.1 10 -5; 0.1 -10 10; 0.1 -10 -5];
% % IDXk = kmeans(X_accel,k,'start', seeds);
% 
% % no seeds
% k = 6
% IDXk = kmeans(X_accel,k);
% 
% % subsets
% for i = 1:size(A,2)
%     IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
% end

%% manual clustering
% At_hor_thresh =2.5;
% At_hor_thresh_min =-3;
% An_hor_thresh =2.5;
% 
% k=9;
% 
% for i=1:length(An_hork)
%     if isnan(An_hork(i)) == 1
%         IDXk(i,1) = nan;
%     elseif An_hork(i) < -An_hor_thresh && At_hork(i) < At_hor_thresh_min
%         IDXk(i,1) = 1;
%     elseif An_hork(i) < -An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 3;
%     elseif An_hork(i) < -An_hor_thresh
%         IDXk(i,1) = 2;
%     elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) < At_hor_thresh_min
%         IDXk(i,1) = 4;
%     elseif abs(An_hork(i)) < An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 6;
%     elseif abs(An_hork(i)) < An_hor_thresh
%         IDXk(i,1) = 5;
%     elseif An_hork(i) > An_hor_thresh && At_hork(i) < At_hor_thresh_min
%         IDXk(i,1) = 7;
%     elseif An_hork(i) > An_hor_thresh && At_hork(i) > At_hor_thresh
%         IDXk(i,1) = 9;
%     elseif An_hork(i) > An_hor_thresh
%         IDXk(i,1) = 8;
%     end
% end
% 
% % subsets
% for i = 1:size(A,2)
%     IDX(:,i) = IDXk((i-1)*size(A,1)+1:i*size(A,1));
% end

%% plot data
% cmap = cmap(1:end/k:end,:);

% % % 9 cluster cmap
% % cmap_k = [0 0 0;  0 0 0; 0 0 0;... % black
% %          0 0 1; .5 .5 .5; 1 0 0;... % blue ; grey; red;
% %          0 1 0;  1 1 0; 1 .5  0];   % green; yellow; orange
% %      
% cmap_k = [.5 0 1;  0  0  1; 0 1 1;... % dark purple; blue; light blue
%          1 0 1; .5 .5 .5; 0 1 0;... % light purple; grey; green;
%          1 0 0;  1 .5  0; 1 1 0];   % red; orange; yellow
% %      
% % cmap_k = [0 1 0;  1  1  0; 1 .5 0;... % green; yellow; orange
% %           0 1 1; .5 .5 .5; 1  0 0;... % cyan; grey; red;
% %           0 0 1; .5  0  1; 1  0 1];   % blue; purple; magenta
% % 
% % cmap_k = [0 1 0;  1  1  0; 1 .5 0;... % green; yellow; orange
% %               0 0 1; .5 .5 .5; 1  0 0]; % blue; grey; red;
% % 
% % cmap_k = [0 0 1; .5 .5 .5; 1  0 0;... % blue; grey; red;
% %           0 1 0;  1  1  0; 1 .5 0]; % green; yellow; orange
% 
% % 6 clusters
% cmap_k = [.5 .5 .5; 1  1  0; 0 1 0;... % grey; yellow; green
%            1 .5  0; 0 0 1; 1  0 0]; % orange; blue; red
% % % 7 clusters
% % cmap_k = [0  0  1; 1 0 1; 0 1 0;... % blue; purple; green
% %          1 .5 0; .5 .5 .5; 1 1 0;... % orange; grey; yellow; 
% %          1 0 0];   % red
%      
% % circular cmap
% cmap_360r = [zeros(45,1); [0:1/(45-1):1]'; ones(3*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(45,1)];
% cmap_360g = [ones(2*45,1); [1:-.5/(45-1):.5]'; [.5:-.5/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]'];
% cmap_360b = [[1:-1/(45-1):0]'; zeros(3*45,1);[0:1/(45-1):1]';ones(3*45,1)];
% 
% cmap_360 = [cmap_360r cmap_360g cmap_360b];
     

%% state clusters
% figure
% skip = 100
% for i = 1:skip:size(IDXk)
%     length(IDXk)-i
%     hold on
%     if isnan(IDXk(i)) == 0
% %         plot3(X_accel(i,1),X_accel(i,2),X_accel(i,3),'.','color',cmap(IDXk(i),:))
%         plot(X_accel(i,2),X_accel(i,3),'ok','markerfacecolor',cmap_k(IDXk(i),:),'markersize',5)
% %         plot(X_accel(i,1),X_accel(i,2),'ok','markerfacecolor',cmap_k(IDXk(i),:),'markersize',5)
%     end
% end
% grid on
% axis equal
% xlabel('An','fontsize',18)
% ylabel('At','fontsize',18)
% set(gca,'xlim',[-20 20],'ylim',[-20 20])
% set(gca,'XTick',[-20:10:20],'fontsize',12)
% set(gca,'YTick',[-20:10:20],'fontsize',12)
% 
% 
% saveas(gca,['kmean_',num2str(k),'clusters_skip',num2str(skip),'.fig'])
% saveas(gca,['kmean_',num2str(k),'clusters_skip',num2str(skip),'.png'])
% plot2svg(['kmean_',num2str(k),'clusters_skip',num2str(skip),'.svg'])
% 

cd ..