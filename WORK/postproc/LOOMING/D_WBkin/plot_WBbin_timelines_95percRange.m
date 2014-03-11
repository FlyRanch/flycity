% calc angle_post

% circ mean & 95% range
X = t_bin;
Y = var;

Y_sort = sort(Y')';

N=length(Y);
N95 = N-.95*N;
N95_high = round(N - N95/2);
N95_low = round(N95/2);

Y_mean = nanmean(Y')';
Y_high = Y_sort(:,N95_high);
Y_low = Y_sort(:,N95_low);


hold on
ciplot(Y_low,Y_high,X,plotcolor)
% alpha(.25)
plot(X,Y_mean,'-','linewidth',1,'color','k')








