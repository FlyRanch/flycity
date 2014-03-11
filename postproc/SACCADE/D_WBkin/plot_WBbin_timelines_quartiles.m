% calc angle_post

% circ mean & 95% range
X = t_bin;
Y = var;

Y_sort = sort(Y')';

N=length(Y);
Nquart = .25*N;
Nquart_high = round(N - Nquart);
Nquart_low = round(Nquart);

Y_mean = nanmean(Y')';
Y_high = Y_sort(:,Nquart_high);
Y_low = Y_sort(:,Nquart_low);


hold on
ciplot(Y_low,Y_high,X,plotcolor)
% alpha(.25)
plot(X,Y_mean,'-','linewidth',1,'color','k')








