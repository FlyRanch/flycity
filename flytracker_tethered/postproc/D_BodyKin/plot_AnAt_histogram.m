%% plot histogram An vs At

AnAt_hist = hist3([At_hist,An_hist], {biny binx});
AnAt_hist_log = log(AnAt_hist);

% imagesc(binx,-biny,AnAt_hist)
imagesc(binx,-biny,AnAt_hist_log)
axis equal
axis([0 20 -20 20])
set(gca,'XTick',0:10:20) 
set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',12)
xlabel('An','fontsize',18) 
ylabel('At','fontsize',18) 

