figure(S)
hold on
axis(AX_theta);
x = linspace(0,0.5,length(dev_ds_R_steady_bins(:,1)));
x2 = linspace(0.5,1,length(dev_us_R_steady_bins(:,1)));

for i = 1:length(dev_wb_R_steady_bins);
y = dev_ds_R_steady_bins(:,i); % downstroke
y2 = dev_us_R_steady_bins(:,i); % upstroke

plot(x,y,'color',[.5 .5 .5])
plot(x2,y2,'color',[.5 .5 .5])


end

hold off