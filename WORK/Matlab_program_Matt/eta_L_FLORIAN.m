% figure(6)
% hold on
% x = linspace(0,1,length(pitch_wb_L_steady_bins(:,1)));
% for i = 1:length(pitch_wb_L_steady_bins);
% y = pitch_wb_L_steady_bins(:,i);
% plot(x,y,'color',[.5 .5 .5])
% end
% hold off

figure(V)
hold on
axis(AX_eta);

x = linspace(0,0.5,length(pitch_ds_L_steady_bins(:,1)));
x2 = linspace(0.5,1,length(pitch_us_L_steady_bins(:,1)));

for i = 1:length(pitch_wb_L_steady_bins);
y = pitch_ds_L_steady_bins(:,i); % downstroke
y2 = pitch_us_L_steady_bins(:,i); % upstroke

plot(x,y,'color',[.5 .5 .5])
plot(x2,y2,'color',[.5 .5 .5])

end

hold off