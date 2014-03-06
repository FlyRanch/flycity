figure(Q)
hold on
axis(AX_phi);
x = linspace(0,1,length(stroke_wb_R_steady_bins(:,1)));
for i = 1:length(stroke_wb_R_steady_bins);
y = stroke_wb_R_steady_bins(:,i);
plot(x,y,'color',[.5 .5 .5])

end
hold off