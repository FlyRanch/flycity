figure
maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
polar(0, maxHistogramValue,'-k')
hold on

for i=1:length(v_pre)
    plot([x_origin(i) u_pre(i)],[y_origin(i) v_pre(i)],'k','linewidth',2)
end
compass_zeroup(0, 0.01,'-k');

figure
maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
polar(0, maxHistogramValue,'-k')
hold on

for i=1:length(v_pre)
    plot([x_origin(i) u_pre(i)],[y_origin(i) v_pre(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
end
compass_zeroup(0, 0.01,'-k');

figure
maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
polar(0, maxHistogramValue,'-k')
hold on

for i=1:length(v_post)
    plot([x_origin(i) u_post(i)],[y_origin(i) v_post(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
    hold on
end
compass_zeroup(0, 0.01,'-k');




