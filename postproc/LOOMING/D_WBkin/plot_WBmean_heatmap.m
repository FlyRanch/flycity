var_hist = hist3([-var(:),t_hist(:)], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

imagesc(binx,-biny,var_hist)
colormap(colmap)
axis([0 1 -y_max -y_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',[x_min:(x_max-x_min)/2:x_max]) 
set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'YTicklabel',-[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)  
% set(gca,'YTick',[-y_max:(y_max-y_min)/2:-y_min],'fontsize',8)

