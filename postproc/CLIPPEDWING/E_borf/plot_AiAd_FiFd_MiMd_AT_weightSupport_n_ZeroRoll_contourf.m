subplot(2,3,1)
ezcontourf(solAd(sol_nr_d),[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([A_min A_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'damaged wing amplitude ratio')
axis equal
axis tight

subplot(2,3,4)
ezcontourf(solAi(sol_nr_i),[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([A_min A_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'intact wing amplitude ratio')
axis equal
axis tight

subplot(2,3,2)
ezcontourf(solFdA,[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([F_min F_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Fz of damaged wing')
axis equal
axis tight

subplot(2,3,5)
ezcontourf(solFiA,[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([F_min F_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Fz of intact wing')
axis equal
axis tight

subplot(2,3,3)
ezcontourf(solMdA,[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([M_min M_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Mx damaged wing - Mx steady')
axis equal
axis tight

subplot(2,3,6)
ezcontourf(solMiA,[S2_min,S2_max,S3_min,S3_max])
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
% caxis([M_min M_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Mx intact wing - Mx steady')
axis equal
axis tight
