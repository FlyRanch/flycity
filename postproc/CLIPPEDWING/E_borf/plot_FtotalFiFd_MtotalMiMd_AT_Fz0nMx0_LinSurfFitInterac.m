subplot(2,3,1)
ezsurf(solFtot,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([-2 0])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',-2:1:0)
title(h,'Fz total')
grid off

subplot(2,3,2)
ezsurf(solFd,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([F_min F_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',F_min:(F_max-F_min)/2:F_max)
title(h,'Fz of damaged wing')
grid off

subplot(2,3,3)
ezsurf(solFi,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([F_min F_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',F_min:(F_max-F_min)/2:F_max)
title(h,'Fz of intact wing')
grid off

subplot(2,3,4)
ezsurf(solMtot,[S2_min,S2_max,S3_min,S3_max])
title([])
view(2)
shading interp
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([M_min M_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',M_min:(M_max-M_min)/2:M_max)
title(h,'Mx total')
axis equal
axis tight
axis([S2_min S2_max S3_min S3_max M_min M_max])
view(2)
grid off

subplot(2,3,5)
ezsurf(solMd,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([M_min M_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',M_min:(M_max-M_min)/2:M_max)
title(h,'Mx damaged wing')
grid off

subplot(2,3,6)
ezsurf(solMi,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([M_min M_max])
set(h,'xtick',M_min:(M_max-M_min)/2:M_max)
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Mx intact wing')
grid off
