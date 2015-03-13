subplot(2,3,1)
ezsurf(solFtotA,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Ftot_min Ftot_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Fz total')
axis equal
axis tight
axis([S2_min S2_max S3_min S3_max Ftot_min Ftot_max])
view(2)
grid off

subplot(2,3,2)
ezsurf(solFdA,[S2_min,S2_max,S3_min,S3_max])
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
title(h,'Fz of damaged wing')
grid off

subplot(2,3,3)
ezsurf(solFiA,[S2_min,S2_max,S3_min,S3_max])
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
title(h,'Fz of intact wing')
grid off

subplot(2,3,4)
ezsurf(solMtotA,[S2_min,S2_max,S3_min,S3_max])
title([])
view(2)
shading interp
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Mtot_min Mtot_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Mx total')
axis equal
axis tight
axis([S2_min S2_max S3_min S3_max Mtot_min Mtot_max])
view(2)
grid off

subplot(2,3,5)
ezsurf(solMdA,[S2_min,S2_max,S3_min,S3_max])
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
title(h,'Mx damaged wing - Mx steady')
grid off

subplot(2,3,6)
ezsurf(solMiA,[S2_min,S2_max,S3_min,S3_max])
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
title(h,'Mx intact wing - Mx steady')
grid off
