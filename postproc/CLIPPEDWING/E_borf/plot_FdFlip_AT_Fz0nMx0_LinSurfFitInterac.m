
subplot(1,2,1)
ezsurf(solFd,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Fd_min Fd_max])
colormap(cmap_surf_neg)
h = colorbar('location','northoutside');
set(h,'xtick',Fd_min:(Fd_max-Fd_min)/2:Fd_max)
title(h,'Fz of damaged wing')
grid off

