subplot(2,3,1)
ezsurf(solFtot,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Ftot_min Ftot_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',Ftot_min:(Ftot_max-Ftot_min)/2:Ftot_max)
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
caxis([Fd_min Fd_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',Fd_min:(Fd_max-Fd_min)/2:Fd_max)
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
caxis([Fi_min Fi_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',Fi_min:(Fi_max-Fi_min)/2:Fi_max)
title(h,'Fz of intact wing')
grid off

subplot(2,3,4)
ezsurf(solMtot,[S2_min,S2_max,S3_min,S3_max])
title([])
view(2)
shading interp
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Mtot_min Mtot_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',Mtot_min:(Mtot_max-Mtot_min)/2:Mtot_max)
title(h,'Mx total')
axis equal
axis tight
axis([S2_min S2_max S3_min S3_max Mtot_min Mtot_max])
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
caxis([Md_min Md_max])
colormap(cmap_surf)
h = colorbar('location','northoutside');
set(h,'xtick',Md_min:(Md_max-Md_min)/2:Md_max)
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
caxis([Mi_min Mi_max])
set(h,'xtick',Mi_min:(Mi_max-Mi_min)/2:Mi_max)
colormap(cmap_surf)
h = colorbar('location','northoutside');
title(h,'Mx intact wing')
grid off
