subplot(2,2,3)
ezsurf(solAd,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([A_min A_max])
colormap(cmap_Aratio)
% h = colorbar('location','northoutside'); 
title('Aratio of damaged wing @ clipped fly freq')

subplot(2,2,4)
ezsurf(solAi,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([A_min A_max])
colormap(cmap_Aratio)
% h = colorbar('location','northoutside'); 
title('Aratio of intact wing @ clipped fly freq') 

