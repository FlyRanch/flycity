subplot(2,3,1)
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
h = colorbar('location','northoutside'); 
title(h,'Aratio damaged wing @ steady freq')

subplot(2,3,2)
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
h = colorbar('location','northoutside'); 
title(h,'Aratio intact wing @ steady freq') 

subplot(2,3,3)
ezsurf(solAiAdRatio,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([AiAdRatio_min AiAdRatio_max])
colormap(cmap_Aratio)
h = colorbar('location','northoutside'); 
title(h,'Adamaged/Aintact @ steady freq') 
