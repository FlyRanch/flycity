

% subplot(1,2,2)
% title('stroke amplitude ratios @ clipped fly wingbeat freq')
% axis off

figure(4)
subplot(1,2,1)
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
title(h,'Aratio damaged wing')
set(h,'xtick',A_min:(A_max-A_min)/2:A_max)

subplot(1,2,2)
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
title(h,'Aratio intact wing') 
set(h,'xtick',A_min:(A_max-A_min)/2:A_max)

figure(7)
subplot(1,2,2)
ezsurf(solAdAiRatio,[S2_min,S2_max,S3_min,S3_max])
view(2)
shading interp
axis equal
axis tight
title([])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([AdAiRatio_min AdAiRatio_max])
colormap(cmap_surf)
h = colorbar('location','northoutside'); 
title(h,'Adamaged/Aintact') 
set(h,'xtick',AdAiRatio_min:(AdAiRatio_max-AdAiRatio_min)/2:AdAiRatio_max)

