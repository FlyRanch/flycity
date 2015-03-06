Amp_min = 1;
Amp_max = 1.25;

S_min = 0;
S_max = 1;

Fx_min = -.5;
Fx_max = .5;
Fy_min = -.5;
Fy_max = .5;
Fz_min = -1.5;
Fz_max = -.5;

Mx_min = -.3;
Mx_max = .3;
My_min = -.3;
My_max = .3;
Mz_min = -.3;
Mz_max = .3;

% Forces
figure(1)
hold off

subplot(1,3,1)
plot(Fx_Amp_S2_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Fx_min Fx_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Fx/mg')
xlabel('Stroke Amplitude ratio')
ylabel('S2 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

subplot(1,3,2)
plot(Fy_Amp_S2_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Fy_min Fy_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Fy/mg')
xlabel('Stroke Amplitude ratio')
ylabel('S2 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

subplot(1,3,3)
plot(Fz_Amp_S2_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Fz_min Fz_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Fz/mg')
xlabel('Stroke Amplitude ratio')
ylabel('S2 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

% Torques
figure(2)
hold off

subplot(1,3,1)
plot(Mx_MinSteady_Amp_S3_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Mx_min Mx_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Mx/mgl')
xlabel('Stroke Amplitude ratio')
ylabel('S3 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

subplot(1,3,2)
plot(My_CoM_Amp_S3_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([My_min My_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'My/mgl')
xlabel('Stroke Amplitude ratio')
ylabel('S3 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

subplot(1,3,3)
plot(Mz_MinSteady_Amp_S3_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Mz_min Mz_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Mz/mgl')
xlabel('Stroke Amplitude ratio')
ylabel('S3 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

% Fz & Mx
figure(3)

subplot(1,2,1)
plot(Fz_Amp_S2_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Fz_min Fz_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Fz/mg')
xlabel('Stroke Amplitude ratio')
ylabel('S2 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

subplot(1,2,2)
plot(Mx_MinSteady_Amp_S3_SurfFit,'XLim',[Amp_min Amp_max],'YLim',[S_min S_max])
view(2)
hold on
shading interp
axis auto
caxis([Mx_min Mx_max])
colormap(jet)
h = colorbar('location','northoutside');
title(h,'Mx/mgl')
xlabel('Stroke Amplitude ratio')
ylabel('S3 ratio')
axis([Amp_min Amp_max S_min S_max])
set(gca,'xtick',Amp_min:(Amp_max-Amp_min):Amp_max)
set(gca,'ytick',S_min:(S_max-S_min):S_max)
legend off

% plot datapoints
for i = 1:length(Fx_norm)
    
    Fx_norm_now = Fx_norm(i);
    Fy_norm_now = Fy_norm(i);
    Fz_norm_now = Fz_norm(i);
    
    Mx_norm_MinSteady_now = Mx_norm_MinSteady(i);
    My_norm_CoM_now = My_norm_CoM(i);
    Mz_norm_MinSteady_now = Mz_norm_MinSteady(i);
    
    Amp_ratio_now = Amp_ratio(i);
    S2_ratio_now = S2_ratio(i);
    S3_ratio_now = S3_ratio(i);
    cut_type_now = cut_type(i);

    if cut_type_now == 0
        symbol_now = 's';
    elseif cut_type_now == 1
        symbol_now = 'o';
    elseif cut_type_now == 2
        symbol_now = 'd';
    end

    % forces
    figure(1)
    
    subplot(1,3,1)
    color_val = (Fx_norm_now - Fx_min) / (Fx_max - Fx_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S2_ratio_now,Fx_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    subplot(1,3,2)
    color_val = (Fy_norm_now - Fy_min) / (Fy_max - Fy_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S2_ratio_now,Fy_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    subplot(1,3,3)
    color_val = (Fz_norm_now - Fz_min) / (Fz_max - Fz_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S2_ratio_now,Fz_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    % torques
    figure(2)
    
    subplot(1,3,1)
    color_val = (Mx_norm_MinSteady_now - Mx_min) / (Mx_max - Mx_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S3_ratio_now,Mx_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    subplot(1,3,2)
    color_val = (My_norm_CoM_now - My_min) / (My_max - My_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S3_ratio_now,My_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    subplot(1,3,3)
    color_val = (Mz_norm_MinSteady_now - Mz_min) / (Mz_max - Mz_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S3_ratio_now,Mz_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    % Fz & Mx
    figure(3)
    
    subplot(1,2,1)
    color_val = (Fz_norm_now - Fz_min) / (Fz_max - Fz_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S2_ratio_now,Fz_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
    
    subplot(1,2,2)
    color_val = (Mx_norm_MinSteady_now - Mx_min) / (Mx_max - Mx_min) * (length(cmap)-1) +1;
    color_now = cmap_surf(round(color_val),:);
    plot3(Amp_ratio_now,S3_ratio_now,Mx_max+.1,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10,'linewidth',2)
end



