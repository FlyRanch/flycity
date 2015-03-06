
for i = 1:length(S2_ratio_fit)
    
    cut_type_now = cut_type_fit(i);
    cut_ratio_now = cut_ratio_fit(i);
    
    S1_ratio_now = S1_ratio_fit(i);
    S2_ratio_now = S2_ratio_fit(i);
    S3_ratio_now = S3_ratio_fit(i);
    
    % force fits
    Fx_fit_now = Fx_Amp_parabfit{i};
    Fy_fit_now = Fy_Amp_parabfit{i};
    Fz_fit_now = Fz_Amp_parabfit{i};
    
    Fy_MinSteady_fit_now = Fy_MinSteady_Amp_parabfit{i};
    
    % force fits symmetric
    Fx_fitsym_now = Fx_Amp_parabfitsym{i};
    Fy_fitsym_now = Fy_Amp_parabfitsym{i};
    Fz_fitsym_now = Fz_Amp_parabfitsym{i};
    
    Fy_MinSteady_fitsym_now = Fy_MinSteady_Amp_parabfitsym{i};
    
    Fx_fit0_now = Fx_Amp_parabfit0{i};
    Fz_fit0_now = Fz_Amp_parabfit0{i};
    
    % moment fits
    Mx_fit_now = Mx_Amp_parabfit{i};
    My_fit_now = My_Amp_parabfit{i};
    Mz_fit_now = Mz_Amp_parabfit{i};
    
    Mx_MinSteady_fit_now = Mx_MinSteady_Amp_parabfit{i};
    My_MinSteady_fit_now = My_MinSteady_Amp_parabfit{i};
    Mz_MinSteady_fit_now = Mz_MinSteady_Amp_parabfit{i};

    My_CoM_fit_now = My_CoM_Amp_parabfit{i};

    % moment fits symmetric
    Mx_fitsym_now = Mx_Amp_parabfitsym{i};
    My_fitsym_now = My_Amp_parabfitsym{i};
    Mz_fitsym_now = Mz_Amp_parabfitsym{i};
    
    Mx_MinSteady_fitsym_now = Mx_MinSteady_Amp_parabfitsym{i};
    My_MinSteady_fitsym_now = My_MinSteady_Amp_parabfitsym{i};
    Mz_MinSteady_fitsym_now = Mz_MinSteady_Amp_parabfitsym{i};

    My_CoM_fitsym_now = My_CoM_Amp_parabfitsym{i};

    % forces cut wing fits
    Fx_cut_fit_now = Fx_cut_Amp_parabfit{i};
    Fy_cut_fit_now = Fy_cut_Amp_parabfit{i};
    Fz_cut_fit_now = Fz_cut_Amp_parabfit{i};
    
    Fy_cut_MinSteady_fit_now = Fy_cut_MinSteady_Amp_parabfit{i};
    Fz_cut_fit0_now = Fz_cut_Amp_parabfit0{i};
    
    % datapoints
    n_now = find(cut_type==cut_type_now & cut_ratio==cut_ratio_now);
    Amp_ratio_now = Amp_ratio(n_now);
    
    % forces
    Fx_now = Fx_norm(n_now);
    Fy_now = Fy_norm(n_now);
    Fz_now = Fz_norm(n_now);
    
    Fy_MinSteady_now = Fy_norm_MinSteady(n_now);

    % torques
    Mx_now = Mx_norm(n_now);
    My_now = My_norm(n_now);
    Mz_now = Mz_norm(n_now);
    
    Mx_MinSteady_now = Mx_norm_MinSteady(n_now);
    My_MinSteady_now = My_norm_MinSteady(n_now);
    Mz_MinSteady_now = Mz_norm_MinSteady(n_now);

    My_CoM_now = My_norm_CoM(n_now);
    
    % forces cut wing
    Fx_cut_now = Fx_norm_cut(n_now);
    Fy_cut_now = Fy_norm_cut(n_now);
    Fz_cut_now = Fz_norm_cut(n_now);
    
    Fy_cut_MinSteady_now = Fy_norm_cut_MinSteady(n_now);
        
    
    % plot symbol
    if cut_type_now == 0
        symbol_now = 's';
    elseif cut_type_now == 1
        symbol_now = 'o';
    elseif cut_type_now == 2
        symbol_now = 'd';
    end
    
    %% plot trends
    % forces
    figure(1)
    color_now = cmap(round(100*S2_ratio_now),:);
    
    subplot(1,3,1)
    plot(Amp_ratio_now,Fx_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fx_fit_now,'k')
    plot(Fx_fitsym_now,'r')
    subplot(1,3,2)
    plot(Amp_ratio_now,Fy_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fy_fit_now,'k')
    plot(Fy_fitsym_now,'r')
    subplot(1,3,3)
    plot(Amp_ratio_now,Fz_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fz_fit_now,'k')
    plot(Fz_fitsym_now,'r')
    plot(Fz_fit0_now,'g')
    
    % Fy symmetric
    figure(2)
    color_now = cmap(round(100*S2_ratio_now),:);

    subplot(1,3,1)
    plot(Amp_ratio_now,Fx_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fx_fit_now,'k')
    plot(Fx_fitsym_now,'r')
    plot(Fx_fit0_now,'g')
    subplot(1,3,2)
    plot(Amp_ratio_now,Fy_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fy_MinSteady_fit_now,'k')
    plot(Fy_MinSteady_fitsym_now,'r')
    subplot(1,3,3)
    plot(Amp_ratio_now,Fz_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Fz_fit_now,'k')
    plot(Fz_fitsym_now,'r')
    plot(Fz_fit0_now,'g')
    
    % Torques
    figure(3)
    color_now = cmap(round(100*S3_ratio_now),:);
    
    subplot(1,3,1)
    plot(Amp_ratio_now,Mx_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mx_fit_now,'k')
    plot(Mx_fitsym_now,'r')
    subplot(1,3,2)
    plot(Amp_ratio_now,My_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(My_fit_now,'k')
    plot(My_fitsym_now,'r')
    subplot(1,3,3)
    plot(Amp_ratio_now,Mz_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mz_fit_now,'k')
    plot(Mz_fitsym_now,'r')
    
    % Torques MINUS steady
    figure(4)
    color_now = cmap(round(100*S3_ratio_now),:);
    
    subplot(1,3,1)
    plot(Amp_ratio_now,Mx_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mx_MinSteady_fit_now,'k')
    plot(Mx_MinSteady_fitsym_now,'r')
    subplot(1,3,2)
    plot(Amp_ratio_now,My_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(My_MinSteady_fit_now,'k')
    plot(My_MinSteady_fitsym_now,'r')
    subplot(1,3,3)
    plot(Amp_ratio_now,Mz_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mz_MinSteady_fit_now,'k')
    plot(Mz_MinSteady_fitsym_now,'r')
    
    % Torques CoM & MINUS steady
    figure(5)
    color_now = cmap(round(100*S3_ratio_now),:);
    
    subplot(1,3,1)
    plot(Amp_ratio_now,Mx_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mx_MinSteady_fit_now,'k')
    plot(Mx_MinSteady_fitsym_now,'r')
    subplot(1,3,2)
    plot(Amp_ratio_now,My_CoM_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(My_CoM_fit_now,'k')
    plot(My_CoM_fitsym_now,'r')
    subplot(1,3,3)
    plot(Amp_ratio_now,Mz_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now)
    hold on
    plot(Mz_MinSteady_fit_now,'k')
    plot(Mz_MinSteady_fitsym_now,'r')
    
end

% F
figure(1)
subplot(1,3,1)
colormap(cmap)
h = colorbar('location','north')
title(h,'S2 perc')

subplot(1,3,1)
xlabel('Amp ratio')
ylabel('Fx/mg')
axis([1 1.25 -.5 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,2)
xlabel('Amp ratio')
ylabel('Fy/mg')
axis([1 1.25 -.5 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,3)
xlabel('Amp ratio')
ylabel('Fz/mg')
axis([1 1.25 -1.5 -.5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

saveas(gcf,'ForceVsStrokeAmpVsS2_allTrends_raw.fig')
saveas(gcf,'ForceVsStrokeAmpVsS2_allTrends_raw.png')
plot2svg('ForceVsStrokeAmpVsS2_allTrends_raw.svg')

% F MINUS steady
figure(2)
subplot(1,3,1)
colormap(cmap)
h = colorbar('location','north')
title(h,'S2 perc')

subplot(1,3,1)
xlabel('Amp ratio')
ylabel('Fx/mg')
axis([1 1.25 -.5 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,2)
xlabel('Amp ratio')
ylabel('Fy/mg - steady')
axis([1 1.25 -.5 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,3)
xlabel('Amp ratio')
ylabel('Fz/mg')
axis([1 1.25 -1.5 -.5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

saveas(gcf,'ForceVsStrokeAmpVsS2_allTrends_minFySteady.fig')
saveas(gcf,'ForceVsStrokeAmpVsS2_allTrends_minFySteady.png')
plot2svg('ForceVsStrokeAmpVsS2_allTrends_minFySteady.svg')

% M
figure(3)
subplot(1,3,3)
colormap(cmap)
h = colorbar('location','north')
title(h,'S3 perc')

subplot(1,3,1)
xlabel('Amp ratio')
ylabel('Mx/mgl')
axis([1 1.25 -.5 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,2)
xlabel('Amp ratio')
ylabel('My/mgl')
axis([1 1.25 -.5 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,3)
xlabel('Amp ratio')
ylabel('Mz/mgl')
axis([1 1.25 -.5 .25])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_raw.fig')
saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_raw.png')
plot2svg('TorqueVsStrokeAmpVsS2_allTrends_raw.svg')

% M MINUS steady
figure(4)
subplot(1,3,3)
colormap(cmap)
h = colorbar('location','north')
title(h,'S3 perc')

subplot(1,3,1)
xlabel('Amp ratio')
ylabel('Mx/mgl - steady')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,2)
xlabel('Amp ratio')
ylabel('My/mgl - steady')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,3)
xlabel('Amp ratio')
ylabel('Mz/mgl - steady')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_minSteady.fig')
saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_minSteady.png')
plot2svg('TorqueVsStrokeAmpVsS2_allTrends_minSteady.svg')

% My@CoM
figure(5)
subplot(1,3,3)
colormap(cmap)
h = colorbar('location','north')
title(h,'S3 perc')

subplot(1,3,1)
xlabel('Amp ratio')
ylabel('Mx/mgl - steady')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,2)
xlabel('Amp ratio')
ylabel('My/mgl @ CoM')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

subplot(1,3,3)
xlabel('Amp ratio')
ylabel('Mz/mgl - steady')
axis([1 1.25 -.25 .5])
set(gca,'xtick',0.75:.25:1.25)
set(gca,'ytick',-1.5:.25:1.5)
legend off

saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_minSteady_MyAtCoM.fig')
saveas(gcf,'TorqueVsStrokeAmpVsS2_allTrends_minSteady_MyAtCoM.png')
plot2svg('TorqueVsStrokeAmpVsS2_allTrends_minSteady_MyAtCoM.svg')

