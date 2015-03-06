
for i = 1:length(S2_ratio_fit)
    
    cut_type_now = cut_type_fit(i);
    cut_ratio_now = cut_ratio_fit(i);
    
    S1_ratio_now = S1_ratio_fit(i);
    S2_ratio_now = S2_ratio_fit(i);
    S3_ratio_now = S3_ratio_fit(i);
    
    % force fits
    Fx_fit_now = Fx_Amp_fit(i,:);
    Fy_fit_now = Fy_Amp_fit(i,:);
    Fz_fit_now = Fz_Amp_fit(i,:);
    
%     Fy_MinSteady_fit_now = Fy_MinSteady_Amp_fit(i,:);
    
    % moment fits
    Mx_MinSteady_fit_now = Mx_MinSteady_Amp_fit(i,:);
%     My_MinSteady_fit_now = My_MinSteady_Amp_fit(i,:);
    Mz_MinSteady_fit_now = Mz_MinSteady_Amp_fit(i,:);

    My_CoM_fit_now = My_CoM_Amp_fit(i,:);

    % datapoints
    n_now = find(cut_type==cut_type_now & cut_ratio==cut_ratio_now);
    Amp_ratio_now = Amp_ratio(n_now);
    
    % forces
    Fx_now = Fx_norm(n_now);
    Fy_now = Fy_norm(n_now);
    Fz_now = Fz_norm(n_now);
    
%     Fy_MinSteady_now = Fy_norm_MinSteady(n_now);

    % torques
    Mx_MinSteady_now = Mx_norm_MinSteady(n_now);
%     My_MinSteady_now = My_norm_MinSteady(n_now);
    Mz_MinSteady_now = Mz_norm_MinSteady(n_now);

    My_CoM_now = My_norm_CoM(n_now);
    
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
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(Fx_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,Fx_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    subplot(1,3,2)
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(Fy_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,Fy_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    subplot(1,3,3)
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(Fz_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,Fz_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    
    % Torques CoM & MINUS steady
    figure(2)
    color_now = cmap(round(100*S3_ratio_now),:);
    
    subplot(1,3,1)
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(Mx_MinSteady_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,Mx_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    subplot(1,3,2)
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(My_CoM_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,My_CoM_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    subplot(1,3,3)
    plot([min(Amp_ratio_now) max(Amp_ratio_now)],polyval(Mz_MinSteady_fit_now,[min(Amp_ratio_now) max(Amp_ratio_now)]),'color',color_now,'LineWidth',2);
    hold on
    plot(Amp_ratio_now,Mz_MinSteady_now,symbol_now,'markeredgecolor','k','markerfacecolor',color_now,'markersize',10)
    
end

% F
figure(1)
subplot(1,3,1)
colormap(cmap)
h = colorbar('location','north');
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

saveas(gcf,'ForceVsStrokeAmpVsS2_LinTrend.fig')
saveas(gcf,'ForceVsStrokeAmpVsS2_LinTrend.png')
plot2svg('ForceVsStrokeAmpVsS2_LinTrend.svg')

% My@CoM &MxMz MINUS steady
figure(2)
subplot(1,3,3)
colormap(cmap)
h = colorbar('location','north');
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

saveas(gcf,'TorqueVsStrokeAmpVsS3_LinTrend_minSteady_MyAtCoM.fig')
saveas(gcf,'TorqueVsStrokeAmpVsS3_LinTrend_minSteady_MyAtCoM.png')
plot2svg('TorqueVsStrokeAmpVsS3_LinTrend_minSteady_MyAtCoM.svg')

