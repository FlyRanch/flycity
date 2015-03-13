
close all
clc
clear

%% constants
DFz0 = -.5;          % Vertical force enhancement from freq modulation !!! read from data !!!
Mx0 = .3;           % Roll torque from a single wing !!! read from data !!!

DMx0 = Mx0 * DFz0;    % Roll torque enhancement from freq modulation, for a single wing !!! read from data !!!

Fz_min = 1.5*DFz0;
Fz_max = 0.5*DFz0;

Mx_min = DMx0;
Mx_max = -DMx0;

fit_type = 'poly22';
plot_fit = 1;

%% intact wing: Fz&M vs Amp ratio reduction
% create linear fit for vertical force change from frequency variation with reducing Amplitude
Ai_DFzi  = [1 0]';
DFzi_DFzi = [DFz0 DFz0/2]';
[DFzi_Ai_fit, DFzi_Ai_fit_error] = polyfit(Ai_DFzi,DFzi_DFzi,1);

% create linear fit for roll torque change from frequency variation with reducing Amplitude
Ai_DMxi  = [1 0]';
DMxi_DMxi = [0 DMx0]';
[DMxi_Ai_fit, DMxi_Ai_fit_error] = polyfit(Ai_DMxi,DMxi_DMxi,1);

% plot results of intact wing
figure(1)
subplot(2,2,1)
plot(Ai_DFzi,DFzi_DFzi,'o-k')
xlabel('amplitude ratio')
ylabel('DFz freq mods')
ylim([Fz_min Fz_max])
set(gca,'ytick',Fz_min:(Fz_max-Fz_min)/2:Fz_max)

subplot(2,2,2)
plot(Ai_DMxi,DMxi_DMxi,'o-k')
xlabel('amplitude ratio')
ylabel('DMx freq mods')
ylim([Mx_min Mx_max])
set(gca,'ytick',Mx_min:(Mx_max-Mx_min)/2:Mx_max)


%% damaged wing: Fz&M vs Amp ratio & S2&S3
% constants

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);
% % black at edges
% cmap_surf(1,:) = [0 0 0];
% cmap_surf(end,:) = [0 0 0];

%% create surface fit for vertical force change from frequency variation in Amplitude & S2 space

% DFzd = DFz0/2 along S2=0 line
Ad_S2 = [0:.1:1.5]';
S2_S2 = zeros(size(Ad_S2));
DFzd_S2 = DFz0/2 * ones(size(Ad_S2));

% DFzd = DFz0/2 along Ad=0 line
S2_Ad = [0:.1:1.5]';
Ad_Ad = zeros(size(S2_Ad));
DFzd_Ad = DFz0/2 * ones(size(S2_Ad));

% DFzd = DFz0 at Ad=1 & S2=1
S2_AdS2 = ones(size(S2_Ad));
Ad_AdS2 = ones(size(S2_Ad));
DFzd_AdS2 = DFz0 * ones(size(S2_Ad));

Ad_DFzd = [Ad_S2;Ad_Ad;Ad_AdS2];
S2_DFzd = [S2_S2;S2_Ad;S2_AdS2];
DFzd_DFzd = [DFzd_S2;DFzd_Ad;DFzd_AdS2];

% DFzd-Aratio-S2ratio SurfFit
figure(1)
subplot(2,2,3)
[DFzd_Amp_S2_SurfFit, DFzd_Amp_S2_SurfFit_error] = createSurfaceFit(Ad_DFzd, S2_DFzd ,DFzd_DFzd ,fit_type ,plot_fit);
view(2)
shading interp
axis equal
axis tight
xlabel('amplitude ratio')
ylabel('S2 ratio')
title([])
% set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
% set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Fz_min Fz_max])
colormap(cmap_surf)
h = colorbar
title('DFz due to frequency modulation')
set(h,'ytick',Fz_min:(Fz_max-Fz_min)/2:Fz_max)


%% create surface fit for roll torque change from frequency variation in Amplitude & S3 space

% DMxd = DMx0/2 along S3=0 line
Ad_S3  = [0:.1:1.5]';
S3_S3  = zeros(size(Ad_S3));
DMxd_S3 = DMx0 * ones(size(Ad_S3));

% DMxd = DMx0/2 along Ad=0 line
S3_Ad  = [0:.1:1.5]';
Ad_Ad  = zeros(size(S3_Ad));
DMxd_Ad = DMx0 * ones(size(S3_Ad));

% DMxd = DMx0 at Ad=1 & S3=1
S3_AdS3  = ones(size(S3_Ad));
Ad_AdS3  = ones(size(S3_Ad));
DMxd_AdS3 = zeros(size(S3_Ad));

Ad_DMxd = [Ad_S3;Ad_Ad;Ad_AdS3];
S3_DMxd = [S3_S3;S3_Ad;S3_AdS3];
DMxd_DMxd = [DMxd_S3;DMxd_Ad;DMxd_AdS3];

% DMxd-Aratio-S3ratio SurfFit & plot
figure(1)
subplot(2,2,4)
[DMxd_Amp_S3_SurfFit, DMxd_Amp_S3_SurfFit_error] = createSurfaceFit(Ad_DMxd, S3_DMxd ,DMxd_DMxd ,fit_type ,plot_fit);
view(2)
shading interp
axis equal
axis tight
xlabel('amplitude ratio')
ylabel('S3 ratio')
title([])
% set(gca,'xtick',S3_min:(S3_max-S3_min)/2:S3_max)
% set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([Mx_min Mx_max])
colormap(cmap_surf)
h = colorbar
title('DMx due to frequency modulation')
set(h,'ytick',Mx_min:(Mx_max-Mx_min)/2:Mx_max)


%% save figs
mkdir('FreqMods_FnMenhancement')
cd('FreqMods_FnMenhancement')

figure(1)
saveas(gcf,'FreqMods_FnM_vs_Amplitude_vs_S2nS3_DamagedSurfFits_IntactLinFits.fig')
saveas(gcf,'FreqMods_FnM_vs_Amplitude_vs_S2nS3_DamagedSurfFits_IntactLinFits.png')
plot2svg('FreqMods_FnM_vs_Amplitude_vs_S2nS3_DamagedSurfFits_IntactLinFits.svg')

cd ..

%% save data
save('roboflyDB_FreqMods_FnM_vs_Amplitude_vs_S2nS3_DamagedSurfFits_IntactLinFits.mat')


