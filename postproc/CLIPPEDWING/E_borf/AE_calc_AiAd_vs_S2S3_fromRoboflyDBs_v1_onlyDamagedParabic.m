
clear
clc
close all

%% const

% colormap

cmap_surf = jet(100);
% cmap_surf = flipud(cmap_surf);
cmap_surf(1,:) = [0 0 0];
cmap_surf(end,:) = [0 0 0];

plot_on = 0;
plot_on = 1;

% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

A_min = 0;
A_max = 2;
S2_min = .5;
S2_max = 1;
S3_min = .4;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -1;
M_max = 1;
sol_nr = 1;

freq_ratio = 220/188;
%% load fit parameters of Non-cut wing (AmpReduce)
load('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF.mat')

% !!! linear fit !!!
% intact wing vertical force vs Amp (linear fit)
pFi1 = Fz_Amp_fit(1);
pFi0 = Fz_Amp_fit(2);

% intact wing roll torque vs Amp (linear fit)
pMi1 = Mx_Amp_fit(1);
pMi0 = Mx_Amp_fit(2);

% !!! parabolic fit !!!
% % intact wing vertical force vs Amp (parabolic fit)
% pFi2 = Fz_Amp_fit2(1);
% pFi1 = Fz_Amp_fit2(2);
% pFi0 = Fz_Amp_fit2(3);
% 
% % intact wing roll torque vs Amp (parabolic fit)
% pMi2 = Mx_Amp_fit2(1);
% pMi1 = Mx_Amp_fit2(2);
% pMi0 = Mx_Amp_fit2(3);

%% load fit parameters of Cut wing FnM vs Amp vs S2&S3
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF.mat')

% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit);

pFd00 = Fz_Amp_S2_SurfFit_coeffs(1);
pFd10 = Fz_Amp_S2_SurfFit_coeffs(2);
pFd01 = Fz_Amp_S2_SurfFit_coeffs(3);
pFd11 = Fz_Amp_S2_SurfFit_coeffs(4);
pFd02 = Fz_Amp_S2_SurfFit_coeffs(5);

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_Amp_S3_SurfFit);

pMd00 = Mx_Amp_S3_SurfFit_coeffs(1);
pMd10 = Mx_Amp_S3_SurfFit_coeffs(2);
pMd01 = Mx_Amp_S3_SurfFit_coeffs(3);
pMd11 = Mx_Amp_S3_SurfFit_coeffs(4);
pMd02 = Mx_Amp_S3_SurfFit_coeffs(5);

%% calc Amp's from vertical force and roll torque balance
syms Ai Ad S2 S3

%% linear fit wrt Ai, quadratic wrt Ad
% vertical force balace
eqnFA = pFi1*Ai + pFi0 + pFd02*Ad^2 + pFd01*Ad + pFd11*Ad*S2 + pFd10*S2 + pFd00 == -2;
% roll torque balance
eqnMA = pMi1*Ai + pMi0 + pMd02*Ad^2 + pMd01*Ad + pMd11*Ad*S3 + pMd10*S3 + pMd00 == 0;

%% test eqns
eqnFA_S2 = subs(eqnFA, S2, 1);
eqnFA_S2Ai = subs(eqnFA_S2, Ai, 1);
Ad_Fsteady = solve(eqnFA_S2Ai,Ad);
Ad_Fsteady = eval(Ad_Fsteady)

eqnMA_S3 = subs(eqnMA, S3, 1);
eqnMA_S3Ai = subs(eqnMA_S3, Ai, 1);
Ad_Msteady = solve(eqnMA_S3Ai,Ad);
Ad_Msteady = eval(Ad_Msteady)

%% parabolic fit wrt ALL Amps
% % vertical force balace
% eqnFA = pFi2*Ai^2 + pFi1*Ai + pFi0 + pFd02*Ad^2 + pFd01*Ad + pFd11*Ad*S2 + pFd10*S2 + pFd00 == -2;
% % roll torque balance
% eqnMA = -pMi2*Ai^2 -pMi1*Ai - pMi0 + pMd02*Ad^2 + pMd01*Ad + pMd11*Ad*S3 + pMd10*S3 + pMd00 == 0;

%% solve Amp as function of S2 & S3 
[solAi, solAd] = solve([eqnFA,eqnMA],[Ai,Ad]);

%% Fz & Mx of intact and damaged wing @ weight support & zero roll torque
solFiA = pFi1*solAi(sol_nr) + pFi0 + .5;
solFdA = pFd02*solAd(sol_nr)^2 + pFd01*solAd(sol_nr) + pFd11*solAd(sol_nr)*S2 + pFd10*S2 + pFd00 + .5;

solMiA = pMi1*solAi(sol_nr) + pMi0;
solMdA = pMd02*solAd(sol_nr)^2 + pMd01*solAd(sol_nr) + pMd11*solAd(sol_nr)*S3 + pMd10*S3 + pMd00;

%% plot Ai & Ad @ weight support & zero roll torque

close all
figure
subplot(2,3,1)
ezsurf(solAd(sol_nr),[S2_min,S2_max,S3_min,S3_max])
view(2)
title('damaged wing amplitude ratio')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([A_min A_max])
colormap(cmap_surf)
% h = colorbar('location','northoutside');

subplot(2,3,4)
ezsurf(solAi(sol_nr),[S2_min,S2_max,S3_min,S3_max])
view(2)
title('intact wing amplitude ratio')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([A_min A_max])
colormap(cmap_surf)

subplot(2,3,2)
ezsurf(solFdA,[S2_min,S2_max,S3_min,S3_max])
view(2)
title('Fz of damaged wing')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([F_min F_max])
colormap(cmap_surf)
% h = colorbar('location','northoutside');

subplot(2,3,5)
ezsurf(solFiA,[S2_min,S2_max,S3_min,S3_max])
view(2)
title('Fz of intact wing')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([F_min F_max])
colormap(cmap_surf)
% h = colorbar('location','northoutside');

subplot(2,3,3)
ezsurf(solMdA,[S2_min,S2_max,S3_min,S3_max])
view(2)
title('Mx damaged wing - Mx steady')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([M_min M_max])
colormap(cmap_surf)
% h = colorbar('location','northoutside');

subplot(2,3,6)
ezsurf(solMiA,[S2_min,S2_max,S3_min,S3_max])
view(2)
title('Mx intact wing - Mx steady')
shading interp
axis equal
axis tight
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
caxis([M_min M_max])
colormap(cmap_surf)
% h = colorbar('location','northoutside');

mkdir('figures_AnFnM_vs_S2nS3_CutnIntactWing')
cd('figures_AnFnM_vs_S2nS3_CutnIntactWing')

saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing.fig')
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing.png')
plot2svg('AnFnM_vs_S2nS3_CutnIntactWing.svg')

cd ..


