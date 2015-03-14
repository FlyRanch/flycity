
clear
clc
close all

%% const

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at edges
cmap_surf(1,:) = [0 0 0];
cmap_surf(end,:) = [0 0 0];

plot_on = 0;
plot_on = 1;

% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

A_min = .5;
A_max = 1.5;
S2_min = .5;
S2_max = 1;
S3_min = .4;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -1;
M_max = 1;

Ftot_min = -1.5;
Ftot_max = -0.5;

Mtot_min = -1;
Mtot_max =  1;

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


%% parabolic fit wrt ALL Amps
% % vertical force balace
% eqnFA = pFi2*Ai^2 + pFi1*Ai + pFi0 + pFd02*Ad^2 + pFd01*Ad + pFd11*Ad*S2 + pFd10*S2 + pFd00 == -2;
% % roll torque balance
% eqnMA = -pMi2*Ai^2 -pMi1*Ai - pMi0 + pMd02*Ad^2 + pMd01*Ad + pMd11*Ad*S3 + pMd10*S3 + pMd00 == 0;

%% solve Amp as function of S2 & S3 
[solAi, solAd] = solve(eqnFA,eqnMA,Ai,Ad);

%% test eqns

% amplitude of uncut wing at weight support
eqnFA_S2 = subs(eqnFA, S2, 1);
eqnFA_S2Ai = subs(eqnFA_S2, Ai, 1);
Ad_Fsteady = solve(eqnFA_S2Ai,Ad);
Ad_Fsteady = eval(Ad_Fsteady)

% amplitude of uncut wing at roll equilibrium
eqnMA_S3 = subs(eqnMA, S3, 1);
eqnMA_S3Ai = subs(eqnMA_S3, Ai, 1);
Ad_Msteady = solve(eqnMA_S3Ai,Ad);
Ad_Msteady = eval(Ad_Msteady)

% amplitude of uncut wing at weight support & roll equilibrium
solAi_S2 = subs(solAi, S2, 1);
solAi_S2S3 = subs(solAi_S2, S3, 1)
Ai_equilibrium = eval(solAi_S2S3)

% amplitude of CUT wing at weight support & roll equilibrium
solAd_S2 = subs(solAd, S2, 1);
solAd_S2S3 = subs(solAd_S2, S3, 1)
Ad_equilibrium = eval(solAd_S2S3)

%% determine equation number for Amplitude ~1 at steady uncut flight
sol_nr_i = find(abs(Ai_equilibrium-1) == min(abs(Ai_equilibrium-1)));
sol_nr_d = find(abs(Ad_equilibrium-1) == min(abs(Ad_equilibrium-1)));

%% Fz & Mx of intact and damaged wing @ weight support & zero roll torque
solFiA = pFi1*solAi(sol_nr_i) + pFi0 + .5;
solFdA = pFd02*solAd(sol_nr_d)^2 + pFd01*solAd(sol_nr_d) + pFd11*solAd(sol_nr_d)*S2 + pFd10*S2 + pFd00 + .5;
solFtotA = solFiA + solFdA;

solMiA = pMi1*solAi(sol_nr_i) + pMi0;
solMdA = pMd02*solAd(sol_nr_d)^2 + pMd01*solAd(sol_nr_d) + pMd11*solAd(sol_nr_d)*S3 + pMd10*S3 + pMd00;
solMtotA = solMiA + solMdA;

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
mkdir('figures_AnFnM_vs_S2nS3_CutnIntactWing')
cd('figures_AnFnM_vs_S2nS3_CutnIntactWing')

close all
figure
plot_AiAd_FiFd_MiMd_AT_weightSupport_n_RollTorqueEquilibrium

saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing.fig')
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing.png')
% plot2svg('AnFnM_vs_S2nS3_CutnIntactWing.svg')

figure
plot_FtotalFiFd_MtotalMiMd_AT_weightSupport_n_RollEquilibrium

saveas(gcf,'FnMtotal_vs_S2nS3.fig')
saveas(gcf,'FnMtotal_vs_S2nS3.png')
% plot2svg('FnMtotal_vs_S2nS3.svg')

cd ..

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3.mat')
