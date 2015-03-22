%% fit param's
% intact wing
pFi1 = Fz_Amp_fit_coeffs(1);
pFi0 = Fz_Amp_fit_coeffs(2);

pMi1 = MxMinSteady_Amp_fit_coeffs(1);
pMi0 = MxMinSteady_Amp_fit_coeffs(2);

% damaged wing
pFd00 = Fz_Amp_S2_SurfFit_coeffs(1);
pFd10 = Fz_Amp_S2_SurfFit_coeffs(2);
pFd01 = Fz_Amp_S2_SurfFit_coeffs(3);
pFd11 = Fz_Amp_S2_SurfFit_coeffs(4);

pMd00 = Mx_Amp_S3_SurfFit_coeffs(1);
pMd10 = Mx_Amp_S3_SurfFit_coeffs(2);
pMd01 = Mx_Amp_S3_SurfFit_coeffs(3);
pMd11 = Mx_Amp_S3_SurfFit_coeffs(4);

%% force & torque equations: Fi = Fd(S2=1) & Mi = -Md(S3=1)
syms Ai Ad S2 S3

% vertical force balace
eqnFA = pFd10*Ad + pFd11*Ad*S2 + pFd01*S2 + pFd00...
      + pFi1*Ai  + pFi0 == -2;
% roll torque balance
eqnMA = - pMd10*Ad - pMd11*Ad*S3 - pMd01*S3 - pMd00...
        + pMi1*Ai  + pMi0 == 0;
  
% solve Amp as function of S2 & S3 
[solAi, solAd] = solve(eqnFA,eqnMA,Ai,Ad);

solAiAdRatio = solAi/solAd;

solAi_intact = subs(solAi,S2,1);
solAi_intact = subs(solAi_intact,S3,1);
solAi_intact = eval(solAi_intact);

solAd_intact = subs(solAd,S2,1);
solAd_intact = subs(solAd_intact,S3,1);
solAd_intact = eval(solAd_intact);

Aratio_intact = mean([solAi_intact solAd_intact]);

%% Fz & Mx of intact and damaged wing @ weight support & zero roll torque
% vertical force balace
solFtot = pFd10*solAd + pFd11*solAd*S2 + pFd01*S2 + pFd00...
        + pFi1*solAi  + pFi0 + 1;
solFd   = pFd10*solAd + pFd11*solAd*S2 + pFd01*S2 + pFd00 + .5;
solFi   = pFi1*solAi  + pFi0 + .5;

% roll torque balance
solMtot = - pMd10*solAd - pMd11*solAd*S3 - pMd01*S3 - pMd00...
          + pMi1*solAi  + pMi0;
solMd   = - pMd10*solAd - pMd11*solAd*S3 - pMd01*S3 - pMd00;
solMi   =   pMi1*solAi  + pMi0;

% % plot
% figure
% subplot(2,3,1)
% ezsurf(solFtot,[0 1],[0 1])
% view(2)
% title('Fz total')
% 
% subplot(2,3,2)
% ezsurf(solFd,[0 1],[0 1])
% view(2)
% title('Fz damaged wing')
% 
% subplot(2,3,3)
% ezsurf(solFi,[0 1],[0 1])
% view(2)
% title('Fz intact wing')
% 
% subplot(2,3,4)
% ezsurf(solM,[0 1],[0 1])
% view(2)
% title('Mx total')
% 
% subplot(2,3,5)
% ezsurf(solMd,[0 1],[0 1])
% view(2)
% title('Mx damaged wing')
% 
% subplot(2,3,6)
% ezsurf(solMi,[0 1],[0 1])
% view(2)
% title('Mx intact wing')

%% test eqns
% % amplitude of uncut wing @ S2 = S3 = 1;
% solAi_S2 = subs(solAi, S2, 1);
% solAi_S2S3 = subs(solAi_S2, S3, 1);
% Ai_equilibrium = eval(solAi_S2S3)
% 
% % amplitude of CUT wing @ S2 = S3 = 1;
% solAd_S2 = subs(solAd, S2, 1);
% solAd_S2S3 = subs(solAd_S2, S3, 1);
% Ad_equilibrium = eval(solAd_S2S3)
% 
% % Fz of uncut wing @ S2 = S3 = 1;
% solFi_S2 = subs(solFi, S2, 1);
% solFi_S2S3 = subs(solFi_S2, S3, 1);
% Fi_equilibrium = eval(solFi_S2S3)
% 
% % Fz of CUT wing @ S2 = S3 = 1;
% solFd_S2 = subs(solFd, S2, 1);
% solFd_S2S3 = subs(solFd_S2, S3, 1);
% Fd_equilibrium = eval(solFd_S2S3)
% 
% % Fz of both wings @ S2 = S3 = 1;
% solFtot_S2 = subs(solFtot, S2, 1);
% solFtot_S2S3 = subs(solFtot_S2, S3, 1);
% Ftot_equilibrium = eval(solFtot_S2S3)
% 
% % Mx of uncut wing @ S2 = S3 = 1;
% solMi_S2 = subs(solMi, S2, 1);
% solMi_S2S3 = subs(solMi_S2, S3, 1);
% Mi_equilibrium = eval(solMi_S2S3)
% 
% % Mx of CUT wing @ S2 = S3 = 1;
% solMd_S2 = subs(solMd, S2, 1);
% solMd_S2S3 = subs(solMd_S2, S3, 1);
% Md_equilibrium = eval(solMd_S2S3)
% 
% % Mx of both wings @ S2 = S3 = 1;
% solMtot_S2 = subs(solMtot, S2, 1);
% solMtot_S2S3 = subs(solMtot_S2, S3, 1);
% Mtot_equilibrium = eval(solMtot_S2S3)
% 











