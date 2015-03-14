%% fit param's

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
eqnFA = pFd10*Ad + pFd11*Ad*S2 + pFd01*S2 + pFd00 ...
      + pFd10*Ai + pFd11*Ai*1  + pFd01*1  + pFd00 == -2;
% roll torque balance
eqnMA = pMd10*Ad + pMd11*Ad*S3 + pMd01*S3 + pMd00 ...
      + pMd10*Ai + pMd11*Ai*1  + pMd01*1  + pMd00 == 0;

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

%% Fz & Mx of intact and damaged wing @ weight support & zero roll torque
solFiA = pFd10*Ai + pFd11*Ai*1  + pFd01*1  + pFd00 + .5;
solFdA = pFd10*Ad + pFd11*Ad*S2 + pFd01*S2 + pFd00 + .5;
solFtotA = solFiA + solFdA;

solMiA = - pMd01*Ai - pMd11*Ai*1  - pMd10*1  - pMd00;
solMdA =   pMd01*Ad + pMd11*Ad*S3 + pMd10*S3 + pMd00;
solMtotA = solMiA + solMdA;
