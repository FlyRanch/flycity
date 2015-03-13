syms Ai Ad S2 S3

%% linear fit wrt Ai, quadratic wrt Ad
% vertical force balace
eqnFA = pFi1*Ai + pFi0 + pFd02*Ad^2 + pFd01*Ad + pFd11*Ad*S2 + pFd10*S2 + pFd00 == -2;
% roll torque balance
eqnMA = pMi1*Ai + pMi0 + pMd02*Ad^2 + pMd01*Ad + pMd11*Ad*S3 + pMd10*S3 + pMd00 == 0;

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