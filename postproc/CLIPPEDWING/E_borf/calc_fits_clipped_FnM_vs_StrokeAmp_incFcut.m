        cut_type_fit(p,1) = cut_type_now;
        cut_ratio_fit(p,1) = cut_ratio_now;

        l_ratio_fit(p,1) = BorfMorphCutData.WingLength_ratio(n_geom);
        CoA_ratio_fit(p,1) = BorfMorphCutData.CoA_ratio(n_geom);
        A_ratio_fit(p,1) = BorfMorphCutData.WingArea_ratio(n_geom);
        S1_ratio_fit(p,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
        S2_ratio_fit(p,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
        S3_ratio_fit(p,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);

        % F-Aratio
        [Fx_Amp_fit(p,:), Fx_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Fx_norm(n_now),1);
        [Fy_Amp_fit(p,:), Fy_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Fy_norm(n_now),1);
        [Fz_Amp_fit(p,:), Fz_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Fz_norm(n_now),1);
        
        [Fy_MinSteady_Amp_fit(p,:), Fy_MinSteady_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Fy_norm_MinSteady(n_now),1);
        
        % M-Aratio
        [Mx_Amp_fit(p,:), Mx_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Mx_norm(n_now),1);
        [My_Amp_fit(p,:), My_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),My_norm(n_now),1);
        [Mz_Amp_fit(p,:), Mz_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Mz_norm(n_now),1);
        
        [Mx_MinSteady_Amp_fit(p,:), Mx_MinSteady_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Mx_norm_MinSteady(n_now),1);
        [My_MinSteady_Amp_fit(p,:), My_MinSteady_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),My_norm_MinSteady(n_now),1);
        [Mz_MinSteady_Amp_fit(p,:), Mz_MinSteady_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),Mz_norm_MinSteady(n_now),1);
        
        [My_CoM_Amp_fit(p,:), My_CoM_Amp_fit_error(p,:)] = polyfit(Amp_ratio(n_now),My_norm_CoM(n_now),1);
        
        % My second order fit
        [My_MinSteady_Amp_fit2(p,:), My_MinSteady_Amp_fit2_error(p,:)] = polyfit(Amp_ratio(n_now),My_norm_MinSteady(n_now),2);
        [My_CoM_Amp_fit2(p,:), My_CoM_Amp_fit2_error(p,:)] = polyfit(Amp_ratio(n_now),My_norm_CoM(n_now),2);
        
        %% parabola fits
        
        %% !!! Fx & Fz: FORCE FUNCTION THROUGH -1/2, Fy: FORCE SYMMETRIC !!!
        Fx_norm0 = Fx_norm_steady/2;
        [Fx_Amp_parabfit0{p}, Fx_Amp_parabfit0_error{p}] = ParabolaFit_FixedConst_symmetric(Amp_ratio(n_now),Fx_norm(n_now),Fx_norm0);
        
        Fz_norm0 = -1/2;
        [Fz_Amp_parabfit0{p}, Fz_Amp_parabfit0_error{p}] = ParabolaFit_FixedConst_symmetric(Amp_ratio(n_now),Fz_norm(n_now),Fz_norm0);
        
        [Fx_Amp_parabfit{p}, Fx_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fx_norm(n_now));
        [Fy_Amp_parabfit{p}, Fy_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fy_norm(n_now));
        [Fz_Amp_parabfit{p}, Fz_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fz_norm(n_now));
        
        [Fx_Amp_parabfitsym{p}, Fx_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fx_norm(n_now));
        [Fy_Amp_parabfitsym{p}, Fy_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fy_norm(n_now));
        [Fz_Amp_parabfitsym{p}, Fz_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fz_norm(n_now));
        
        [Fy_MinSteady_Amp_parabfit{p}, Fy_MinSteady_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fy_norm_MinSteady(n_now));
        [Fy_MinSteady_Amp_parabfitsym{p}, Fy_MinSteady_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fy_norm_MinSteady(n_now));
        
        %% !!! Moments: FORCE SYMMETRIC !!!
        [Mx_Amp_parabfit{p}, Mx_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Mx_norm(n_now));
        [My_Amp_parabfit{p}, My_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),My_norm(n_now));
        [Mz_Amp_parabfit{p}, Mz_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Mz_norm(n_now));
        
        [Mx_MinSteady_Amp_parabfit{p}, Mx_MinSteady_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Mx_norm_MinSteady(n_now));
        [My_MinSteady_Amp_parabfit{p}, My_MinSteady_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),My_norm_MinSteady(n_now));
        [Mz_MinSteady_Amp_parabfit{p}, Mz_MinSteady_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Mz_norm_MinSteady(n_now));
        
        [My_CoM_Amp_parabfit{p}, My_CoM_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),My_norm_CoM(n_now));
        
        % symmetric
        [Mx_Amp_parabfitsym{p}, Mx_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Mx_norm(n_now));
        [My_Amp_parabfitsym{p}, My_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),My_norm(n_now));
        [Mz_Amp_parabfitsym{p}, Mz_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Mz_norm(n_now));
        
        [Mx_MinSteady_Amp_parabfitsym{p}, Mx_MinSteady_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Mx_norm_MinSteady(n_now));
        [My_MinSteady_Amp_parabfitsym{p}, My_MinSteady_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),My_norm_MinSteady(n_now));
        [Mz_MinSteady_Amp_parabfitsym{p}, Mz_MinSteady_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Mz_norm_MinSteady(n_now));
        
        [My_CoM_Amp_parabfitsym{p}, My_CoM_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),My_norm_CoM(n_now));
        
        %% Fcut vs Aratio !!! FORCE Fz FUNCTION THROUGH ZERO !!!
        Fz_cut0 = 0;
        [Fz_cut_Amp_parabfit0{p}, Fz_cut_Amp_parabfit0_error{p}] = ParabolaFit_FixedConst_symmetric(Amp_ratio(n_now),Fz_norm_cut(n_now),Fz_cut0);
        
        [Fx_cut_Amp_parabfit{p}, Fx_cut_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fx_norm_cut(n_now));
        [Fy_cut_Amp_parabfit{p}, Fy_cut_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fy_norm_cut(n_now));
        [Fz_cut_Amp_parabfit{p}, Fz_cut_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fz_norm_cut(n_now));
        
        [Fy_cut_MinSteady_Amp_parabfit{p}, Fy_cut_MinSteady_Amp_parabfit_error{p}] = ParabolaFit(Amp_ratio(n_now),Fy_norm_cut_MinSteady(n_now));
        
        % symmetric
        [Fx_cut_Amp_parabfitsym{p}, Fx_cut_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fx_norm_cut(n_now));
        [Fy_cut_Amp_parabfitsym{p}, Fy_cut_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fy_norm_cut(n_now));
        [Fz_cut_Amp_parabfitsym{p}, Fz_cut_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fz_norm_cut(n_now));
        
        [Fy_cut_MinSteady_Amp_parabfitsym{p}, Fy_cut_MinSteady_Amp_parabfitsym_error{p}] = ParabolaFit_symmetric(Amp_ratio(n_now),Fy_norm_cut_MinSteady(n_now));
        
        