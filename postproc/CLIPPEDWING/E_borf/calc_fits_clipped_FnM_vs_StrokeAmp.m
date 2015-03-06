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
