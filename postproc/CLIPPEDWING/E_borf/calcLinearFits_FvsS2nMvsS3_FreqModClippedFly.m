% F-S
[Fx_S2_fit_freqMod, Fx_S2_fit_error_freqMod,] = polyfit(S2_ratio,Fx_norm_freqMod,1);
[Fy_S2_fit_freqMod, Fy_S2_fit_error_freqMod,] = polyfit(S2_ratio,Fy_norm_freqMod,1);
[Fz_S2_fit_freqMod, Fz_S2_fit_error_freqMod,] = polyfit(S2_ratio,Fz_norm_freqMod,1);

% M-S
[Mx_S3_fit_freqMod, Mx_S3_fit_error_freqMod,] = polyfit(S3_ratio,Mx_norm_freqMod,1);
[My_S3_fit_freqMod, My_S3_fit_error_freqMod,] = polyfit(S3_ratio,My_norm_freqMod,1);
[Mz_S3_fit_freqMod, Mz_S3_fit_error_freqMod,] = polyfit(S3_ratio,Mz_norm_freqMod,1);


[MxMinSteady_S3_fit_freqMod, MxMinSteady_S3_fit_error_freqMod,] = polyfit(S3_ratio,Mx_norm_freqMod-Mx_norm_freqMod(cut_type==0),1);
[MyMinSteady_S3_fit_freqMod, MyMinSteady_S3_fit_error_freqMod,] = polyfit(S3_ratio,My_norm_freqMod-My_norm_freqMod(cut_type==0),1);
[MzMinSteady_S3_fit_freqMod, MzMinSteady_S3_fit_error_freqMod,] = polyfit(S3_ratio,Mz_norm_freqMod-Mz_norm_freqMod(cut_type==0),1);

[My_CoM_S3_fit_freqMod, My_CoM_S3_fit_error_freqMod,] = polyfit(S3_ratio,My_norm_freqMod_CoM,1);
[MyMinSteady_CoM_S3_fit_freqMod, MyMinSteady_CoM_S3_fit_error_freqMod,] = polyfit(S3_ratio,My_norm_freqMod_CoM-My_norm_freqMod_CoM(cut_type==0),1);
