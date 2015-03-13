% F-S
[Fx_S2_fit, Fx_S2_fit_error] = polyfit(S2_ratio,Fx_norm,1);
[Fy_S2_fit, Fy_S2_fit_error] = polyfit(S2_ratio,Fy_norm,1);
[Fz_S2_fit, Fz_S2_fit_error] = polyfit(S2_ratio,Fz_norm,1);

% M-S
[Mx_S3_fit, Mx_S3_fit_error] = polyfit(S3_ratio,Mx_norm,1);
[My_S3_fit, My_S3_fit_error] = polyfit(S3_ratio,My_norm,1);
[Mz_S3_fit, Mz_S3_fit_error] = polyfit(S3_ratio,Mz_norm,1);

[MxMinSteady_S3_fit, MxMinSteady_S3_fit_error] = polyfit(S3_ratio,Mx_norm-Mx_norm(cut_type==0),1);
[MyMinSteady_S3_fit, MyMinSteady_S3_fit_error] = polyfit(S3_ratio,My_norm-My_norm(cut_type==0),1);
[MzMinSteady_S3_fit, MzMinSteady_S3_fit_error] = polyfit(S3_ratio,Mz_norm-Mz_norm(cut_type==0),1);

[My_CoM_S3_fit, My_CoM_S3_fit_error] = polyfit(S3_ratio,My_norm_CoM,1);
[MyMinSteady_CoM_S3_fit, MyMinSteady_CoM_S3_fit_error] = polyfit(S3_ratio,My_norm_CoM-My_norm_CoM(cut_type==0),1);
