% calib sensor allignment

Fx_steady = [];
Fy_steady = [];
Fz_steady = [];
Mx_steady = [];
My_steady = [];
Mz_steady = [];

Fx_steady(end+1,1) = Fx_mean_all(find(mod_value_all == 0));
Fx_steady(end+1,1) = Fx_mean_freq(find(mod_value_freq == 0));
Fx_steady(end+1,1) = Fx_mean_stroke(find(mod_value_stroke == 0));
Fx_steady(end+1,1) = Fx_mean_pitch(find(mod_value_pitch == 0));
Fx_steady(end+1,1) = Fx_mean_dev(find(mod_value_dev == 0));

Fy_steady(end+1,1) = Fy_mean_all(find(mod_value_all == 0));
Fy_steady(end+1,1) = Fy_mean_freq(find(mod_value_freq == 0));
Fy_steady(end+1,1) = Fy_mean_stroke(find(mod_value_stroke == 0));
Fy_steady(end+1,1) = Fy_mean_pitch(find(mod_value_pitch == 0));
Fy_steady(end+1,1) = Fy_mean_dev(find(mod_value_dev == 0));

Fz_steady(end+1,1) = Fz_mean_all(find(mod_value_all == 0));
Fz_steady(end+1,1) = Fz_mean_freq(find(mod_value_freq == 0));
Fz_steady(end+1,1) = Fz_mean_stroke(find(mod_value_stroke == 0));
Fz_steady(end+1,1) = Fz_mean_pitch(find(mod_value_pitch == 0));
Fz_steady(end+1,1) = Fz_mean_dev(find(mod_value_dev == 0));

Mx_steady(end+1,1) = Mx_mean_all(find(mod_value_all == 0));
Mx_steady(end+1,1) = Mx_mean_freq(find(mod_value_freq == 0));
Mx_steady(end+1,1) = Mx_mean_stroke(find(mod_value_stroke == 0));
Mx_steady(end+1,1) = Mx_mean_pitch(find(mod_value_pitch == 0));
Mx_steady(end+1,1) = Mx_mean_dev(find(mod_value_dev == 0));

My_steady(end+1,1) = My_mean_all(find(mod_value_all == 0));
My_steady(end+1,1) = My_mean_freq(find(mod_value_freq == 0));
My_steady(end+1,1) = My_mean_stroke(find(mod_value_stroke == 0));
My_steady(end+1,1) = My_mean_pitch(find(mod_value_pitch == 0));
My_steady(end+1,1) = My_mean_dev(find(mod_value_dev == 0));

Mz_steady(end+1,1) = Mz_mean_all(find(mod_value_all == 0));
Mz_steady(end+1,1) = Mz_mean_freq(find(mod_value_freq == 0));
Mz_steady(end+1,1) = Mz_mean_stroke(find(mod_value_stroke == 0));
Mz_steady(end+1,1) = Mz_mean_pitch(find(mod_value_pitch == 0));
Mz_steady(end+1,1) = Mz_mean_dev(find(mod_value_dev == 0));

Fx_steady_mean = nanmean(Fx_steady);
Fy_steady_mean = nanmean(Fy_steady);
Fz_steady_mean = nanmean(Fz_steady);

Mx_steady_mean = nanmean(Mx_steady);
My_steady_mean = nanmean(My_steady);
Mz_steady_mean = nanmean(Mz_steady);

% rotation matrix
angz_F = atan(Fy_steady_mean/Fx_steady_mean)
Tz_F = [cos(angz_F) sin(angz_F) 0; -sin(angz_F) cos(angz_F) 0; 0 0 1];

% switch My n Mz
angz_M = atan(-Mx_steady_mean/Mz_steady_mean)
Tz_M = [cos(angz_M) sin(angz_M) 0; -sin(angz_M) cos(angz_M) 0; 0 0 1];

angz_mean = mean([angz_F,angz_M]);
Tz_mean = [cos(angz_mean) sin(angz_mean) 0; -sin(angz_mean) cos(angz_mean) 0; 0 0 1];

Fpre = [Fx_steady_mean;Fy_steady_mean;Fz_steady_mean]
Fpost_F = Tz_F * Fpre
Fpost_M = Tz_M * Fpre
Fpost_mean = Tz_mean * Fpre

Mpre = [Mx_steady_mean;Mz_steady_mean;My_steady_mean]
Mpost_F = Tz_F * Mpre
Mpost_M = Tz_M * Mpre
Mpost_mean = Tz_mean * Mpre


save('FM_steady.mat','Fx_steady','Fy_steady','Fz_steady','Mx_steady','My_steady','Mz_steady'...
    ,'Fx_steady_mean','Fy_steady_mean','Fz_steady_mean','Mx_steady_mean','My_steady_mean','Mz_steady_mean'...
    ,'angz_F','Tz_F','angz_M','Tz_M','angz_mean','Tz_mean')



