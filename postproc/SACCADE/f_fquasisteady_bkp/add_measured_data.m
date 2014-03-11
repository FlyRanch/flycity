load('borf_db_Fenhance_NOcali_alldata.mat')

loc = find(mod_value_all==0)

subplot(3,2,1)
Fx = Fx_all(:,3,loc)./Mg_fly;
Fx = Fx(isnan(Fx)==0);
t = [0:1/(length(Fx)-1):1];
plot(t,Fx,'c')

subplot(3,2,3)
Fy = Fy_all(:,3,loc)./Mg_fly;
Fy = Fy(isnan(Fy)==0);
t = [0:1/(length(Fy)-1):1];
plot(t,Fy,'c')

subplot(3,2,5)
Fz = Fz_all(:,3,loc)./Mg_fly;
Fz = Fz(isnan(Fz)==0);
t = [0:1/(length(Fz)-1):1];
plot(t,-Fz,'c')

subplot(3,2,2)
Mx = Mx_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;
Mx = Mx(isnan(Mx)==0);
t = [0:1/(length(Mx)-1):1];
plot(t,Mx,'c')

subplot(3,2,4)
My = My_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;
My = My(isnan(My)==0);
t = [0:1/(length(My)-1):1];
plot(t,My,'c')

subplot(3,2,6)
Mz = Mz_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;
Mz = Mz(isnan(Mz)==0);
t = [0:1/(length(Mz)-1):1];
plot(t,-Mz,'c')



F = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
plot(t,F,'c')

