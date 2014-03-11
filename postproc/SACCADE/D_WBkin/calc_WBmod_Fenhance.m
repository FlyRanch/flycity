
%% Force
t_bins = t_F;
binx = t_bins(:,1);

var = Dstroke_F;
for i = 1:length(binx)

    var_now = var(t_bins == binx(i));
    [mu ul ll] = circ_mean_deg_nonan(var_now);
    [std std0] = circ_std_deg_nonan(var_now);

    var_meanCIstd(i,:) = [mu ul ll std];
end
Dstroke_F_meanCIstd = var_meanCIstd;

var = Dpitch_F;
for i = 1:length(binx)

    var_now = var(t_bins == binx(i));
    [mu ul ll] = circ_mean_deg_nonan(var_now);
    [std std0] = circ_std_deg_nonan(var_now);

    var_meanCIstd(i,:) = [mu ul ll std];
end
Dpitch_F_meanCIstd = var_meanCIstd;

var = Ddev_F;
for i = 1:length(binx)

    var_now = var(t_bins == binx(i));
    [mu ul ll] = circ_mean_deg_nonan(var_now);
    [std std0] = circ_std_deg_nonan(var_now);

    var_meanCIstd(i,:) = [mu ul ll std];
end
Ddev_F_meanCIstd = var_meanCIstd;

% legendre polys
n_pol_WBmod = 10; % Order of used polynomials 
    t_loc = binx;
    Dstroke_loc = Dstroke_F_meanCIstd(:,1);
    Dpitch_loc = Dpitch_F_meanCIstd(:,1);
    Ddev_loc = Ddev_F_meanCIstd(:,1);

    [Dstroke_Fenhance_fit_binmean, Dstroke_Fenhance_fit_binmean_periodic] = WBfitting_singlevar(n_pol_WBmod,t_loc,Dstroke_loc);
    [Dpitch_Fenhance_fit_binmean, Dpitch_Fenhance_fit_binmean_periodic] = WBfitting_singlevar(n_pol_WBmod,t_loc,Dpitch_loc);
    [Ddev_Fenhance_fit_binmean, Ddev_Fenhance_fit_binmean_periodic] = WBfitting_singlevar(n_pol_WBmod,t_loc,Ddev_loc);

plot_WBmod_heatmap_Fenhance
saveas(gca,'WBmod_Fenhance.fig')
saveas(gca,'WBmod_Fenhance.png')
plot2svg('WBmod_Fenhance.svg')

% save data
cd ..
save('WBmod_Fenhance_data.mat',...
    't_F','Dstroke_F','Dpitch_F','Ddev_F',...
    'Dstroke_F_meanCIstd','Dpitch_F_meanCIstd','Ddev_F_meanCIstd',...
    'Dstroke_Fenhance_fit_binmean_periodic','Dpitch_Fenhance_fit_binmean_periodic','Ddev_Fenhance_fit_binmean_periodic',...
    'Fenhance_3std')
cd('WBmean_figs')


