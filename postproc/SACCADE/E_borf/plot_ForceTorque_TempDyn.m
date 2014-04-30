% pitch roll
figure
hold on

subplot(4,1,1)
plot(t_mean,-Fz_mean_mean_norm,'-g.')
ylabel('-Fz')
grid on
ylim([1,1.5])

subplot(4,1,2)
plot(t_mean,Mx_mean_mean_norm,'-g.')
ylabel('Troll')
grid on
ylim([-.03,.03])

subplot(4,1,3)
plot(t_mean,My_mean_mean_minWBstart_norm,'-g.')
ylabel('Tpitch')
grid on
ylim([-.03,.03])

subplot(4,1,4)
plot(t_mean,Mz_mean_mean_norm,'-g.')
ylabel('Tyaw')
xlabel('t')
grid on
ylim([0,.06])

saveas(gca, [save_name,'_roll_pitch.png'])

%% mu1 mu2
figure
hold on

subplot(4,1,1)
plot(t_mean,-Fz_mean_mean_norm,'-g.')
ylabel('-Fz')
grid on
ylim([1,1.5])

subplot(4,1,2)
plot(t_mean,M_R_mean_mean_minWBstart_norm,'-g.')
ylabel('Tm1')
grid on
ylim([-.03,.03])

subplot(4,1,3)
plot(t_mean,M_L_mean_mean_minWBstart_norm,'-g.')
ylabel('Tm2')
grid on
ylim([-.03,.03])

subplot(4,1,4)
plot(t_mean,Mz_mean_mean_norm,'-g.')
ylabel('Tyaw')
xlabel('t')
grid on
ylim([0,.06])

saveas(gca, [save_name,'_mu1mu2.png'])
