% make borf plots F&M timelines

clc
clear

load('borf_db_Fenhance_steady_NOcali_alldata.mat')
load('dt.mat')

for j = 10:5:20
    
close all
    
    
% butters
butter_cutoff = j;
butter_order = 5;

F_all = sqrt(Fx_all.^2 + Fz_all.^2 + Fz_all.^2);
M_all = sqrt(Mx_all.^2 + Mz_all.^2 + Mz_all.^2);

mod_values = mod_value_all;
color_map = jet(length(mod_values));
for i = 1:length(mod_values)
    
    figure(1)
    mod_now = mod_values(i);
    n = i;
    
    % Fx ALL
    val_now(:,:) = Fx_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now(:,:) = Fy_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now(:,:) = Fz_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now(:,:) = Mx_all(:,:,n)/Mg_fly/Lwing;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    
    % My ALL
    val_now(:,:) = My_all(:,:,n)/Mg_fly/Lwing;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
    
    % Mz ALL
    val_now(:,:) = Mz_all(:,:,n)/Mg_fly/Lwing;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
    val_now(:,:) = F_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
    % Mnorm ALL
    val_now(:,:) = M_all(:,:,n)/Mg_fly/Lwing;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M')
    
    % -Fz ALL
    val_now(:,:) = -Fz_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    val_mean_butter = filter_borf_butter_nWB(val_mean,butter_cutoff,butter_order,n_wb);
%     val_mean_butter = filter_borf_butter(val_mean,butter_cutoff,butter_order);
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
%     plot(t_now,val_mean_butter,'--','color',color_map(i,:),'linewidth',1)
    plot(t_now,val_mean_butter,'-k','linewidth',2)
%     plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_mean,.9999)
%     fnplt(val_pp);
%         ylim([-1 4])
    ylabel('-Fz/Mg')
    
    figure(2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
%     plot(t_now,val_mean_butter,'--','color',color_map(i,:),'linewidth',1)
    plot(t_now,val_mean_butter,'-k','linewidth',2)
%     plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_mean,.9999)
%     fnplt(val_pp);
%         ylim([-1 4])
    ylabel('-Fz/Mg')
    
    figure(3)
    plot(val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
%     plot(t_now,val_mean_butter,'--','color',color_map(i,:),'linewidth',1)
    plot(val_mean_butter,'-k','linewidth',2)
%     plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_mean,.9999)
%     fnplt(val_pp);
%         ylim([-1 4])
    ylabel('-Fz/Mg')
end

mkdir('steady_lowfreq_figs')
cd('steady_lowfreq_figs')
figure(1)
saveas(gca, ['steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.fig'])
saveas(gca, ['steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.png'])
figure(2)
saveas(gca, ['Fz_tnorm_steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.fig'])
saveas(gca, ['Fz_tnorm_steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.png'])
figure(3)
saveas(gca, ['Fz_t_steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.fig'])
saveas(gca, ['Fz_t_steady_lowfreq_timelines_butter',num2str(butter_cutoff),'.png'])
cd ..

end
