% make borf plots F&M timelines

clc
close all


%% Myaw
clear
load('MOD_norm_data.mat')

if exist(yawDB_all) == 2
load(yawDB_all)
load('MOD_norm_data.mat')

% steady Mz
mod_now = 0;
mod_diff = abs(mod_value_allNOfreq - mod_now);
n = find(mod_diff==min(mod_diff));
    clear val_now
val_now(:,:) = Mz_allNOfreq(:,:,n)/Mg_fly/Lwing;
val_steady = nanmean(val_now,2);

% mod plots
mod_values = YawMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(YawMods) max(YawMods)];
color_map = [.5 .5 .5; 1 .5 0];

F_allNOfreq = sqrt(Fx_allNOfreq.^2 + Fy_allNOfreq.^2 + Fz_allNOfreq.^2);
M_allNOfreq = sqrt(Mx_allNOfreq.^2 + My_allNOfreq.^2 + Mz_allNOfreq.^2);
for i = 1:length(mod_values)
figure(1)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
        clear val_now
val_now(:,:) = Fx_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
        clear val_now
val_now(:,:) = Fy_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
        clear val_now
val_now(:,:) = Fz_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % M_R ALL
        clear val_now
val_now(:,:) = M_R_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R')
    
    % M_L ALL
        clear val_now
val_now(:,:) = M_L_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L')
    
    % Mz ALL
       clear val_now
 val_now(:,:) = Mz_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
        clear val_now
val_now(:,:) = F_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
 %     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
   subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
    % Mnorm ALL
        clear val_now
val_now(:,:) = M_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M')
    
    % Mz ALL
        clear val_now
val_now(:,:) = Mz_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('Mz')
    ylim([-.2 .2])
    
    figure(3)
    subplot(2,2,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('Mz')
    ylim([-.1 .2])
end
end

%% M_R
clear
load('MOD_norm_data.mat')

if exist(RaxisDB_all) == 2
load(RaxisDB_all)
load('MOD_norm_data.mat')

% steady M_L
mod_now = 0;
mod_diff = abs(mod_value_allNOfreq - mod_now);
n = find(mod_diff==min(mod_diff));
    clear val_now
val_now(:,:) = M_L_allNOfreq(:,:,n)/Mg_fly/Lwing;
% val_steady = nanmean(val_now,2);

    j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_steady_L = nanmean(val_now(:));

% steady M_R
mod_now = 0;
mod_diff = abs(mod_value_allNOfreq - mod_now);
n = find(mod_diff==min(mod_diff));
    clear val_now
val_now(:,:) = M_R_allNOfreq(:,:,n)/Mg_fly/Lwing;
% val_steady = nanmean(val_now,2);

    j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_steady_R = nanmean(val_now(:));

% mod plots
mod_values = RaxisMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
mod_values = [-max(RaxisMods) min(RaxisMods) max(RaxisMods)];
color_map = [0 .5 1; .5 .5 .5; 1 .5 0];

F_allNOfreq = sqrt(Fx_allNOfreq.^2 + Fy_allNOfreq.^2 + Fz_allNOfreq.^2);
M_allNOfreq = sqrt(Mx_allNOfreq.^2 + My_allNOfreq.^2 + Mz_allNOfreq.^2);
for i = 1:length(mod_values)
figure(2)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
        clear val_now
val_now(:,:) = Fx_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
        clear val_now
val_now(:,:) = Fy_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
        clear val_now
val_now(:,:) = Fz_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % M_R ALL
        clear val_now
val_now(:,:) = M_R_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R')
    
    % M_L ALL
        clear val_now
val_now(:,:) = M_L_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L')
    
    % Mz ALL
       clear val_now
 val_now(:,:) = Mz_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
        clear val_now
val_now(:,:) = F_allNOfreq(:,:,n)/Mg_fly;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
 %     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
   subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
%     % Mnorm ALL
%         clear val_now
% val_now(:,:) = M_allNOfreq(:,:,n)/Mg_fly/Lwing;
%         j_wb = size(val_now,2);
% 
%     val_now=val_now(:);
%     val_now(val_now==0)=[];
%     val_now(isnan(val_now))=[];
%     val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);
% 
%     clear val_now
%     J_wb = length(val_now_butter)/j_wb;
%     for j=1:j_wb
%         val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
%     end
%     
% val_mean = nanmean(val_now,2);
%     val_std = nanstd(val_now')';
%     
%     val_mean = val_mean(isnan(val_mean)==0);
%     val_std = val_std(isnan(val_std)==0);
%     t_now = [0:1/(length(val_mean)-1):1];
%     
% %     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
%     
%     subplot(3,3,6)
%     plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
%     hold on
%     ylabel('M')
    
    % M_R-M_R_steady_mean ALL
        clear val_now
val_now(:,:) = M_R_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2) - val_steady_R;

    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('M_R - M_Rsteady')
    ylim([-.5 .5])
    
    figure(3)
    subplot(2,2,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('M_R - M_Rsteady')
    ylim([-.5 .5])
    
    
    % M_L-M_L_steady_mean ALL
        clear val_now
val_now(:,:) = M_L_allNOfreq(:,:,n)/Mg_fly/Lwing;
        j_wb = size(val_now,2);

    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
val_mean = nanmean(val_now,2) - val_steady_L;

    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     val_mean = filter_borf_butter_nWB(val_mean,butter_cut,butter_n,n_wb); 
    
    figure(2)
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('M_L - M_Lsteady')
    ylim([-.5 .5])
    
    figure(3)
    subplot(2,2,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    % fnplt(val_pp); 
    ylabel('M_L - M_Lsteady')
    ylim([-.5 .5])
    
end
end



mkdir(['MSfigs_butter_string',num2str(butter_cut)])
cd(['MSfigs_butter_string',num2str(butter_cut)])

figure(1)
saveas(gca, 'MSfig_Myaw_timelines_butter.fig')
saveas(gca, 'MSfig_Myaw_timelines_butter.png')
plot2svg(['MSfig_Myaw_timelines_butter.svg'])

figure(1)
saveas(gca, 'MSfig_M_R_timelines_butter.fig')
saveas(gca, 'MSfig_M_R_timelines_butter.png')
plot2svg(['MSfig_M_R_timelines_butter.svg'])

figure(3)
saveas(gca, 'MSfig_FnM_timelines_butter.fig')
saveas(gca, 'MSfig_FnM_timelines_butter.png')
plot2svg(['MSfig_FnM_timelines_butter.svg'])    
cd ..

