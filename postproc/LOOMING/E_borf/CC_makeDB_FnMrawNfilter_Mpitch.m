% make borf plots F&M timelines

clc
close all

%% Mpitch
clear
load('MOD_norm_data.mat')
load(pitchDB_all)
load('MOD_norm_data.mat')
load('DBkinematics_forces_torques.mat')

step = 10;

% Fx_allNOfreq_raw = nan(size(t_allNOfreq));
% Fy_allNOfreq_raw = nan(size(t_allNOfreq));
% Fz_allNOfreq_raw = nan(size(t_allNOfreq));
% Mx_allNOfreq_raw = nan(size(t_allNOfreq));
% My_allNOfreq_raw = nan(size(t_allNOfreq));
% Mz_allNOfreq_raw = nan(size(t_allNOfreq));
% 
% Fx_allNOfreq_butter = nan(size(t_allNOfreq));
% Fy_allNOfreq_butter = nan(size(t_allNOfreq));
% Fz_allNOfreq_butter = nan(size(t_allNOfreq));
% Mx_allNOfreq_butter = nan(size(t_allNOfreq));
% My_allNOfreq_butter = nan(size(t_allNOfreq));
% Mz_allNOfreq_butter = nan(size(t_allNOfreq));
% 
% Fx_freq_raw = nan(size(t_freq));
% Fy_freq_raw = nan(size(t_freq));
% Fz_freq_raw = nan(size(t_freq));
% Mx_freq_raw = nan(size(t_freq));
% My_freq_raw = nan(size(t_freq));
% Mz_freq_raw = nan(size(t_freq));
% 
% Fx_freq_butter = nan(size(t_freq));
% Fy_freq_butter = nan(size(t_freq));
% Fz_freq_butter = nan(size(t_freq));
% Mx_freq_butter = nan(size(t_freq));
% My_freq_butter = nan(size(t_freq));
% Mz_freq_butter = nan(size(t_freq));
% 
% Fx_stroke_raw = nan(size(t_stroke));
% Fy_stroke_raw = nan(size(t_stroke));
% Fz_stroke_raw = nan(size(t_stroke));
% Mx_stroke_raw = nan(size(t_stroke));
% My_stroke_raw = nan(size(t_stroke));
% Mz_stroke_raw = nan(size(t_stroke));
% 
% Fx_stroke_butter = nan(size(t_stroke));
% Fy_stroke_butter = nan(size(t_stroke));
% Fz_stroke_butter = nan(size(t_stroke));
% Mx_stroke_butter = nan(size(t_stroke));
% My_stroke_butter = nan(size(t_stroke));
% Mz_stroke_butter = nan(size(t_stroke));
% 
% Fx_dev_raw = nan(size(t_dev));
% Fy_dev_raw = nan(size(t_dev));
% Fz_dev_raw = nan(size(t_dev));
% Mx_dev_raw = nan(size(t_dev));
% My_dev_raw = nan(size(t_dev));
% Mz_dev_raw = nan(size(t_dev));
% 
% Fx_dev_butter = nan(size(t_dev));
% Fy_dev_butter = nan(size(t_dev));
% Fz_dev_butter = nan(size(t_dev));
% Mx_dev_butter = nan(size(t_dev));
% My_dev_butter = nan(size(t_dev));
% Mz_dev_butter = nan(size(t_dev));
% 
% Fx_pitch_raw = nan(size(t_pitch));
% Fy_pitch_raw = nan(size(t_pitch));
% Fz_pitch_raw = nan(size(t_pitch));
% Mx_pitch_raw = nan(size(t_pitch));
% My_pitch_raw = nan(size(t_pitch));
% Mz_pitch_raw = nan(size(t_pitch));
% 
% Fx_pitch_butter = nan(size(t_pitch));
% Fy_pitch_butter = nan(size(t_pitch));
% Fz_pitch_butter = nan(size(t_pitch));
% Mx_pitch_butter = nan(size(t_pitch));
% My_pitch_butter = nan(size(t_pitch));
% Mz_pitch_butter = nan(size(t_pitch));

for n = 1:size(t_allNOfreq,2)
    %% Fx allNOfreq
    clear val_now
    val_now(:,:) = Fx_allNOfreq(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[];
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fx_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    Fx_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    %% Fy allNOfreq
    clear val_now
    val_now(:,:) = Fy_allNOfreq(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fy_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    Fy_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Fz allNOfreq
    clear val_now
    val_now(:,:) = Fz_allNOfreq(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fz_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    Fz_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    
    
    %% Mx allNOfreq
    clear val_now
    val_now(:,:) = Mx_allNOfreq(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mx_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    Mx_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    %% My allNOfreq
    clear val_now
    val_now(:,:) = My_allNOfreq(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    My_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    My_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Mz allNOfreq
    clear val_now
    val_now(:,:) = Mz_allNOfreq(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mz_allNOfreq_raw(1:length(val_raw),n) = val_raw;
    Mz_allNOfreq_butter(1:length(val_butter),n) = val_butter;
    
    
    
    
    
    
    
    
%     %% Fx freq
%     clear val_now
%     val_now(:,:) = Fx_freq(:,:,n)/Mg_fly;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     Fx_freq_raw(1:length(val_raw),n) = val_raw;
%     Fx_freq_butter(1:length(val_butter),n) = val_butter;
%     
%     %% Fy freq
%     clear val_now
%     val_now(:,:) = Fy_freq(:,:,n)/Mg_fly;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     Fy_freq_raw(1:length(val_raw),n) = val_raw;
%     Fy_freq_butter(1:length(val_butter),n) = val_butter;
%     
%     
%     %% Fz freq
%     clear val_now
%     val_now(:,:) = Fz_freq(:,:,n)/Mg_fly;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     Fz_freq_raw(1:length(val_raw),n) = val_raw;
%     Fz_freq_butter(1:length(val_butter),n) = val_butter;
%     
%     
%     
%     %% Mx freq
%     clear val_now
%     val_now(:,:) = Mx_freq(:,:,n)/Mg_fly/Lwing;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     Mx_freq_raw(1:length(val_raw),n) = val_raw;
%     Mx_freq_butter(1:length(val_butter),n) = val_butter;
%     
%     %% My freq
%     clear val_now
%     val_now(:,:) = My_freq(:,:,n)/Mg_fly/Lwing;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     My_freq_raw(1:length(val_raw),n) = val_raw;
%     My_freq_butter(1:length(val_butter),n) = val_butter;
%     
%     
%     %% Mz freq
%     clear val_now
%     val_now(:,:) = Mz_freq(:,:,n)/Mg_fly/Lwing;
%     j_wb = size(val_now,2);
%     
%     % raw data
%     val_raw = nanmean(val_now,2);
%      val_raw(val_raw==0)=[];
%     val_raw(isnan(val_raw))=[]; 
% %     plot(val_raw)
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
%     val_butter = nanmean(val_now,2);
% %     hold on
% %     plot(val_butter,'r')
% 
%     Mz_freq_raw(1:length(val_raw),n) = val_raw;
%     Mz_freq_butter(1:length(val_butter),n) = val_butter;
    
    
    
    
    
    
    
    
    %% Fx stroke
    clear val_now
    val_now(:,:) = Fx_stroke(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
    val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[];
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fx_stroke_raw(1:length(val_raw),n) = val_raw;
    Fx_stroke_butter(1:length(val_butter),n) = val_butter;
    
    %% Fy stroke
    clear val_now
    val_now(:,:) = Fy_stroke(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
    val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[];
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fy_stroke_raw(1:length(val_raw),n) = val_raw;
    Fy_stroke_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Fz stroke
    clear val_now
    val_now(:,:) = Fz_stroke(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
    val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[];
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fz_stroke_raw(1:length(val_raw),n) = val_raw;
    Fz_stroke_butter(1:length(val_butter),n) = val_butter;
    
    
    
    %% Mx stroke
    clear val_now
    val_now(:,:) = Mx_stroke(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mx_stroke_raw(1:length(val_raw),n) = val_raw;
    Mx_stroke_butter(1:length(val_butter),n) = val_butter;
    
    %% My stroke
    clear val_now
    val_now(:,:) = My_stroke(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    My_stroke_raw(1:length(val_raw),n) = val_raw;
    My_stroke_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Mz stroke
    clear val_now
    val_now(:,:) = Mz_stroke(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mz_stroke_raw(1:length(val_raw),n) = val_raw;
    Mz_stroke_butter(1:length(val_butter),n) = val_butter;
    
    
    
    
    
    
    
    
    %% Fx pitch
    clear val_now
    val_now(:,:) = Fx_pitch(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fx_pitch_raw(1:length(val_raw),n) = val_raw;
    Fx_pitch_butter(1:length(val_butter),n) = val_butter;
    
    %% Fy pitch
    clear val_now
    val_now(:,:) = Fy_pitch(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fy_pitch_raw(1:length(val_raw),n) = val_raw;
    Fy_pitch_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Fz pitch
    clear val_now
    val_now(:,:) = Fz_pitch(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fz_pitch_raw(1:length(val_raw),n) = val_raw;
    Fz_pitch_butter(1:length(val_butter),n) = val_butter;
    
    
    
    %% Mx pitch
    clear val_now
    val_now(:,:) = Mx_pitch(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mx_pitch_raw(1:length(val_raw),n) = val_raw;
    Mx_pitch_butter(1:length(val_butter),n) = val_butter;
    
    %% My pitch
    clear val_now
    val_now(:,:) = My_pitch(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    My_pitch_raw(1:length(val_raw),n) = val_raw;
    My_pitch_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Mz pitch
    clear val_now
    val_now(:,:) = Mz_pitch(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mz_pitch_raw(1:length(val_raw),n) = val_raw;
    Mz_pitch_butter(1:length(val_butter),n) = val_butter;
    
    
    
    
    
    
    
    
    %% Fx dev
    clear val_now
    val_now(:,:) = Fx_dev(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fx_dev_raw(1:length(val_raw),n) = val_raw;
    Fx_dev_butter(1:length(val_butter),n) = val_butter;
    
    %% Fy dev
    clear val_now
    val_now(:,:) = Fy_dev(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fy_dev_raw(1:length(val_raw),n) = val_raw;
    Fy_dev_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Fz dev
    clear val_now
    val_now(:,:) = Fz_dev(:,:,n)/Mg_fly;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Fz_dev_raw(1:length(val_raw),n) = val_raw;
    Fz_dev_butter(1:length(val_butter),n) = val_butter;
    
    
    
    %% Mx dev
    clear val_now
    val_now(:,:) = Mx_dev(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mx_dev_raw(1:length(val_raw),n) = val_raw;
    Mx_dev_butter(1:length(val_butter),n) = val_butter;
    
    %% My dev
    clear val_now
    val_now(:,:) = My_dev(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    My_dev_raw(1:length(val_raw),n) = val_raw;
    My_dev_butter(1:length(val_butter),n) = val_butter;
    
    
    %% Mz dev
    clear val_now
    val_now(:,:) = Mz_dev(:,:,n)/Mg_fly/Lwing;
    j_wb = size(val_now,2);
    
    % raw data
    val_raw = nanmean(val_now,2);
     val_raw(val_raw==0)=[];
    val_raw(isnan(val_raw))=[]; 
%     plot(val_raw)
    
    val_now=val_now(:);
    val_now(val_now==0)=[];
    val_now(isnan(val_now))=[];
    val_now_butter = filter_borf_butter_nWB(val_now,butter_cut,butter_n,n_wb);

    clear val_now
    J_wb = length(val_now_butter)/j_wb;
    for j=1:j_wb
        val_now(:,j) =  val_now_butter((j-1)*J_wb+1:j*J_wb,1);
    end
    
    val_butter = nanmean(val_now,2);
%     hold on
%     plot(val_butter,'r')

    Mz_dev_raw(1:length(val_raw),n) = val_raw;
    Mz_dev_butter(1:length(val_butter),n) = val_butter;
    
end

% remove zeros
Fx_allNOfreq_raw(Fx_allNOfreq_raw==0)=nan;
Fy_allNOfreq_raw(Fy_allNOfreq_raw==0)=nan;
Fz_allNOfreq_raw(Fz_allNOfreq_raw==0)=nan;
Mx_allNOfreq_raw(Mx_allNOfreq_raw==0)=nan;
My_allNOfreq_raw(My_allNOfreq_raw==0)=nan;
Mz_allNOfreq_raw(Mz_allNOfreq_raw==0)=nan;

Fx_allNOfreq_butter(Fx_allNOfreq_butter==0)=nan;
Fy_allNOfreq_butter(Fy_allNOfreq_butter==0)=nan;
Fz_allNOfreq_butter(Fz_allNOfreq_butter==0)=nan;
Mx_allNOfreq_butter(Mx_allNOfreq_butter==0)=nan;
My_allNOfreq_butter(My_allNOfreq_butter==0)=nan;
Mz_allNOfreq_butter(Mz_allNOfreq_butter==0)=nan;
% 
% Fx_freq_raw(Fx_freq_raw==0)=nan;
% Fy_freq_raw(Fy_freq_raw==0)=nan;
% Fz_freq_raw(Fz_freq_raw==0)=nan;
% Mx_freq_raw(Mx_freq_raw==0)=nan;
% My_freq_raw(My_freq_raw==0)=nan;
% Mz_freq_raw(Mz_freq_raw==0)=nan;
% 
% Fx_freq_butter(Fx_freq_butter==0)=nan;
% Fy_freq_butter(Fy_freq_butter==0)=nan;
% Fz_freq_butter(Fz_freq_butter==0)=nan;
% Mx_freq_butter(Mx_freq_butter==0)=nan;
% My_freq_butter(My_freq_butter==0)=nan;
% Mz_freq_butter(Mz_freq_butter==0)=nan;
% 
Fx_stroke_raw(Fx_stroke_raw==0)=nan;
Fy_stroke_raw(Fy_stroke_raw==0)=nan;
Fz_stroke_raw(Fz_stroke_raw==0)=nan;
Mx_stroke_raw(Mx_stroke_raw==0)=nan;
My_stroke_raw(My_stroke_raw==0)=nan;
Mz_stroke_raw(Mz_stroke_raw==0)=nan;

Fx_stroke_butter(Fx_stroke_butter==0)=nan;
Fy_stroke_butter(Fy_stroke_butter==0)=nan;
Fz_stroke_butter(Fz_stroke_butter==0)=nan;
Mx_stroke_butter(Mx_stroke_butter==0)=nan;
My_stroke_butter(My_stroke_butter==0)=nan;
Mz_stroke_butter(Mz_stroke_butter==0)=nan;

Fx_dev_raw(Fx_dev_raw==0)=nan;
Fy_dev_raw(Fy_dev_raw==0)=nan;
Fz_dev_raw(Fz_dev_raw==0)=nan;
Mx_dev_raw(Mx_dev_raw==0)=nan;
My_dev_raw(My_dev_raw==0)=nan;
Mz_dev_raw(Mz_dev_raw==0)=nan;

Fx_dev_butter(Fx_dev_butter==0)=nan;
Fy_dev_butter(Fy_dev_butter==0)=nan;
Fz_dev_butter(Fz_dev_butter==0)=nan;
Mx_dev_butter(Mx_dev_butter==0)=nan;
My_dev_butter(My_dev_butter==0)=nan;
Mz_dev_butter(Mz_dev_butter==0)=nan;

Fx_pitch_raw(Fx_pitch_raw==0)=nan;
Fy_pitch_raw(Fy_pitch_raw==0)=nan;
Fz_pitch_raw(Fz_pitch_raw==0)=nan;
Mx_pitch_raw(Mx_pitch_raw==0)=nan;
My_pitch_raw(My_pitch_raw==0)=nan;
Mz_pitch_raw(Mz_pitch_raw==0)=nan;

Fx_pitch_butter(Fx_pitch_butter==0)=nan;
Fy_pitch_butter(Fy_pitch_butter==0)=nan;
Fz_pitch_butter(Fz_pitch_butter==0)=nan;
Mx_pitch_butter(Mx_pitch_butter==0)=nan;
My_pitch_butter(My_pitch_butter==0)=nan;
Mz_pitch_butter(Mz_pitch_butter==0)=nan;

% subset
Fx_allNOfreq_raw_sub = Fx_allNOfreq_raw(1:step:end,:);
Fy_allNOfreq_raw_sub = Fy_allNOfreq_raw(1:step:end,:);
Fz_allNOfreq_raw_sub = Fz_allNOfreq_raw(1:step:end,:);
Mx_allNOfreq_raw_sub = Mx_allNOfreq_raw(1:step:end,:);
My_allNOfreq_raw_sub = My_allNOfreq_raw(1:step:end,:);
Mz_allNOfreq_raw_sub = Mz_allNOfreq_raw(1:step:end,:);

Fx_allNOfreq_butter_sub = Fx_allNOfreq_butter(1:step:end,:);
Fy_allNOfreq_butter_sub = Fy_allNOfreq_butter(1:step:end,:);
Fz_allNOfreq_butter_sub = Fz_allNOfreq_butter(1:step:end,:);
Mx_allNOfreq_butter_sub = Mx_allNOfreq_butter(1:step:end,:);
My_allNOfreq_butter_sub = My_allNOfreq_butter(1:step:end,:);
Mz_allNOfreq_butter_sub = Mz_allNOfreq_butter(1:step:end,:);

% Fx_freq_raw_sub = Fx_freq_raw(1:step:end,:);
% Fy_freq_raw_sub = Fy_freq_raw(1:step:end,:);
% Fz_freq_raw_sub = Fz_freq_raw(1:step:end,:);
% Mx_freq_raw_sub = Mx_freq_raw(1:step:end,:);
% My_freq_raw_sub = My_freq_raw(1:step:end,:);
% Mz_freq_raw_sub = Mz_freq_raw(1:step:end,:);
% 
% Fx_freq_butter_sub = Fx_freq_butter(1:step:end,:);
% Fy_freq_butter_sub = Fy_freq_butter(1:step:end,:);
% Fz_freq_butter_sub = Fz_freq_butter(1:step:end,:);
% Mx_freq_butter_sub = Mx_freq_butter(1:step:end,:);
% My_freq_butter_sub = My_freq_butter(1:step:end,:);
% Mz_freq_butter_sub = Mz_freq_butter(1:step:end,:);

Fx_stroke_raw_sub = Fx_stroke_raw(1:step:end,:);
Fy_stroke_raw_sub = Fy_stroke_raw(1:step:end,:);
Fz_stroke_raw_sub = Fz_stroke_raw(1:step:end,:);
Mx_stroke_raw_sub = Mx_stroke_raw(1:step:end,:);
My_stroke_raw_sub = My_stroke_raw(1:step:end,:);
Mz_stroke_raw_sub = Mz_stroke_raw(1:step:end,:);

Fx_stroke_butter_sub = Fx_stroke_butter(1:step:end,:);
Fy_stroke_butter_sub = Fy_stroke_butter(1:step:end,:);
Fz_stroke_butter_sub = Fz_stroke_butter(1:step:end,:);
Mx_stroke_butter_sub = Mx_stroke_butter(1:step:end,:);
My_stroke_butter_sub = My_stroke_butter(1:step:end,:);
Mz_stroke_butter_sub = Mz_stroke_butter(1:step:end,:);

Fx_dev_raw_sub = Fx_dev_raw(1:step:end,:);
Fy_dev_raw_sub = Fy_dev_raw(1:step:end,:);
Fz_dev_raw_sub = Fz_dev_raw(1:step:end,:);
Mx_dev_raw_sub = Mx_dev_raw(1:step:end,:);
My_dev_raw_sub = My_dev_raw(1:step:end,:);
Mz_dev_raw_sub = Mz_dev_raw(1:step:end,:);

Fx_dev_butter_sub = Fx_dev_butter(1:step:end,:);
Fy_dev_butter_sub = Fy_dev_butter(1:step:end,:);
Fz_dev_butter_sub = Fz_dev_butter(1:step:end,:);
Mx_dev_butter_sub = Mx_dev_butter(1:step:end,:);
My_dev_butter_sub = My_dev_butter(1:step:end,:);
Mz_dev_butter_sub = Mz_dev_butter(1:step:end,:);

Fx_pitch_raw_sub = Fx_pitch_raw(1:step:end,:);
Fy_pitch_raw_sub = Fy_pitch_raw(1:step:end,:);
Fz_pitch_raw_sub = Fz_pitch_raw(1:step:end,:);
Mx_pitch_raw_sub = Mx_pitch_raw(1:step:end,:);
My_pitch_raw_sub = My_pitch_raw(1:step:end,:);
Mz_pitch_raw_sub = Mz_pitch_raw(1:step:end,:);

Fx_pitch_butter_sub = Fx_pitch_butter(1:step:end,:);
Fy_pitch_butter_sub = Fy_pitch_butter(1:step:end,:);
Fz_pitch_butter_sub = Fz_pitch_butter(1:step:end,:);
Mx_pitch_butter_sub = Mx_pitch_butter(1:step:end,:);
My_pitch_butter_sub = My_pitch_butter(1:step:end,:);
Mz_pitch_butter_sub = Mz_pitch_butter(1:step:end,:);

figure
subplot(2,3,1)
plot(pitch.t,Fx_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,Fx_allNOfreq_butter_sub,'r')
    
subplot(2,3,2)
plot(pitch.t,Fy_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,Fy_allNOfreq_butter_sub,'r')
    
subplot(2,3,3)
plot(pitch.t,Fz_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,Fz_allNOfreq_butter_sub,'r')
    
subplot(2,3,4)
plot(pitch.t,Mx_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,Mx_allNOfreq_butter_sub,'r')
    
subplot(2,3,5)
plot(pitch.t,My_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,My_allNOfreq_butter_sub,'r')
    
subplot(2,3,6)
plot(pitch.t,Mz_allNOfreq_raw_sub,'b')
hold on
plot(pitch.t,Mz_allNOfreq_butter_sub,'r')

% figure
% subplot(2,3,1)
% plot(pitch.t_INCfreq,Fx_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,Fx_freq_butter_sub,'r')
%     
% subplot(2,3,2)
% plot(pitch.t_INCfreq,Fy_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,Fy_freq_butter_sub,'r')
%     
% subplot(2,3,3)
% plot(pitch.t_INCfreq,Fz_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,Fz_freq_butter_sub,'r')
%     
% subplot(2,3,4)
% plot(pitch.t_INCfreq,Mx_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,Mx_freq_butter_sub,'r')
%     
% subplot(2,3,5)
% plot(pitch.t_INCfreq,My_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,My_freq_butter_sub,'r')
%     
% subplot(2,3,6)
% plot(pitch.t_INCfreq,Mz_freq_raw_sub,'b')
% hold on
% plot(pitch.t_INCfreq,Mz_freq_butter_sub,'r')
% 
figure
subplot(2,3,1)
plot(pitch.t,Fx_stroke_raw_sub,'b')
hold on
plot(pitch.t,Fx_stroke_butter_sub,'r')
    
subplot(2,3,2)
plot(pitch.t,Fy_stroke_raw_sub,'b')
hold on
plot(pitch.t,Fy_stroke_butter_sub,'r')
    
subplot(2,3,3)
plot(pitch.t,Fz_stroke_raw_sub,'b')
hold on
plot(pitch.t,Fz_stroke_butter_sub,'r')
    
subplot(2,3,4)
plot(pitch.t,Mx_stroke_raw_sub,'b')
hold on
plot(pitch.t,Mx_stroke_butter_sub,'r')
    
subplot(2,3,5)
plot(pitch.t,My_stroke_raw_sub,'b')
hold on
plot(pitch.t,My_stroke_butter_sub,'r')
    
subplot(2,3,6)
plot(pitch.t,Mz_stroke_raw_sub,'b')
hold on
plot(pitch.t,Mz_stroke_butter_sub,'r')

figure
subplot(2,3,1)
plot(pitch.t,Fx_dev_raw_sub,'b')
hold on
plot(pitch.t,Fx_dev_butter_sub,'r')
    
subplot(2,3,2)
plot(pitch.t,Fy_dev_raw_sub,'b')
hold on
plot(pitch.t,Fy_dev_butter_sub,'r')
    
subplot(2,3,3)
plot(pitch.t,Fz_dev_raw_sub,'b')
hold on
plot(pitch.t,Fz_dev_butter_sub,'r')
    
subplot(2,3,4)
plot(pitch.t,Mx_dev_raw_sub,'b')
hold on
plot(pitch.t,Mx_dev_butter_sub,'r')
    
subplot(2,3,5)
plot(pitch.t,My_dev_raw_sub,'b')
hold on
plot(pitch.t,My_dev_butter_sub,'r')
    
subplot(2,3,6)
plot(pitch.t,Mz_dev_raw_sub,'b')
hold on
plot(pitch.t,Mz_dev_butter_sub,'r')

figure
subplot(2,3,1)
plot(pitch.t,Fx_pitch_raw_sub,'b')
hold on
plot(pitch.t,Fx_pitch_butter_sub,'r')
    
subplot(2,3,2)
plot(pitch.t,Fy_pitch_raw_sub,'b')
hold on
plot(pitch.t,Fy_pitch_butter_sub,'r')
    
subplot(2,3,3)
plot(pitch.t,Fz_pitch_raw_sub,'b')
hold on
plot(pitch.t,Fz_pitch_butter_sub,'r')
    
subplot(2,3,4)
plot(pitch.t,Mx_pitch_raw_sub,'b')
hold on
plot(pitch.t,Mx_pitch_butter_sub,'r')
    
subplot(2,3,5)
plot(pitch.t,My_pitch_raw_sub,'b')
hold on
plot(pitch.t,My_pitch_butter_sub,'r')
    
subplot(2,3,6)
plot(pitch.t,Mz_pitch_raw_sub,'b')
hold on
plot(pitch.t,Mz_pitch_butter_sub,'r')
    
%% add to DB
pitch.Fx_norm_all = Fx_allNOfreq_raw_sub;
pitch.Fy_norm_all = Fy_allNOfreq_raw_sub;
pitch.Fz_norm_all = Fz_allNOfreq_raw_sub;

pitch.Mx_norm_all = Mx_allNOfreq_raw_sub;
pitch.My_norm_all = My_allNOfreq_raw_sub;
pitch.Mz_norm_all = Mz_allNOfreq_raw_sub;

pitch.Fx_norm_all_butterfilt = Fx_allNOfreq_butter_sub;
pitch.Fy_norm_all_butterfilt = Fy_allNOfreq_butter_sub;
pitch.Fz_norm_all_butterfilt = Fz_allNOfreq_butter_sub;

pitch.Mx_norm_all_butterfilt = Mx_allNOfreq_butter_sub;
pitch.My_norm_all_butterfilt = My_allNOfreq_butter_sub;
pitch.Mz_norm_all_butterfilt = Mz_allNOfreq_butter_sub;

% pitch.Fx_norm_freq = Fx_freq_raw_sub;
% pitch.Fy_norm_freq = Fy_freq_raw_sub;
% pitch.Fz_norm_freq = Fz_freq_raw_sub;
% 
% pitch.Mx_norm_freq = Mx_freq_raw_sub;
% pitch.My_norm_freq = My_freq_raw_sub;
% pitch.Mz_norm_freq = Mz_freq_raw_sub;
% 
% pitch.Fx_norm_freq_butterfilt = Fx_freq_butter_sub;
% pitch.Fy_norm_freq_butterfilt = Fy_freq_butter_sub;
% pitch.Fz_norm_freq_butterfilt = Fz_freq_butter_sub;
% 
% pitch.Mx_norm_freq_butterfilt = Mx_freq_butter_sub;
% pitch.My_norm_freq_butterfilt = My_freq_butter_sub;
% pitch.Mz_norm_freq_butterfilt = Mz_freq_butter_sub;

pitch.Fx_norm_stroke = Fx_stroke_raw_sub;
pitch.Fy_norm_stroke = Fy_stroke_raw_sub;
pitch.Fz_norm_stroke = Fz_stroke_raw_sub;

pitch.Mx_norm_stroke = Mx_stroke_raw_sub;
pitch.My_norm_stroke = My_stroke_raw_sub;
pitch.Mz_norm_stroke = Mz_stroke_raw_sub;

pitch.Fx_norm_stroke_butterfilt = Fx_stroke_butter_sub;
pitch.Fy_norm_stroke_butterfilt = Fy_stroke_butter_sub;
pitch.Fz_norm_stroke_butterfilt = Fz_stroke_butter_sub;

pitch.Mx_norm_stroke_butterfilt = Mx_stroke_butter_sub;
pitch.My_norm_stroke_butterfilt = My_stroke_butter_sub;
pitch.Mz_norm_stroke_butterfilt = Mz_stroke_butter_sub;

pitch.Fx_norm_rotation = Fx_pitch_raw_sub;
pitch.Fy_norm_rotation = Fy_pitch_raw_sub;
pitch.Fz_norm_rotation = Fz_pitch_raw_sub;

pitch.Mx_norm_rotation = Mx_pitch_raw_sub;
pitch.My_norm_rotation = My_pitch_raw_sub;
pitch.Mz_norm_rotation = Mz_pitch_raw_sub;

pitch.Fx_norm_rotation_butterfilt = Fx_pitch_butter_sub;
pitch.Fy_norm_rotation_butterfilt = Fy_pitch_butter_sub;
pitch.Fz_norm_rotation_butterfilt = Fz_pitch_butter_sub;

pitch.Mx_norm_rotation_butterfilt = Mx_pitch_butter_sub;
pitch.My_norm_rotation_butterfilt = My_pitch_butter_sub;
pitch.Mz_norm_rotation_butterfilt = Mz_pitch_butter_sub;

pitch.Fx_norm_deviation = Fx_dev_raw_sub;
pitch.Fy_norm_deviation = Fy_dev_raw_sub;
pitch.Fz_norm_deviation = Fz_dev_raw_sub;

pitch.Mx_norm_deviation = Mx_dev_raw_sub;
pitch.My_norm_deviation = My_dev_raw_sub;
pitch.Mz_norm_deviation = Mz_dev_raw_sub;

pitch.Fx_norm_deviation_butterfilt = Fx_dev_butter_sub;
pitch.Fy_norm_deviation_butterfilt = Fy_dev_butter_sub;
pitch.Fz_norm_deviation_butterfilt = Fz_dev_butter_sub;

pitch.Mx_norm_deviation_butterfilt = Mx_dev_butter_sub;
pitch.My_norm_deviation_butterfilt = My_dev_butter_sub;
pitch.Mz_norm_deviation_butterfilt = Mz_dev_butter_sub;



save('DBkinematics_pitchs_torques_NEWpitch.mat','pitch')
