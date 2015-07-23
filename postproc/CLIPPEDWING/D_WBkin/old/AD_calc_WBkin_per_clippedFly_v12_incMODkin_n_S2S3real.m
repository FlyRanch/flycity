clc
clear
close all

Eqname=dir('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq*')
Eqname=Eqname.name;
load(Eqname)

MODname=dir('WBdataset_steadyNclipMods_S2S3AmpRatioFunc.mat')
MODname=MODname.name;
load(MODname)

loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

plot_on = 1;
plot_on = 0;

%% plot dir
if plot_on == 1
    mkdir('steadyWBkinNwbMODkin_seqs_figs_S2S3AmpRatioFunc')
    cd('steadyWBkinNwbMODkin_seqs_figs_S2S3AmpRatioFunc')
end

%% steady wb
freq_steady = f_wb_steady_meanCIstd;

stroke_steady = stroke_wb_steady_bins_meanCIstd;
rot_steady = pitch_wb_steady_bins_meanCIstd;
dev_steady = dev_wb_steady_bins_meanCIstd;

Astroke_steady = max(stroke_steady(:,1)) - min(stroke_steady(:,1));

%% wb MODS
strokeMODclipped_S2S3AmpRatioFunc = strokeMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd;
rotMODclipped_S2S3AmpRatioFunc = pitchMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd;
devMODclipped_S2S3AmpRatioFunc = devMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd;

strokeMODintact_S2S3AmpRatioFunc = strokeMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd;
rotMODintact_S2S3AmpRatioFunc = pitchMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd;
devMODintact_S2S3AmpRatioFunc = devMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd;
            
%% loop through SecondMomentRatio's etc
SecondMomentRatio_list      = unique(SecondMomentRatio);

% make nan DB arrays
L_ratio_all = nan(150,length(SecondMomentRatio_list));
S2_ratio_all = nan(150,length(SecondMomentRatio_list));
S3_ratio_all = nan(150,length(SecondMomentRatio_list));

Area_normClipped_all  = nan(150,length(SecondMomentRatio_list));
S2_normClipped_all  = nan(150,length(SecondMomentRatio_list));
S3_normClipped_all  = nan(150,length(SecondMomentRatio_list));

Area_normIntact_all  = nan(150,length(SecondMomentRatio_list));
S2_normIntact_all  = nan(150,length(SecondMomentRatio_list));
S3_normIntact_all  = nan(150,length(SecondMomentRatio_list));

L_all = nan(150,length(SecondMomentRatio_list));
freq_all = nan(150,length(SecondMomentRatio_list));

Astroke_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_intact_all = nan(150,length(SecondMomentRatio_list));

Astroke_ratio_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_intact_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_clip_intact_all = nan(150,length(SecondMomentRatio_list));

clip_type_all = nan(150,length(SecondMomentRatio_list));

for seq_now = 1:length(SecondMomentRatio_list)
    
    S2_now = SecondMomentRatio_list(seq_now)
    wb = find(SecondMomentRatio==S2_now & steady_nr_mean_wb==1);
    
                % current wb
            seq_nr_now = seq_nr(wb);
            wb_nr_now = wb_nr(wb);
    
%             if mean(seq_nr_now) ~= min(seq_nr_now);
%                 ERROR = 1
%                 break
%             end
            
            %% clip data
            clip_side_now = clip_side(wb);
            clip_type_now = clip_type(wb);

            FirstMomentNormClipped_now = FirstMomentNormClipped(wb);
            FirstMomentNormIntact_now = FirstMomentNormIntact(wb);
            FirstMomentRatio_now = FirstMomentRatio(wb);

            SecondMomentNormClipped_now = SecondMomentNormClipped(wb);
            SecondMomentNormIntact_now = SecondMomentNormIntact(wb);
            SecondMomentRatio_now = SecondMomentRatio(wb);

            ThirdMomentNormClipped_now = ThirdMomentNormClipped(wb);
            ThirdMomentNormIntact_now = ThirdMomentNormIntact(wb);
            ThirdMomentRatio_now = ThirdMomentRatio(wb);

            AreaNormClipped_now = AreaNormClipped(wb);
            AreaNormIntact_now = AreaNormIntact(wb);
            AreaRatio_now = AreaRatio(wb);

            LengthClipped_pixels_now = LengthClipped_pixels(wb);
            LengthIntact_pixels_now = LengthIntact_pixels(wb);
            LengthRatio_now = LengthRatio(wb);
            
            %% S2S3AmpRatioFunc
            sol = subs(solAdAiRatio,S2,mean(SecondMomentRatio_now));
            sol = subs(sol,S3,mean(ThirdMomentRatio_now));
            S2S3AmpRatioFunc_now = eval(sol);
            
            sol = subs(solAd,S2,mean(SecondMomentRatio_now));
            sol = subs(sol,S3,mean(ThirdMomentRatio_now));
            S2S3AmpFuncClipped_now = eval(sol);
            
            sol = subs(solAi,S2,mean(SecondMomentRatio_now));
            sol = subs(sol,S3,mean(ThirdMomentRatio_now));
            S2S3AmpFuncIntact_now = eval(sol);
            
            %% WB kin based on S2S3AmpRatioFunc MODs
            strokeClipped_S2S3AmpRatioFunc_now   = stroke_steady(:,1)    + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * strokeMODclipped_S2S3AmpRatioFunc(:,1);
            rotClipped_S2S3AmpRatioFunc_now      = rot_steady(:,1)       + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * rotMODclipped_S2S3AmpRatioFunc(:,1);
            devClipped_S2S3AmpRatioFunc_now      = dev_steady(:,1)       + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * devMODclipped_S2S3AmpRatioFunc(:,1);
            
            strokeIntact_S2S3AmpRatioFunc_now   = stroke_steady(:,1)    + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * strokeMODintact_S2S3AmpRatioFunc(:,1);
            rotIntact_S2S3AmpRatioFunc_now      = rot_steady(:,1)       + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * rotMODintact_S2S3AmpRatioFunc(:,1);
            devIntact_S2S3AmpRatioFunc_now      = dev_steady(:,1)       + (nanmean(S2S3AmpRatioFunc_now)-S2S3AmpRatioFunc_NONclipped) * devMODintact_S2S3AmpRatioFunc(:,1);
            
            % Stroke amplitudes from MODs
%             AstrokeSteady = max(stroke_steady(:,1)) - min(stroke_steady(:,1));
            
            AstrokeIntact_S2S3AmpRatioFunc_now = max(strokeIntact_S2S3AmpRatioFunc_now) - min(strokeIntact_S2S3AmpRatioFunc_now);
            AstrokeClipped_S2S3AmpRatioFunc_now = max(strokeClipped_S2S3AmpRatioFunc_now) - min(strokeClipped_S2S3AmpRatioFunc_now);
            AstrokeRatio_S2S3AmpRatioFunc_now = AstrokeClipped_S2S3AmpRatioFunc_now / AstrokeIntact_S2S3AmpRatioFunc_now;

            AstrokeRatioClipped_S2S3AmpRatioFunc_now = AstrokeClipped_S2S3AmpRatioFunc_now / Astroke_steady;
            AstrokeRatioIntact_S2S3AmpRatioFunc_now = AstrokeIntact_S2S3AmpRatioFunc_now / Astroke_steady;

            %% body kin data
            vel_now = V_mean_wb(wb);
            
            slip_now = slip_mean_wb(wb);
            roll_now = roll_mean_wb(wb);
            pitch_now = pitch_mean_wb(wb);
            
            Fsp_pitch_now = Fsp_pitch_mean_wb(wb);
            Fsp_roll_now = Fsp_roll_mean_wb(wb);
            
            % wingbeat kin data
            l_wb_now = l_wing_mean_wb(wb);
            
            f_wb_now = mean([f_wb_L(wb),f_wb_R(wb)],2);
            
            stroke_clip_now = stroke_wb_R_bins(:,wb);
            stroke_intact_now = stroke_wb_L_bins(:,wb);

            pitch_wb_clip_now = pitch_wb_R_bins(:,wb);
            pitch_wb_intact_now = pitch_wb_L_bins(:,wb);

            dev_wb_clip_now = dev_wb_R_bins(:,wb);
            dev_wb_intact_now = dev_wb_L_bins(:,wb);

            U_wb_clip_now = U_wb_R_bins(:,wb);
            U_wb_intact_now = U_wb_L_bins(:,wb);

            aoa_wb_clip_now = aoa_wb_R_bins(:,wb);
            aoa_wb_intact_now = aoa_wb_L_bins(:,wb);

            Dstroke = Dstroke_wb_bins(:,wb);
            Dpitch = Dpitch_wb_bins(:,wb);
            Ddev = Ddev_wb_bins(:,wb);
            DU = DU_wb_bins(:,wb);
            Daoa = Daoa_wb_bins(:,wb);
            
            %% calc amplitudes
            Astroke_clip_now = (max(stroke_clip_now)-min(stroke_clip_now))';
            Astroke_intact_now = (max(stroke_intact_now)-min(stroke_intact_now))';
            
            Astroke_ratio_clip_now = Astroke_clip_now ./ Astroke_steady;
            Astroke_ratio_intact_now = Astroke_intact_now ./ Astroke_steady;
            Astroke_ratio_clip_intact_now = Astroke_clip_now ./ Astroke_intact_now;

            %% store data in DB
            
            % model data
            S2S3AmpFuncClipped_all(seq_now,1) = S2S3AmpFuncClipped_now;
            S2S3AmpFuncIntact_all(seq_now,1) = S2S3AmpFuncIntact_now;
            S2S3AmpRatioFunc_all(seq_now,1) = S2S3AmpRatioFunc_now;
            
            % MOD data
            AstrokeRatioClipped_S2S3AmpRatioFunc_all(seq_now,1) = AstrokeRatioClipped_S2S3AmpRatioFunc_now;
            AstrokeRatioIntact_S2S3AmpRatioFunc_all(seq_now,1) = AstrokeRatioIntact_S2S3AmpRatioFunc_now;
            AstrokeRatio_S2S3AmpRatioFunc_all(seq_now,1) = AstrokeRatio_S2S3AmpRatioFunc_now;
            
            % real fly data
            L_ratio_all(1:length(LengthRatio_now),seq_now) = LengthRatio_now;
            S2_ratio_all(1:length(SecondMomentRatio_now),seq_now) = SecondMomentRatio_now;
            S3_ratio_all(1:length(ThirdMomentRatio_now),seq_now) = ThirdMomentRatio_now;
            
            Area_normClipped_all(1:length(AreaNormClipped_now),seq_now) = AreaNormClipped_now;
            S2_normClipped_all(1:length(SecondMomentNormClipped_now),seq_now) = SecondMomentNormClipped_now;
            S3_normClipped_all(1:length(ThirdMomentNormClipped_now),seq_now) = ThirdMomentNormClipped_now;

            Area_normIntact_all(1:length(AreaNormIntact_now),seq_now) = AreaNormIntact_now;
            S2_normIntact_all(1:length(SecondMomentNormIntact_now),seq_now) = SecondMomentNormIntact_now;
            S3_normIntact_all(1:length(ThirdMomentNormIntact_now),seq_now) = ThirdMomentNormIntact_now;

            L_all(1:length(l_wb_now),seq_now) = l_wb_now;
            freq_all(1:length(f_wb_now),seq_now) = f_wb_now;
   
            Astroke_clip_all(1:length(Astroke_clip_now),seq_now) = Astroke_clip_now;
            Astroke_intact_all(1:length(Astroke_intact_now),seq_now) = Astroke_intact_now;
   
            Astroke_ratio_clip_all(1:length(Astroke_ratio_clip_now),seq_now) = Astroke_ratio_clip_now;
            Astroke_ratio_intact_all(1:length(Astroke_ratio_intact_now),seq_now) = Astroke_ratio_intact_now;
            Astroke_ratio_clip_intact_all(1:length(Astroke_ratio_clip_intact_now),seq_now) = Astroke_ratio_clip_intact_now;
            
            clip_type_all(1:length(clip_type_now),seq_now) = clip_type_now;
   
            
            %% plot
            if plot_on == 1
                
                % all WBs
                subplot(2,2,1)
                plot(stroke_clip_now,'r')
                hold on
                plot(stroke_intact_now,'c')
                plot(strokeIntact_S2S3AmpRatioFunc_now(:,1),'b','linewidth',2)
                plot(strokeClipped_S2S3AmpRatioFunc_now(:,1),'color',[1 .5 0],'linewidth',2)
%                 plot(stroke_steady(:,1),'g','linewidth',2)
                ylim([-90 90])
                hold off

                subplot(2,2,2)
                plot(pitch_wb_clip_now-90,'r')
                hold on
                plot(pitch_wb_intact_now-90,'c')
                plot(rotIntact_S2S3AmpRatioFunc_now(:,1)-90,'b','linewidth',2)
                plot(rotClipped_S2S3AmpRatioFunc_now(:,1)-90,'color',[1 .5 0],'linewidth',2)
%                 plot(rot_steady(:,1)-90,'g','linewidth',2)
                ylim([-90 90])
                hold off

                subplot(2,2,3)
                plot(dev_wb_clip_now,'r')
                hold on
                plot(dev_wb_intact_now,'c')
                plot(devIntact_S2S3AmpRatioFunc_now(:,1),'b','linewidth',2)
                plot(devClipped_S2S3AmpRatioFunc_now(:,1),'color',[1 .5 0],'linewidth',2)
%                 plot(dev_steady(:,1),'g','linewidth',2)
                ylim([-90 90])
                hold off
                
                subplot(2,2,4)
                plot(0,0,'-r','linewidth',2)
                hold on
                plot(0,0,'-','color',[1 .5 0],'linewidth',2)
                plot(0,0,'-c','linewidth',2)
                plot(0,0,'-b','linewidth',2)
%                 plot(0,0,'-g','linewidth',2)
                hold off
%                 legend('clipped wing','intact wing','steady wingbeat','clipped MOD','intact MOD')
                legend('clipped wing','clipped MOD','intact wing','intact MOD')
                axis off

                saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.fig'])
                saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.png'])
%                 saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.svg'])
                plot2svg(['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.svg'])

                
                %% MEAN WBs
                ny = 200;
                biny_min = -90;
                biny_max = 90;
                biny = biny_min: (biny_max - biny_min)/ny :biny_max;

                % stroke
                var_clip = [stroke_clip_now];
                var_intact = [stroke_intact_now];
                
                t_hist_clip = t_wb_steady_bins(:,1:size(var_clip,2));
                var_clip = -var_clip(:);
                t_hist_clip = t_hist_clip(:);
                var_hist_clip = hist3([var_clip,t_hist_clip], {biny binx});

                t_hist_intact = t_wb_steady_bins(:,1:size(var_intact,2));
                var_intact = -var_intact(:);
                t_hist_intact = t_hist_intact(:);
                var_hist_intact = hist3([var_intact,t_hist_intact], {biny binx});

                for i = 1:length(binx)
                    var_hist_clip(:,i) = var_hist_clip(:,i) / max(var_hist_clip(:,i)); % normalize per time bin
                    var_hist_intact(:,i) = var_hist_intact(:,i) / max(var_hist_intact(:,i)); % normalize per time bin
                end

                % inverse
                var_hist_clip = 1 - var_hist_clip;
                var_hist_intact = 1 - var_hist_intact;

                clear var_hist
                var_hist(:,:) = var_hist_intact;
                var_hist(:,:,2) = var_hist_clip;
                var_hist(:,:,3) = var_hist_clip;


                subplot(2,2,1)
                imagesc(binx,biny,var_hist)
                axis([0 1 biny_min biny_max])
                %     set(gca,'XTick',[],'XTickLabel',[]) 
                % set(gca,'YTick',[],'fontsize',8) 
                % axis off

                hold on
%                 plot(t_wb_bin,-nanmean(stroke_clip_now,2),'--k','linewidth',2)
%                 plot(t_wb_bin,-nanmean(stroke_intact_now,2),'--k','linewidth',2)
                plot(t_wb_bin,-strokeIntact_S2S3AmpRatioFunc_now(:,1),'b','linewidth',2)
                plot(t_wb_bin,-strokeClipped_S2S3AmpRatioFunc_now(:,1),'color',[1 .5 0],'linewidth',2)
%                 plot(stroke_steady(:,1),'g','linewidth',2)
%                 ylim([-90 90])
                hold off
                
                % ROTATION
                var_clip = [pitch_wb_clip_now]-90;
                var_intact = [pitch_wb_intact_now]-90;

                t_hist_clip = t_wb_steady_bins(:,1:size(var_clip,2));
                var_clip = -var_clip(:);
                t_hist_clip = t_hist_clip(:);
                var_hist_clip = hist3([var_clip,t_hist_clip], {biny binx});

                t_hist_intact = t_wb_steady_bins(:,1:size(var_intact,2));
                var_intact = -var_intact(:);
                t_hist_intact = t_hist_intact(:);
                var_hist_intact = hist3([var_intact,t_hist_intact], {biny binx});

                for i = 1:length(binx)
                    var_hist_clip(:,i) = var_hist_clip(:,i) / max(var_hist_clip(:,i)); % normalize per time bin
                    var_hist_intact(:,i) = var_hist_intact(:,i) / max(var_hist_intact(:,i)); % normalize per time bin
                end

                % inverse
                var_hist_clip = 1 - var_hist_clip;
                var_hist_intact = 1 - var_hist_intact;

                clear var_hist
                var_hist(:,:) = var_hist_intact;
                var_hist(:,:,2) = var_hist_clip;
                var_hist(:,:,3) = var_hist_clip;


                subplot(2,2,2)
                imagesc(binx,biny,var_hist)
                axis([0 1 biny_min biny_max])
                %     set(gca,'XTick',[],'XTickLabel',[]) 
                % set(gca,'YTick',[],'fontsize',8) 
                % axis off

                hold on
%                 plot(t_wb_bin,-nanmean(pitch_wb_clip_now,2),'--k','linewidth',2)
%                 plot(t_wb_bin,-nanmean(pitch_wb_intact_now,2),'--k','linewidth',2)
                plot(t_wb_bin,-rotIntact_S2S3AmpRatioFunc_now(:,1)+90,'b','linewidth',2)
                plot(t_wb_bin,-rotClipped_S2S3AmpRatioFunc_now(:,1)+90,'color',[1 .5 0],'linewidth',2)
%                 plot(stroke_steady(:,1),'g','linewidth',2)
%                 ylim([-90 90])
                hold off
                
                % dev
                var_clip = [dev_wb_clip_now];
                var_intact = [dev_wb_intact_now];

                t_hist_clip = t_wb_steady_bins(:,1:size(var_clip,2));
                var_clip = -var_clip(:);
                t_hist_clip = t_hist_clip(:);
                var_hist_clip = hist3([var_clip,t_hist_clip], {biny binx});

                t_hist_intact = t_wb_steady_bins(:,1:size(var_intact,2));
                var_intact = -var_intact(:);
                t_hist_intact = t_hist_intact(:);
                var_hist_intact = hist3([var_intact,t_hist_intact], {biny binx});

                for i = 1:length(binx)
                    var_hist_clip(:,i) = var_hist_clip(:,i) / max(var_hist_clip(:,i)); % normalize per time bin
                    var_hist_intact(:,i) = var_hist_intact(:,i) / max(var_hist_intact(:,i)); % normalize per time bin
                end

                % inverse
                var_hist_clip = 1 - var_hist_clip;
                var_hist_intact = 1 - var_hist_intact;

                clear var_hist
                var_hist(:,:) = var_hist_intact;
                var_hist(:,:,2) = var_hist_clip;
                var_hist(:,:,3) = var_hist_clip;


                subplot(2,2,3)
                imagesc(binx,biny,var_hist)
                axis([0 1 biny_min biny_max])
                %     set(gca,'XTick',[],'XTickLabel',[]) 
                % set(gca,'YTick',[],'fontsize',8) 
                % axis off

                hold on
%                 plot(t_wb_bin,-nanmean(dev_wb_clip_now,2),'--k','linewidth',2)
%                 plot(t_wb_bin,-nanmean(dev_wb_intact_now,2),'--k','linewidth',2)
                plot(t_wb_bin,-devIntact_S2S3AmpRatioFunc_now(:,1),'b','linewidth',2)
                plot(t_wb_bin,-devClipped_S2S3AmpRatioFunc_now(:,1),'color',[1 .5 0],'linewidth',2)
%                 plot(stroke_steady(:,1),'g','linewidth',2)
%                 ylim([-90 90])
                hold off
                
                subplot(2,2,4)
                plot(0,0,'-r','linewidth',2)
                hold on
                plot(0,0,'-','color',[1 .5 0],'linewidth',2)
                plot(0,0,'-c','linewidth',2)
%                 plot(0,0,'-g','linewidth',2)
                plot(0,0,'-b','linewidth',2)
                hold off
%                 legend('clipped wing','intact wing','steady wingbeat','clipped MOD','intact MOD')
                legend('clipped wing','clipped MOD','intact wing','intact MOD')
                axis off
                
                saveas(gcf,['HEATMAPsteadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.fig'])
                saveas(gcf,['HEATMAPsteadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.png'])
%                 saveas(gcf,['HEATMAPsteadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.svg'])
                plot2svg(['HEATMAPsteadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_cliptype=',num2str(clip_type_now(1)),'_n=',num2str(length(seq_nr_now)),'.svg'])

            end
end

if plot_on == 1
    cd ..
end

L_ratio_mean = nanmean(L_ratio_all)';
S2_ratio_mean = nanmean(S2_ratio_all)';
S3_ratio_mean = nanmean(S3_ratio_all)';

Area_normClipped_mean = nanmean(Area_normClipped_all)';
S2_normClipped_mean = nanmean(S2_normClipped_all)';
S3_normClipped_mean = nanmean(S3_normClipped_all)';

Area_normIntact_mean = nanmean(Area_normIntact_all)';
S2_normIntact_mean = nanmean(S2_normIntact_all)';
S3_normIntact_mean = nanmean(S3_normIntact_all)';

L_mean = nanmean(L_all)';
freq_mean = nanmean(freq_all)';
freqRatio_mean = freq_mean/freq_steady(1);

Astroke_clip_mean = nanmean(Astroke_clip_all)';
Astroke_intact_mean = nanmean(Astroke_intact_all)';

Astroke_ratio_clip_mean = nanmean(Astroke_ratio_clip_all)';
Astroke_ratio_clip_std = nanstd(Astroke_ratio_clip_all)';
Astroke_ratio_clip_n = nansum(0*(Astroke_ratio_clip_all)+1)';
Astroke_ratio_clip_ste = Astroke_ratio_clip_std ./ sqrt(Astroke_ratio_clip_n);

Astroke_ratio_intact_mean = nanmean(Astroke_ratio_intact_all)';
Astroke_ratio_intact_std = nanstd(Astroke_ratio_intact_all)';
Astroke_ratio_intact_n = nansum(0*(Astroke_ratio_intact_all)+1)';
Astroke_ratio_intact_ste = Astroke_ratio_intact_std ./ sqrt(Astroke_ratio_intact_n);

Astroke_ratio_clip_intact_mean = nanmean(Astroke_ratio_clip_intact_all)';
Astroke_ratio_clip_intact_std = nanstd(Astroke_ratio_clip_intact_all)';
Astroke_ratio_clip_intact_n = nansum(0*(Astroke_ratio_clip_intact_all)+1)';
Astroke_ratio_clip_intact_ste = Astroke_ratio_clip_intact_std ./ sqrt(Astroke_ratio_clip_intact_n);

clip_type_mean = nanmean(clip_type_all)';

%% calc real length & S2 & S3 of iontact and clipped wings

% wing length in [m], VARIABLE INTACT WING SIZE
L_Intact_mean = 2.*L_mean ./ (L_ratio_mean + 1) ./1000;
L_Clipped_mean = L_Intact_mean .* L_ratio_mean;

% % wing length in [m], CONSTANT INTACT WING SIZE
% L_Intact_mean = 2.*mean(L_mean) ./ (L_ratio_mean + 1) ./1000;
% L_Clipped_mean = L_Intact_mean .* L_ratio_mean;

Area_Clipped_mean = Area_normClipped_mean .* (L_Clipped_mean).^2;
S2_Clipped_mean = S2_normClipped_mean .* Area_normClipped_mean .* (L_Clipped_mean).^4;
S3_Clipped_mean = S3_normClipped_mean .* Area_normClipped_mean .* (L_Clipped_mean).^5;

AR_Clipped_mean = 1./Area_normClipped_mean;
r2_Clipped_mean = S2_Clipped_mean.^(1/4);
r3_Clipped_mean = S3_Clipped_mean.^(1/5);

Area_Intact_mean = Area_normIntact_mean .* (L_Intact_mean).^2;
S2_Intact_mean = S2_normIntact_mean .* Area_normIntact_mean .* (L_Intact_mean).^4;
S3_Intact_mean = S3_normIntact_mean .* Area_normIntact_mean .* (L_Intact_mean).^5;

AR_Intact_mean = 1./Area_normIntact_mean;
r2_Intact_mean = S2_Intact_mean.^(1/4);
r3_Intact_mean = S3_Intact_mean.^(1/5);

%% plot
mkdir('clippedfly_steadyWBkin_param_figs')
cd('clippedfly_steadyWBkin_param_figs')

% datapoints with color
% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);
cmap_Aratio = cmap_surf;

S2norm_min = .5;
S2norm_max = 1;
S3norm_min = .5;
S3norm_max = 1;

freqRatio_max = 1.25;
freqRatio_min = .75;

S2_min = 4e-12;
S2_max = 11e-12;
S3_min = .5e-14;
S3_max = 2.5e-14;
Area_min = 1.5e-6;
Area_max = 3.5e-6;

r2_min = 1.4e-3;
r2_max = 2e-3;
r3_min = 1.4e-3;
r3_max = 2e-3;
AR_min = 2;
AR_max = 6;

%% r2 & r3 NON-normailzed VS freq_clippedfly & AR
figure

% freq
subplot(2,2,3)
color_nr_steady = round(99/(freqRatio_max-freqRatio_min)*(1-freqRatio_min)+1);

for i = 1:length(freqRatio_mean)
    color_nr = round(99/(freqRatio_max-freqRatio_min)*(freqRatio_mean(i)-freqRatio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    plot(r2_Intact_mean(i),r3_Intact_mean(i),'sk','markerfacecolor',cmap_Aratio(color_nr_steady,:),'markersize',6)
    if clip_type_mean(i) > 1.5
        plot(r2_Clipped_mean(i),r3_Clipped_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
        plot(r2_Clipped_mean(i),r3_Clipped_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
    hold on
end
axis equal
% axis square
axis tight
axis([r2_min,r2_max,r3_min,r3_max])
set(gca,'xtick',r2_min:(r2_max-r2_min):r2_max)
set(gca,'ytick',r3_min:(r3_max-r3_min):r3_max)
xlabel('r_2')
ylabel('r_3')
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])

subplot(2,2,1)
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])
h = colorbar('location','northoutside'); 
title(h,'freq/freq_s_t_e_a_d_y')
set(h,'xtick',freqRatio_min:(freqRatio_max-freqRatio_min)/2:freqRatio_max)

% wing AR
subplot(2,2,4)

for i = 1:length(AR_Clipped_mean)
    color_nr_clipped = round(99/(AR_max-AR_min)*(AR_Clipped_mean(i)-AR_min)+1);
    color_nr_intact = round(99/(AR_max-AR_min)*(AR_Intact_mean(i)-AR_min)+1);
    
    if color_nr_clipped<1
        color_nr_clipped=1
    elseif color_nr_clipped>size(cmap_Aratio,1)
        color_nr_clipped=size(cmap_Aratio,1)
    end
    
    if color_nr_intact<1
        color_nr_intact=1
    elseif color_nr_intact>size(cmap_Aratio,1)
        color_nr_intact=size(cmap_Aratio,1)
    end
    
    plot(r2_Intact_mean(i),r3_Intact_mean(i),'sk','markerfacecolor',cmap_Aratio(color_nr_intact,:),'markersize',6)
    if clip_type_mean(i) > 1.5
        plot(r2_Clipped_mean(i),r3_Clipped_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr_clipped,:),'markersize',8)
    else
        plot(r2_Clipped_mean(i),r3_Clipped_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr_clipped,:),'markersize',8)
    end
    hold on
end
axis equal
% axis square
axis tight
axis([r2_min,r2_max,r3_min,r3_max])
set(gca,'xtick',r2_min:(r2_max-r2_min):r2_max)
set(gca,'ytick',r3_min:(r3_max-r3_min):r3_max)
xlabel('r_2')
ylabel('r_3')
colormap(cmap_Aratio)
caxis([AR_min AR_max])

subplot(2,2,2)
colormap(cmap_Aratio)
caxis([AR_min AR_max])
h = colorbar('location','northoutside'); 
title(h,'Wing AR')
set(h,'xtick',AR_min:(AR_max-AR_min):AR_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_r2nr3_vs_WBfreqnAR.fig'])
saveas(gcf,['clippedfly_steadyWBkin_r2nr3_vs_WBfreqnAR.png'])
% saveas(gcf,['clippedfly_steadyWBkin_r2nr3_vs_WBfreqnAR.svg'])
plot2svg(['clippedfly_steadyWBkin_r2nr3_vs_WBfreqnAR.svg'])

%% S2 & S3 NON-normailzed VS freq_clippedfly & WingArea
figure

% freq
subplot(2,2,3)
color_nr_steady = round(99/(freqRatio_max-freqRatio_min)*(1-freqRatio_min)+1);

for i = 1:length(freqRatio_mean)
    color_nr = round(99/(freqRatio_max-freqRatio_min)*(freqRatio_mean(i)-freqRatio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    plot(S2_Intact_mean(i),S3_Intact_mean(i),'sk','markerfacecolor',cmap_Aratio(color_nr_steady,:),'markersize',3)
    if clip_type_mean(i) > 1.5
        plot(S2_Clipped_mean(i),S3_Clipped_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2_Clipped_mean(i),S3_Clipped_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    hold on
end
% axis equal
axis square
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min):S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min):S3_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])

subplot(2,2,1)
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])
h = colorbar('location','northoutside'); 
title(h,'freq/freq_s_t_e_a_d_y')
set(h,'xtick',freqRatio_min:(freqRatio_max-freqRatio_min)/2:freqRatio_max)

% wing area
subplot(2,2,4)

for i = 1:length(Area_Clipped_mean)
    color_nr_clipped = round(99/(Area_max-Area_min)*(Area_Clipped_mean(i)-Area_min)+1);
    color_nr_intact = round(99/(Area_max-Area_min)*(Area_Intact_mean(i)-Area_min)+1);
    
    if color_nr_clipped<1
        color_nr_clipped=1
    elseif color_nr_clipped>size(cmap_Aratio,1)
        color_nr_clipped=size(cmap_Aratio,1)
    end
    
    if color_nr_intact<1
        color_nr_intact=1
    elseif color_nr_intact>size(cmap_Aratio,1)
        color_nr_intact=size(cmap_Aratio,1)
    end
    
    plot(S2_Intact_mean(i),S3_Intact_mean(i),'sk','markerfacecolor',cmap_Aratio(color_nr_intact,:),'markersize',3)
    if clip_type_mean(i) > 1.5
        plot(S2_Clipped_mean(i),S3_Clipped_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr_clipped,:),'markersize',5)
    else
        plot(S2_Clipped_mean(i),S3_Clipped_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr_clipped,:),'markersize',5)
    end
    hold on
end
% axis equal
axis square
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min):S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min):S3_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([Area_min Area_max])

subplot(2,2,2)
colormap(cmap_Aratio)
caxis([Area_min Area_max])
h = colorbar('location','northoutside'); 
title(h,'Wing area')
set(h,'xtick',Area_min:(Area_max-Area_min):Area_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_S2nS3_vs_WBfreqnArea.fig'])
saveas(gcf,['clippedfly_steadyWBkin_S2nS3_vs_WBfreqnArea.png'])
% saveas(gcf,['clippedfly_steadyWBkin_S2nS3_vs_WBfreqnArea.svg'])
plot2svg(['clippedfly_steadyWBkin_S2nS3_vs_WBfreqnArea.svg'])


%% stroke amplitude ratios & Freq: measured, model & MODvals
Aratio_max = 1.3;
Aratio_min = 1;

% Ad/Ai model(x) VS measurement(y) VS MODvals(color)
figure
% subplot(2,2,1)
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(AstrokeRatio_S2S3AmpRatioFunc_all(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
%         errorbar(S2S3AmpRatioFunc_all(i),Astroke_ratio_clip_intact_mean(i),Astroke_ratio_clip_intact_ste(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
        plot(S2S3AmpRatioFunc_all(i),Astroke_ratio_clip_intact_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
%         errorbar(S2S3AmpRatioFunc_all(i),Astroke_ratio_clip_intact_mean(i),Astroke_ratio_clip_intact_ste(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
        plot(S2S3AmpRatioFunc_all(i),Astroke_ratio_clip_intact_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ad/Ai model')
ylabel('Ad/Ai real fly')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Ai MOD')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_realFly_MODfunc.fig'])
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_realFly_MODfunc.png'])
% saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_realFly_MODfunc2.svg'])
plot2svg(['clippedfly_steadyWBkin_strokeAmpRatio_model_realFly_MODfunc.svg'])

% Ad/Ai model(x) VS MODvals(y) VS measurement(color)
figure
% subplot(2,2,2)
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_clip_intact_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2S3AmpRatioFunc_all(i),AstrokeRatio_S2S3AmpRatioFunc_all(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2S3AmpRatioFunc_all(i),AstrokeRatio_S2S3AmpRatioFunc_all(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ad/Ai model')
ylabel('Ad/Ai MOD')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Ai real fly')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_MODfunc_realFly.fig'])
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_MODfunc_realFly.png'])
% saveas(gcf,['clippedfly_steadyWBkin_strokeAmpRatio_model_MODfunc_realFly.svg'])
plot2svg(['clippedfly_steadyWBkin_strokeAmpRatio_model_MODfunc_realFly.svg'])

%% stroke amplitude clipped & intact wing: measured, model & MODvals
Aratio_max = 1.2;
Aratio_min = .8;
figure

% Ad/Asteady clipped model(x) VS measurement(y) VS MODvals(color)
subplot(2,2,1)
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(AstrokeRatioClipped_S2S3AmpRatioFunc_all(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        errorbar(S2S3AmpFuncClipped_all(i),Astroke_ratio_clip_mean(i),Astroke_ratio_clip_ste(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        errorbar(S2S3AmpFuncClipped_all(i),Astroke_ratio_clip_mean(i),Astroke_ratio_clip_ste(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ad/Asteady model')
ylabel('Ad/Asteady real fly')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Asteady MOD')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% Ad/Asteady clipped model(x) VS MODvals(y) VS measurement(color)
subplot(2,2,2)
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_clip_intact_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2S3AmpFuncClipped_all(i),AstrokeRatioClipped_S2S3AmpRatioFunc_all(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2S3AmpFuncClipped_all(i),AstrokeRatioClipped_S2S3AmpRatioFunc_all(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ad/Asteady model')
ylabel('Ad/Asteady MOD')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Asteady real fly')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)


% Ai/Asteady clipped model(x) VS measurement(y) VS MODvals(color)
Aratio_max = 1.2;
Aratio_min = .8;

subplot(2,2,3)
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(AstrokeRatioIntact_S2S3AmpRatioFunc_all(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        errorbar(S2S3AmpFuncIntact_all(i),Astroke_ratio_intact_mean(i),Astroke_ratio_intact_ste(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        errorbar(S2S3AmpFuncIntact_all(i),Astroke_ratio_intact_mean(i),Astroke_ratio_intact_ste(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ai/Asteady model')
ylabel('Ai/Asteady real fly')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ai/Asteady MOD')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% Ai/Asteady clipped model(x) VS MODvals(y) VS measurement(color)
subplot(2,2,4)
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_intact_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2S3AmpFuncIntact_all(i),AstrokeRatioIntact_S2S3AmpRatioFunc_all(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2S3AmpFuncIntact_all(i),AstrokeRatioIntact_S2S3AmpRatioFunc_all(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
plot([Aratio_min,Aratio_max],[Aratio_min,Aratio_max],'k')
axis equal
axis tight
axis([Aratio_min,Aratio_max,Aratio_min,Aratio_max])
set(gca,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
set(gca,'ytick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)
xlabel('Ai/Asteady model')
ylabel('Ai/Asteady MOD')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ai/Asteady real fly')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpNorm_model_MODfunc_realFly.fig'])
saveas(gcf,['clippedfly_steadyWBkin_strokeAmpNorm_model_MODfunc_realFly.png'])
% saveas(gcf,['clippedfly_steadyWBkin_strokeAmpNorm_model_MODfunc_realFly.svg'])
plot2svg(['clippedfly_steadyWBkin_strokeAmpNorm_model_MODfunc_realFly.svg'])


%% stroke amplitude ratios & wingbeat freq
Aratio_max = 1.2;
Aratio_min = .8;

% Ad/Asteady
figure
subplot(2,2,1)
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_clip_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    hold on
end
axis equal
axis tight
axis([S2norm_min,S2norm_max,S3norm_min,S3norm_max])
set(gca,'xtick',S2norm_min:(S2norm_max-S2norm_min)/2:S2norm_max)
set(gca,'ytick',S3norm_min:(S3norm_max-S3norm_min)/2:S3norm_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Asteady')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% Ai/Asteady
subplot(2,2,2)
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_intact_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    hold on
end
axis equal
axis tight
axis([S2norm_min,S2norm_max,S3norm_min,S3norm_max])
set(gca,'xtick',S2norm_min:(S2norm_max-S2norm_min)/2:S2norm_max)
set(gca,'ytick',S3norm_min:(S3norm_max-S3norm_min)/2:S3norm_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ai/Asteady')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% Ad/Ai
subplot(2,2,3)
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(Aratio_max-Aratio_min)*(Astroke_ratio_clip_intact_mean(i)-Aratio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    
    hold on
end
axis equal
axis tight
axis([S2norm_min,S2norm_max,S3norm_min,S3norm_max])
set(gca,'xtick',S2norm_min:(S2norm_max-S2norm_min)/2:S2norm_max)
set(gca,'ytick',S3norm_min:(S3norm_max-S3norm_min)/2:S3norm_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([Aratio_min Aratio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Ai')
set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min)/2:Aratio_max)

% freq_clippedfly/freq_steady
subplot(2,2,4)
for i = 1:length(freqRatio_mean)
    color_nr = round(99/(freqRatio_max-freqRatio_min)*(freqRatio_mean(i)-freqRatio_min)+1);
    if color_nr<1
        color_nr=1
    elseif color_nr>size(cmap_Aratio,1)
        color_nr=size(cmap_Aratio,1)
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    end
    hold on
end
axis equal
axis tight
axis([S2norm_min,S2norm_max,S3norm_min,S3norm_max])
set(gca,'xtick',S2norm_min:(S2norm_max-S2norm_min)/2:S2norm_max)
set(gca,'ytick',S3norm_min:(S3norm_max-S3norm_min)/2:S3norm_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])
h = colorbar('location','northoutside'); 
title(h,'freq/freq_s_t_e_a_d_y')
set(h,'xtick',freqRatio_min:(freqRatio_max-freqRatio_min)/2:freqRatio_max)

% save fig
saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.fig'])
saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.png'])
% saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.svg'])
plot2svg(['clippedfly_steadyWBkin_strokeAmp_WBfreq.svg'])


cd ..

%% save data
save('WBdataset_ClipNintact_wingbeat_kin_S2S3AmpRatioFunc.mat','S2_ratio_mean','S3_ratio_mean',...
    'freq_mean','freqRatio_mean','clip_type_mean',...
    'Astroke_clip_mean','Astroke_intact_mean',...
    'Astroke_ratio_clip_mean','Astroke_ratio_intact_mean','Astroke_ratio_clip_intact_mean',...
    'Astroke_ratio_clip_n','Astroke_ratio_intact_n','Astroke_ratio_clip_intact_n',...
    'Astroke_ratio_clip_ste','Astroke_ratio_intact_ste','Astroke_ratio_clip_intact_ste',...
    'S2S3AmpFuncClipped_all',...
    'S2S3AmpFuncIntact_all',...
    'S2S3AmpRatioFunc_all',...
    'AstrokeRatioClipped_S2S3AmpRatioFunc_all',...
    'AstrokeRatioIntact_S2S3AmpRatioFunc_all',...
    'AstrokeRatio_S2S3AmpRatioFunc_all');

