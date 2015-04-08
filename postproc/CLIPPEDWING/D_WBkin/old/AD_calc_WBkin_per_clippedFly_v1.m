clc
clear
close all

plot_on = 1;
plot_on = 0;

loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

%% plot dir
if plot_on == 1
    mkdir('steadyWBkin_seqs_figs')
    cd('steadyWBkin_seqs_figs')
end

%% steady wb
freq_steady = f_wb_steady_meanCIstd;
stroke_steady = stroke_wb_steady_bins_meanCIstd;
Astroke_steady = max(stroke_steady(:,1)) - min(stroke_steady(:,1));

            
%% loop through SecondMomentRatio's
SecondMomentRatio_list = unique(SecondMomentRatio);


%% clipped fly WBs
% make nan DB arrays
S2_ratio_all = nan(150,length(SecondMomentRatio_list));
S3_ratio_all = nan(150,length(SecondMomentRatio_list));

freq_all = nan(150,length(SecondMomentRatio_list));

Astroke_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_intact_all = nan(150,length(SecondMomentRatio_list));

Astroke_ratio_clip_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_intact_all = nan(150,length(SecondMomentRatio_list));
Astroke_ratio_clip_intact_all = nan(150,length(SecondMomentRatio_list));


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

            % body kin data
            vel_now = V_mean_wb(wb);
            
            slip_now = slip_mean_wb(wb);
            roll_now = roll_mean_wb(wb);
            pitch_now = pitch_mean_wb(wb);
            
            Fsp_pitch_now = Fsp_pitch_mean_wb(wb);
            Fsp_roll_now = Fsp_roll_mean_wb(wb);
            
            % wingbeat kin data
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
            S2_ratio_all(1:length(SecondMomentRatio_now),seq_now) = SecondMomentRatio_now;
            S3_ratio_all(1:length(ThirdMomentRatio_now),seq_now) = ThirdMomentRatio_now;
            
            freq_all(1:length(f_wb_now),seq_now) = f_wb_now;
   
            Astroke_clip_all(1:length(Astroke_clip_now),seq_now) = Astroke_clip_now;
            Astroke_intact_all(1:length(Astroke_intact_now),seq_now) = Astroke_intact_now;
   
            Astroke_ratio_clip_all(1:length(Astroke_ratio_clip_now),seq_now) = Astroke_ratio_clip_now;
            Astroke_ratio_intact_all(1:length(Astroke_ratio_intact_now),seq_now) = Astroke_ratio_intact_now;
            Astroke_ratio_clip_intact_all(1:length(Astroke_ratio_clip_intact_now),seq_now) = Astroke_ratio_clip_intact_now;
            
            
            %% plot
            if plot_on == 1
                subplot(2,2,1)
                plot(stroke_clip_now,'r')
                hold on
                plot(stroke_intact_now,'b')
                plot(stroke_wb_steady_bins_meanCIstd(:,1),'g','linewidth',2)
                ylim([-90 90])
                hold off

                subplot(2,2,2)
                plot(pitch_wb_clip_now-90,'r')
                hold on
                plot(pitch_wb_intact_now-90,'b')
                plot(pitch_wb_steady_bins_meanCIstd(:,1)-90,'g','linewidth',2)
                ylim([-90 90])
                hold off

                subplot(2,2,3)
                plot(0,0,'-r')
                hold on
                plot(0,0,'-b')
                plot(0,0,'-g')
                plot(dev_wb_clip_now,'r')
                plot(dev_wb_intact_now,'b')
                plot(dev_wb_steady_bins_meanCIstd(:,1),'g','linewidth',2)
                ylim([-90 90])
                legend('clipped wing','intact wing','steady wingbeat')
                hold off

                saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_n=',num2str(length(seq_nr_now)),'.fig'])
                saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_n=',num2str(length(seq_nr_now)),'.png'])
%                 saveas(gcf,['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_n=',num2str(length(seq_nr_now)),'.svg'])
                plot2svg(['steadyWBkin_S2ratio=',num2str(round(100*S2_now)/100),'_S3ratio=',num2str(round(100*mean(ThirdMomentRatio_now))/100),'_Aratio=',num2str(round(100*mean(Astroke_ratio_clip_intact_now))/100),'_n=',num2str(length(seq_nr_now)),'.svg'])

%                 seq_nr_now
%                 seq_now
%                 mean(Astroke_ratio_clip_intact_now)
%                 S2_now
%                 pause
            end
end

if plot_on == 1
    cd ..
end

S2_ratio_mean = nanmean(S2_ratio_all)';
S3_ratio_mean = nanmean(S3_ratio_all)';

freq_mean = nanmean(freq_all)';
freqRatio_mean = freq_mean/freq_steady(1);

Astroke_clip_mean = nanmean(Astroke_clip_all)';
Astroke_intact_mean = nanmean(Astroke_intact_all)';

Astroke_ratio_clip_mean = nanmean(Astroke_ratio_clip_all)';
Astroke_ratio_intact_mean = nanmean(Astroke_ratio_intact_all)';
Astroke_ratio_clip_intact_mean = nanmean(Astroke_ratio_clip_intact_all)';


%% plot
% datapoints with color
% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);
cmap_Aratio = cmap_surf;

S2_min = .5;
S2_max = 1;
S3_min = .5;
S3_max = 1;
Aratio_max = 1.2;
Aratio_min = .8;
freqRatio_max = 1.25;
freqRatio_min = .75;

%% stroke amplitude ratios & wingbeat freq

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
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    hold on
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
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
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    hold on
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
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
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    hold on
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
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
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',5)
    hold on
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')
colormap(cmap_Aratio)
caxis([freqRatio_min freqRatio_max])
h = colorbar('location','northoutside'); 
title(h,'Ad/Ai')
set(h,'xtick',freqRatio_min:(freqRatio_max-freqRatio_min)/2:freqRatio_max)

% save fig
mkdir('FIG_clippedfly_steadyWBkin_strokeAmp_WBfreq')
cd('FIG_clippedfly_steadyWBkin_strokeAmp_WBfreq')

saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.fig'])
saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.png'])
% saveas(gcf,['clippedfly_steadyWBkin_strokeAmp_WBfreq.svg'])
plot2svg(['clippedfly_steadyWBkin_strokeAmp_WBfreq.svg'])

cd ..

%% save data
save('WBdataset_ClipNintact_wingbeat_kin.mat','S2_ratio_mean','S3_ratio_mean',...
    'freq_mean','freqRatio_mean',...
    'Astroke_clip_mean','Astroke_intact_mean',...
    'Astroke_ratio_clip_mean','Astroke_ratio_intact_mean','Astroke_ratio_clip_intact_mean');





