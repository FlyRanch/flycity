clear;
clc;
close all
warning off


freq_asymFitNr = 10;
peakloc = 3;

%% colormap for A+: reds
% cmap_Reds =cbrewer('seq','Reds',100);
% cmap_AdAi = cmap_Reds;
% 
% cmap_YlOrRd =cbrewer('seq','YlOrRd',100);
% cmap_AdAi = cmap_YlOrRd;

cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_AdAi = cmap_hot;

if peakloc == 0
    Arot_min = 90;
    Arot_max = 135;

    dM_min = .1;
    dM_max = .6;

elseif peakloc ==2
    Arot_min = 80;
    Arot_max = 105;

    dM_min = .05;
    dM_max = .2;
    
elseif peakloc ==3
    Arot_min = 70;
    Arot_max = 95;

    dM_min = -.15;
    dM_max = .15;
end

Aplus_min = 1.0;
Aplus_max = 1.3;

        %% plot TIP cut wing data
    load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.mat'])

    S2ratios_TEnTip = S2ratios;
    S3ratios_TEnTip = S3ratios;
    S2S3AmpRatioFuncs_TEnTip = S2S3AmpRatioFuncs;
    
    Arot_damagedwing_all_TEnTip = Arot_damagedwing_all;
    dMwing_total_damagedwing_all_TEnTip = dMwing_total_damagedwing_all;
    
    Arot_intactwing_all_TEnTip = Arot_intactwing_all;
    dMwing_total_intactwing_all_TEnTip = dMwing_total_intactwing_all;
    
    
    % rotation angle amplitude VS change in spanwise torque
%     figure
%     subplot(2,2,1)
    hold on

    plot(Arot_intactwing_all,dMwing_total_intactwing_all,'-','color',[.5 .5 .5])
    plot(Arot_damagedwing_all,dMwing_total_damagedwing_all,'-','color',[.5 .5 .5])

    for i = 1:length(S2S3AmpRatioFuncs)
        color_nr = round(99/(Aplus_max-Aplus_min)*(S2S3AmpRatioFuncs(i)-Aplus_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_AdAi,1)
            color_nr=size(cmap_AdAi,1)
        end

        plot(Arot_intactwing_all(i),dMwing_total_intactwing_all(i),'sk-','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
        plot(Arot_damagedwing_all(i),dMwing_total_damagedwing_all(i),'ok-','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',10)
    end


%     axis tight
%     axis square
%     axis([Arot_min,Arot_max,dM_min,dM_max])
%     set(gca,'xtick',Arot_min:(Arot_max-Arot_min)/2:Arot_max)
%     set(gca,'ytick',dM_min:(dM_max-dM_min)/2:dM_max)
%     xlabel('rotation angle amplitude [deg]')
%     ylabel('normalized spanwise torque change')
%     colormap(cmap_AdAi)
%     caxis([Aplus_min Aplus_max])
%     h = colorbar('location','northoutside'); 
%     title(h,'Ad/Ai MOD')
%     set(h,'xtick',Aplus_min:(Aplus_max-Aplus_min)/2:Aplus_max)
    
    %% plot TE cut wing data
    load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.mat'])

    S2ratios_TEnTip = [S2ratios_TEnTip;S2ratios];
    S3ratios_TEnTip = [S3ratios_TEnTip;S3ratios];
    S2S3AmpRatioFuncs_TEnTip = [S2S3AmpRatioFuncs_TEnTip;S2S3AmpRatioFuncs];
    
    Arot_damagedwing_all_TEnTip = [Arot_damagedwing_all_TEnTip;Arot_damagedwing_all];
    dMwing_total_damagedwing_all_TEnTip = [dMwing_total_damagedwing_all_TEnTip;dMwing_total_damagedwing_all];
    
    Arot_intactwing_all_TEnTip = [Arot_intactwing_all_TEnTip;Arot_intactwing_all];
    dMwing_total_intactwing_all_TEnTip = [dMwing_total_intactwing_all_TEnTip;dMwing_total_intactwing_all];
    
    
    % rotation angle amplitude VS change in spanwise torque
%     figure
%     subplot(2,2,1)
    hold on

        plot(Arot_intactwing_all,dMwing_total_intactwing_all,'k-')
        plot(Arot_damagedwing_all,dMwing_total_damagedwing_all,'k-')

    for i = 1:length(S2S3AmpRatioFuncs)
        color_nr = round(99/(Aplus_max-Aplus_min)*(S2S3AmpRatioFuncs(i)-Aplus_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_AdAi,1)
            color_nr=size(cmap_AdAi,1)
        end

        plot(Arot_intactwing_all(i),dMwing_total_intactwing_all(i),'sk-','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
        plot(Arot_damagedwing_all(i),dMwing_total_damagedwing_all(i),'dk-','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',10)
    end
    
    axis tight
    axis square
    axis([Arot_min,Arot_max,dM_min,dM_max])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min)/2:Arot_max)
    set(gca,'ytick',dM_min:(dM_max-dM_min)/2:dM_max)
    xlabel('rotation angle amplitude [deg]')
    ylabel('normalized spanwise torque change')
    colormap(cmap_AdAi)
    caxis([Aplus_min Aplus_max])
    h = colorbar('location','northoutside'); 
    title(h,'Ad/Ai MOD')
    set(h,'xtick',Aplus_min:(Aplus_max-Aplus_min)/2:Aplus_max)
    
    
%     linfit_intactwing = polyfit(Arot_intactwing_all_TEnTip,dMwing_total_intactwing_all_TEnTip,1);
%     linfit_damwing = polyfit(Arot_intactwing_all_TEnTip,dMwing_total_intactwing_all_TEnTip,1);

% save plot
mkdir('qsModel_FnM_TEnTipCut')
cd('qsModel_FnM_TEnTipCut')

saveas(gca,['Arot_vs_DMwing_TEnTipClip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.fig'])
saveas(gca,['Arot_vs_DMwing_TEnTipClip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.png'])
plot2svg(['Arot_vs_DMwing_TEnTipClip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.svg'])

cd .. 
    
    
    