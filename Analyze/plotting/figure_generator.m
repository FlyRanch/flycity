% %% Construct Bins by Individual
% 
% stim_exp_90 = [30:32,48:50,51:53,64:66,74:76];
% stim_exp_270 = [36:41,58:60,70:73];
% stim_exp_0 = [11:29,42:44,54:56,61:63,77:79];
% stim_exp_180 = [33:35,45:47,57,67:70];
% 
% stim_drum_1hz = [85,93];
% stim_drum_5hz = [87,95,97];
% stim_drum_10hz = [89,91,99,102,103,104,111,112];
% 
% 
% %% Bin data and interpolate
% 
% bins = 200;
% % seqs_max = 83;
% 
% % Hack into bins of length 'bins' and interpolate
% for i = [stim_exp_90,stim_exp_270];
%     stroke_wb_L_MATT_bins(:,:,i) = bin(stroke_wb_L_MATT(:,:,i),bins); % input, number of bins
%     stroke_wb_R_MATT_bins(:,:,i) = bin(stroke_wb_R_MATT(:,:,i),bins);
%     dev_wb_R_MATT_bins(:,:,i) = bin(dev_wb_R_MATT(:,:,i),bins);
%     dev_wb_L_MATT_bins(:,:,i) = bin(dev_wb_L_MATT(:,:,i),bins);
%     pitch_wb_R_MATT_bins(:,:,i) = bin(pitch_wb_R_MATT(:,:,i),bins);
%     pitch_wb_L_MATT_bins(:,:,i) = bin(pitch_wb_L_MATT(:,:,i),bins);
% end
% 
% stroke_wb_L_MATT_bins(stroke_wb_L_MATT_bins == 0) = nan;
% stroke_wb_R_MATT_bins(stroke_wb_R_MATT_bins == 0) = nan;
% dev_wb_R_MATT_bins(dev_wb_R_MATT_bins == 0) = nan;
% dev_wb_L_MATT_bins(dev_wb_L_MATT_bins == 0) = nan;
% pitch_wb_R_MATT_bins(pitch_wb_R_MATT_bins == 0) = nan;
% pitch_wb_L_MATT_bins(pitch_wb_L_MATT_bins == 0) = nan;
% 
% for i = [stim_exp_90,stim_exp_270];
%     for j = 1:length(stroke_wb_L_MATT_bins(1,:,i))
%         if turning_def(j,7,i) > 0.15
%             
%                 if strcmp(wing{j,1,i},'right') == 1
%                     stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
%                 
%                     stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_L_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_L_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_L_MATT_bins(:,j,i);
%                 
%                 elseif strcmp(wing{j,1,i},'left') == 1
%                     stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_L_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_L_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_L_MATT_bins(:,j,i);
%                 
%                     stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
%                 end
%                 
%         elseif turning_def(j,7,i) < 0.15
%             
%                     stroke_wb_MATT_bins_up_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_up_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_up_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
%                     
%                     stroke_wb_MATT_bins_dn_turn(:,j,i) = stroke_wb_R_MATT_bins(:,j,i);
%                     dev_wb_MATT_bins_dn_turn(:,j,i) = dev_wb_R_MATT_bins(:,j,i);
%                     pitch_wb_MATT_bins_dn_turn(:,j,i) = pitch_wb_R_MATT_bins(:,j,i);
%                     
%         end
%     end
% end
% 
% stroke_wb_MATT_bins_up_turn(stroke_wb_MATT_bins_up_turn == 0) = nan;
% dev_wb_MATT_bins_up_turn(dev_wb_MATT_bins_up_turn == 0) = nan;
% pitch_wb_MATT_bins_up_turn(pitch_wb_MATT_bins_up_turn == 0) = nan;
% stroke_wb_MATT_bins_dn_turn(stroke_wb_MATT_bins_dn_turn == 0) = nan;
% dev_wb_MATT_bins_dn_turn(dev_wb_MATT_bins_dn_turn == 0) = nan;
% pitch_wb_MATT_bins_dn_turn(pitch_wb_MATT_bins_dn_turn == 0) = nan;
% 
% 
% 
% %% plot all parameters in seqs
% clear A B C D E F
% 
% 
% 
% for i = 1:length(stroke_wb_MATT_bins_up_turn(1,:,1)) % wingbeats
% 
% for j = 1:length(stroke_wb_MATT_bins_up_turn(1,1,:)) % individual
%                     
%     A(:,j) = stroke_wb_MATT_bins_up_turn(:,i,j);
%     B(:,j) = dev_wb_MATT_bins_up_turn(:,i,j);
%     C(:,j) = pitch_wb_MATT_bins_up_turn(:,i,j);
%     
%     D(:,j) = stroke_wb_MATT_bins_dn_turn(:,i,j);
%     E(:,j) = dev_wb_MATT_bins_dn_turn(:,i,j);
%     F(:,j) = pitch_wb_MATT_bins_dn_turn(:,i,j);
%     
% end
% 
% count = 0;
% for k = 1:length(A(:,1));
%     count = count+1;
% 
%     stroke_up(count,1) = nanmean(A(k,:));
%     stroke_up(count,2) = (stroke_up(count,1)+((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
%     stroke_up(count,3) = (stroke_up(count,1)-((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
%     
%     dev_up(count,1) = nanmean(B(k,:));
%     dev_up(count,2) = (dev_up(count,1)+((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
%     dev_up(count,3) = (dev_up(count,1)-((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
%     
%     pitch_up(count,1) = nanmean(C(k,:));
%     pitch_up(count,2) = (pitch_up(count,1)+((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
%     pitch_up(count,3) = (pitch_up(count,1)-((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
%     
%     stroke_dn(count,1) = nanmean(D(k,:));
%     stroke_dn(count,2) = (stroke_dn(count,1)+((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
%     stroke_dn(count,3) = (stroke_dn(count,1)-((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
%     
%     dev_dn(count,1) = nanmean(E(k,:));
%     dev_dn(count,2) = (dev_dn(count,1)+((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
%     dev_dn(count,3) = (dev_dn(count,1)-((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
%     
%     pitch_dn(count,1) = nanmean(F(k,:));
%     pitch_dn(count,2) = (pitch_dn(count,1)+((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
%     pitch_dn(count,3) = (pitch_dn(count,1)-((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
%     
% end
% 
% stroke_up_wingbeats(:,1:3,i) = stroke_up;
% stroke_dn_wingbeats(:,1:3,i) = stroke_dn;
% dev_up_wingbeats(:,1:3,i) = dev_up;
% dev_dn_wingbeats(:,1:3,i) = dev_dn;
% pitch_up_wingbeats(:,1:3,i) = pitch_up;
% pitch_dn_wingbeats(:,1:3,i) = pitch_dn;
% 
% end



%% stroke
figure(1)
title('Stroke Angle')
box on; grid off;
fr = 1/7500;
xstep = [400:200:4800];
xlabel = fr.*(xstep-200);
set(gca,'XTick',xstep,'XTickLabel',xlabel,'TickLength',[0 0],'LineWidth', 1.2);
xlim([400 4800]); ylim([-60 95]);
x = [1:200]';

for i = 1:length(stroke_dn_wingbeats(1,1,:))
    hold on
    z = x';

    y1 = stroke_up_wingbeats(:,3,i)';
    y2 = stroke_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, stroke_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, stroke_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = stroke_dn_wingbeats(:,1)';
    y1 = stroke_dn_wingbeats(:,3,i)';
    y2 = stroke_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, stroke_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, stroke_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

    title('stroke - up(red) and down(blue)')
    hold off
    
    x = [1:200]' + max(x);
end

    
%% dev
     z = x';
    figure(2)
    hold on
    set(gca,'LineWidth', 1.2);
%     mid = stroke_up_wingbeats(:,1)';
    y1 = dev_up_wingbeats(:,3,i)';
    y2 = dev_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, dev_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, dev_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = dev_dn_wingbeats(:,1)';
    y1 = dev_dn_wingbeats(:,3,i)';
    y2 = dev_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, dev_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, dev_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

    title('dev - up(red) and down(blue)')
    hold off
    
%% pitch
     z = x';
    figure(3)
    hold on
    set(gca,'LineWidth', 1.2);
%     mid = stroke_up(:,1)';
    y1 = pitch_up_wingbeats(:,3,i)';
    y2 = pitch_up_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h1 = fill(X,Y,1/255*[250,107,107],'EdgeColor','none');
     
    plot(x, pitch_up_wingbeats(:,2,i),'color',1/255*[250,107,107])
    plot(x, pitch_up_wingbeats(:,3,i),'color',1/255*[250,107,107])
    
%     mid = pitch_dn_wingbeats(:,1)';
    y1 = pitch_dn_wingbeats(:,3,i)';
    y2 = pitch_dn_wingbeats(:,2,i)';
    X = [z,fliplr(z)];
    Y = [y1,fliplr(y2)];

    h2 = fill(X,Y,1/255*[107,145,250],'EdgeColor','none');
    alpha(0.8)
    plot(x, pitch_dn_wingbeats(:,2,i),'color',1/255*[107,145,250])
    plot(x, pitch_dn_wingbeats(:,3,i),'color',1/255*[107,145,250])

    title('pitch - up(red) and down(blue)')
    hold off
    
%     figure(1)
%     subplot(2,3,1)
%     hold on
%     plot(x, stroke_up(:,1),'-r')
%     plot(x, stroke_up(:,2),'-b')
%     plot(x, stroke_up(:,3), '-g')
%     title('stroke - up wing')
% %     axis([0 4000 -80 100])
%     hold off
%     
%     subplot(2,3,2)
%     hold on
%     plot(x, dev_up(:,1),'-r')
%     plot(x, dev_up(:,2),'-b')
%     plot(x, dev_up(:,3), '-g')
%     title('dev - up wing')
% %     axis([0 4000 -80 100])
%     hold off
% %     
%     subplot(2,3,3)
%     hold on
%     plot(x, pitch_up(:,1),'-r')
%     plot(x, pitch_up(:,2),'-b')
%     plot(x, pitch_up(:,3), '-g')
% %     plot(0, 0, '*')
% %     plot(x, dev_wb_steady_bins_meanCIstd(:,1), '-', 'color', 'black')
%     title('pitch - up wing')
% %     axis([0 4000 -25 25])
%     hold off
%     
%     subplot(2,3,4)
%     hold on
%     plot(x, stroke_dn(:,1),'-r')
%     plot(x, stroke_dn(:,2),'-b')
%     plot(x, stroke_dn(:,3), '-g')
%     title('stroke - down wing')
% %     axis([0 4000 -25 25])
%     hold off
%    
%     subplot(2,3,5)
%     hold on
%     plot(x, dev_dn(:,1),'-r')
%     plot(x, dev_dn(:,2),'-b')
%     plot(x, dev_dn(:,3), '-g')
%     title('dev - down wing')
% %     axis([0 4000 20 160])
%     hold off
%     
%     subplot(2,3,6)
%     hold on
%     plot(x, pitch_dn(:,1),'-r')
%     plot(x, pitch_dn(:,2),'-b')
%     plot(x, pitch_dn(:,3), '-g')
% %     plot(x(1),0,'*')
%     title('pitch - down wing')
% %     plot(x(1),159,'*','color','red','linewidth',4)
% %     axis([0 4000 20 160])
%     hold off
    
    

end

% figure(1)
% 
% subplot(2,3,1)
% hold on
% axis([0 200 -70 100]);
% A = plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% % plot(stroke_wb_steady_bins_meanCIstd(:,1:3), '-k')
% B = plot(stroke_wb_L_MATT_bins_mean_steady,'-','color','blue','linewidth', 1.2);
% C = plot(stroke_wb_L_MATT_bins_mean_turning_right,'-','color','red','linewidth', 1.2);
% D = plot(stroke_wb_L_MATT_bins_mean_turning_left,'-','color','yellow','linewidth', 1.2);
% title('stroke - left wing')
% legend([A(1) B(1) C(1) D(1)],'free flight','tethered flight','tethered - turning right','tethered - turning left')
% hold off
% 
% subplot(2,3,4)
% hold on
% axis([0 200 -70 100]);
% plot(stroke_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% plot(stroke_wb_R_MATT_bins_mean_turning_right,'-r', 'linewidth', 1.2)
% plot(stroke_wb_R_MATT_bins_mean_turning_left,'-y', 'linewidth', 1.2)
% plot(stroke_wb_R_MATT_bins_mean_steady,'-b', 'linewidth', 1.2)
% title('stroke - right wing')
% hold off
% 
% subplot(2,3,2)
% hold on
% axis([0 200 -20 25]);
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% plot(dev_wb_L_MATT_bins_mean_turning_right,'-r')
% plot(dev_wb_L_MATT_bins_mean_turning_left,'-y')
% plot(dev_wb_L_MATT_bins_mean_steady,'-b')
% title('deviation - left wing')
% hold off
% 
% subplot(2,3,5)
% hold on
% axis([0 200 -20 25]);
% plot(dev_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% plot(dev_wb_R_MATT_bins_mean_turning_right,'-r')
% plot(dev_wb_R_MATT_bins_mean_turning_left,'-y')
% plot(dev_wb_R_MATT_bins_mean_steady,'-b')
% title('deviation - right wing')
% hold off
% % 
% subplot(2,3,3)
% hold on
% axis([0 200 10 160]);
% plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% plot(pitch_wb_L_MATT_bins_mean_turning_right,'-r')
% plot(pitch_wb_L_MATT_bins_mean_turning_left,'-y')
% plot(pitch_wb_L_MATT_bins_mean_steady,'-b')
% title('pitch - left wing')
% hold off
% % 
% subplot(2,3,6)
% hold on
% axis([0 200 10 160]);
% plot(pitch_wb_steady_bins_meanCIstd(:,1:3),'-', 'color','green','linewidth', 1.2);
% plot(pitch_wb_R_MATT_bins_mean_turning_right,'-r')
% plot(pitch_wb_R_MATT_bins_mean_turning_left,'-y')
% plot(pitch_wb_R_MATT_bins_mean_steady,'-b')
% title('pitch - right wing')
% hold off