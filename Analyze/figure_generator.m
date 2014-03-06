%% plot all parameters in seqs
clear A B C D E F

x = [1:200]';

for i = 1:38 % wingbeats

for j = [86,87,89,91,93,95,97,99] % individual
% for j = [1:29,42:44,54:56,61:63,77:79] % zero
    
%     if max(abs(stroke_wb_L_MATT_bins(:,i,j)-stroke_wb_R_MATT_bins(:,i,j))) > 15
    
    A(:,j) = stroke_wb_L_MATT_bins(:,i,j);
    B(:,j) = stroke_wb_R_MATT_bins(:,i,j);
    C(:,j) = dev_wb_L_MATT_bins(:,i,j);
    D(:,j) = dev_wb_R_MATT_bins(:,i,j);
    E(:,j) = pitch_wb_L_MATT_bins(:,i,j);
    F(:,j) = pitch_wb_R_MATT_bins(:,i,j);
    
%     end
    
end

A(A==0) = nan;
B(B==0) = nan;
C(C==0) = nan;
D(D==0) = nan;
E(E==0) = nan;
F(F==0) = nan;


count = 0;
for k = 1:length(A(:,1));
    count = count+1;

    stroke_L(count,1) = nanmean(A(k,:));
    stroke_L(count,2) = (stroke_L(count,1)+((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
    stroke_L(count,3) = (stroke_L(count,1)-((nanstd(A(k,:))/sqrt(length(A(k,:))))*1.96));
    
    stroke_R(count,1) = nanmean(B(k,:));
    stroke_R(count,2) = (stroke_R(count,1)+((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
    stroke_R(count,3) = (stroke_R(count,1)-((nanstd(B(k,:))/sqrt(length(B(k,:))))*1.96));
    
    dev_L(count,1) = nanmean(C(k,:));
    dev_L(count,2) = (dev_L(count,1)+((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
    dev_L(count,3) = (dev_L(count,1)-((nanstd(C(k,:))/sqrt(length(C(k,:))))*1.96));
    
    dev_R(count,1) = nanmean(D(k,:));
    dev_R(count,2) = (dev_R(count,1)+((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
    dev_R(count,3) = (dev_R(count,1)-((nanstd(D(k,:))/sqrt(length(D(k,:))))*1.96));
    
    pitch_L(count,1) = nanmean(E(k,:));
    pitch_L(count,2) = (pitch_L(count,1)+((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
    pitch_L(count,3) = (pitch_L(count,1)-((nanstd(E(k,:))/sqrt(length(E(k,:))))*1.96));
    
    pitch_R(count,1) = nanmean(F(k,:));
    pitch_R(count,2) = (pitch_R(count,1)+((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
    pitch_R(count,3) = (pitch_R(count,1)-((nanstd(F(k,:))/sqrt(length(F(k,:))))*1.96));
    
end

    figure(1)
    subplot(2,3,1)
    hold on
    plot(x, stroke_L(:,1),'-r')
    plot(x, stroke_L(:,2),'-b')
    plot(x, stroke_L(:,3), '-g')
    title('stroke - left')
    axis([0 4000 -80 100])
    hold off
    
    subplot(2,3,4)
    hold on
    plot(x, stroke_R(:,1),'-r')
    plot(x, stroke_R(:,2),'-b')
    plot(x, stroke_R(:,3), '-g')
    title('stroke - right')
    axis([0 4000 -80 100])
    hold off
%     
    subplot(2,3,2)
    hold on
    plot(x, dev_L(:,1),'-r')
    plot(x, dev_L(:,2),'-b')
    plot(x, dev_L(:,3), '-g')
%     plot(0, 0, '*')
%     plot(x, dev_wb_steady_bins_meanCIstd(:,1), '-', 'color', 'black')
    title('dev - left')
    axis([0 4000 -25 25])
    hold off
    
    subplot(2,3,5)
    hold on
    plot(x, dev_R(:,1),'-r')
    plot(x, dev_R(:,2),'-b')
    plot(x, dev_R(:,3), '-g')
    title('dev - right')
    axis([0 4000 -25 25])
    hold off
   
    subplot(2,3,3)
    hold on
    plot(x, pitch_L(:,1),'-r')
    plot(x, pitch_L(:,2),'-b')
    plot(x, pitch_L(:,3), '-g')
    title('pitch - left')
    axis([0 4000 20 160])
    hold off
    
    subplot(2,3,6)
    hold on
    plot(x, pitch_R(:,1),'-r')
    plot(x, pitch_R(:,2),'-b')
    plot(x, pitch_R(:,3), '-g')
    plot(x(1),0,'*')
    title('pitch - right')
    plot(x(1),159,'*','color','red','linewidth',4)
    axis([0 4000 20 160])
    hold off
    
    x = [1:200]' + max(x);

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