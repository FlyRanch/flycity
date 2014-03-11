
addpath('/home/matt/Desktop/MOD');
addpath('/home/matt/Dropbox/Analyze/');

load('WBdataset_steady_1807WBs.mat','stroke_wb*','pitch_wb*','dev_wb*');
load('WBdataset_steady_1807WBs.mat');
load('WBmod_PitchAccel_641WBs.mat');
load('WBmod_RollAccel_625WBs.mat');
load('WBmod_YawAccel_672WBs.mat');


MODval_90deg_rollmax = 0.5;
MODval_90deg_yawmax = -0.5;
MODval_90deg_pitchmax = 1;
MODval_90deg_forcemax = 1;


stroke_wb_90deg_rollmax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax * DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1);
stroke_wb_90deg_yawmax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax * DstrokeMOD_wb_YawAccel_bins_meanCIstd(:,1);
% stroke_up_90deg_pitchmax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax * DstrokeMOD_wb_PitchAccel_bins_meanCIstd(:,1);
% stroke_up_90deg_forcemax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax * strokeMOD_us_PitchAccel_bins_meanCIstd(:,1);

dev_wb_90deg_rollmax = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax * DdevMOD_wb_RollAccel_bins_meanCIstd(:,1);
dev_wb_90deg_yawmax = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax * DdevMOD_wb_YawAccel_bins_meanCIstd(:,1);
% dev_up_90deg_pitchmax = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax * DdevMOD_wb_PitchAccel_bins_meanCIstd(:,1);
% stroke_up_90deg_forcemax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax * strokeMOD_us_PitchAccel_bins_meanCIstd(:,1);

pitch_wb_90deg_rollmax = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax * DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1);
pitch_wb_90deg_yawmax = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax * DpitchMOD_wb_YawAccel_bins_meanCIstd(:,1);
% pitch_up_90deg_pitchmax = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax * DpitchMOD_wb_PitchAccel_bins_meanCIstd(:,1);
% stroke_up_90deg_forcemax = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax * strokeMOD_us_PitchAccel_bins_meanCIstd(:,1);



% figure(1)
% hold on
% plot(stroke_wb_steady_bins_meanCIstd(:,1),'-r')
% 
% plot(stroke_up_90deg_rollmax,'-b')
% plot(stroke_up_90deg_yawmax,'-g')
% % plot(stroke_up_90deg_pitchmax,'-y')
% 
% hold off
% 
% figure(2)
% hold on
% plot(dev_wb_steady_bins_meanCIstd(:,1),'-r')
% 
% plot(dev_up_90deg_rollmax,'-b')
% plot(dev_up_90deg_yawmax,'-g')
% % plot(dev_up_90deg_pitchmax,'-y')
% 
% hold off
% 
% figure(3)
% hold on
% plot(pitch_wb_steady_bins_meanCIstd(:,1),'-r')
% 
% plot(pitch_up_90deg_rollmax,'-b')
% plot(pitch_up_90deg_yawmax,'-g')
% % plot(pitch_up_90deg_pitchmax,'-y')
% 
% hold off

clearvars -except stroke_wb_90deg_rollmax stroke_wb_90deg_yawmax dev_wb_90deg_rollmax dev_wb_90deg_yawmax pitch_wb_90deg_rollmax pitch_wb_90deg_yawmax
