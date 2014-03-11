% cd('/home/matt/Desktop/MOD')
addpath('/home/matt/Desktop/MOD');
addpath('D:\Dropbox\Michael\WBmod_data\looming_resp');

load('WBdataset_steady_1807WBs.mat');
load('WBmod_PitchAccel_641WBs.mat');
load('WBmod_RollAccel_625WBs.mat');
load('WBmod_YawAccel_672WBs.mat');
load('WBmod_Fenhanced_806WBs.mat');


MODval_90deg_rollmax = [-5:0.1:5]; % 0.5
MODval_90deg_yawmax = [-5:0.1:5]; % -0.5
MODval_90deg_pitchmax = [-5:0.1:5]; % 1
MODval_90deg_forcemax = [-5:0.1:5]; % 1
MODval = [-5:0.1:5];

for i = 1:length(MODval)
stroke_wb_down_90deg_rollmax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * strokeMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
stroke_wb_up_90deg_rollmax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * strokeMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
stroke_wb_fwd_90deg_yawmax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
stroke_wb_rwd_90deg_yawmax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
stroke_wb_90deg_pitchmax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax(i) * strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1);
stroke_wb_90deg_forcemax(:,i) = stroke_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax(i) * strokeMOD_wb_Fenhance_bins_meanCIstd(:,1);

dev_wb_down_90deg_rollmax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * devMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
dev_wb_up_90deg_rollmax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * devMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
dev_wb_fwd_90deg_yawmax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * devMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
dev_wb_rwd_90deg_yawmax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * devMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
dev_wb_90deg_pitchmax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax(i) * devMOD_wb_PitchAccel_bins_meanCIstd(:,1);
dev_wb_90deg_forcemax(:,i) = dev_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax(i) * devMOD_wb_Fenhance_bins_meanCIstd(:,1);

pitch_wb_down_90deg_rollmax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * pitchMOD_wb_down_RollAccel_bins_meanCIstd(:,1);
pitch_wb_up_90deg_rollmax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_rollmax(i) * pitchMOD_wb_up_RollAccel_bins_meanCIstd(:,1);
pitch_wb_fwd_90deg_yawmax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * pitchMOD_wb_fwd_YawAccel_bins_meanCIstd(:,1);
pitch_wb_rwd_90deg_yawmax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_yawmax(i) * pitchMOD_wb_rwd_YawAccel_bins_meanCIstd(:,1);
pitch_wb_90deg_pitchmax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_pitchmax(i) * pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1);
pitch_wb_90deg_forcemax(:,i) = pitch_wb_steady_bins_meanCIstd(:,1) + MODval_90deg_forcemax(i) * pitchMOD_wb_Fenhance_bins_meanCIstd(:,1);
end

MOD = {};
MOD.stroke_wb_down_90deg_rollmax = stroke_wb_down_90deg_rollmax;
MOD.stroke_wb_up_90deg_rollmax = stroke_wb_up_90deg_rollmax;
MOD.stroke_wb_fwd_90deg_yawmax = stroke_wb_fwd_90deg_yawmax;
MOD.stroke_wb_rwd_90deg_yawmax = stroke_wb_rwd_90deg_yawmax;
MOD.stroke_wb_90deg_pitchmax = stroke_wb_90deg_pitchmax;
MOD.stroke_wb_90deg_forcemax = stroke_wb_90deg_forcemax;
MOD.dev_wb_down_90deg_rollmax = dev_wb_down_90deg_rollmax;
MOD.dev_wb_up_90deg_rollmax = dev_wb_up_90deg_rollmax;
MOD.dev_wb_fwd_90deg_yawmax = dev_wb_fwd_90deg_yawmax;
MOD.dev_wb_rwd_90deg_yawmax = dev_wb_rwd_90deg_yawmax;
MOD.dev_wb_90deg_pitchmax = dev_wb_90deg_pitchmax;
MOD.dev_wb_90deg_forcemax = dev_wb_90deg_forcemax;
MOD.pitch_wb_down_90deg_rollmax = pitch_wb_down_90deg_rollmax;
MOD.pitch_wb_up_90deg_rollmax = pitch_wb_up_90deg_rollmax;
MOD.pitch_wb_fwd_90deg_yawmax = pitch_wb_fwd_90deg_yawmax;
MOD.pitch_wb_rwd_90deg_yawmax = pitch_wb_rwd_90deg_yawmax;
MOD.pitch_wb_90deg_pitchmax = pitch_wb_90deg_pitchmax;
MOD.pitch_wb_90deg_forcemax = pitch_wb_90deg_forcemax;
MOD.MODval = MODval;

clearvars -except MOD
