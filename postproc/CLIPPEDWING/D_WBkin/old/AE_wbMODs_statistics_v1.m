clear
clc

load('WBdataset_all_1252WBs_S2nS3.mat')

f_wb = mean([f_wb_L,f_wb_R]')';

RS2_steady = SecondMomentRatio(steady_nr_mean_wb==1);
RS3_steady = ThirdMomentRatio(steady_nr_mean_wb==1);
f_steady = f_wb(steady_nr_mean_wb==1);
pitch_steady = pitch_mean_wb(steady_nr_mean_wb==1);
roll_steady = roll_mean_wb(steady_nr_mean_wb==1);
slip_steady = slip_mean_wb(steady_nr_mean_wb==1);


RS2_steady_unique = unique(RS2_steady);
for i = 1:length(RS2_steady_unique)
    n_now = find(RS2_steady==RS2_steady_unique(i));
    
    RS2_mean(i,1) = mean(RS2_steady(n_now));
    RS3_mean(i,1) = mean(RS3_steady(n_now));
    
    f_mean(i,1) = mean(f_steady(n_now));
    pitch_mean(i,1) = mean(pitch_steady(n_now));
    roll_mean(i,1) = mean(roll_steady(n_now));
    slip_mean(i,1) = mean(slip_steady(n_now));
end