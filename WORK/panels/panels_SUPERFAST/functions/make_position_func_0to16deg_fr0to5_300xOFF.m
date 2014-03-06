% make_func_script.m
% make and save position function files

func_path = 'E:\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

% 9 steps up then 10 steps constant
func = [[0:5] zeros(1,300)];
save([ 'position_function_0to16deg_fr0to5_300xOFF.mat'], 'func');
