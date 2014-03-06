% make_func_script.m
% make and save position function files

func_path = 'E:\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

func = [15:-1:0];
save([ 'position_function_drumrotation_fr15to0_REVERSE.mat'], 'func');
