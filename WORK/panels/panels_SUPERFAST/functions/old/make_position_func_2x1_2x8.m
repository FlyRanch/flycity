% make_func_script.m
% make and save position function files

func_path = 'C:\Users\florian\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

% 9 steps up then 10 steps constant
func = [1 1 8 8];
save([ 'position_function_2x1_2x8.mat'], 'func');
