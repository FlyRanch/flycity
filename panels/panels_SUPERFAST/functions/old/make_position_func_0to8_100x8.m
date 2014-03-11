% make_func_script.m
% make and save position function files

func_path = 'C:\Users\florian\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

% 9 steps up then 10 steps constant
func = [[0:8] 8*ones(1,100)];
save([ 'position_function_0to8_100x8.mat'], 'func');
