% make_func_script.m
% make and save position function files

func_path = 'E:\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

% 9 steps up then 10 steps constant
func = [1 17*ones(1,300)];
save([ 'position_function_const64deg_fr1_300xfr17.mat'], 'func');
