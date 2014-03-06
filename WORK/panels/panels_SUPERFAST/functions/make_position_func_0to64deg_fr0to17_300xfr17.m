% make_func_script.m
% make and save position function files

func_path = 'C:\Users\florian\Dropbox\WORK\panels\panels_SUPERFAST\functions';
cd(func_path);

% 9 steps up then 10 steps constant
func = [[0:17] 17*ones(1,300)];
save([ 'position_function_0to64deg_fr0to17_300xfr17.mat'], 'func');
