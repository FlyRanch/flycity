clc
clear
close all
warning off

% list = dir('*S0*');
% for ls = 1:length(list)

list = dir;
for ls = 3:length(list)
    
    
    if list(ls).isdir == 1
        
        length(list) - ls
        dir_now = list(ls).name
        cd(dir_now)
        
        calc_BodyWing_kinematics_byJohan
    end
end
