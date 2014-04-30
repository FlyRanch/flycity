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
        
            A_BodyWing_kinematics_v1_checkFlyData
%             
%             A_Wing_kinematics_v2_checkFlyData_rootsave
%             A_Wing_kinematics_v2_checkFlyData_overwrite_rootsave
% 
%             A_Wing_kinematics_v2_checkFlyData_overwrite
    end
end
