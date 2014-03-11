clear all;
close all;
clc;

% Program to analyze Robofly data:


add_paths = {'W:/RoboFly_data_30_05_2013'};


addpath(char(add_paths(1)))

folder_names = [];

temp_dir = dir;

dir_names = {temp_dir.name};

temp_isdir = [temp_dir.isdir];

for i = 1:length(dir_names)
           
    if temp_isdir(i) == 0
       
        folder_names = [folder_names; dir_names(i)];
        
    end
    
end

Test = [];

for i = 1:length(folder_names)
    
    temp_file = load((char(folder_names(i))));
    
    temp_name = char(folder_names(i));
    
    pos1 = findstr(temp_name,'_');
    
    pos2 = findstr(temp_name,'.');
    
    temp_nr = str2num(temp_name((pos1(end)+1):(pos2-1)));
    
    Test.([char(temp_file.exp_type) '_' int2str(temp_nr+1) ]) = temp_file;
    
end

Test_names = fieldnames(Test);

for i = 1:length(Test_names)
    
    Kine_temp = Test.(char(Test_names(i))).kine;
    
    F_temp = Test.(char(Test_names(i))).ft;
    
    t_temp = Test.(char(Test_names(i))).t;
    
    figure()
    plot(t_temp,F_temp)
    
    figure()
    plot(t_temp,Kine_temp)
    
end