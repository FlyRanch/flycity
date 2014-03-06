% set working dir to /media/""/solutions/

clear all; close all; clc;

addpath('/home/matt/Dropbox/WORK/flytracker/');

% homedirectory = '/media/Matt01/solutions/';
homedirectory = cd;
cd(homedirectory)

list = dir('20*');

% matlabpool(8)

% for i = [95,101,104,107,110,113,116]
for i = 1:length(list)
    if list(i).isdir == 1
        cd(list(i).name)
        currentfolder = pwd % display folder change
    
    list2 = dir();
    list2(1:2) = [];

        if exist('flytracks') == 7
            cd flytracks
            currentfolder = pwd
            
                ls = dir('fly0*.mat');
                if isempty(ls) == 0
%                 if exist('fly000800.mat') == 2
                    paste_top_projection_makemovie_BATCH
                end
        end
    end
    cd(homedirectory)
 end

% matlabpool close
    