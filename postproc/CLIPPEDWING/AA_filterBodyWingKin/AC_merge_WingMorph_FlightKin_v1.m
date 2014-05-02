clc
clear
close all
warning off

plot_on = 0;

load('WingClipDatabase.mat')

%% TrailingEdgeClip
cd('TrailingEdgeClip')
clip = 1

list = dir('*.mat');
for i = 1:length(list)
        
        length(list) - i
        
        % load data
        pathDB_now = load(name_now);
        
        % seq info
        name_now = list(i).name;
        date_now = str2num(name_now(1:8));
        seq_now = str2num(name_now(11:14));
        
        sub_now = str2num(name_now(16));
        if isempty(sub_now)
            sub_now = 1;
        end
        
        PathDB.(i,:)
        
        
        
        
        
        
        
        cd ..
    end
end

cd ..
save('WingClipDatabase.mat','WingClipData')