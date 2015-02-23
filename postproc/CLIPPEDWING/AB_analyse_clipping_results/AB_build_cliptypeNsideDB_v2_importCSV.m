clear
clc

% info
clip_typeNside.info.type{1,1}='distal_tip';
clip_typeNside.info.type{2,1}='trailing_edge';

clip_typeNside.info.side{1,1}='left';
clip_typeNside.info.side{2,1}='right';

% T=readtable('clip_typeNside.csv','ReadVariableNames',false)
T=importdata('clip_typeNside_v2.csv');

for i=1:size(T.textdata,1)
    
    date_now = T.textdata{i,1};
    date_now = str2double([date_now(1:4) date_now(6:7) date_now(9:10)]);
    clip_typeNside.date(i,1) = date_now;
end

clip_typeNside.seq = T.data(:,1);
clip_typeNside.clip_type = T.data(:,2);
clip_typeNside.clip_side = T.data(:,3);
    
save('clip_typeNsideDB.mat','clip_typeNside')